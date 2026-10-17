nextflow.enable.dsl = 2

process PREREGISTER {
  publishDir "${params.outdir}/provenance", mode:'copy', overwrite:true
  input: path expected; path sources; path coords
  output: path 'provenance', emit: provenance
  script:
  """
  mkdir provenance
  cp $expected provenance/expected_results.tsv
  cp $sources provenance/sources.tsv
  cp $coords provenance/coordinates.tsv
  sha256sum provenance/* > provenance/sha256sums.txt
  """
}

process SLICE_DNA {
  tag "$sample"
  publishDir "${params.outdir}/dna_slices", mode:'copy', overwrite:true
  input: tuple val(sample), val(url), val(index_url), val(kind)
  output: tuple val(sample), path("${sample}.bam"), emit: bams
  script:
  def idx = kind == 'cram' ? "${sample}.cram.crai" : "${sample}.bam.bai"
  def refarg = kind == 'cram' ? "-T ${params.grch38_reference}" : ""
  """
  curl -L --fail --retry 5 '$index_url' -o '$idx'
  samtools view -@ ${task.cpus} -P -X -b $refarg -o ${sample}.bam '$url' '$idx' chr13:28033000-28035000
  samtools quickcheck -v ${sample}.bam
  """
}

process CALL_ITD {
  tag "$sample"
  publishDir "${params.outdir}/itd", mode:'copy', overwrite:true
  input: tuple val(sample), path(bam); path ref; path anno
  output: tuple val(sample), path("${sample}_getitd"), emit: calls
  script:
  """
  samtools sort -n -@ ${task.cpus} -o namesort.bam $bam
  samtools fastq namesort.bam | gzip -c > ${sample}.fastq.gz
  getitd $sample ${sample}.fastq.gz -reference $ref -anno $anno \
    -require_indel_free_primers False -infer_sense_from_alignment True \
    -min_read_copies 1 -nkern ${task.cpus} -plot_coverage False
  """
}

process SLICE_RNA {
  tag "$sample"
  publishDir "${params.outdir}/rna_slices", mode:'copy', overwrite:true
  input: tuple val(sample), val(url), val(index_url), val(regions)
  output: tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), emit: reads
  script:
  """
  curl -L --fail --retry 5 '$index_url' -o source.bam.bai
  samtools view -@ ${task.cpus} -P -X -b -o partners.bam '$url' source.bam.bai $regions
  samtools sort -n -@ ${task.cpus} partners.bam -o namesort.bam
  samtools fastq -1 ${sample}_R1.fastq -2 ${sample}_R2.fastq \
    -0 /dev/null -s /dev/null -n namesort.bam
  gzip ${sample}_R1.fastq ${sample}_R2.fastq
  """
}

process FUSION_CANDIDATES {
  tag "$sample"
  publishDir "${params.outdir}/fusions", mode:'copy', overwrite:true
  input: tuple val(sample), path(r1), path(r2); path index
  output: tuple val(sample), path("${sample}.fusion.txt"), emit: candidates
  script:
  """
  kallisto quant --fusion -t ${task.cpus} -i $index -o ${sample}.kallisto $r1 $r2
  cp ${sample}.kallisto/fusion.txt ${sample}.fusion.txt
  """
}

process PIZZLY_FILTER {
  tag "$sample"
  publishDir "${params.outdir}/fusions", mode:'copy', overwrite:true
  input: tuple val(sample), path(fusion); path gtf; path fasta
  output: tuple val(sample), path("${sample}.fusion.txt"), path("${sample}.pizzly.json"), emit: filtered
  script:
  """
  first_header=\$(grep -m 1 '^>' $fasta)
  first_identifier=\${first_header#>}
  first_identifier=\${first_identifier%% *}
  case "\$first_identifier" in
    *'|'*) echo 'ERROR: pizzly FASTA identifiers still contain pipe-delimited GENCODE metadata' >&2; exit 1 ;;
  esac
  pizzly -k 31 --gtf $gtf --cache gencode.cache --align-score 2 \
    --insert-size 400 --fasta $fasta --output ${sample}.pizzly ${sample}.fusion.txt
  """
}

process REPORT_TRACK2 {
  publishDir "${params.outdir}", mode:'copy', overwrite:true
  input: path itd_dirs; path fusions; path pizzly_jsons; path summary_script
  output: path 'known_positive_recovery.tsv'; path 'summary.txt'
  script:
  """
  mv=\$(find -L . -path '*MV411_getitd/*_hc.tsv' | head -1)
  ka=\$(find -L . -path '*KASUMI1_getitd/*_hc.tsv' | head -1)
  test -n "\$mv"
  # getITD emits no hc file when a sample has no calls.
  test -n "\$ka" || ka=KASUMI1.no_calls.tsv
  touch "\$ka"
  python $summary_script \
    --mv-itd "\$mv" --kasumi-itd "\$ka" \
    --kasumi-fusion \$(find . -name 'KASUMI1.fusion.txt' | head -1) \
    --nb4-fusion \$(find . -name 'NB4.fusion.txt' | head -1) \
    --kasumi-pizzly \$(find . -name 'KASUMI1.pizzly.json' | head -1) \
    --nb4-pizzly \$(find . -name 'NB4.pizzly.json' | head -1) --outdir .
  """
}

workflow {
  PREREGISTER(file(params.expected), file(params.sources), file(params.coords))
  dna = Channel.of(
    ['MV411', params.mv_dna, params.mv_dna_index, 'bam'],
    ['KASUMI1', params.kasumi_dna, params.kasumi_dna_index, params.kasumi_kind])
  sliced_dna = SLICE_DNA(dna)
  itd = CALL_ITD(sliced_dna, file(params.flt3_ref), file(params.flt3_anno))
  rna = Channel.of(
    ['KASUMI1', params.kasumi_rna, params.kasumi_rna_index, '8:92967203-93115514 21:36160098-37376965'],
    ['NB4', params.nb4_rna, params.nb4_rna_index, '15:74287014-74340153 17:38465444-38513094'])
  sliced_rna = SLICE_RNA(rna)
  candidates = FUSION_CANDIDATES(sliced_rna, file(params.kallisto_index))
  filtered = PIZZLY_FILTER(candidates.candidates, file(params.gencode_gtf), file(params.gencode_fasta))
  REPORT_TRACK2(
    itd.calls.map{ it[1] }.collect(),
    filtered.filtered.map{ it[1] }.collect(),
    filtered.filtered.map{ it[2] }.collect(),
    file("${projectDir}/scripts/summarize_track2.py"))
}
