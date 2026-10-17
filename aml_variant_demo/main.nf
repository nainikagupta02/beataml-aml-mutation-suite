nextflow.enable.dsl = 2

process FETCH_RESOURCES {
    tag 'SEQC2 + GENCODE v50 + IMPACT468 + GRCh38.d1.vd1'
    publishDir "${params.outdir}/resources", mode: 'symlink', overwrite: true

    output:
    path 'resources', emit: resources

    script:
    """
    mkdir -p resources

    curl -L --fail --retry 5 --retry-delay 5 \
      https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz \
      -o resources/sSNV.truth.vcf.gz
    curl -L --fail --retry 5 --retry-delay 5 \
      https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz \
      -o resources/sINDEL.truth.vcf.gz
    curl -L --fail --retry 5 --retry-delay 5 \
      https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/High-Confidence_Regions_v1.2.bed \
      -o resources/High-Confidence_Regions.bed
    curl -L --fail --retry 5 --retry-delay 5 \
      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz \
      -o resources/gencode.v50.basic.annotation.gtf.gz
    curl -L --fail --retry 5 --retry-delay 5 \
      https://www.cbioportal.org/api/gene-panels/IMPACT468 \
      -o resources/IMPACT468.json
    curl -L --fail --retry 5 --retry-delay 5 \
      https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 \
      -o GRCh38.d1.vd1.fa.tar.gz

    echo '3ffbcfe2d05d43206f57f81ebb251dc9  GRCh38.d1.vd1.fa.tar.gz' | md5sum -c -
    gzip -t resources/sSNV.truth.vcf.gz
    gzip -t resources/sINDEL.truth.vcf.gz
    gzip -t resources/gencode.v50.basic.annotation.gtf.gz
    tar -xzf GRCh38.d1.vd1.fa.tar.gz -C resources
    test -s resources/GRCh38.d1.vd1.fa

    samtools faidx resources/GRCh38.d1.vd1.fa
    gatk --java-options '-Xmx2g' CreateSequenceDictionary \
      -R resources/GRCh38.d1.vd1.fa \
      -O resources/GRCh38.d1.vd1.dict
    tabix -f -p vcf resources/sSNV.truth.vcf.gz
    tabix -f -p vcf resources/sINDEL.truth.vcf.gz

    python ${projectDir}/scripts/parse_impact_panel.py \
      resources/IMPACT468.json resources/impact_genes.txt resources/impact_panel_report.txt

    cat > resources/source_manifest.tsv <<'EOF'
item	build	source
tumor_bam	GRCh38.d1.vd1	https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_LL_T_1.bwa.dedup.bam
normal_bam	GRCh38.d1.vd1	https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_LL_N_1.bwa.dedup.bam
truth_snv	GRCh38.d1.vd1	https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
truth_indel	GRCh38.d1.vd1	https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz
truth_regions	GRCh38.d1.vd1	https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/High-Confidence_Regions_v1.2.bed
annotation	GRCh38.p14	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
gene_panel	HGNC symbols	https://www.cbioportal.org/api/gene-panels/IMPACT468
reference	GRCh38.d1.vd1	https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834
EOF
    """
}

process BUILD_VALIDATION_REGIONS {
    tag 'truth-count gate'
    publishDir "${params.outdir}/regions", mode: 'copy', overwrite: true

    input:
    path resources

    output:
    path 'regions', emit: regions

    script:
    """
    # fallback_scope=GENCODE_v50_basic_protein_coding_exons
    mkdir -p regions
    python ${projectDir}/scripts/build_validation_regions.py \
      --gtf resources/gencode.v50.basic.annotation.gtf.gz \
      --genes resources/impact_genes.txt \
      --truth-bed resources/High-Confidence_Regions.bed \
      --snv resources/sSNV.truth.vcf.gz \
      --indel resources/sINDEL.truth.vcf.gz \
      --fai resources/GRCh38.d1.vd1.fa.fai \
      --min-truth ${params.min_truth_variants} \
      --outdir regions
    """
}

process SLICE_REMOTE_BAMS {
    tag 'WES_LL_T_1 + WES_LL_N_1'
    publishDir "${params.outdir}/sliced_bams", mode: 'symlink', overwrite: true

    input:
    path regions

    output:
    path 'HCC1395.tumor.bam', emit: tumor_bam
    path 'HCC1395.tumor.bam.bai', emit: tumor_bai
    path 'HCC1395BL.normal.bam', emit: normal_bam
    path 'HCC1395BL.normal.bam.bai', emit: normal_bai
    path 'slice_report.txt', emit: slice_report

    script:
    """
    samtools view -@ ${task.cpus} -b -h -L regions/track1_slice.padded200.bed \
      '${params.tumor_url}' -o HCC1395.tumor.bam
    samtools index -@ ${task.cpus} HCC1395.tumor.bam
    samtools view -@ ${task.cpus} -b -h -L regions/track1_slice.padded200.bed \
      '${params.normal_url}' -o HCC1395BL.normal.bam
    samtools index -@ ${task.cpus} HCC1395BL.normal.bam
    samtools quickcheck -v HCC1395.tumor.bam HCC1395BL.normal.bam

    {
      echo -e 'sample\tread_group_sample\treads\tbytes'
      echo -e 'tumor\t'\$(samtools samples HCC1395.tumor.bam | cut -f1 | tail -n 1)'\t'\$(samtools view -c HCC1395.tumor.bam)'\t'\$(/usr/bin/stat -f %z HCC1395.tumor.bam)
      echo -e 'normal\t'\$(samtools samples HCC1395BL.normal.bam | cut -f1 | tail -n 1)'\t'\$(samtools view -c HCC1395BL.normal.bam)'\t'\$(/usr/bin/stat -f %z HCC1395BL.normal.bam)
    } > slice_report.txt
    """
}

process MUTECT2_CALL {
    tag 'HCC1395 tumor-normal'
    publishDir "${params.outdir}/calls", mode: 'symlink', overwrite: true

    input:
    path resources
    path regions
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai

    output:
    path 'HCC1395.unfiltered.vcf.gz', emit: unfiltered_vcf
    path 'HCC1395.unfiltered.vcf.gz.tbi', emit: unfiltered_tbi
    path 'HCC1395.unfiltered.vcf.gz.stats', emit: unfiltered_stats
    path 'HCC1395.f1r2.tar.gz', emit: f1r2

    script:
    """
    gatk --java-options '-Xmx8g -XX:ParallelGCThreads=2' Mutect2 \
      -R resources/GRCh38.d1.vd1.fa \
      -I ${tumor_bam} -tumor '${params.tumor_sample}' \
      -I ${normal_bam} -normal '${params.normal_sample}' \
      -L regions/track1_slice.padded200.bed \
      --native-pair-hmm-threads ${task.cpus} \
      --f1r2-tar-gz HCC1395.f1r2.tar.gz \
      -O HCC1395.unfiltered.vcf.gz
    """
}

process FILTER_CALLS {
    tag 'FilterMutectCalls + exact BED'
    publishDir "${params.outdir}/calls", mode: 'copy', overwrite: true

    input:
    path resources
    path regions
    path unfiltered_vcf
    path unfiltered_tbi
    path unfiltered_stats

    output:
    path 'HCC1395.pass.validation.vcf.gz', emit: pass_vcf
    path 'HCC1395.pass.validation.vcf.gz.tbi', emit: pass_tbi
    path 'call_counts.txt', emit: call_counts

    script:
    """
    gatk --java-options '-Xmx4g' FilterMutectCalls \
      -R resources/GRCh38.d1.vd1.fa \
      -V ${unfiltered_vcf} \
      -O HCC1395.filtered.vcf.gz

    bcftools view -f PASS -R regions/track1_validation.bed \
      HCC1395.filtered.vcf.gz -Ou \
      | bcftools norm -f resources/GRCh38.d1.vd1.fa -m -any -Oz \
      -o HCC1395.pass.validation.vcf.gz
    tabix -f -p vcf HCC1395.pass.validation.vcf.gz

    {
      echo -e 'type\tcount'
      echo -e 'SNV\t'\$(bcftools view -H -v snps HCC1395.pass.validation.vcf.gz | wc -l | tr -d ' ')
      echo -e 'INDEL\t'\$(bcftools view -H -v indels HCC1395.pass.validation.vcf.gz | wc -l | tr -d ' ')
    } > call_counts.txt
    """
}

process BENCHMARK_CALLS {
    tag 'bcftools isec SNV + indel'
    publishDir "${params.outdir}/benchmark", mode: 'copy', overwrite: true

    input:
    path resources
    path regions
    path pass_vcf
    path pass_tbi

    output:
    path 'benchmark', emit: benchmark

    script:
    """
    mkdir -p benchmark/isec_snv benchmark/isec_indel

    bcftools view -R regions/track1_validation.bed resources/sSNV.truth.vcf.gz -Ou \
      | bcftools norm -f resources/GRCh38.d1.vd1.fa -m -any -Oz \
      -o benchmark/truth.snv.vcf.gz
    bcftools view -R regions/track1_validation.bed resources/sINDEL.truth.vcf.gz -Ou \
      | bcftools norm -f resources/GRCh38.d1.vd1.fa -m -any -Oz \
      -o benchmark/truth.indel.vcf.gz
    bcftools view -v snps ${pass_vcf} -Oz -o benchmark/calls.snv.vcf.gz
    bcftools view -v indels ${pass_vcf} -Oz -o benchmark/calls.indel.vcf.gz
    for vcf in benchmark/*.vcf.gz; do tabix -f -p vcf "\$vcf"; done

    bcftools isec -p benchmark/isec_snv benchmark/truth.snv.vcf.gz benchmark/calls.snv.vcf.gz
    bcftools isec -p benchmark/isec_indel benchmark/truth.indel.vcf.gz benchmark/calls.indel.vcf.gz

    python ${projectDir}/scripts/calculate_metrics.py \
      --snv-dir benchmark/isec_snv \
      --indel-dir benchmark/isec_indel \
      --output-tsv benchmark/metrics.tsv \
      --output-txt benchmark/metrics.txt
    cp regions/gate_report.txt benchmark/region_gate.txt
    cp regions/region_stats.tsv benchmark/region_stats.tsv
    """
}

workflow {
    fetched = FETCH_RESOURCES()
    built = BUILD_VALIDATION_REGIONS(fetched.resources)
    sliced = SLICE_REMOTE_BAMS(built.regions)
    called = MUTECT2_CALL(
        fetched.resources,
        built.regions,
        sliced.tumor_bam,
        sliced.tumor_bai,
        sliced.normal_bam,
        sliced.normal_bai
    )
    filtered = FILTER_CALLS(
        fetched.resources,
        built.regions,
        called.unfiltered_vcf,
        called.unfiltered_tbi,
        called.unfiltered_stats
    )
    BENCHMARK_CALLS(
        fetched.resources,
        built.regions,
        filtered.pass_vcf,
        filtered.pass_tbi
    )
}
