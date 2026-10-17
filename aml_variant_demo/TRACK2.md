# Track 2 — AML cell-line known-positive recovery

## Result first

This track is a **known-positive recovery exercise**, not a sensitivity or
specificity validation. All four preregistered checks passed. MV4-11 carried
the expected FLT3-ITD, the required KASUMI-1 ITD negative control stayed
negative, and pizzly recovered both expected RNA fusions with the corrected
full GENCODE v50 reference.

| Cell line | Event | Expected | Observed | Outcome |
|---|---|---|---|---|
| MV4-11 | FLT3-ITD | positive | 30-bp ITD, 6 unique reads, 23.077% VAF | ✅ |
| KASUMI-1 | FLT3-ITD | negative | no ITD detected; low residual depth | ✅ |
| KASUMI-1 | RUNX1–RUNX1T1 | positive | pizzly call: 199 split + 1 paired supporting records | ✅ |
| NB-4 | PML–RARA | positive | pizzly call: 65 split supporting records | ✅ |

The fusion evidence is now caller-backed. It supersedes the earlier zero-call
run, which used unmodified GENCODE FASTA headers. That failure was a reference
formatting error: kallisto retained pipe-delimited metadata in transcript names,
while pizzly's GTF cache used bare transcript IDs. Pizzly's documented GENCODE
preparation requires replacing pipes in FASTA headers before building the
kallisto index. `scripts/prepare_pizzly_fasta.py` performs and audits that step,
and `PIZZLY_FILTER` now fails early if a pipe remains in the first identifier.

A small four-gene diagnostic first confirmed the cause. The official results
above come from a clean rerun against the complete 670,670-transcript index,
not from the diagnostic reference. Pizzly retained 209 of 386 raw KASUMI-1
candidate records and 74 of 131 NB-4 records. Its GTF parser warned that 26,378
of the FASTA transcripts lacked annotations, but the expected partner genes
were annotated and successfully called.

## Material source deviations

The current DepMap AWS inventory does not expose WES for MV4-11 or KASUMI-1.
It exposes GRCh38 WGS for both, so the DNA arm used only remotely sliced WGS
reads over the FLT3 target. It did not download either full alignment and must
not be called WES.

The open RNA alignments available for KASUMI-1 and NB-4 are legacy GRCh37
TopHat BAMs. No current GRCh38 RNA BAM or open FASTQ was found for these exact
samples. Those BAMs were sliced at verified GRCh37 partner loci, converted to
FASTQ, and pseudoaligned to GENCODE v50 GRCh38 transcripts. The fusions were
recovered, but this still fails the strict “GRCh38 for everything” condition
at the source-alignment step.

## Verified identities and coordinates

| Name | DepMap model | Omics sample/profile | Source build |
|---|---|---|---|
| MV4-11 | ACH-000045 | CDS-XBgVc6 / PR-46SoC2 (WGS) | hg38 |
| KASUMI-1 | ACH-000263 | CDS-pouIEM / PR-V5gLbW (WGS); CDS-KnmxYu / PR-ahX9nb (RNA) | hg38 DNA; hg19 RNA |
| NB-4 | ACH-000294 | CDS-zQSwaP / PR-A2ryTk (RNA) | hg19 RNA |

FLT3 is minus-strand on chromosome 13. GENCODE v50 canonical/basic transcript
`ENST00001131816.1` places exon 14 at `chr13:28034082-28034214` and exon 15 at
`chr13:28033887-28033991` (GRCh38, one-based closed). The getITD reference is
the reverse-complemented 330-bp span `chr13:28033886-28034215`; its annotation
was shifted from getITD's supplied hg19 annotation and checked against the GDC
GRCh38.d1.vd1 FASTA.

Legacy RNA slice coordinates were derived from official GENCODE v19 gene
bounds: RUNX1T1 `8:92967203-93115514`, RUNX1
`21:36160098-37376965`, PML `15:74287014-74340153`, and RARA
`17:38465444-38513094`.

## Reference preparation and reproduction

Install and activate the environment as described in `INSTALL.md`. Create the
pizzly-compatible FASTA and index once:

```bash
source activate_env.sh
python scripts/prepare_pizzly_fasta.py \
  --input track2_discovery/pizzly/gencode.v50.transcripts.fa \
  --output track2_discovery/pizzly/gencode.v50.transcripts.pizzly.fa \
  --report track2_discovery/pizzly/pizzly_fasta_conversion.tsv
kallisto index -k 31 \
  -i track2_discovery/pizzly/gencode.v50.pizzly.full.kallisto.idx \
  track2_discovery/pizzly/gencode.v50.transcripts.pizzly.fa
```

The audit recorded 670,670 headers changed and a first identifier of
`ENST00000832824.1`. The corrected fully remote run used:

```bash
nextflow -C nextflow_track2.config run track2.nf \
  -profile replay --outdir results/track2_corrected
```

The run used macOS 26.6.2 on Apple silicon with the x86_64 Conda environment
under Rosetta 2. The seed is `20260917`; none of the called tools is stochastic
under these settings. Authoritative source URLs are in `track2_sources.tsv`.

## Resource evidence

The corrected remote run completed in about 4 minutes 3 seconds from the first
slice through the report and retained 15 MiB of results. It downloaded only
indexed target regions: the DNA BAMs total 272 KiB and the RNA FASTQs total
about 10.7 MiB.

Nextflow 26.04.6 does not report task RSS on macOS without a container engine,
so representative commands were measured separately. getITD peaked at 74 MiB,
full-reference kallisto quant at 7.44 GiB, and pizzly at 2.65 GiB. The corrected
full index build took 681 seconds and produced a 4.84-GiB index; its peak RSS
was not captured and remains `NA`. The active GTF, corrected FASTA, index, and
cache use about 10.9 GiB. The directory is larger because the incompatible
original FASTA/index were retained as provenance. See
`results/track2_corrected/resource_report.tsv`.

## Limitations

There are only three cell lines and no matched normals. This demonstrates
recovery of preregistered known events; it does not estimate population
sensitivity, specificity, or false-positive rate. KASUMI-1 had only five
getITD reads left after quality/alignment filtering, so its negative ITD result
has limited coverage support. RNA began as legacy GRCh37 TopHat alignments and
only partner-locus reads were retained, so the fusion result is a focused
laptop-scale demonstration rather than an unbiased transcriptome-wide fusion
survey. KASUMI-1 also produced additional pizzly gene-pair records; because the
track lacks a comprehensive truth set, their precision was not assessed.

Exact inputs are in `track2_sources.tsv`, preregistered outcomes in
`track2_expected_results.tsv`, coordinates in `track2_coordinates.tsv`, and
the corrected evidence in `results/track2_corrected/`. FASTA, GTF, and index
checksums are in `results/track2_corrected/reference_correction_manifest.tsv`.
