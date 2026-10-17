# Laptop-scale AML variant workflow demo

This repository demonstrates reproducible somatic-variant workflows that can
run locally on a laptop. Track 1 is complete: it validates a matched
tumor-normal **whole-exome sequencing (WES)** SNV/indel workflow against the
SEQC2 HCC1395 truth set. HCC1395 is a breast-cancer reference cell line, so the
result validates the workflow and its implementation; it is not an estimate of
variant-calling performance in AML samples.

## Track 1 result

The final exact-allele benchmark used PASS Mutect2 calls and SEQC2 v1.2.1
high-confidence truth variants, restricted on both sides to the same unpadded
validation BED and normalized against the same GRCh38.d1.vd1 reference.

| Type | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| SNV | 705 | 283 | 700 | 0.7136 | 0.5018 | 0.5892 |
| Indel | 32 | 16 | 19 | 0.6667 | 0.6275 | 0.6465 |

This is one WES replicate pair (`WES_LL_T_1` and `WES_LL_N_1`). Recall is
therefore conditional on that replicate's coverage and should not be read as a
multi-replicate estimate.

## Region gate and final scope

MSK-IMPACT was used only to define the first, small validation-region candidate.
Its high-confidence exon intersection contained 44 truth SNVs and one truth
indel, below the predeclared minimum of 200 variants. The workflow therefore
initially widened to all GENCODE basic exons. That trial spanned 151,135,889
bases, which was broader than a defensible coding-exome scope, so the final
rule was tightened to GENCODE v50 basic **protein-coding** exons intersected
with the SEQC2 high-confidence BED: 225,074 intervals, 96,549,992 bases, 1,405
truth SNVs, and 51 truth indels. The reported result is consequently a
high-confidence exonic WES benchmark, not a targeted-panel benchmark. The
151 Mb trial was abandoned before the final read slicing and benchmark.

The current cBioPortal `IMPACT468` endpoint describes a 468-gene panel but
returns 469 unique HGNC symbols, including `GTF2I`. The source list was
preserved, and the discrepancy is reported rather than silently editing it.

## Evaluation safeguards

- Mutect2 assembled reads in the final regions plus 200 bp of padding.
- `FilterMutectCalls` PASS variants were restricted to the exact unpadded BED.
- Both truth VCFs were independently restricted to that same unpadded BED.
- Calls and truth were decomposed and normalized with `bcftools norm -f` against
  the same GRCh38.d1.vd1 FASTA before `bcftools isec`.
- A closeout audit found zero call or truth records outside the BED. Re-running
  normalization found zero additional splits or realignments in all four
  benchmark VCFs.

The complete evidence is in
[`results/track1/closeout_audit.txt`](results/track1/closeout_audit.txt). The
detailed workflow report and reproduction command are in
[`TRACK1.md`](TRACK1.md).

## Laptop feasibility

The run used at most four CPU cores. The heaviest process, Mutect2, took 1 h 51
min and peaked at 1,019.3 MiB RSS. Remote BAM slicing took 2 min 33 s and kept
only region-overlapping reads, so no full WES BAM was downloaded. Final-run work
directories occupied about 8.3 GiB; the environment and package cache occupied
5.4 GiB.

The simple matched-pair workflow omits a panel of normals, gnomAD germline
prior, contamination estimation, and orientation-bias filtering. Those choices
keep the demonstration compact and reproducible, but they likely reduce
precision. Exact matching can also count differently represented complex
alleles as errors even after standard normalization.

## Track 2 result

Track 2 tests preregistered AML cell-line events with region-only downloads.
All four checks passed. MV4-11 yielded a 30-bp FLT3-ITD at 23.077% VAF with
six unique supporting reads, while the required KASUMI-1 negative control
stayed negative. Pizzly recovered RUNX1–RUNX1T1 in KASUMI-1 with 199 split and
one paired supporting record, and PML–RARA in NB-4 with 65 split records.

An earlier zero-call fusion run was traced to unmodified pipe-delimited GENCODE
FASTA headers. The corrected workflow converts and audits the headers before
building the kallisto index, and fails early if incompatible identifiers are
supplied. The official calls come from the complete 670,670-transcript index.

The current open inputs still force two deviations: DNA is GRCh38 WGS rather
than WES, and RNA begins as legacy GRCh37 TopHat BAMs before partner-locus reads
are pseudoaligned to GENCODE v50 transcripts. Track 2 is known-positive
recovery, not a cohort sensitivity/specificity estimate.

See [`TRACK2.md`](TRACK2.md) for sample IDs, coordinates, reference preparation,
resources, and limitations. The corrected table is in
[`results/track2_corrected/known_positive_recovery.tsv`](results/track2_corrected/known_positive_recovery.tsv).
