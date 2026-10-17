# Laptop-scale AML variant workflow demo

This repository demonstrates a reproducible somatic-variant workflow that can
run locally on a laptop: a preregistered recovery of known AML cell-line
events using region-only downloads, with samtools, getITD, kallisto and
Pizzly.

## Result

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
are pseudoaligned to GENCODE v50 transcripts. This is known-positive recovery,
not a cohort sensitivity/specificity estimate.

See [`TRACK2.md`](TRACK2.md) for sample IDs, coordinates, reference preparation,
resources, and limitations. The corrected table is in
[`results/track2_corrected/known_positive_recovery.tsv`](results/track2_corrected/known_positive_recovery.tsv).
