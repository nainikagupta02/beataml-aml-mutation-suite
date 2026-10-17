# Track 1: HCC1395 WES somatic SNV/indel validation

Track 1 completed successfully. It validates a matched tumor-normal WES
workflow against one SEQC2 HCC1395 replicate pair. HCC1395 is a breast-cancer
reference cell line; this tests the workflow, not AML-specific biology or
clinical performance.

## Results

| Type | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| SNV | 705 | 283 | 700 | 0.7136 | 0.5018 | 0.5892 |
| Indel | 32 | 16 | 19 | 0.6667 | 0.6275 | 0.6465 |

These are exact-allele results for one WES replicate pair. They should not be
generalized to all SEQC2 replicates or to AML specimens.

## Region selection

The first candidate was the intersection of cBioPortal's `IMPACT468` genes,
GENCODE v50 basic exons, and the SEQC2 high-confidence BED. It contained only
44 truth SNVs and one truth indel, below the configured minimum of 200. The
pipeline followed its go/no-go rule and first widened to all GENCODE basic
exons. That trial covered 151,135,889 bases, broader than a defensible coding
exome, so the final fallback was tightened to GENCODE v50 basic protein-coding
exons intersected with the high-confidence BED. The final scope has 225,074
intervals covering 96,549,992 bases and contains 1,405 truth SNVs and 51 truth
indels. This makes the result a high-confidence exonic WES benchmark, not a
targeted-panel benchmark. The 151 Mb trial did not supply the final read slices
or metrics.

The cBioPortal endpoint currently returns 469 unique symbols despite the
`IMPACT468` name and description. `GTF2I` is the extra returned symbol; it is
not a duplicate. The downloaded source list was preserved.

## Evaluation checks

Mutect2 used `track1_slice.padded200.bed` for local assembly. The PASS calls and
both truth VCFs were then restricted independently to the exact unpadded
`track1_validation.bed`. Calls and truth were decomposed and normalized against
the same GRCh38.d1.vd1 FASTA before `bcftools isec`.

The closeout audit reports zero records outside the unpadded BED in the final
callset and all four type-specific benchmark VCFs. A second normalization pass
reported zero additional splits and zero realignments for calls and truth. See
`results/track1/closeout_audit.txt` for the recorded checks.

## Software and data versions

- Python 3.10.14
- OpenJDK 17.0.11
- Nextflow 26.04.6
- samtools, bcftools, and htslib 1.21
- bedtools 2.31.1
- GATK 4.6.2.0
- GENCODE human release 50, GRCh38.p14 basic annotation
- GDC GRCh38.d1.vd1 analysis reference
- SEQC2 truth calls v1.2.1 and high-confidence regions v1.2

The run was performed on macOS 26.6.2 on Apple silicon using an isolated
`osx-64` Conda environment through Rosetta 2.

## Reproduce locally

```bash
source activate_env.sh
python scripts/monitor_pipeline.py &
nextflow run main.nf -profile local -resume -ansi-log false
```

The pipeline has six processes: resource acquisition, validation-region
construction and truth-count gating, indexed remote BAM slicing, Mutect2,
FilterMutectCalls, and exact-allele benchmarking with `bcftools isec`.

Final results are under `results/track1`. The concise report is
`results/track1/summary.txt`; metrics are in
`results/track1/benchmark/benchmark/metrics.tsv`; per-step resource use is in
`results/track1/resource_usage.tsv`.

The local profile uses the project-local Conda environment, at most four CPUs,
and at most 10 GB for Mutect2. `PYTHONHASHSEED=20260916` is set for scripted
steps. Source datasets are read-only; generated files go to Nextflow `work/`
and `results/track1/`.

No full WES BAM is downloaded. `samtools view -L` streams reads overlapping the
final validation intervals plus 200 bp from the indexed remote tumor and normal
BAMs.

## Measured local resources

| Process | Runtime | Peak RAM | Work-directory disk |
|---|---:|---:|---:|
| Fetch resources | 2 min 40 s | 777.6 MiB proxy | 3.9 GiB |
| Build validation regions | 39.5 s | 836.9 MiB | 20 MiB |
| Slice remote BAMs | 2 min 33 s | 32.0 MiB | 4.4 GiB |
| Mutect2 | 1 h 51 min | 1,019.3 MiB | 2.0 MiB excluding linked BAMs |
| Filter calls | 11.4 s | 657.3 MiB | 848 KiB |
| Benchmark | 1.7 s | 26.2 MiB | 3.5 MiB |

The resource-fetch RAM value is the highest-memory reference-preparation
command, because macOS Nextflow does not expose task RSS without a container
engine. The other stated peak values are sampled or deterministic replay
measurements described in `results/track1/resource_usage.tsv`.

## Limitations

Mutect2 was deliberately run without a panel of normals, gnomAD germline prior,
contamination estimation, or orientation-bias filtering. These omissions likely
reduce precision. The SEQC2 truth set draws on richer multi-platform evidence
than this single WES replicate, so missed variants can reflect coverage as well
as caller behavior. Exact-allele comparison uses no fuzzy matching.
