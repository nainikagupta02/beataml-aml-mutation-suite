# Frozen BEAT-AML models tested on TCGA-LAML

Read `results/summary.txt` first. The two elastic-net classifiers, selected
genes, training means/SDs, coefficients, intercepts and threshold are frozen.
No training, tuning, calibration fitting, batch correction or feature selection
is performed on TCGA. No BEAT-AML artifacts are modified.

## Data

`data/gdc_files.json` is the saved GDC manifest for TCGA-LAML, access=open,
data_type=Gene Expression Quantification, analysis.workflow_type=STAR - Counts.
There are 151 files/patients, about 640 MB. Every file size and MD5 is checked.
The expression column is **unstranded integer counts**, not a normalized column.
Ignore N_unmapped/N_multimapping/N_noFeature/N_ambiguous. Each file contains
60,660 gene rows with GENCODE v36 annotation. All 151 specimens are primary
blood-derived cancers; there is one RNA-seq file per patient, so no replicate
selection is needed. Patient barcodes, sample barcodes and file UUIDs are saved.

Labels come from **TCGA NEJM 2013 Table S1**, available on GDC's marker-paper page:
https://gdc.cancer.gov/about-data/publications/laml_2012

Direct supplement: https://api.gdc.cancer.gov/data/b9196563-a05d-40b8-80dc-640ec712eb06

The 200 sequenced patients have a complete author-curated mutation matrix.
`TCGA Patient ID` (e.g. 2803) maps to TCGA-AB-2803. The final aggregate/totals
row is excluded. In this particular complete matrix, blank entries denote no
reported mutation; a missing patient row would remain missing, not negative.

- FLT3-ITD positives are the 39 original `in_frame_ins` calls. FLT3 D835/D839
  point mutations are ITD-negative. This case set agrees exactly with the
  author-released TCGA reannotation at
  https://github.com/gerstung-lab/AML-multistage/tree/master/data.
  The latter contains TCGA cases only and is used only as a label cross-check.
  Current cBioPortal consequence annotations can label some of these insertions
  as splice-region, nonsense or even in-frame deletion; they are not used to
  overwrite the original curated calls. Open masked-MAF is not used.
- NPM1 positives are the 54 reported NPM1-mutant cases, with exact positive-case
  agreement in cBioPortal `laml_tcga`. The RNA-seq intersection includes one
  noncanonical missense p.K263R (TCGA-AB-2915). Primary analysis retains the
  curated label; a predeclared sensitivity analysis excludes this patient using
  the same saved predictions, with no model or threshold changes.
- The 151-case intersection has **FLT3-ITD 34+/117-**, **NPM1 38+/113-**.
  No target labels are missing. Another 10 cases have FLT3 point mutations.
  ~30% is closer to overall FLT3 mutation prevalence, not ITD-only prevalence.

## Gene harmonization

Strip Ensembl version suffixes and preserve PAR_Y suffixes as distinct identities.
Match only stable Ensembl IDs; do not guess equivalences using gene symbols or
apply coordinate liftover to gene-level quantifications. Preserve the exact
gene order in BEAT-AML `gene_annotation.csv` (51,016 pipeline input genes).

46,973 input genes are present; 4,043 have no direct ID match. Build the frozen
input matrix from matched raw counts, using zero placeholders for unavailable
input columns solely to execute sample-local CPM. Its denominator is the sum
over the original input-gene universe, before the frozen 5,000-gene filter.
Additional TCGA genes are excluded, as requested; retained-library fractions
are reported for every patient. Annotation differences make this denominator
an approximation and are a remaining transfer limitation.

Apply the **frozen** LogCPM transform, filter indices and StandardScaler. Then
set every unavailable selected gene to **standardized zero**, as instructed,
so it contributes exactly zero to the linear predictor. Do not assign a
negative standardized value as though the gene had measured expression zero.
Present genes with real zero counts are processed normally.

Coverage:

| Target | Selected present | Nonzero predictors present | Absolute coefficient weight retained |
|---|---:|---:|---:|
| FLT3-ITD | 4871/5000 (97.42%) | 89/93 | 96.78% |
| NPM1 | 4872/5000 (97.44%) | 86/88 | 97.32% |

The 129/128 absent selected genes are neutralized; 4/2 of these have nonzero
coefficients. Coverage below 90% coefficient weight triggers a prominent warning.
Every selected gene's availability and coefficient are in the coverage CSVs.

## Evaluation and reproduction

Seed 20260916. Threshold 0.50 is carried over unchanged. Point estimates and
95% percentile confidence intervals are reported for AUROC, AUPRC (average
precision), sensitivity, specificity, PPV, NPV and Brier score. The 5,000 ordinary
patient-bootstrap resamples vary prevalence; here each patient has one sample.
Single-class resamples are skipped for discrimination metrics, and the number
of defined replicates is saved. No models are fitted during bootstrapping.
Calibration uses five quantile bins and the original probabilities.

Use the original environment recorded in ../beataml_modeling/requirements.txt.
From this directory:

```bash
python3 -B acquire.py --verify-only
python3 -B -m unittest test_validation.py
# Existing results already contain the completed one-time inference.
python3 -B validate.py report
python3 -B audit.py
```

`acquire.py` without `--verify-only` downloads missing public files with curl,
using the saved GDC manifest and checking checksums. Label data hashes prevent
silent changes when an upstream source is revised. Cached source files and
the saved manifest define the original reproducible data snapshot.

To reproduce the same fixed protocol in a separate output folder:

```bash
python3 -B validate.py prepare --out results_reproduction
python3 -B validate.py evaluate --out results_reproduction
python3 -B validate.py report --out results_reproduction
```

Prepare records hashes and locks the protocol before external inference.
Evaluate refuses a second inference in the same output folder. Report uses
saved predictions only. Every phase verifies all original BEAT-AML files
remain unchanged; inference additionally checks the fitted pipeline's in-memory
hash before/after transforms. Do not change the protocol or select variants
based on these external results.

The four tests verify full-coverage equivalence to the original pipeline,
zero contribution of absent genes without fitting, genuine raw-count parsing,
and unchanged BEAT-AML files. Toy numerical inputs are only software tests,
not an additional dataset or evidence of predictive performance.

## Outputs

- `results/summary.txt`: full plain-text report and limitations.
- `results/metrics.csv`: point estimates, intervals, denominators and confusion matrices.
- `results/external_diagnostics.png`: reliability curves, ROC and precision-recall plots.
- `results/predictions.csv`: immutable external predictions with patient/sample IDs.
- `results/cohort.csv`, `curated_labels_all200.csv`: exact joins and original label strings.
- `results/coverage_summary.csv`, `*_gene_coverage.csv`: gene and coefficient coverage.
- `results/protocol.json`, source/prepared/model hashes: provenance and frozen settings.

This is one external-cohort feasibility test with only 34/38 positive patients.
Matching raw-count units does not remove counting-pipeline, annotation,
cell-composition, specimen-source or population shifts. Bootstrap intervals
condition on the already-trained classifier and do not establish clinical validity.
