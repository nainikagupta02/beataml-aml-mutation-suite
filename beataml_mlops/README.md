# MLOps additions for the BEAT-AML / TCGA-LAML study

Read `../beataml_modeling/README.md` and `../tcga_laml_validation/README.md` first.
This folder adds operational tooling around those two frozen, already-audited
projects. Nothing here fits, tunes, recalibrates or modifies either project; it
only reads their saved `results/` files (hashed source data, frozen `.joblib`
pipelines, and already-computed CV/held-out/external metrics).

## Scope

Built now:

- **MLflow logging and registry** (`log_models.py`) — logs each target's frozen
  elastic-net model (the one actually used downstream) with its hyperparameters,
  development CV metrics, BEAT-AML held-out metrics, TCGA-LAML external metrics,
  and coefficient table, then registers the already-fitted pipeline in the MLflow
  Model Registry under alias `frozen`.
- **Drift detection** (`drift.py`) — compares BEAT-AML (training, n=445) and
  TCGA-LAML (external, n=151) gene-expression distributions in the exact frozen
  LogCPM space, and decomposes the shift into each gene's contribution to the
  frozen linear predictor, tying the drift directly to the AUROC/Brier change
  already observed between the two evaluations.
- **Monitoring check** (`monitor.py`) — a per-sample out-of-distribution flag,
  generalizing the aggregate QC check already computed in
  `tcga_laml_validation/validate.py`, runnable on new data before trusting a
  prediction from it.
- **PyTorch comparison** (`torch_model.py`) — a small, heavily regularized MLP
  trained per target on the exact same BEAT-AML development rows and CV folds
  as the frozen elastic-net study, then compared to it in three places: BEAT-AML
  CV, the BEAT-AML held-out partition, and the TCGA-LAML external cohort (the
  TCGA comparison is the cleanest one — see caveats below).

- **Containerized API** (`api/`) — a FastAPI service wrapping the frozen
  elastic-net pipelines (`POST /predict/{target}`), with the same
  out-of-distribution heuristic as `monitor.py` built in. Built and run
  locally as a Docker image; verified end-to-end against real TCGA-LAML
  samples through the running container (see `api/README.md`).
- **AWS deployment config** (`deploy/aws/`) — Terraform for an ECR repo + App
  Runner service to run that image. **Not applied** — see
  `deploy/aws/README.md` for the manual build/push/apply steps and cost
  before doing so.

## Run

From this directory, with the frozen `develop`+`finalize` phases of
`beataml_modeling` and the `prepare`/`evaluate`/`report` phases of
`tcga_laml_validation` already completed (they are, under each project's
`results/`):

```bash
python3 log_models.py
python3 drift.py
python3 monitor.py
python3 torch_model.py   # ~1-2 minutes on CPU; refuses to rerun its locked held-out evaluation
```

`log_models.py` writes to a local file-based MLflow store at `./mlruns`; browse
it with `mlflow ui --backend-store-uri file://$(pwd)/mlruns`. The rest write to
`./results/`. For the API, see `api/README.md`; for the (unapplied) AWS config,
see `deploy/aws/README.md`.

## Notes and limitations

- MLflow logging only covers the elastic-net family (the selected model for
  both targets). The other three families (majority, marker logistic, random
  forest) are not registered; their CV numbers are already in
  `beataml_modeling/results/cv_summary.csv` if needed for comparison.
- The drift decomposition assumes the frozen coefficients remain valid
  out-of-cohort; it explains average linear-predictor bias, not per-patient
  calibration, and standardized mean difference / KS statistics are
  descriptive, not causal — they cannot separate sequencing/annotation
  differences from true population shift.
- Genes absent from TCGA (by Ensembl ID) are zeroed under the same frozen rule
  used for external inference. `drift.py` reports this explicitly: it hides
  those genes' real drift rather than resolving it.
- The monitoring check reuses the same 5/10 standardized-z thresholds already
  validated as a cohort-level QC check in `tcga_laml_validation`; it is a
  heuristic screen, not a calibrated anomaly detector, and a flagged sample is
  not automatically a bad prediction.
- No model here was retrained, retuned or recalibrated on TCGA or any other
  external data, consistent with both parent projects' audited protocols.
- `torch_model.py` uses one fixed, prespecified architecture (MLP 64-16, dropout
  0.3, weight decay 1e-3, 60 epochs) with no hyperparameter search, unlike the
  elastic-net's small grid — results are not tuned to win the comparison.
  At n~350-360 training specimens and up to 5,000 features this is a severe
  small-n/high-p regime for a neural network; across both targets the frozen
  elastic-net matched or beat the MLP on CV and TCGA AUPRC/Brier, so this is a
  feasibility comparison, not a claim that deep learning helps here. Reusing
  the BEAT-AML held-out partition for the MLP's single evaluation is a
  secondary analysis (see `results/pytorch_vs_elastic_net.txt`); the TCGA
  comparison is the cleaner one, since both models get an identical frozen
  external cohort and gene zero-imputation rule.
