# BEAT-AML mutation-probability API

A FastAPI service wrapping the two frozen elastic-net pipelines from
`beataml_modeling`. It never fits, retrains, tunes or recalibrates anything;
it loads the already-saved `.joblib` pipelines and applies them.

Read `../../beataml_modeling/README.md` and `../../tcga_laml_validation/README.md`
first: this is a feasibility-study model, not a validated clinical diagnostic.

## Run locally

```bash
python3 prepare_artifacts.py   # copies frozen files into ./artifacts (read-only)
docker build -t beataml-mutation-api:latest .
docker run -d -p 8000:8000 --name beataml-api beataml-mutation-api:latest
curl http://localhost:8000/health
```

Or without Docker: `pip install -r requirements.txt && python3 prepare_artifacts.py
&& uvicorn main:app --reload` from this directory.

## Endpoints

- `GET /health` — status and the two target names.
- `GET /targets/{target}/genes` — the exact selected genes (Ensembl stable_id,
  symbol, frozen coefficient) this target's model weights, so a caller knows
  what it needs to supply.
- `POST /predict/{target}` — body `{"gene_counts": {"ENSG...": <raw count>, ...}}`.
  Returns the predicted probability, the 0.5 threshold decision, how many of
  the model's selected genes were actually supplied, and an out-of-distribution
  flag (same heuristic as `../monitor.py`, generalized to one sample).

### Important input contract detail

The frozen `LogCPM` step needs a sample's **total** library size, computed
over the full ~51,016-gene training input universe, not just the ~5,000
selected genes a given model weights. If you supply fewer genes than that,
the library-size denominator is approximate — the response's
`coverage_warning` field says so explicitly, mirroring the same approximation
`tcga_laml_validation` documents for its own external cohort. Genes with
unrecognized Ensembl IDs are ignored and listed back in
`unknown_gene_ids_ignored`.

Verified against real data: fed real TCGA-LAML counts through this API and
confirmed the returned probabilities match `tcga_laml_validation/results/predictions.csv`
(the frozen pipeline's own saved predictions) to within the expected
difference from a sparser gene payload.

## Limitations

- No authentication. Do not expose this publicly without adding one (see
  `../deploy/aws/README.md`).
- Serves only the elastic-net models (the ones actually selected and
  externally validated), not the PyTorch comparison model or the other CV
  families.
- The out-of-distribution flag is a heuristic screen (same thresholds as
  `../monitor.py`), not a calibrated anomaly detector.
