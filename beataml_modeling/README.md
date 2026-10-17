# BEAT-AML expression-to-mutation feasibility study

Two independent binary tasks: FLT3-ITD and NPM1. Only the supplied BEAT-AML
counts, clinical workbook and mapping workbook are loaded. No external cohort
or normalized-expression matrix is used. Source files are never written.

## Run

From this directory, with the package versions in `requirements.txt`:

```bash
python3 -m unittest test_workflow.py
python3 -u run.py --phase develop
python3 audit_results.py
python3 -u run.py --phase finalize
python3 audit_results.py
python3 report.py
```

Development writes to `results/` and never predicts on the test partition.
It locks model selection for both targets before finalization. Finalization
refuses to run twice; it computes all test metrics and plots from a single
prediction call per target. To reproduce the *same prespecified analysis*, use
`--out results_reproduction` for both phases; this is not a license to search
test seeds or tune after examining test results. `--data-dir /path/to/files`
overrides the default parent directory. Keep `workflow.py` on the Python import
path when loading the saved joblib pipelines. Use `run.py` as the entry point.

## Cohort and features

- Match counts sample columns to mapping and clinical `dbgap_rnaseq_sample`.
- Exclude every mapping `rna_control` value other than `No`.
- Keep exact broad AML category and `Initial Diagnosis` disease stage.
- Drop missing target labels separately. Preserve specimen-specific outcomes.
- Group on `dbgap_subject_id`; no one-specimen selection or patient-label pooling.
- Use Ensembl `stable_id` as the predictor identity. Annotation symbols are
  only for literature-panel lookup and output interpretation. All corresponding
  Ensembl rows are retained if a panel symbol has multiple mappings.
- Drop genes zero across all 707 original counts samples, before splitting,
  as expressly allowed. No other global normalization or feature filtering.

Every nontrivial estimator is an sklearn Pipeline:

1. Per-specimen CPM over all remaining input genes, then natural `log1p`.
   This normalization is sample-local: it estimates no reference factor across
   samples and never uses another held-out specimen. Filtering follows CPM, so
   marker models do not normalize just to their small gene panels.
2. Fold-fitted expression/variance filter: CPM >=1 in at least 10% of training
   specimens and variance of log1p CPM >0.01. For tractability in this small-n
   feasibility run, retain at most 5,000 highest-training-variance eligible
   genes (stable tie order). This additional, prespecified cap is a limitation:
   low-variance predictive genes may be lost. No outcome-informed gene ranking.
3. StandardScaler fitted only on training observations.
4. Classifier, without resampling or class weights; prevalence is preserved.

## Prespecified model progression

- Majority floor: DummyClassifier(strategy='prior'). At threshold 0.50 it
  predicts the majority class, with the training prevalence as its constant
  probability. This gives meaningful AUROC/AP/Brier floor values rather than
  a hard zero-probability forecast.
- Marker logistic: L2 logistic, C in [0.1, 1, 10].
  NPM1: HOXA9, HOXA10, HOXB2, HOXB3, HOXB6, HOXD4, MEIS1.
  FLT3-ITD: PRDM16, IGFBP2, PTP4A3, CCNA1, IL2RA, SOCS2, DUSP6.
  These small literature-informed panels are not a claimed reproduction of
  a published classifier. They were fixed before fitting.
- Broad elastic-net logistic: saga, l1_ratio=0.5, C in [0.01, 0.05, 0.2],
  tolerance 0.001, max_iter 5,000. A convergence warning aborts the run.
- Random forest: 300 trees, sqrt features per split, minimum leaf size 2 or 5.
  Included instead of gradient boosting under the requested nonlinear option.
- No optional MLP: the study deliberately avoids unnecessary deep models at n~445.

Literature supporting the panels:

- NPM1 HOX/MEIS: https://pubmed.ncbi.nlm.nih.gov/17597811/
- Adult NPM1 expression: https://ashpublications.org/blood/article/106/12/3747/109726/Mutations-in-nucleophosmin-NPM1-in-acute-myeloid
- FLT3-ITD target genes: https://pmc.ncbi.nlm.nih.gov/articles/PMC4726916/
- SOCS2/IL2RA: https://www.nature.com/articles/s41598-019-45579-0
- DUSP6: https://pmc.ncbi.nlm.nih.gov/articles/PMC3464674/

## Evaluation

Global seed: 20260916. Each target uses the first split returned by a shuffled
5-fold StratifiedGroupKFold for the approximately 20% holdout. No favorable-fold
selection. Inside development, use 5 folds x 3 repetitions, seeds 20261016,
20261017, 20261018, identical folds for every candidate within a target.
Assertions check both classes and zero patient overlap in every split.

Choose each family's parameters, then the winning family, by mean fold
average precision; ties use AUROC then deterministic candidate/family order.
CV mean +/- sample SD uses 15 folds, not independent observations or a confidence
interval. Hyperparameter selection uses those same CV results, so these are
development/tuning estimates and potentially optimistic. Grouped nested CV is
the more rigorous internal alternative but was not implemented in this run;
the untouched holdout evaluates the complete selected procedure once.

Report AP as AUPRC (sklearn average_precision_score, noninterpolated, not
trapezoidal integration), AUROC, sensitivity, specificity, PPV, NPV and Brier.
Threshold is fixed in advance at 0.50, not optimized on validation or test labels.
PPV is undefined when there are no predicted positives; report NA, not zero.
All measures are specimen-weighted, with patient grouping for splitting.

Calibration plots use five quantile bins and uncalibrated model probabilities.
CV calibration plots use repeat 1 only, once per development specimen, for the
selected hyperparameters of each family; they are descriptive tuning plots.
Held-out plots show only the chosen family. Small test bins are noisy. No
calibration model is fitted to held-out predictions.

Only the winning family receives held-out metrics for each target. All families
receive CV metrics. This follows the one-selected-model holdout requirement;
there is intentionally no test-set leaderboard.

## Interpretation and limitations

Feature lists and plots use development-only fitted models. Elastic-net weights
are log-odds coefficients per one training SD of log1p CPM, conditional on other
genes. Forest importances are training impurity importances, unsigned, not
held-out permutation effects or causal evidence. Correlated genes can exchange
weights/importance; no stable or externally validated biomarker claim is made.
Panel membership is evidence of prior biological support, not proof that other
features are novel discoveries. Mutation co-occurrence, lineage, blast content,
sample source and cohort batches may explain part of the signal. Random
patient-grouped splits do not demonstrate transportability across these shifts.

See `results/summary.txt` for the narrative and complete CV/holdout table;
`cv_all_candidates.csv` contains every candidate and fold, and the split/CV
manifests permit independent patient-overlap checks. The small toy tests only
verify code invariants; no synthetic observations are used in the study.
