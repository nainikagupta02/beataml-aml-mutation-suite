"""Out-of-distribution check: flag a sample before trusting its prediction, using
only the frozen pipeline's own standardized features (no model is fit, retrained
or changed here). Reusable as a library (standardize/ood_flags) for a future
serving layer, and runnable as a CLI report over the two known cohorts."""
import os
import sys
sys.dont_write_bytecode = True
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BEAT = ROOT.parent / 'beataml_modeling'
TCGA = ROOT.parent / 'tcga_laml_validation'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/beataml-mpl')
sys.path.insert(0, str(BEAT))
import argparse
import json
import joblib
import numpy as np
import pandas as pd
import workflow

OUT = ROOT / 'results'
Z_WARN, Z_CRITICAL, WARN_FRACTION = 5.0, 10.0, 0.01


def standardize(pipe, X, available=None):
    """Apply the frozen log_cpm/filter/scale steps; optionally zero out genes
    unavailable in an external cohort, matching tcga_laml_validation's predict_frozen."""
    log_values = pipe.named_steps['log_cpm'].transform(X)
    selected = pipe.named_steps['filter'].transform(log_values)
    z = pipe.named_steps['scale'].transform(selected)
    if available is not None:
        missing = ~available[pipe.named_steps['filter'].indices_]
        z[:, missing] = 0.0
    return z


def ood_flags(z, coefficients):
    """Flag rows whose standardized values, on the genes this model actually weights,
    look unlike the BEAT-AML training distribution. Thresholds match the diagnostic
    already used in tcga_laml_validation/validate.py (inference_qc.csv), generalized
    to a per-sample flag instead of a cohort-level aggregate."""
    active = coefficients != 0
    za = np.abs(z[:, active])
    warn_fraction = (za > Z_WARN).mean(axis=1)
    max_z = za.max(axis=1)
    flag = (warn_fraction > WARN_FRACTION) | (max_z > Z_CRITICAL)
    return pd.DataFrame({'warn_fraction_abs_z_gt_5': warn_fraction, 'max_abs_z': max_z, 'flag_ood': flag})


def score(target, X, ids, available=None, label=''):
    pipe = joblib.load(BEAT / f'results/{target}_elastic_net.joblib')
    z = standardize(pipe, X, available)
    coefficients = pipe.named_steps['model'].coef_[0]
    frame = ood_flags(z, coefficients)
    frame.insert(0, 'sample', ids)
    frame.insert(0, 'target', target)
    frame.insert(0, 'cohort', label)
    return frame


def beat_heldout(target):
    X = np.load(BEAT / 'results/counts_input.npy', allow_pickle=False)
    cohort = pd.read_csv(BEAT / 'results/cohort.csv')
    rows = json.loads((BEAT / f'results/{target}_indices.json').read_text())['test_rows']
    return X[rows], cohort.iloc[rows].dbgap_rnaseq_sample.to_numpy()


def tcga_external(target):
    X = np.load(TCGA / 'results/aligned_raw_counts.npy', allow_pickle=False)
    available = np.load(TCGA / 'results/available_genes.npy', allow_pickle=False)
    cohort = pd.read_csv(TCGA / 'results/cohort.csv')
    valid = cohort[target].isin(['positive', 'negative']).to_numpy()
    return X[valid], cohort.loc[valid, 'patient'].to_numpy(), available


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', choices=['beataml_heldout', 'tcga', 'both'], default='both')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    frames = []
    for target in workflow.TARGETS:
        if args.cohort in ('beataml_heldout', 'both'):
            X, ids = beat_heldout(target)
            frames.append(score(target, X, ids, available=None, label='beataml_heldout'))
        if args.cohort in ('tcga', 'both'):
            X, ids, available = tcga_external(target)
            frames.append(score(target, X, ids, available=available, label='tcga_external'))
    report = pd.concat(frames, ignore_index=True)
    report.to_csv(OUT / 'monitor_flags.csv', index=False)
    summary = report.groupby(['target', 'cohort']).agg(
        n=('flag_ood', 'size'), flagged=('flag_ood', 'sum'),
        flagged_fraction=('flag_ood', 'mean'), max_abs_z=('max_abs_z', 'max')).reset_index()
    summary.to_csv(OUT / 'monitor_summary.csv', index=False)
    print(summary.to_string(index=False))
    print(
        '\nPer-sample flags: results/monitor_flags.csv. A flagged sample is not automatically '
        "wrong; it means its gene-expression profile is unlike BEAT-AML training data on the "
        'genes this frozen model actually uses, so its predicted probability should be treated '
        'with more caution before acting on it. The BEAT-AML held-out rows are run as a sanity '
        'baseline: an in-distribution cohort should flag rarely.'
    )


if __name__ == '__main__':
    main()
