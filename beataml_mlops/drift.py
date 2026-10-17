"""Compare BEAT-AML (training) vs TCGA-LAML (external) gene-expression distributions
in the exact frozen LogCPM space, and translate the shift into its effect on each
target's frozen linear predictor. No model is fit, retrained or modified here;
inputs are the already-saved counts/coverage arrays from both sibling projects."""
import os
import sys
sys.dont_write_bytecode = True
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BEAT = ROOT.parent / 'beataml_modeling'
TCGA = ROOT.parent / 'tcga_laml_validation'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/beataml-mpl')
sys.path.insert(0, str(BEAT))
import joblib
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
import workflow

OUT = ROOT / 'results'


def load_logcpm_matrices():
    X_beat = np.load(BEAT / 'results/counts_input.npy', allow_pickle=False)
    X_tcga = np.load(TCGA / 'results/aligned_raw_counts.npy', allow_pickle=False)
    available = np.load(TCGA / 'results/available_genes.npy', allow_pickle=False)
    annotation = pd.read_csv(BEAT / 'results/gene_annotation.csv')
    assert X_beat.shape[1] == X_tcga.shape[1] == len(annotation) == len(available)
    log_cpm = workflow.LogCPM().fit(X_beat)  # stateless beyond gene count; same transform as the frozen pipelines
    return log_cpm.transform(X_beat), log_cpm.transform(X_tcga), annotation, available


def gene_drift_table(L_beat, L_tcga, annotation, available):
    mean_b, sd_b = L_beat.mean(axis=0), L_beat.std(axis=0, ddof=1)
    mean_t, sd_t = L_tcga.mean(axis=0), L_tcga.std(axis=0, ddof=1)
    pooled_sd = np.sqrt((sd_b.astype(np.float64) ** 2 + sd_t.astype(np.float64) ** 2) / 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        smd = np.where(pooled_sd > 0, (mean_t - mean_b) / pooled_sd, np.nan)
    table = annotation.copy()
    table['present_in_TCGA'] = available
    table['beat_mean_log1p_cpm'] = mean_b
    table['tcga_mean_log1p_cpm'] = mean_t
    table['standardized_mean_difference'] = smd
    table['large_shift'] = table.present_in_TCGA & (np.abs(table.standardized_mean_difference) > 1.0)
    return table


def predictor_shift(target, table, annotation, L_beat, L_tcga):
    """Decompose the average change in the frozen logit into a per-gene contribution:
    coefficient * (TCGA mean, standardized by the frozen scaler fit on BEAT-AML dev data).
    Because the BEAT-AML standardized mean is 0 by construction, this sum IS the average
    shift in log-odds caused purely by the gene-expression distribution shift, holding the
    frozen coefficients and threshold fixed. Absent genes are zeroed, matching the frozen
    external-inference rule (tcga_laml_validation/validate.py:predict_frozen)."""
    pipe = joblib.load(BEAT / f'results/{target}_elastic_net.joblib')
    scaler = pipe.named_steps['scale']
    indices = pipe.named_steps['filter'].indices_
    coefficients = pipe.named_steps['model'].coef_[0]
    present = table.present_in_TCGA.to_numpy()[indices]
    tcga_mean = table.tcga_mean_log1p_cpm.to_numpy()[indices]
    standardized_tcga_mean = (tcga_mean - scaler.mean_) / scaler.scale_
    standardized_tcga_mean = np.where(present, standardized_tcga_mean, 0.0)
    contribution = coefficients * standardized_tcga_mean
    selected = annotation.iloc[indices].copy()
    selected['coefficient'] = coefficients
    selected['present_in_TCGA'] = present
    selected['standardized_tcga_mean'] = standardized_tcga_mean
    selected['logit_contribution'] = contribution
    ks_stat = np.full(len(indices), np.nan)
    ks_p = np.full(len(indices), np.nan)
    for row, gi in enumerate(indices):
        if present[row]:
            ks_stat[row], ks_p[row] = ks_2samp(L_beat[:, gi], L_tcga[:, gi])
    selected['ks_statistic'] = ks_stat
    selected['ks_pvalue'] = ks_p
    selected = selected.reindex(selected.logit_contribution.abs().sort_values(ascending=False).index)
    return float(contribution.sum()), selected


def main():
    OUT.mkdir(exist_ok=True)
    L_beat, L_tcga, annotation, available = load_logcpm_matrices()
    table = gene_drift_table(L_beat, L_tcga, annotation, available)
    table.to_csv(OUT / 'gene_drift_all.csv', index=False)

    heldout = pd.read_csv(BEAT / 'results/heldout_metrics.csv')
    tcga_metrics = pd.read_csv(TCGA / 'results/metrics.csv')
    present_count = int(available.sum())
    large = int(table.large_shift.sum())

    lines = [
        'BEAT-AML (training cohort, n=445) vs TCGA-LAML (external, n=151) gene-expression drift',
        'Natural log1p(CPM) space; same frozen LogCPM transform used for training and external inference.',
        '',
        f'{present_count}/{len(annotation)} genes present in both cohorts by Ensembl stable ID; '
        f'{large} of those ({large / present_count:.1%}) have |standardized mean difference| > 1.0.',
        '',
    ]
    for target in workflow.TARGETS:
        total_shift, selected = predictor_shift(target, table, annotation, L_beat, L_tcga)
        selected.to_csv(OUT / f'{target}_drift_contribution.csv', index=False)
        held = heldout[heldout.target == target].iloc[0]
        ext = tcga_metrics[tcga_metrics.target == target].iloc[0]
        lines.append(
            f'{target}: average TCGA-driven shift in the frozen logit = {total_shift:+.3f} '
            '(from gene-distribution shift alone, frozen coefficients/threshold unchanged).'
        )
        lines.append(
            f'  BEAT-AML held-out AUROC {held.AUROC:.3f} -> TCGA AUROC {ext.AUROC:.3f}; '
            f'Brier {held.Brier:.3f} -> {ext.Brier:.3f}.'
        )
        top = selected.head(5)[['stable_id', 'display_label', 'coefficient', 'standardized_tcga_mean', 'logit_contribution']]
        lines.append('  Top 5 selected genes by |contribution to the shift|:')
        lines.append('  ' + top.to_string(index=False).replace('\n', '\n  '))
        lines.append('')

    lines += [
        'Interpretation: standardized mean difference and the KS test are descriptive, not causal;',
        'they cannot separate sequencing/annotation-pipeline differences from real population shift.',
        'The logit-shift decomposition assumes the frozen coefficients remain valid out-of-cohort;',
        'it explains average bias in the linear predictor, not per-patient calibration.',
        'Genes absent from TCGA contribute exactly 0 here, matching the frozen external-inference rule;',
        'that rule hides their drift rather than resolving it (see tcga_laml_validation coverage CSVs).',
        'This is descriptive monitoring, not a recalibration, domain-adaptation or retraining procedure.',
    ]
    (OUT / 'drift_summary.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
