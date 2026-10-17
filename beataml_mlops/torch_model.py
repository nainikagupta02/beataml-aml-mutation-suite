"""Small PyTorch MLP for FLT3-ITD and NPM1, trained on the exact BEAT-AML
development rows/CV folds used for the frozen elastic-net study, then compared
to it on the BEAT-AML held-out test partition and on the external TCGA-LAML
cohort. Reuses workflow.py's LogCPM/FoldGeneFilter classes (fold-fit, no
leakage) so the feature pipeline is identical; the elastic-net models and
all of beataml_modeling/tcga_laml_validation's saved results are untouched.

This is a feasibility comparison at n~350 training specimens with up to 5,000
input features: a small, heavily regularized network, not a claim that deep
learning outperforms the frozen elastic-net here. See README.md for caveats,
in particular that evaluating a second model on the SAME locked held-out
partition is a secondary analysis, not an independent held-out test.
"""
import os
import sys
sys.dont_write_bytecode = True
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BEAT = ROOT.parent / 'beataml_modeling'
TCGA = ROOT.parent / 'tcga_laml_validation'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/beataml-mpl')
sys.path.insert(0, str(BEAT))
import json
import random
import time
import joblib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, precision_recall_curve
from sklearn.calibration import calibration_curve
import workflow  # original frozen LogCPM/FoldGeneFilter/metrics; never retrains the elastic-net models

OUT = ROOT / 'results'
SEED = workflow.SEED
TARGETS = workflow.TARGETS
METRICS = workflow.METRICS
MAX_GENES = 5000           # same cap as the broad elastic-net/random-forest filter, for a fair comparison
HIDDEN = (64, 16)
DROPOUT = 0.3
WEIGHT_DECAY = 1e-3
LR = 1e-3
EPOCHS = 60
BATCH_SIZE = 32


class MLP(nn.Module):
    def __init__(self, n_features):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, HIDDEN[0]), nn.ReLU(), nn.Dropout(DROPOUT),
            nn.Linear(HIDDEN[0], HIDDEN[1]), nn.ReLU(), nn.Dropout(DROPOUT),
            nn.Linear(HIDDEN[1], 1),
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def train_mlp(X, y, seed):
    set_seed(seed)
    model = MLP(X.shape[1])
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    loss_fn = nn.BCEWithLogitsLoss()
    Xt = torch.tensor(X, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32)
    n = len(y)
    model.train()
    for _ in range(EPOCHS):
        perm = torch.randperm(n)
        for i in range(0, n, BATCH_SIZE):
            idx = perm[i:i + BATCH_SIZE]
            opt.zero_grad()
            loss = loss_fn(model(Xt[idx]), yt[idx])
            loss.backward()
            opt.step()
    model.eval()
    return model


def predict_proba(model, X):
    with torch.no_grad():
        return torch.sigmoid(model(torch.tensor(X, dtype=torch.float32))).numpy()


def fit_transform(Xtr_raw, Xva_raw=None):
    """Fold-fit (or final-fit) feature pipeline: LogCPM -> FoldGeneFilter -> StandardScaler,
    fit on Xtr_raw only. Mirrors workflow.py's broad-model pipeline exactly."""
    log_cpm = workflow.LogCPM().fit(Xtr_raw)
    Ltr = log_cpm.transform(Xtr_raw)
    filt = workflow.FoldGeneFilter(max_genes=MAX_GENES).fit(Ltr)
    Ftr = filt.transform(Ltr)
    scaler = StandardScaler().fit(Ftr)
    Str = scaler.transform(Ftr).astype(np.float32)
    if Xva_raw is None:
        return Str, log_cpm, filt, scaler
    Lva = log_cpm.transform(Xva_raw)
    Sva = scaler.transform(filt.transform(Lva)).astype(np.float32)
    return Str, Sva, log_cpm, filt, scaler


def load_target(target):
    X = np.load(BEAT / 'results/counts_input.npy', allow_pickle=False)
    cohort = pd.read_csv(BEAT / 'results/cohort.csv')
    idx = json.loads((BEAT / f'results/{target}_indices.json').read_text())
    dev, test = np.array(idx['development_rows']), np.array(idx['test_rows'])
    y_dev = cohort.loc[dev, target].eq('positive').to_numpy(dtype=int)
    y_test = cohort.loc[test, target].eq('positive').to_numpy(dtype=int)
    groups_dev = cohort.loc[dev, 'dbgap_subject_id'].to_numpy()
    return X[dev], y_dev, groups_dev, X[test], y_test, dev, test, cohort


def cross_validate(target, xd, yd, gd):
    rows = []
    for repeat in range(3):
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED + 100 + repeat)
        for fold, (tr, va) in enumerate(cv.split(xd, yd, gd)):
            workflow.assert_group_split(tr, va, gd, yd)
            start = time.time()
            Str, Sva, *_ = fit_transform(xd[tr], xd[va])
            model = train_mlp(Str, yd[tr], seed=SEED + 100 + repeat)
            p = predict_proba(model, Sva)
            row = {'target': target, 'model': 'pytorch_mlp', 'repeat': repeat, 'fold': fold,
                   'n_train': len(tr), 'n_valid': len(va), **workflow.metrics(yd[va], p)}
            rows.append(row)
            print(f'{target} pytorch_mlp repeat {repeat + 1} fold {fold + 1}: '
                  f'AUROC={row["AUROC"]:.3f} AUPRC={row["AUPRC"]:.3f} ({time.time() - start:.1f}s)', flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / f'{target}_pytorch_cv_folds.csv', index=False)
    summary = {'target': target, 'model': 'pytorch_mlp'}
    for metric in METRICS:
        summary[f'{metric}_mean'] = float(frame[metric].mean())
        summary[f'{metric}_sd'] = float(frame[metric].std(ddof=1))
    return summary


def finalize_and_evaluate(target, xd, yd, xtest, ytest):
    heldout_path = OUT / f'{target}_pytorch_heldout_metrics.json'
    if heldout_path.exists():
        raise FileExistsError(f'{heldout_path} exists; PyTorch held-out evaluation for {target} already ran once.')
    Str, log_cpm, filt, scaler = fit_transform(xd)
    model = train_mlp(Str, yd, seed=SEED)
    joblib.dump({'filter': filt, 'scaler': scaler, 'log_cpm': log_cpm}, OUT / f'{target}_pytorch_preprocess.joblib')
    torch.save(model.state_dict(), OUT / f'{target}_pytorch_model.pt')
    Stest = scaler.transform(filt.transform(log_cpm.transform(xtest))).astype(np.float32)
    p_test = predict_proba(model, Stest)
    held = {'target': target, 'model': 'pytorch_mlp', 'n': len(ytest), **workflow.metrics(ytest, p_test)}
    heldout_path.write_text(json.dumps(held, indent=2) + '\n')
    return held, model, filt, scaler, log_cpm


def evaluate_tcga(target, model, filt, scaler, log_cpm):
    tcga_path = TCGA / 'results'
    X = np.load(tcga_path / 'aligned_raw_counts.npy', allow_pickle=False)
    available = np.load(tcga_path / 'available_genes.npy', allow_pickle=False)
    cohort = pd.read_csv(tcga_path / 'cohort.csv')
    valid = cohort[target].isin(['positive', 'negative']).to_numpy()
    y = cohort.loc[valid, target].eq('positive').to_numpy(dtype=int)
    Z = scaler.transform(filt.transform(log_cpm.transform(X[valid]))).astype(np.float32)
    missing = ~available[filt.indices_]
    Z[:, missing] = 0.0  # same frozen zero-impute rule used for the elastic-net external evaluation
    p = predict_proba(model, Z)
    result = {'target': target, 'model': 'pytorch_mlp', 'n': len(y),
              'present_selected_genes': int((~missing).sum()), 'selected_genes': len(filt.indices_),
              **workflow.metrics(y, p)}
    return result, y, p


def comparison_plot(all_tcga_preds):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for row, target in enumerate(TARGETS):
        y, p = all_tcga_preds[target]
        actual, pred = calibration_curve(y, p, n_bins=5, strategy='quantile')
        axes[row, 0].plot(pred, actual, 'o-', label='pytorch_mlp')
        axes[row, 0].plot([0, 1], [0, 1], 'k--', alpha=0.4)
        axes[row, 0].set(xlim=(0, 1), ylim=(0, 1), xlabel='Mean predicted probability',
                          ylabel='Observed positive fraction', title=f'{target}: TCGA calibration')
        axes[row, 0].legend(fontsize=8)
        fpr, tpr, _ = roc_curve(y, p)
        axes[row, 1].plot(fpr, tpr)
        axes[row, 1].plot([0, 1], [0, 1], 'k--', alpha=0.4)
        axes[row, 1].set(xlabel='False positive rate', ylabel='Sensitivity', title=f'{target}: TCGA ROC')
    fig.suptitle('PyTorch MLP on external TCGA-LAML (same frozen evaluation as the elastic-net models)')
    fig.tight_layout()
    fig.savefig(OUT / 'pytorch_tcga_diagnostics.png', dpi=170)
    plt.close(fig)


def write_comparison(cv_summaries, heldout_results, tcga_results):
    beat_cv = pd.read_csv(BEAT / 'results/cv_summary.csv')
    beat_heldout = pd.read_csv(BEAT / 'results/heldout_metrics.csv')
    tcga_elastic = pd.read_csv(TCGA / 'results/metrics.csv')
    lines = ['PyTorch MLP vs frozen elastic-net: BEAT-AML CV, BEAT-AML held-out, TCGA-LAML external', '']
    for target in TARGETS:
        en_cv = beat_cv[(beat_cv.target == target) & (beat_cv.model == 'elastic_net')].iloc[0]
        nn_cv = cv_summaries[target]
        en_held = beat_heldout[beat_heldout.target == target].iloc[0]
        nn_held = heldout_results[target]
        en_tcga = tcga_elastic[tcga_elastic.target == target].iloc[0]
        nn_tcga = tcga_results[target]
        lines.append(f'{target}:')
        lines.append(f'  CV AUPRC:      elastic-net {en_cv.AUPRC_mean:.3f} +/- {en_cv.AUPRC_sd:.3f}   '
                     f'pytorch_mlp {nn_cv["AUPRC_mean"]:.3f} +/- {nn_cv["AUPRC_sd"]:.3f}')
        lines.append(f'  CV AUROC:      elastic-net {en_cv.AUROC_mean:.3f} +/- {en_cv.AUROC_sd:.3f}   '
                     f'pytorch_mlp {nn_cv["AUROC_mean"]:.3f} +/- {nn_cv["AUROC_sd"]:.3f}')
        lines.append(f'  BEAT held-out AUROC: elastic-net {en_held.AUROC:.3f}   pytorch_mlp {nn_held["AUROC"]:.3f}')
        lines.append(f'  TCGA AUROC:          elastic-net {en_tcga.AUROC:.3f}   pytorch_mlp {nn_tcga["AUROC"]:.3f}')
        lines.append(f'  TCGA AUPRC:          elastic-net {en_tcga.AUPRC:.3f}   pytorch_mlp {nn_tcga["AUPRC"]:.3f}')
        lines.append(f'  TCGA Brier:          elastic-net {en_tcga.Brier:.3f}   pytorch_mlp {nn_tcga["Brier"]:.3f}')
        delta = nn_tcga['AUPRC'] - en_tcga.AUPRC
        lines.append(f'  pytorch_mlp {"exceeded" if delta > 0 else "did not beat"} elastic-net on TCGA AUPRC '
                     f'(difference {delta:+.3f}); descriptive only, not a significance test.')
        lines.append('')
    lines += [
        f'Architecture: MLP({HIDDEN[0]}, {HIDDEN[1]}), ReLU, dropout {DROPOUT}, up to {MAX_GENES} genes '
        f'(same fold-fit filter as the broad elastic-net/random-forest models), Adam lr={LR}, '
        f'weight_decay={WEIGHT_DECAY}, {EPOCHS} epochs, batch size {BATCH_SIZE}. One fixed configuration,',
        'prespecified before fitting; no architecture/hyperparameter search was run, unlike the elastic-net grid.',
        'Caveats:',
        '- n~350-360 training specimens with up to 5,000 features is a severe small-n/high-p regime for a',
        '  neural network; this is a feasibility comparison, not a claim that deep learning helps here.',
        '- The BEAT-AML held-out numbers reuse the SAME locked test partition the elastic-net study already',
        '  evaluated once. Evaluating a second model on that partition is a secondary analysis: it is still a',
        '  held-out comparison for THIS model (never used for its own fitting/filtering/tuning), but it is not',
        '  an independent test of the overall two-model comparison procedure.',
        '- The TCGA comparison is the primary, cleaner comparison: both models receive the identical frozen',
        '  external cohort and zero-imputation rule for genes absent from TCGA, with no refitting either way.',
        '- No calibration, threshold or architecture tuning was performed after seeing TCGA or held-out results.',
    ]
    (OUT / 'pytorch_vs_elastic_net.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    OUT.mkdir(exist_ok=True)
    cv_summaries, heldout_results, tcga_results, tcga_preds = {}, {}, {}, {}
    for target in TARGETS:
        xd, yd, gd, xtest, ytest, dev_rows, test_rows, cohort = load_target(target)
        cv_summaries[target] = cross_validate(target, xd, yd, gd)
        pd.DataFrame([cv_summaries[target]]).to_csv(OUT / f'{target}_pytorch_cv_summary.csv', index=False)
        held, model, filt, scaler, log_cpm = finalize_and_evaluate(target, xd, yd, xtest, ytest)
        heldout_results[target] = held
        tcga_result, y_tcga, p_tcga = evaluate_tcga(target, model, filt, scaler, log_cpm)
        tcga_results[target] = tcga_result
        tcga_preds[target] = (y_tcga, p_tcga)
        pd.DataFrame([tcga_result]).to_csv(OUT / f'{target}_pytorch_tcga_metrics.csv', index=False)
        print(f'{target}: held-out AUROC={held["AUROC"]:.3f}, TCGA AUROC={tcga_result["AUROC"]:.3f}', flush=True)
    comparison_plot(tcga_preds)
    write_comparison(cv_summaries, heldout_results, tcga_results)


if __name__ == '__main__':
    main()
