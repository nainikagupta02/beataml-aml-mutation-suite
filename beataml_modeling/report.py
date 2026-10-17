"""Render existing saved metrics as plain text. Does not fit or predict."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results')
    out=parser.parse_args().out
    cv=pd.read_csv(out/'cv_summary.csv')
    test=pd.read_csv(out/'heldout_metrics.csv')
    flow=json.loads((out/'cohort_flow.json').read_text())
    splits=pd.read_csv(out/'split_manifest.csv')
    lines=[]
    for _,r in test.iterrows():
        lines.append(f'{r.target}: selected {r.model}; held-out AUROC {r.AUROC:.3f}, AUPRC {r.AUPRC:.3f}, '
                     f'sensitivity {r.sensitivity:.3f}, specificity {r.specificity:.3f}, Brier {r.Brier:.3f}.')
    lines += ['', 'All thresholds were fixed at 0.50 before fitting. AUPRC means average precision.',
              'Test results are from one locked evaluation per target. Other models were not tested.',
              'CV entries below are mean +/- sample SD over 15 grouped folds (5 folds x 3 repeats).',
              'CV selected hyperparameters: these are tuning estimates, not nested-CV estimates or confidence intervals.', '']
    columns=['AUROC','AUPRC','sensitivity','specificity','PPV','NPV','Brier']
    display=[]
    for _,r in cv.iterrows():
        d={'target':r.target,'model':r.model,'evaluation':'CV'}
        for k in columns:
            d[k]='NA' if pd.isna(r[k+'_mean']) else f'{r[k+"_mean"]:.3f} +/- {r[k+"_sd"]:.3f}'
        display.append(d)
    for _,r in test.iterrows():
        display.append({'target':r.target,'model':r.model,'evaluation':'holdout',
                        **{k:'NA' if pd.isna(r[k]) else f'{r[k]:.3f}' for k in columns}})
    table=pd.DataFrame(display)
    table.to_csv(out/'results_table.csv',index=False)
    lines.extend([table.to_string(index=False), '',
                  'NA PPV means there were no positive predictions. Undefined-fold counts are in cv_summary.csv.',
                  'Majority probabilities equal training prevalence; classifications are all negative at 0.50.', '',
                  'Cohort: '+json.dumps(flow), ''])
    for (target,partition),g in splits.groupby(['target','partition'],sort=False):
        lines.append(f'{target} {partition}: {len(g)} specimens, {g.patient.nunique()} patients; '
                     f'{g.y.sum()} positive / {len(g)-g.y.sum()} negative.')
    lines += ['', 'Held-out confusion matrices at 0.50:']
    for _,r in test.iterrows():
        lines.append(f'{r.target}: TP={r.TP}, FN={r.FN}, TN={r.TN}, FP={r.FP}; '
                     f'PPV={r.PPV:.3f}, NPV={r.NPV:.3f}.')
    lines += ['', 'Model comparison (development CV, not a test leaderboard):']
    for target in ['FLT3-ITD','NPM1']:
        v=cv[cv.target==target].set_index('model')
        delta=v.loc['random_forest','AUPRC_mean']-v.loc['elastic_net','AUPRC_mean']
        direction='exceeded' if delta>0 else 'did not beat'
        lines.append(f'{target}: random forest {direction} elastic-net mean CV AUPRC '
                     f'({v.loc["random_forest","AUPRC_mean"]:.3f} vs '
                     f'{v.loc["elastic_net","AUPRC_mean"]:.3f}; difference {delta:+.3f}). '
                     'This is a descriptive comparison, not evidence of statistical significance.')
    lines += ['The optional MLP was not run; no deep-learning superiority claim is made.', '',
              'Top genes below come from development-only fits. All full lists retain Ensembl IDs.',
              'Elastic-net coefficients: log odds per training SD; signs are conditional associations.',
              'Forest importances: unsigned training impurity importance; biased and unstable with correlated genes.',
              'Prior-panel membership is annotated in feature CSVs; unlisted genes are not proven novel biomarkers.', '']
    for target in ['FLT3-ITD','NPM1']:
        for family in ['marker_logistic','elastic_net','random_forest']:
            features=pd.read_csv(out/f'{target}_{family}_features.csv')
            lines.append(f'{target} / {family}:')
            lines.append(features[['stable_id','display_label','weight','literature_panel_member']].head(15).to_string(index=False))
            lines.append('')
    notes=out.parent/'biology_notes.txt'
    if notes.exists(): lines += [notes.read_text(), '']
    lines += ['Limitations:',
              '- One cohort, 445 specimens from 443 patients, no external validation; no diagnostic-use claim.',
              '- Small holdouts and few positive cases make sensitivity and calibration imprecise.',
              '- Repeated CV folds overlap; fold SD is not a confidence interval. Tuning CV can be optimistic.',
              '- Grouped nested CV would better estimate the selection procedure internally; not run due to added fitting cost.',
              '- Fixed 0.50 threshold is not a clinically optimized operating point; PPV/NPV depend on prevalence.',
              '- Cohort wave, marrow/blood source, cell composition and co-mutations can confound gene associations.',
              '- Prespecified 5,000-gene variance cap may lose useful low-variance genes; CPM is not composition-robust normalization.',
              '- Importance lists are exploratory, not causal, stable or externally validated biomarker panels.',
              '- No model was refit on the held-out specimens, and no clinical/mutation-derived field entered X.', '',
              'Reproduction: see ../README.md, ../requirements.txt, protocol.json, data_hashes.json,',
              'split_manifest.csv, cv_manifest.csv and cv_all_candidates.csv.',
              'Plots: calibration_cv.png, heldout_diagnostics.png, and *_features.png.']
    report='\n'.join(lines)+'\n'
    (out/'summary.txt').write_text(report)
    print(report)

if __name__=='__main__': main()
