"""Audit saved outputs without refitting or making new test predictions."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from workflow import TARGETS, METRICS, metrics, digest

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results')
    parser.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parent.parent)
    args=parser.parse_args(); out=args.out
    splits=pd.read_csv(out/'split_manifest.csv')
    folds=pd.read_csv(out/'cv_manifest.csv')
    cv=pd.read_csv(out/'cv_all_candidates.csv')
    summary=pd.read_csv(out/'cv_summary.csv')
    oof=pd.read_csv(out/'oof_predictions.csv')
    selection=json.loads((out/'selection.json').read_text())
    cohort=pd.read_csv(out/'cohort.csv')
    assert cohort.dbgap_rnaseq_sample.is_unique
    checks=[]
    for target in TARGETS:
        s=splits[splits.target==target]
        dev=s[s.partition=='development']; test=s[s.partition=='test']
        assert set(dev.patient).isdisjoint(test.patient)
        assert len(s)==len(cohort[cohort[target].notna()])
        assert set(s['sample'])==set(cohort.loc[cohort[target].notna(),'dbgap_rnaseq_sample'])
        f=folds[folds.target==target]
        for repeat in range(3):
            r=f[f.repeat==repeat]
            assert r['sample'].is_unique and set(r['sample'])==set(dev['sample'])
            assert (r.groupby('patient').fold.nunique()==1).all()
            assert set(r.fold)==set(range(5))
        for _,chosen in summary[summary.target==target].iterrows():
            scores=cv[(cv.target==target)&(cv.model==chosen.model)&(cv.candidate==chosen.candidate)]
            assert len(scores)==15 and (scores.group_overlap==0).all()
            predictions=oof[(oof.target==target)&(oof.model==chosen.model)]
            assert len(predictions)==3*len(dev)
            assert set(predictions['sample'])==set(dev['sample'])
            for (repeat,fold),assignment in f.groupby(['repeat','fold']):
                p=predictions[(predictions.repeat==repeat)&predictions['sample'].isin(assignment['sample'])]
                assert len(p)==len(assignment)
                calculated=metrics(p.y,p.probability)
                saved=scores[(scores.repeat==repeat)&(scores.fold==fold)].iloc[0]
                for metric in METRICS:
                    np.testing.assert_allclose(calculated[metric],saved[metric],rtol=1e-10,atol=1e-10,equal_nan=True)
            for metric in METRICS:
                np.testing.assert_allclose(scores[metric].mean(),chosen[metric+'_mean'],atol=1e-10,equal_nan=True)
                np.testing.assert_allclose(scores[metric].std(ddof=1),chosen[metric+'_sd'],atol=1e-10,equal_nan=True)
        winner=summary[summary.target==target].sort_values(['AUPRC_mean','AUROC_mean'],ascending=False).iloc[0]
        assert selection[target]['model']==winner.model
        checks.append(f'{target}: no patient overlap; all 15 folds partition development correctly; OOF metrics reconciled; selection matches CV.')
    for name,expected in json.loads((out/'data_hashes.json').read_text()).items():
        assert digest(args.data_dir/name)==expected['sha256']
    checks.append('All three source files match their pre-analysis SHA256 hashes.')
    if (out/'heldout_metrics.csv').exists():
        test_predictions=pd.read_csv(out/'heldout_predictions.csv')
        test_metrics=pd.read_csv(out/'heldout_metrics.csv')
        for target in TARGETS:
            p=test_predictions[test_predictions.target==target]
            expected=splits[(splits.target==target)&(splits.partition=='test')]
            assert p['sample'].is_unique and set(p['sample'])==set(expected['sample'])
            saved=test_metrics[test_metrics.target==target].iloc[0]
            assert set(p.model)=={selection[target]['model']}
            actual=metrics(p.y,p.probability)
            for key,value in actual.items():
                np.testing.assert_allclose(value,saved[key],atol=1e-10,equal_nan=True)
        assert digest(out/'selection.json')==json.loads((out/'test_evaluation_started.json').read_text())['selection_sha256']
        checks.append('Held-out metrics reconcile from saved predictions; model selections are unchanged since evaluation began.')
    (out/'audit.txt').write_text('\n'.join(checks)+'\n')
    print('\n'.join(checks))

if __name__=='__main__': main()
