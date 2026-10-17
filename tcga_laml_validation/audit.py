"""Reconcile saved external outputs without invoking the classifiers again."""
import sys
sys.dont_write_bytecode=True
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
from validate import OUT,BEAT,TARGETS,sha,check_frozen_files

def main():
    checks=[]
    n=check_frozen_files();checks.append(f'All {n} original BEAT-AML files are unchanged; no files added or removed.')
    protocol=json.loads((OUT/'protocol.json').read_text())
    for target,h in protocol['models'].items():assert sha(BEAT/'results'/f'{target}_elastic_net.joblib')==h
    for name,h in json.loads((OUT/'prepared_hashes.json').read_text()).items():assert sha(OUT/name)==h
    assert sha(OUT/'predictions.csv')==json.loads((OUT/'evaluation_completed.json').read_text())['predictions_sha256']
    predictions=pd.read_csv(OUT/'predictions.csv');cohort=pd.read_csv(OUT/'cohort.csv')
    metrics=pd.read_csv(OUT/'metrics.csv');coverage=pd.read_csv(OUT/'coverage_summary.csv')
    assert cohort.patient.is_unique and len(cohort)==151
    for target in TARGETS:
        p=predictions[predictions.target==target];m=metrics[metrics.target==target].iloc[0]
        assert p.patient.is_unique and set(p.patient)==set(cohort.patient)
        truth=cohort.set_index('patient').loc[p.patient,target].eq('positive').to_numpy(dtype=int)
        np.testing.assert_array_equal(p.y,truth)
        for key,value in {'AUROC':roc_auc_score(truth,p.probability),
                          'AUPRC':average_precision_score(truth,p.probability),
                          'Brier':brier_score_loss(truth,p.probability)}.items():
            np.testing.assert_allclose(m[key],value,atol=1e-12)
        pred=p.probability.to_numpy()>=.5
        for key,value in {'TP':((truth==1)&pred).sum(),'FN':((truth==1)&~pred).sum(),
                          'TN':((truth==0)&~pred).sum(),'FP':((truth==0)&pred).sum()}.items():assert m[key]==value
        np.testing.assert_allclose(m.sensitivity,m.TP/(m.TP+m.FN))
        np.testing.assert_allclose(m.specificity,m.TN/(m.TN+m.FP))
        np.testing.assert_allclose(m.PPV,m.TP/(m.TP+m.FP))
        np.testing.assert_allclose(m.NPV,m.TN/(m.TN+m.FN))
        for metric in ['AUROC','AUPRC','sensitivity','specificity','PPV','NPV','Brier']:
            assert 0<=m[metric+'_CI_low']<=m[metric+'_CI_high']<=1
            assert m[metric+'_bootstrap_defined']==5000
        genes=pd.read_csv(OUT/f'{target}_gene_coverage.csv');cov=coverage[coverage.target==target].iloc[0]
        weight=genes.coefficient.abs(); observed=genes.present_in_TCGA
        np.testing.assert_allclose(cov.absolute_coefficient_weight_coverage,weight[observed].sum()/weight.sum())
        assert cov.selected_absent==int((~observed).sum())
        checks.append(f'{target}: patient labels, prediction threshold, metrics and coefficient coverage reconcile; 5000 valid bootstrap replicates.')
    checks.append('Prepared inputs, protocol and saved predictions retain their original SHA256 hashes.')
    (OUT/'audit.txt').write_text('\n'.join(checks)+'\n');print('\n'.join(checks))

if __name__=='__main__':main()
