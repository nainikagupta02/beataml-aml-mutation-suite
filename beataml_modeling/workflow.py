"""Patient-grouped BEAT-AML feasibility study; never edits source data.

Run development first, inspect diagnostics, then finalize exactly once.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/beataml-mpl')
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import argparse
import hashlib
import json
import platform
import random
import sys
import time
import warnings
from pathlib import Path
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.calibration import calibration_curve
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, average_precision_score, brier_score_loss,
                             confusion_matrix, roc_curve, precision_recall_curve)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import xlsx_read

SEED = 20260916
TARGETS = ('FLT3-ITD', 'NPM1')
THRESHOLD = 0.5
AML = 'ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS'
PANELS = {
    'FLT3-ITD': ['PRDM16', 'IGFBP2', 'PTP4A3', 'CCNA1', 'IL2RA', 'SOCS2', 'DUSP6'],
    'NPM1': ['HOXA9', 'HOXA10', 'HOXB2', 'HOXB3', 'HOXB6', 'HOXD4', 'MEIS1'],
}
SOURCES = {
    'FLT3-ITD': ['https://pmc.ncbi.nlm.nih.gov/articles/PMC4726916/',
                 'https://www.nature.com/articles/s41598-019-45579-0',
                 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3464674/'],
    'NPM1': ['https://pubmed.ncbi.nlm.nih.gov/17597811/',
             'https://ashpublications.org/blood/article/106/12/3747/109726/Mutations-in-nucleophosmin-NPM1-in-acute-myeloid'],
}
METRICS = ['AUROC', 'AUPRC', 'sensitivity', 'specificity', 'PPV', 'NPV', 'Brier']

def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

class LogCPM(BaseEstimator, TransformerMixin):
    """Sample-local library normalization: no cross-sample factor is estimated.

    Denominator includes ALL nonzero input genes, before any learned filtering.
    A test specimen's own total is needed, but other test specimens are not.
    """
    def fit(self, X, y=None):
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X):
        if X.shape[1] != self.n_features_in_: raise ValueError('Gene order/width changed')
        a = np.asarray(X, dtype=np.float64)
        totals = a.sum(axis=1)
        if not np.isfinite(a).all() or (a < 0).any() or (totals <= 0).any():
            raise ValueError('Invalid counts or empty library')
        return np.log1p(a / totals[:, None] * 1e6).astype(np.float32)

class FoldGeneFilter(BaseEstimator, TransformerMixin):
    """Train-only low-expression/variance selection, optionally a fixed panel.

    CPM >= 1 in >=10% of training specimens; log1p CPM variance >0.01.
    The broad models retain at most 5,000 highest-variance eligible genes.
    Panel indices come only from annotation and literature, never outcomes.
    """
    def __init__(self, panel_indices=None, max_genes=5000, min_fraction=0.1, min_variance=0.01):
        self.panel_indices = panel_indices
        self.max_genes = max_genes
        self.min_fraction = min_fraction
        self.min_variance = min_variance

    def fit(self, X, y=None):
        self.variance_ = np.var(X, axis=0, dtype=np.float64)
        self.prevalence_ = np.mean(X >= np.log1p(1.0), axis=0)
        keep = (self.prevalence_ >= self.min_fraction) & (self.variance_ > self.min_variance)
        if self.panel_indices is not None:
            panel = np.zeros(X.shape[1], dtype=bool)
            panel[list(self.panel_indices)] = True
            keep &= panel
        indices = np.flatnonzero(keep)
        self.eligible_count_ = len(indices)
        if self.max_genes is not None and len(indices) > self.max_genes:
            ranked = np.argsort(-self.variance_[indices], kind='stable')[:self.max_genes]
            indices = np.sort(indices[ranked])
        if len(indices) == 0: raise ValueError('No genes passed training filter')
        self.indices_ = indices
        return self

    def transform(self, X):
        return X[:, self.indices_]

def metrics(y, p):
    tn, fp, fn, tp = confusion_matrix(y, p >= THRESHOLD, labels=[0, 1]).ravel()
    ratio = lambda a, b: float(a / b) if b else np.nan
    return dict(AUROC=roc_auc_score(y,p), AUPRC=average_precision_score(y,p),
                sensitivity=ratio(tp,tp+fn), specificity=ratio(tn,tn+fp),
                PPV=ratio(tp,tp+fp), NPV=ratio(tn,tn+fn), Brier=brier_score_loss(y,p),
                TN=int(tn), FP=int(fp), FN=int(fn), TP=int(tp))

def definitions(panel_indices):
    return {
        'majority': ({'model': DummyClassifier(strategy='prior', random_state=SEED)}, [{}]),
        'marker_logistic': ({'model': LogisticRegression(solver='lbfgs', max_iter=3000, random_state=SEED),
                             'panel_indices': panel_indices}, [{'model__C': c} for c in [0.1, 1.0, 10.0]]),
        'elastic_net': ({'model': LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5,
                             max_iter=5000, tol=1e-3, random_state=SEED)},
                        [{'model__C': c} for c in [0.01, 0.05, 0.2]]),
        'random_forest': ({'model': RandomForestClassifier(n_estimators=300, max_features='sqrt',
                              n_jobs=2, random_state=SEED)},
                          [{'model__min_samples_leaf': k} for k in [2, 5]]),
    }

def make_pipeline(definition):
    if isinstance(definition['model'], DummyClassifier):
        return Pipeline([('model', clone(definition['model']))])
    return Pipeline([('log_cpm', LogCPM()),
                     ('filter', FoldGeneFilter(panel_indices=definition.get('panel_indices'))),
                     ('scale', StandardScaler()), ('model', clone(definition['model']))])

def load_cohort(root, out):
    files = ['beataml_wv1to4_clinical.xlsx', 'beataml_waves1to4_sample_mapping.xlsx',
             'beataml_waves1to4_counts_dbgap.txt']
    hashes = {}
    for name in files:
        p = root/name
        with p.open('rb') as f: beginning = f.read(200)
        if p.stat().st_size < 1000 or b'version https://git-lfs.github.com/spec/v1' in beginning:
            raise ValueError(f'Missing or pointer data: {p}')
        hashes[name] = {'bytes':p.stat().st_size, 'sha256':digest(p)}
    xlsx_read.ROOT = root
    c = xlsx_read.workbook(files[0])['summary']
    m = xlsx_read.workbook(files[1])['sample_map']
    d = pd.read_csv(root/files[2], sep='\t', keep_default_na=False)
    assert d.columns[:4].tolist() == ['stable_id','display_label','description','biotype']
    assert d.stable_id.is_unique and d.stable_id.str.startswith('ENSG').all()
    annotation = d.iloc[:, :4].copy()
    counts = d.iloc[:, 4:].to_numpy(dtype=np.float64).T
    assert np.isfinite(counts).all() and (counts >= 0).all() and (counts == np.floor(counts)).all()
    assert (counts.sum(axis=1) > 0).all()
    sample_ids = d.columns[4:].tolist()
    assert len(set(sample_ids)) == len(sample_ids)
    sample = pd.DataFrame({'dbgap_rnaseq_sample':sample_ids, 'count_row':np.arange(len(sample_ids))})
    mapped = sample.merge(m.dropna(subset=['dbgap_rnaseq_sample']), how='left',
                          on='dbgap_rnaseq_sample', validate='one_to_one', indicator=True)
    assert (mapped._merge == 'both').all()
    assert mapped.rna_control.notna().all()
    controls = mapped.rna_control.ne('No')
    joined = mapped.loc[~controls, ['dbgap_rnaseq_sample','count_row','dbgap_subject_id']].rename(
        columns={'dbgap_subject_id':'mapping_subject_id'}).merge(
        c.dropna(subset=['dbgap_rnaseq_sample']), on='dbgap_rnaseq_sample', how='left',
        validate='one_to_one', indicator=True)
    assert (joined._merge == 'both').all()
    assert (joined.mapping_subject_id == joined.dbgap_subject_id).all()
    aml = joined[joined.dxAtSpecimenAcquisition == AML]
    cohort = aml[aml.diseaseStageAtSpecimenCollection.str.strip().str.casefold() == 'initial diagnosis'].copy()
    # Only label-independent global all-zero removal; ALL other selection stays inside Pipeline.
    nonzero = np.any(counts != 0, axis=0)
    annotation = annotation.loc[nonzero].reset_index(drop=True)
    X = counts[cohort.count_row.to_numpy(), :][:, nonzero].astype(np.float32)
    cohort = cohort.reset_index(drop=True)
    # No clinical column, including mutation-related text, is ever used to construct X.
    assert X.shape == (len(cohort), len(annotation))
    manifest_cols = ['dbgap_rnaseq_sample','dbgap_subject_id','cohort','specimenType',
                     'dxAtSpecimenAcquisition','diseaseStageAtSpecimenCollection',*TARGETS]
    cohort[manifest_cols].to_csv(out/'cohort.csv',index=False)
    annotation.to_csv(out/'gene_annotation.csv', index=False)
    flow = {'input_samples':len(sample_ids), 'controls_removed':int(controls.sum()),
            'matched_noncontrol':len(joined), 'broad_AML':len(aml), 'initial_AML':len(cohort),
            'patients':int(cohort.dbgap_subject_id.nunique()),
            'patients_with_multiple_specimens':int((cohort.groupby('dbgap_subject_id').size()>1).sum()),
            'input_genes':len(nonzero), 'globally_all_zero_removed':int((~nonzero).sum()),
            'input_genes_after_zero_removal':int(nonzero.sum())}
    flow['targets'] = {}
    for target in TARGETS:
        labels = cohort[target]
        assert set(labels.dropna()) <= {'positive','negative'}
        flow['targets'][target] = {'positive':int(labels.eq('positive').sum()),
                                   'negative':int(labels.eq('negative').sum()),
                                   'missing':int(labels.isna().sum()),
                                   'patients':int(cohort.loc[labels.notna(),'dbgap_subject_id'].nunique())}
    save_json(out/'cohort_flow.json',flow)
    save_json(out/'data_hashes.json', hashes)
    print('COHORT\n'+json.dumps(flow,indent=2),flush=True)
    return X, cohort, annotation

def fit_checked(pipe, X, y):
    # A nonconverged run fails instead of silently supplying a result.
    with warnings.catch_warnings():
        warnings.simplefilter('error', ConvergenceWarning)
        pipe.fit(X,y)
    return pipe

def assert_group_split(train, test, groups, y):
    assert set(groups[train]).isdisjoint(groups[test])
    assert set(y[train]) == {0,1} and set(y[test]) == {0,1}

def reliability(ax, y, p, label):
    actual, pred = calibration_curve(y,p,n_bins=5,strategy='quantile')
    ax.plot(pred,actual,marker='o',label=label)

def export_features(pipe, annotation, target, family, out):
    filtered = annotation.iloc[pipe.named_steps['filter'].indices_].copy()
    estimator = pipe.named_steps['model']
    if hasattr(estimator,'coef_'):
        filtered['weight'] = estimator.coef_[0]
        filtered['importance_type'] = 'standardized_log_odds_coefficient'
    else:
        filtered['weight'] = estimator.feature_importances_
        filtered['importance_type'] = 'training_impurity_importance'
    filtered['absolute_weight'] = filtered.weight.abs()
    filtered['literature_panel_member'] = filtered.display_label.isin(PANELS[target])
    filtered = filtered.sort_values('absolute_weight',ascending=False,kind='stable')
    filtered.to_csv(out/f'{target}_{family}_features.csv',index=False)
    print(f'{target} {family} top genes:\n'+filtered[['stable_id','display_label','weight']].head(10).to_string(index=False),flush=True)
    fig,ax = plt.subplots(figsize=(8,6))
    top = filtered.head(20).iloc[::-1]
    ax.barh(top.display_label.astype(str)+' | '+top.stable_id,top.weight,
            color=np.where(top.weight>=0,'#287f8e','#bd5b47'))
    ax.set_title(f'{target}: {family} (development data only)')
    ax.set_xlabel(filtered.importance_type.iloc[0].replace('_',' '))
    fig.tight_layout(); fig.savefig(out/f'{target}_{family}_features.png',dpi=160); plt.close(fig)

def develop(root, out):
    if (out/'selection.json').exists() or (out/'heldout_metrics.csv').exists():
        raise FileExistsError('Completed development exists; use a NEW output directory for a fresh study.')
    random.seed(SEED); np.random.seed(SEED)
    X, cohort, annotation = load_cohort(root,out)
    config = {'seed':SEED, 'holdout':'first fold of shuffled 5-fold StratifiedGroupKFold',
              'cv':'5 folds x 3 repetitions of StratifiedGroupKFold',
              'cv_seeds':[SEED+100+r for r in range(3)],
              'selection_metric':'mean validation average precision (AUPRC), ties AUROC then family order',
              'threshold':THRESHOLD, 'normalization':'natural log1p of CPM; own sample library total',
              'filter':'CPM>=1 in >=10% training rows; variance>0.01; cap 5000 by training variance',
              'panels':PANELS, 'panel_sources':SOURCES, 'nested_cv':False,
              'MLP':'omitted: optional at modest sample size; no deep-model claim',
              'versions':{k:v for k,v in [('python',platform.python_version()),('numpy',np.__version__),
                         ('pandas',pd.__version__),('scipy',scipy.__version__),('scikit-learn',sklearn.__version__),
                         ('matplotlib',matplotlib.__version__),('joblib',joblib.__version__)]}}
    save_json(out/'protocol.json',config)
    # Preserve exact numeric inputs and ordering so finalization cannot reconstruct a different cohort.
    np.save(out/'counts_input.npy',X,allow_pickle=False)
    all_fold_metrics=[]; summaries=[]; oof_records=[]; selection={}; split_records=[]; cv_records=[]
    for target in TARGETS:
        valid = cohort[target].notna().to_numpy()
        idx = np.flatnonzero(valid)
        y = cohort.loc[valid,target].eq('positive').to_numpy(dtype=int)
        groups = cohort.loc[valid,'dbgap_subject_id'].to_numpy()
        outer = StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=SEED)
        dev, test = next(outer.split(X[idx],y,groups))
        assert_group_split(dev,test,groups,y)
        save_json(out/f'{target}_indices.json',{'development_rows':idx[dev].tolist(),'test_rows':idx[test].tolist()})
        for partition, positions in [('development',dev),('test',test)]:
            print(f'{target} {partition}: n={len(positions)}, patients={len(set(groups[positions]))}, '+
                  f'positive={int(y[positions].sum())}, negative={int((1-y[positions]).sum())}',flush=True)
            for k in positions:
                split_records.append({'target':target,'partition':partition,'sample':cohort.iloc[idx[k]].dbgap_rnaseq_sample,
                                      'patient':int(groups[k]),'y':int(y[k])})
        xd, yd, gd = X[idx[dev]], y[dev], groups[dev]
        folds=[]
        for repeat in range(3):
            cv = StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=SEED+100+repeat)
            for fold,(tr,va) in enumerate(cv.split(xd,yd,gd)):
                assert_group_split(tr,va,gd,yd)
                folds.append((repeat,fold,tr,va))
                for k in va:
                    cv_records.append({'target':target,'repeat':repeat,'fold':fold,
                        'sample':cohort.iloc[idx[dev[k]]].dbgap_rnaseq_sample,'patient':int(gd[k])})
        panel_indices=tuple(np.flatnonzero(annotation.display_label.isin(PANELS[target])).tolist())
        found=set(annotation.iloc[list(panel_indices)].display_label)
        if found != set(PANELS[target]): raise ValueError(f'Missing panel genes: {set(PANELS[target])-found}')
        families=definitions(panel_indices); family_scores=[]
        for family,(definition,grid) in families.items():
            candidate_scores=[]; candidate_predictions={}
            for candidate,params in enumerate(grid):
                start=time.time(); rows=[]; pred=np.full((3,len(dev)),np.nan)
                for repeat,fold,tr,va in folds:
                    pipe=make_pipeline(definition).set_params(**params)
                    fit_checked(pipe,xd[tr],yd[tr])
                    p=pipe.predict_proba(xd[va])[:,1]
                    pred[repeat,va]=p
                    row={'target':target,'model':family,'candidate':candidate,'params':json.dumps(params,sort_keys=True),
                         'repeat':repeat,'fold':fold,'n_train':len(tr),'n_valid':len(va),
                         'train_patients':len(set(gd[tr])),'valid_patients':len(set(gd[va])),
                         'group_overlap':0,'threshold':THRESHOLD, **metrics(yd[va],p)}
                    if 'filter' in pipe.named_steps:
                        row['genes']=len(pipe.named_steps['filter'].indices_)
                        row['eligible_genes']=pipe.named_steps['filter'].eligible_count_
                    row['iterations']=int(np.max(getattr(pipe.named_steps['model'],'n_iter_', [0])))
                    rows.append(row); all_fold_metrics.append(row)
                    print(f'{target} {family} candidate {candidate+1}/{len(grid)} repeat {repeat+1} fold {fold+1} done',flush=True)
                assert np.isfinite(pred).all()
                frame=pd.DataFrame(rows)
                score={'candidate':candidate,'params':params}
                for metric in METRICS:
                    score[metric+'_mean']=float(frame[metric].mean())
                    score[metric+'_sd']=float(frame[metric].std(ddof=1))
                    score[metric+'_defined_folds']=int(frame[metric].notna().sum())
                candidate_scores.append(score); candidate_predictions[candidate]=pred
                pd.DataFrame(all_fold_metrics).to_csv(out/'cv_all_candidates.csv',index=False)
                print(f'{target} {family} {params}: AP={score["AUPRC_mean"]:.3f}, '+
                      f'AUC={score["AUROC_mean"]:.3f}, elapsed={time.time()-start:.1f}s',flush=True)
            best=sorted(candidate_scores,key=lambda z:(-z['AUPRC_mean'],-z['AUROC_mean'],z['candidate']))[0]
            summaries.append({'target':target,'model':family,'threshold':THRESHOLD,
                              **{k:json.dumps(v,sort_keys=True) if k=='params' else v for k,v in best.items()}})
            family_scores.append((family,best))
            pred=candidate_predictions[best['candidate']]
            for repeat in range(3):
                for k in range(len(dev)):
                    oof_records.append({'target':target,'model':family,'repeat':repeat,
                        'sample':cohort.iloc[idx[dev[k]]].dbgap_rnaseq_sample,'patient':int(gd[k]),
                        'y':int(yd[k]),'probability':float(pred[repeat,k])})
            # Fits use development data only, including models retained for feature interpretation.
            final_pipe=fit_checked(make_pipeline(definition).set_params(**best['params']),xd,yd)
            joblib.dump(final_pipe,out/f'{target}_{family}.joblib')
            if family!='majority': export_features(final_pipe,annotation,target,family,out)
            pd.DataFrame(summaries).to_csv(out/'cv_summary.csv',index=False)
        winner,best=sorted(family_scores,key=lambda z:(-z[1]['AUPRC_mean'],-z[1]['AUROC_mean']))[0]
        selection[target]={'model':winner,'params':best['params'],'threshold':THRESHOLD,
                           'test_n':len(test),'development_n':len(dev),'CV_AUPRC':best['AUPRC_mean']}
        print(f'LOCKED selection for {target}: {winner} {best["params"]}',flush=True)
    pd.DataFrame(split_records).to_csv(out/'split_manifest.csv',index=False)
    pd.DataFrame(cv_records).to_csv(out/'cv_manifest.csv',index=False)
    oof=pd.DataFrame(oof_records); oof.to_csv(out/'oof_predictions.csv',index=False)
    save_json(out/'selection.json',selection)
    fig, axes=plt.subplots(1,2,figsize=(12,5))
    for target,ax in zip(TARGETS,axes):
        for family in definitions(()):
            # Each dev sample appears once in this curve; no pseudo-independent pooling of repeats.
            v=oof[(oof.target==target)&(oof.model==family)&(oof.repeat==0)]
            reliability(ax,v.y,v.probability,family)
        ax.plot([0,1],[0,1],'k--',alpha=.5); ax.set(xlim=(0,1),ylim=(0,1),
            xlabel='Mean predicted probability',ylabel='Observed positive fraction',title=target+' development OOF, repeat 1')
        ax.legend(fontsize=8)
    fig.suptitle('Tuning-selected models; descriptive CV calibration, 5 quantile bins')
    fig.tight_layout(); fig.savefig(out/'calibration_cv.png',dpi=170); plt.close(fig)
    print('Development finished. No held-out predictions have been made.',flush=True)

def finalize(root,out):
    if (out/'heldout_metrics.csv').exists() or (out/'test_evaluation_started.json').exists():
        raise FileExistsError('Test evaluation already started/completed; never rerun for model selection.')
    selection=json.loads((out/'selection.json').read_text())
    for name,expected in json.loads((out/'data_hashes.json').read_text()).items():
        assert digest(root/name)==expected['sha256'],f'Source changed: {name}'
    X=np.load(out/'counts_input.npy',allow_pickle=False)
    cohort=pd.read_csv(out/'cohort.csv')
    save_json(out/'test_evaluation_started.json',{'selection_sha256':digest(out/'selection.json')})
    results=[]; predictions=[]
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    for row,target in enumerate(TARGETS):
        chosen=selection[target]
        test=np.array(json.loads((out/f'{target}_indices.json').read_text())['test_rows'])
        model=joblib.load(out/f'{target}_{chosen["model"]}.joblib')
        p=model.predict_proba(X[test])[:,1]  # One prediction call on the locked test partition.
        y=cohort.iloc[test][target].eq('positive').to_numpy(dtype=int)
        result={'target':target,'model':chosen['model'],'n':len(test),'positive':int(y.sum()),
                'negative':int((1-y).sum()),'threshold':THRESHOLD,**metrics(y,p)}
        results.append(result)
        for k,prob,label in zip(test,p,y):
            predictions.append({'target':target,'model':chosen['model'],
                                'sample':cohort.iloc[k].dbgap_rnaseq_sample,
                                'patient':int(cohort.iloc[k].dbgap_subject_id),'y':int(label),'probability':float(prob)})
        reliability(axes[row,0],y,p,chosen['model'])
        axes[row,0].plot([0,1],[0,1],'k--',alpha=.5)
        axes[row,0].set(xlim=(0,1),ylim=(0,1),xlabel='Predicted probability',ylabel='Observed fraction',
                        title=f'{target}: held-out calibration (n={len(y)})')
        axes[row,0].legend(fontsize=8)
        fpr,tpr,_=roc_curve(y,p); axes[row,1].plot(fpr,tpr)
        axes[row,1].plot([0,1],[0,1],'k--',alpha=.5)
        axes[row,1].set(xlabel='False positive rate',ylabel='Sensitivity',title=f'AUROC {result["AUROC"]:.3f}')
        precision,recall,_=precision_recall_curve(y,p); axes[row,2].plot(recall,precision)
        axes[row,2].axhline(y.mean(),linestyle='--',color='gray',label='Test prevalence')
        axes[row,2].set(xlabel='Recall',ylabel='Precision',title=f'Average precision {result["AUPRC"]:.3f}')
        axes[row,2].legend()
    fig.suptitle('Locked selections: single held-out evaluation; threshold fixed at 0.50')
    fig.tight_layout(); fig.savefig(out/'heldout_diagnostics.png',dpi=180); plt.close(fig)
    pd.DataFrame(results).to_csv(out/'heldout_metrics.csv',index=False)
    pd.DataFrame(predictions).to_csv(out/'heldout_predictions.csv',index=False)
    print(pd.DataFrame(results).to_string(index=False),flush=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results')
    parser.add_argument('--phase',choices=['develop','finalize'],required=True)
    args=parser.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    if args.phase=='develop': develop(args.data_dir.resolve(),args.out.resolve())
    else: finalize(args.data_dir.resolve(),args.out.resolve())

if __name__=='__main__': main()
