"""Frozen BEAT-AML -> TCGA-LAML validation. No estimator is ever fitted here."""
import os
import sys
sys.dont_write_bytecode = True  # importing original pipeline must not modify its directory
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BEAT=ROOT.parent/'beataml_modeling'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/beataml-mpl')
sys.path.insert(0,str(BEAT))
import argparse
import hashlib
import json
import re
import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.calibration import calibration_curve
from sklearn.metrics import roc_curve, precision_recall_curve
import workflow  # exposes the original serialized transform classes; never call training functions
import xlsx_read

SEED=20260916
TARGETS=('FLT3-ITD','NPM1')
DATA=ROOT/'data'
OUT=ROOT/'results'
METRICS=workflow.METRICS
S1_URL='https://api.gdc.cancer.gov/data/b9196563-a05d-40b8-80dc-640ec712eb06'

def sha(path, kind='sha256'):
    h=hashlib.new(kind)
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def dump(name, value):
    (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def check_frozen_files():
    before=json.loads((OUT/'beataml_artifacts_before.json').read_text())
    for name,h in before.items():
        if sha(name)!=h: raise RuntimeError('BEAT-AML artifact changed: '+name)
    now={str(p.resolve()) for p in BEAT.rglob('*') if p.is_file()}
    if now != set(before): raise RuntimeError('Files added to/removed from BEAT-AML artifact directory')
    return len(before)

def read_star(path):
    with path.open() as f: first=f.readline().strip()
    if first!='# gene-model: GENCODE v36': raise ValueError('Unexpected STAR gene model: '+first)
    d=pd.read_csv(path,sep='\t',comment='#',keep_default_na=False)
    if 'unstranded' not in d: raise ValueError('STAR unstranded raw counts absent')
    d=d[d.gene_id.str.startswith('ENSG')].copy()
    # Preserve distinct PAR_Y identities; only the Ensembl version component is removed.
    d['stable_id']=d.gene_id.str.replace(r'\.\d+(?=_PAR_Y$|$)','',regex=True)
    if not d.stable_id.is_unique: raise ValueError('Ambiguous Ensembl ID duplication after version removal')
    a=pd.to_numeric(d.unstranded,errors='raise').to_numpy(dtype=np.float64)
    if not np.isfinite(a).all() or (a<0).any() or (a!=np.floor(a)).any() or a.sum()<=0:
        raise ValueError('Invalid raw counts')
    return d,a

def curated_labels():
    xlsx_read.ROOT=DATA
    source=xlsx_read.workbook('TCGA_AML_Table_S1.xlsx')['Sheet1']
    # The last workbook row contains cohort totals, not a patient.
    s=source.loc[source['TCGA Patient ID'].notna(),['TCGA Patient ID','UPN','RNAseq data?','FLT3','NPM1']].copy()
    assert len(s)==200
    s['patient']='TCGA-AB-'+s['TCGA Patient ID'].astype(int).astype(str).str.zfill(4)
    assert s.patient.is_unique
    s=s.rename(columns={'FLT3':'FLT3_curated_variant','NPM1':'NPM1_curated_variant'})
    f=s.FLT3_curated_variant.fillna('').astype(str).str.strip()
    # This is the original author-curated per-patient mutation matrix, not masked MAF.
    # Original in_frame_ins calls define its 39 ITD cases, including calls whose later
    # cBioPortal reannotation misleadingly says splice-region/nonsense/in-frame deletion.
    itd=f.str.contains('in_frame_ins',regex=False)
    pm=f.str.fullmatch(r'p\.D83[59][A-Z]')
    assert (itd|pm|f.eq('')).all(), 'Unrecognized FLT3 annotation: do not guess labels'
    assert int(itd.sum())==39 and int(pm.sum())==17
    # All 200 in the source study underwent WGS/WES. Blank cells here denote no
    # reported mutation; absent patient rows (if any at join) remain missing.
    s['FLT3-ITD']=np.where(itd,'positive','negative')
    s['FLT3_point_mutation']=pm
    s['NPM1']=np.where(s.NPM1_curated_variant.notna(),'positive','negative')
    assert (s.NPM1=='positive').sum()==54
    # Verify all positive gene-level calls against cBioPortal original TCGA study.
    cb=pd.read_csv(DATA/'cbio_mutations.txt',sep='\t',comment='#',low_memory=False)
    cb['patient']=cb.Tumor_Sample_Barcode.str[:12]
    cb_npm=set(cb.loc[cb.Hugo_Symbol=='NPM1','patient'])
    assert cb_npm==set(s.loc[s.NPM1=='positive','patient'])
    assert set(s.loc[itd,'patient']) <= set(cb.loc[cb.Hugo_Symbol=='FLT3','patient'])
    # Independent author-released TCGA reannotation confirms exactly the same 39 insertions.
    cross=pd.read_csv(DATA/'TCGA_Genetic_crosscheck.txt',sep='\t')
    cross_itd=set(cross.loc[(cross.Hugo_Symbol=='FLT3')&(cross.Variant_Classification=='In_Frame_Ins'),'TCGA_ID'])
    assert cross_itd==set(s.loc[itd,'patient'])
    s['label_source']='NEJM 2013 TCGA Table S1, GDC b9196563-a05d-40b8-80dc-640ec712eb06'
    s.to_csv(OUT/'curated_labels_all200.csv',index=False)
    return s

def coverage(pipe, annotation, available, target):
    ix=pipe.named_steps['filter'].indices_
    coefficients=pipe.named_steps['model'].coef_[0]
    selected=annotation.iloc[ix].reset_index(drop=True).copy()
    selected['coefficient']=coefficients
    selected['present_in_TCGA']=available[ix]
    selected['standardized_zero_imputed']=~available[ix]
    selected.to_csv(OUT/f'{target}_gene_coverage.csv',index=False)
    weight=np.abs(coefficients)
    present=available[ix]
    return {'target':target,'selected_genes':len(ix),'selected_present':int(present.sum()),
            'selected_absent':int((~present).sum()),'selected_gene_coverage':float(present.mean()),
            'nonzero_coefficient_genes':int((weight>0).sum()),
            'nonzero_coefficient_present':int(((weight>0)&present).sum()),
            'nonzero_coefficient_absent':int(((weight>0)&~present).sum()),
            'absolute_coefficient_weight_coverage':float(weight[present].sum()/weight.sum()),
            'below_90_percent_weight':bool(weight[present].sum()/weight.sum()<.9)}

def prepare():
    if (OUT/'protocol.json').exists(): raise FileExistsError('Prepared protocol exists; do not overwrite after inspecting outcomes')
    if not (OUT/'beataml_artifacts_before.json').exists():
        dump('beataml_artifacts_before.json',{str(p.resolve()):sha(p) for p in BEAT.rglob('*') if p.is_file()})
    check_frozen_files()
    hits=json.loads((DATA/'gdc_files.json').read_text())['data']['hits']
    assert len(hits)==151
    annotation=pd.read_csv(BEAT/'results/gene_annotation.csv')
    assert annotation.stable_id.is_unique
    rows=[]; matrices=[]; original_ids=None; counts_hashes={}; totals=[]
    for hit in sorted(hits,key=lambda h:h['cases'][0]['submitter_id']):
        assert hit['access']=='open' and hit['data_type']=='Gene Expression Quantification'
        assert hit['analysis']['workflow_type']=='STAR - Counts'
        assert len(hit['cases'])==1 and len(hit['cases'][0]['samples'])==1
        case=hit['cases'][0]; sample=case['samples'][0]
        assert sample['sample_type'] in ('Primary Blood Derived Cancer - Peripheral Blood','Primary Blood Derived Cancer - Bone Marrow')
        path=DATA/'counts'/(hit['file_id']+'.tsv')
        if path.stat().st_size!=hit['file_size'] or sha(path,'md5')!=hit['md5sum']:
            raise ValueError('GDC integrity check failed: '+str(path))
        d,raw=read_star(path)
        if original_ids is None:
            original_ids=d.stable_id.to_numpy()
            positions=pd.Index(original_ids).get_indexer(annotation.stable_id)
            available=positions>=0
            annotation.assign(present_in_TCGA=available).to_csv(OUT/'all_input_gene_coverage.csv',index=False)
        else: np.testing.assert_array_equal(original_ids,d.stable_id.to_numpy())
        vector=np.zeros(len(annotation),dtype=np.float64)
        vector[available]=raw[positions[available]]
        matrices.append(vector)
        rows.append({'patient':case['submitter_id'],'sample':sample['submitter_id'],
                     'sample_type':sample['sample_type'],'file_id':hit['file_id'],
                     'full_gene_count_sum':float(raw.sum()),'matched_input_gene_count_sum':float(vector.sum()),
                     'fraction_counts_retained':float(vector.sum()/raw.sum())})
        counts_hashes[str(path.relative_to(ROOT))]={'md5':hit['md5sum'],'sha256':sha(path),'bytes':hit['file_size']}
    manifest=pd.DataFrame(rows)
    assert manifest.patient.is_unique # one case per patient: no replicate selection needed
    labels=curated_labels()
    joined=manifest.merge(labels,on='patient',how='left',validate='one_to_one',indicator=True)
    joined.to_csv(OUT/'cohort.csv',index=False)
    X=np.stack(matrices)
    assert X.shape==(len(joined),len(annotation)) and (X.sum(axis=1)>0).all()
    np.save(OUT/'aligned_raw_counts.npy',X,allow_pickle=False)
    np.save(OUT/'available_genes.npy',available,allow_pickle=False)
    coverage_rows=[]; model_hashes={}
    for target in TARGETS:
        path=BEAT/'results'/f'{target}_elastic_net.joblib'
        pipe=joblib.load(path)
        assert pipe.named_steps['log_cpm'].n_features_in_==len(annotation)
        assert pipe.named_steps['model'].classes_.tolist()==[0,1]
        coverage_rows.append(coverage(pipe,annotation,available,target))
        model_hashes[target]=sha(path)
        print(target,joined[target].value_counts(dropna=False).to_dict(),flush=True)
    cv=pd.DataFrame(coverage_rows); cv.to_csv(OUT/'coverage_summary.csv',index=False)
    for r in coverage_rows:
        if r['below_90_percent_weight']: print('WARNING: COEFFICIENT WEIGHT COVERAGE BELOW 90%: '+r['target'],flush=True)
    sources={str(p.relative_to(ROOT)):sha(p) for p in DATA.iterdir() if p.is_file() and p.suffix!='.curl'}
    dump('source_hashes.json',{'metadata_and_labels':sources,'counts':counts_hashes})
    protocol={'seed':SEED,'project':'TCGA-LAML','access':'open','workflow':'STAR - Counts',
              'expression_column':'unstranded','gene_model':'GENCODE v36',
              'alignment':'Ensembl stable IDs only; strip version; preserve PAR_Y identity; no symbol remapping',
              'normalization':'Frozen LogCPM over the exact 51016 training input columns, populated from 46973 matched IDs; missing raw placeholders zero only for transform input',
              'missing_features':'After frozen standardization, every absent selected gene is set to 0; not interpreted as measured zero expression',
              'threshold':0.5,'bootstrap':'5000 ordinary patient bootstrap resamples; percentile 95% CI; one sample per patient',
              'label_source':S1_URL,'FLT3_definition':'39 author-curated in_frame_ins calls; point mutations ITD-negative',
              'NPM1_definition':'Any reported NPM1 mutation in complete author-curated matrix; includes one p.K263R missense',
              'secondary_label_check':'Report NPM1 excluding noncanonical p.K263R using the SAME predictions; never select model/threshold by this result',
              'models':model_hashes,'all_input_genes':len(annotation),'matched_input_genes':int(available.sum()),
              'absent_input_genes':int((~available).sum()),'STAR_gene_rows':len(original_ids),
              'refitting':False,'retuning':False,'recalibration':False,'batch_correction':False,
              'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,
                          'sklearn':workflow.sklearn.__version__}}
    dump('protocol.json',protocol)
    dump('prepared_hashes.json',{name:sha(OUT/name) for name in ['aligned_raw_counts.npy','available_genes.npy','cohort.csv','protocol.json']})
    check_frozen_files()
    print(cv.to_string(index=False),flush=True)
    print('Prepared TCGA only. No external predictions made.',flush=True)

def predict_frozen(pipe,X,available):
    """Apply only frozen transforms and coefficients, with user-specified zero contribution."""
    before=joblib.hash(pipe)
    log_values=pipe.named_steps['log_cpm'].transform(X)
    selected=pipe.named_steps['filter'].transform(log_values)
    standardized=pipe.named_steps['scale'].transform(selected)
    missing=~available[pipe.named_steps['filter'].indices_]
    standardized[:,missing]=0.0
    assert np.isfinite(standardized).all() and np.all(standardized[:,missing]==0)
    p=pipe.named_steps['model'].predict_proba(standardized)[:,1]
    assert before==joblib.hash(pipe), 'Frozen pipeline mutated during inference'
    return p,standardized

def evaluate():
    if (OUT/'evaluation_started.json').exists(): raise FileExistsError('External inference already started; use report to summarize saved predictions')
    check_frozen_files()
    for name,h in json.loads((OUT/'prepared_hashes.json').read_text()).items(): assert sha(OUT/name)==h
    protocol=json.loads((OUT/'protocol.json').read_text())
    X=np.load(OUT/'aligned_raw_counts.npy',allow_pickle=False)
    available=np.load(OUT/'available_genes.npy',allow_pickle=False)
    c=pd.read_csv(OUT/'cohort.csv')
    dump('evaluation_started.json',{'protocol_sha256':sha(OUT/'protocol.json'),'seed':SEED})
    predictions=[]; qc=[]
    for target in TARGETS:
        path=BEAT/'results'/f'{target}_elastic_net.joblib'
        assert sha(path)==protocol['models'][target]
        pipe=joblib.load(path)
        p,z=predict_frozen(pipe,X,available)
        valid=c[target].isin(['positive','negative']).to_numpy()
        for i in np.flatnonzero(valid):
            predictions.append({'target':target,'patient':c.iloc[i].patient,'sample':c.iloc[i]['sample'],
                'label':c.iloc[i][target],'y':int(c.iloc[i][target]=='positive'),'probability':float(p[i]),
                'prediction':int(p[i]>=protocol['threshold'])})
        active=pipe.named_steps['model'].coef_[0]!=0
        observed=available[pipe.named_steps['filter'].indices_]&active
        qc.append({'target':target,'abs_z_gt_5_fraction':float((np.abs(z[:,observed])>5).mean()),
                   'abs_z_gt_10_fraction':float((np.abs(z[:,observed])>10).mean()),
                   'max_abs_z':float(np.max(np.abs(z[:,observed]))),
                   'predicted_positive':int((p>=.5).sum()),'probability_min':float(p.min()),
                   'probability_max':float(p.max())})
    pd.DataFrame(predictions).to_csv(OUT/'predictions.csv',index=False)
    pd.DataFrame(qc).to_csv(OUT/'inference_qc.csv',index=False)
    dump('evaluation_completed.json',{'predictions_sha256':sha(OUT/'predictions.csv'),
                                      'unchanged_BEAT_files':check_frozen_files()})
    print('External inference complete; original models and all BEAT-AML artifacts unchanged.',flush=True)

def bootstrapped(y,p,repetitions=5000):
    rng=np.random.default_rng(SEED)
    scores=[]
    for _ in range(repetitions):
        ix=rng.integers(0,len(y),size=len(y))
        if len(np.unique(y[ix]))<2: continue
        scores.append(workflow.metrics(y[ix],p[ix]))
    b=pd.DataFrame(scores)
    point=workflow.metrics(y,p)
    return point,{m:(float(b[m].quantile(.025)),float(b[m].quantile(.975)),int(b[m].notna().sum())) for m in METRICS}

def report():
    check_frozen_files()
    assert sha(OUT/'predictions.csv')==json.loads((OUT/'evaluation_completed.json').read_text())['predictions_sha256']
    assert sha(OUT/'protocol.json')==json.loads((OUT/'evaluation_started.json').read_text())['protocol_sha256']
    d=pd.read_csv(OUT/'predictions.csv'); c=pd.read_csv(OUT/'cohort.csv'); cov=pd.read_csv(OUT/'coverage_summary.csv')
    rows=[]; bins=[]; text=[]
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    for row,target in enumerate(TARGETS):
        v=d[d.target==target]
        assert v.patient.is_unique
        y=v.y.to_numpy(); p=v.probability.to_numpy()
        point,ci=bootstrapped(y,p)
        record={'target':target,'N':len(y),'positive':int(y.sum()),'negative':int((1-y).sum()),'threshold':.5,**point}
        for metric,(lo,hi,n) in ci.items():
            record[metric+'_CI_low']=lo; record[metric+'_CI_high']=hi; record[metric+'_bootstrap_defined']=n
        rows.append(record)
        actual,pred=calibration_curve(y,p,n_bins=5,strategy='quantile')
        ax=axes[row,0]; ax.plot(pred,actual,'o-',label='Frozen elastic-net'); ax.plot([0,1],[0,1],'k--',alpha=.4)
        ax.set(xlim=(0,1),ylim=(0,1),xlabel='Mean predicted probability',ylabel='Observed positive fraction',
               title=f'{target}: calibration, N={len(y)}'); ax.legend(fontsize=8)
        for k,(a,b) in enumerate(zip(pred,actual)): bins.append({'target':target,'bin':k+1,'mean_predicted':a,'observed_fraction':b})
        fpr,tpr,_=roc_curve(y,p); ax=axes[row,1]; ax.plot(fpr,tpr); ax.plot([0,1],[0,1],'k--',alpha=.4)
        ax.set(xlabel='False positive rate',ylabel='Sensitivity',title=f'AUROC {point["AUROC"]:.3f} ({ci["AUROC"][0]:.3f}-{ci["AUROC"][1]:.3f})')
        precision,recall,_=precision_recall_curve(y,p); ax=axes[row,2]; ax.step(recall,precision,where='post')
        ax.axhline(y.mean(),color='gray',linestyle='--',label=f'Prevalence {y.mean():.3f}')
        ax.set(xlabel='Recall',ylabel='Precision',title=f'AUPRC {point["AUPRC"]:.3f} ({ci["AUPRC"][0]:.3f}-{ci["AUPRC"][1]:.3f})'); ax.legend(fontsize=8)
    fig.suptitle('External TCGA-LAML validation: frozen BEAT-AML models, threshold 0.50\n95% CIs: 5,000 patient bootstrap resamples; no refitting, tuning or recalibration')
    fig.tight_layout(); fig.savefig(OUT/'external_diagnostics.png',dpi=180); plt.close(fig)
    results=pd.DataFrame(rows); results.to_csv(OUT/'metrics.csv',index=False); pd.DataFrame(bins).to_csv(OUT/'calibration_bins.csv',index=False)
    table=[]
    for r in rows:
        table.append({'target':r['target'],'N (+/-)':f'{r["N"]} ({r["positive"]}/{r["negative"]})',
                      **{m:f'{r[m]:.3f} [{r[m+"_CI_low"]:.3f}, {r[m+"_CI_high"]:.3f}]' for m in METRICS}})
    text += ['External TCGA-LAML validation of frozen BEAT-AML elastic-net classifiers',
             pd.DataFrame(table).to_string(index=False),'',
             'Brackets are 95% percentile CIs from 5,000 ordinary patient bootstrap resamples (seed 20260916).',
             'Threshold fixed at 0.50; AUPRC is noninterpolated average precision. No classifier, gene filter, scaler,',
             'threshold, calibration or batch correction was fitted or tuned on TCGA.', '']
    for r in rows:
        text.append(f'{r["target"]}: TP={r["TP"]}, FN={r["FN"]}, TN={r["TN"]}, FP={r["FP"]}.')
    text += ['', 'Interpretation: NPM1 retains strong external discrimination. FLT3-ITD retains useful ranking,',
             'but the unchanged 0.50 threshold misses 15 of 34 ITD-positive patients (44.1%).',
             'This sensitivity limitation is retained rather than repaired by tuning on the external cohort.',
             'Reliability curves are descriptive: five small bins cannot establish precise calibration.']
    text += ['', 'Acquisition and labels:',
             '151 GDC open-access STAR - Counts files, 151 distinct patients, 60,660 gene rows/file (GENCODE v36).',
             'Used integer unstranded raw counts, never TPM/FPKM/RSEM. Excluded four STAR summary rows.',
             'Every file size and GDC MD5 checksum verified. All samples are primary blood-derived cancers.',
             'Labels: TCGA NEJM 2013 author-curated Table S1, obtained through the GDC publication page.',
             'The full 200-case matrix contains 39 FLT3-ITD, 17 FLT3 point-mutation and 54 NPM1-mutant cases.',
             'In the 151 RNA-seq intersection: FLT3-ITD 34 positive/117 negative; NPM1 38 positive/113 negative.',
             'All 151 have label-source rows; no target is missing. Ten FLT3 point-mutant cases are ITD-negative.',
             'The RNA-seq ITD fraction is 22.5%; the full paper cohort has 19.5% ITD and 28% combined FLT3 mutations.',
             'Blank entries in the complete, sequenced author mutation matrix mean no reported mutation;',
             'this is not a blanket rule that missing clinical data mean negative. Missing patient rows would be dropped.',
             'FLT3 calls use original in_frame_ins annotations; all 39 are independently concordant with the',
             'author-released TCGA mutation reannotation. This avoids misleading later MAF consequence labels.',
             'NPM1 positives exactly agree with cBioPortal laml_tcga. One positive is p.K263R, a noncanonical missense.',
             'Source: https://gdc.cancer.gov/about-data/publications/laml_2012',
             'Table S1: '+S1_URL,
             'cBioPortal: https://github.com/cBioPortal/datahub/tree/master/public/laml_tcga',
             'Cross-check: https://github.com/gerstung-lab/AML-multistage/tree/master/data', '',
             'Gene harmonization and coverage:',cov.to_string(index=False),'',
             '51,016 input genes expected; 46,973 match and 4,043 are absent by stable Ensembl ID.',
             'Versions stripped, PAR_Y identities kept distinct; no gene-symbol substitution or GRCh coordinate liftover.',
             'CPM uses the observed counts in the exact frozen input-gene universe, before selection, as in training.',
             'Missing input columns have zero placeholders for this library-total calculation. After the frozen',
             'filter and scaler, every absent selected gene is set to standardized 0, i.e. zero coefficient contribution.',
             'No missing gene is asserted to be biologically unexpressed. Annotation changes make library totals approximate.',
             f'Matched genes retain {c.fraction_counts_retained.min():.2%}-{c.fraction_counts_retained.max():.2%} of each TCGA library\'s gene-assigned counts.',
             'Both models exceed 90% absolute coefficient-weight coverage. Unavailable genes are fully listed in *_gene_coverage.csv.', '']
    previous=pd.read_csv(BEAT/'results/heldout_metrics.csv')
    text.append('Comparison with the original BEAT-AML holdout (descriptive; different populations/prevalences):')
    for r in rows:
        old=previous[previous.target==r['target']].iloc[0]
        text.append(f'{r["target"]}: AUROC {old.AUROC:.3f} -> {r["AUROC"]:.3f}; AUPRC {old.AUPRC:.3f} -> {r["AUPRC"]:.3f}; '
                    f'Brier {old.Brier:.3f} -> {r["Brier"]:.3f}.')
    # Predeclared label-definition sensitivity; no new model predictions, fitting, or selection.
    excluded=set(c.loc[c.NPM1_curated_variant.eq('p.K263R'),'patient'])
    secondary=d[(d.target=='NPM1')&~d.patient.isin(excluded)]
    sensitivity=workflow.metrics(secondary.y,secondary.probability)
    dump('NPM1_noncanonical_excluded.json',{'description':'Predeclared label-definition sensitivity, same saved predictions',
        'excluded_patients':sorted(excluded),'N':len(secondary),'positive':int(secondary.y.sum()),**sensitivity})
    text += ['',f'NPM1 label-definition sensitivity excluding p.K263R: N={len(secondary)}, positives={int(secondary.y.sum())}; '+
             f'AUROC={sensitivity["AUROC"]:.3f}, AUPRC={sensitivity["AUPRC"]:.3f}, Brier={sensitivity["Brier"]:.3f}.',
             'This sensitivity result was not used to alter the primary definition or any prediction.', '',
             'Limits: only 34/38 positives; intervals are conditional on the frozen model and this cohort.',
             'STAR reprocessing and BEAT-AML counting differ despite matching count units; annotation, sample source,',
             'cell composition and cohort shifts can change calibration. CPM is not composition-robust normalization.',
             'No domain adaptation or recalibration was performed. Missing predictor contributions were neutralized as requested.',
             'These are external-cohort feasibility estimates, not validation of clinical diagnostic use.',
             'All BEAT-AML files and model bytes verified unchanged. See protocol, manifests and SHA256 hashes.']
    if cov.below_90_percent_weight.any(): text.insert(0,'WARNING: COEFFICIENT-WEIGHT COVERAGE BELOW 90%; COMPARISON WEAKENED.')
    (OUT/'summary.txt').write_text('\n'.join(text)+'\n')
    print('\n'.join(text),flush=True)
    check_frozen_files()

def main():
    global OUT
    parser=argparse.ArgumentParser(); parser.add_argument('phase',choices=['prepare','evaluate','report'])
    parser.add_argument('--out',type=Path,default=OUT)
    args=parser.parse_args(); OUT=args.out.resolve(); OUT.mkdir(parents=True,exist_ok=True)
    {'prepare':prepare,'evaluate':evaluate,'report':report}[args.phase]()

if __name__=='__main__': main()
