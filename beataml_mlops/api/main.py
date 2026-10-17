"""FastAPI service wrapping the frozen BEAT-AML elastic-net pipelines.

Loads read-only copies of the frozen models plus the original LogCPM/
FoldGeneFilter transformer classes (./artifacts/, produced by
prepare_artifacts.py); this process never fits, retrains, tunes or
otherwise modifies any model. See ../README.md for what this is and is not
validated for before using it in any real decision.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent
ARTIFACTS = ROOT / 'artifacts'
sys.path.insert(0, str(ARTIFACTS))  # exposes workflow.py so joblib can resolve the pipeline's classes

import joblib
import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
import workflow  # noqa: F401  (import kept so LogCPM/FoldGeneFilter are resolvable when unpickling)

TARGETS = ('FLT3-ITD', 'NPM1')
THRESHOLD = 0.5
Z_WARN, Z_CRITICAL, WARN_FRACTION = 5.0, 10.0, 0.01

ANNOTATION = pd.read_csv(ARTIFACTS / 'gene_annotation.csv')
GENE_ORDER = ANNOTATION.stable_id.tolist()
GENE_INDEX = {gene_id: position for position, gene_id in enumerate(GENE_ORDER)}
N_GENES = len(GENE_ORDER)
PIPELINES = {target: joblib.load(ARTIFACTS / f'{target}_elastic_net.joblib') for target in TARGETS}

app = FastAPI(
    title='BEAT-AML mutation-probability API',
    description=(
        'Frozen elastic-net FLT3-ITD/NPM1 classifiers from beataml_modeling, externally '
        'validated once on TCGA-LAML without refitting. This is a feasibility-study model, '
        'not a validated clinical diagnostic; see beataml_modeling/README.md and '
        'tcga_laml_validation/README.md for the full protocol and limitations.'
    ),
    version='1.0.0',
)


class ExpressionRequest(BaseModel):
    gene_counts: Dict[str, float] = Field(..., description=(
        'Raw (unnormalized) RNA-seq counts keyed by Ensembl stable_id, e.g. "ENSG00000171246". '
        'Genes absent from this dict are treated as zero, both for the sample library-size total '
        'and for the model. Supplying only a subset of the training gene universe makes the '
        'CPM library-size denominator approximate; see "coverage_warning" in the response.'
    ))


class SelectedGene(BaseModel):
    stable_id: str
    display_label: Optional[str] = None
    coefficient: float


class GeneListResponse(BaseModel):
    target: str
    n_selected_genes: int
    genes: List[SelectedGene]


def to_vector(gene_counts: Dict[str, float]):
    vector = np.zeros(N_GENES, dtype=np.float64)
    supplied = np.zeros(N_GENES, dtype=bool)
    unknown = []
    for gene_id, value in gene_counts.items():
        position = GENE_INDEX.get(gene_id)
        if position is None:
            unknown.append(gene_id)
            continue
        if value < 0 or not np.isfinite(value):
            raise HTTPException(422, f'Invalid count for {gene_id}: {value}')
        vector[position] = value
        supplied[position] = True
    if vector.sum() <= 0:
        raise HTTPException(422, 'All supplied counts are zero or unrecognized; cannot compute a library size.')
    return vector, supplied, unknown


def standardize(pipe, X):
    log_values = pipe.named_steps['log_cpm'].transform(X)
    selected = pipe.named_steps['filter'].transform(log_values)
    return pipe.named_steps['scale'].transform(selected)


def ood_flags(z, coefficients):
    """Same heuristic (and thresholds) as beataml_mlops/monitor.py, generalized to one sample."""
    active = coefficients != 0
    za = np.abs(z[:, active])
    warn_fraction = float((za > Z_WARN).mean())
    max_z = float(za.max())
    flag = (warn_fraction > WARN_FRACTION) or (max_z > Z_CRITICAL)
    return warn_fraction, max_z, flag


@app.get('/health')
def health():
    return {'status': 'ok', 'targets': list(TARGETS), 'threshold': THRESHOLD, 'n_input_genes': N_GENES}


@app.get('/targets/{target}/genes', response_model=GeneListResponse)
def required_genes(target: str):
    if target not in PIPELINES:
        raise HTTPException(404, f'Unknown target "{target}"; choose one of {TARGETS}')
    pipe = PIPELINES[target]
    indices = pipe.named_steps['filter'].indices_
    coefficients = pipe.named_steps['model'].coef_[0]
    genes = ANNOTATION.iloc[indices][['stable_id', 'display_label']].copy()
    genes['coefficient'] = coefficients
    return {'target': target, 'n_selected_genes': len(indices), 'genes': genes.to_dict(orient='records')}


@app.post('/predict/{target}')
def predict(target: str, request: ExpressionRequest):
    if target not in PIPELINES:
        raise HTTPException(404, f'Unknown target "{target}"; choose one of {TARGETS}')
    pipe = PIPELINES[target]
    X, supplied, unknown = to_vector(request.gene_counts)
    z = standardize(pipe, X[None, :])
    indices = pipe.named_steps['filter'].indices_
    coefficients = pipe.named_steps['model'].coef_[0]
    missing_selected = ~supplied[indices]
    z[:, missing_selected] = 0.0  # same zero-impute rule as the TCGA external validation
    probability = float(pipe.named_steps['model'].predict_proba(z)[0, 1])
    warn_fraction, max_z, flag = ood_flags(z, coefficients)
    return {
        'target': target,
        'threshold': THRESHOLD,
        'probability': probability,
        'prediction': int(probability >= THRESHOLD),
        'selected_gene_coverage': {
            'present': int((~missing_selected).sum()),
            'total': len(indices),
            'fraction': float((~missing_selected).mean()),
        },
        'out_of_distribution': {'flag': flag, 'warn_fraction_abs_z_gt_5': warn_fraction, 'max_abs_z': max_z},
        'unknown_gene_ids_ignored': unknown[:20],
        'coverage_warning': (
            None if supplied.all() else
            'Fewer than all 51,016 training input genes were supplied; the CPM library-size '
            'denominator is approximate, as in the TCGA external validation.'
        ),
    }
