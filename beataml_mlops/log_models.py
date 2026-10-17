"""Log the frozen BEAT-AML elastic-net models (params, CV/held-out/TCGA metrics,
coefficients) to MLflow and register them. No model is fitted, retrained, tuned
or modified here; everything comes from already-saved results/ files."""
import os
import sys
sys.dont_write_bytecode = True  # importing the original pipeline must not touch its directory
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BEAT = ROOT.parent / 'beataml_modeling'
TCGA = ROOT.parent / 'tcga_laml_validation'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/beataml-mpl')
sys.path.insert(0, str(BEAT))
import json
import joblib
import mlflow
import mlflow.sklearn
import pandas as pd
import workflow  # original frozen pipeline classes; never call training functions

FAMILY = 'elastic_net'
MLRUNS = ROOT / 'mlruns'
EXPERIMENT = 'beataml-frozen-models'


def load_frozen_metrics(target):
    selection = json.loads((BEAT / 'results/selection.json').read_text())[target]
    cv = pd.read_csv(BEAT / 'results/cv_summary.csv')
    cv_row = cv[(cv.target == target) & (cv.model == FAMILY)].iloc[0]
    heldout = pd.read_csv(BEAT / 'results/heldout_metrics.csv')
    held = heldout[heldout.target == target].iloc[0]
    tcga_metrics_path = TCGA / 'results/metrics.csv'
    tcga = None
    if tcga_metrics_path.exists():
        t = pd.read_csv(tcga_metrics_path)
        tcga = t[t.target == target].iloc[0]
    return selection, cv_row, held, tcga


def log_target(target):
    pipe = joblib.load(BEAT / f'results/{target}_{FAMILY}.joblib')
    selection, cv_row, held, tcga = load_frozen_metrics(target)
    features = pd.read_csv(BEAT / f'results/{target}_{FAMILY}_features.csv')
    model_name = f'beataml-{target.lower()}-{FAMILY}'
    with mlflow.start_run(run_name=f'{target}-{FAMILY}-frozen') as run:
        mlflow.set_tags({
            'target': target, 'model_family': FAMILY, 'frozen': 'true',
            'source': 'beataml_modeling/results (no refitting)',
            'selection_rule': 'mean CV average precision, ties AUROC then family order',
        })
        mlflow.log_params({
            **selection['params'],
            'threshold': selection['threshold'],
            'development_n': selection['development_n'],
            'test_n': selection['test_n'],
            'seed': workflow.SEED,
        })
        for metric in workflow.METRICS:
            mlflow.log_metric(f'cv_{metric}_mean', cv_row[f'{metric}_mean'])
            mlflow.log_metric(f'cv_{metric}_sd', cv_row[f'{metric}_sd'])
            mlflow.log_metric(f'heldout_{metric}', held[metric])
        if tcga is not None:
            for metric in workflow.METRICS:
                if metric in tcga and pd.notna(tcga[metric]):
                    mlflow.log_metric(f'tcga_external_{metric}', tcga[metric])

        coeff_path = ROOT / f'_tmp_{target}_coefficients.csv'
        features[['stable_id', 'display_label', 'weight', 'literature_panel_member']].to_csv(coeff_path, index=False)
        mlflow.log_artifact(str(coeff_path), artifact_path='coefficients')
        coeff_path.unlink()
        png = BEAT / f'results/{target}_{FAMILY}_features.png'
        if png.exists():
            mlflow.log_artifact(str(png), artifact_path='plots')

        model_info = mlflow.sklearn.log_model(pipe, name='model', registered_model_name=model_name)
        client = mlflow.tracking.MlflowClient()
        version = model_info.registered_model_version
        client.set_registered_model_alias(model_name, 'frozen', version)
        client.set_model_version_tag(model_name, version, 'description', (
            f'Frozen {FAMILY} pipeline for {target} from beataml_modeling; fitted once on '
            'development data and never refit. See the logging run for CV, BEAT-AML held-out '
            'and TCGA-LAML external-validation metrics.'
        ))
        print(f'{target}: run {run.info.run_id}, registered {model_name} v{version} (alias "frozen")')


def main():
    MLRUNS.mkdir(exist_ok=True)
    mlflow.set_tracking_uri(f'file://{MLRUNS}')
    mlflow.set_experiment(EXPERIMENT)
    for target in workflow.TARGETS:
        log_target(target)
    print(f"\nRun 'mlflow ui --backend-store-uri file://{MLRUNS}' from {ROOT} to browse runs and the registry.")


if __name__ == '__main__':
    main()
