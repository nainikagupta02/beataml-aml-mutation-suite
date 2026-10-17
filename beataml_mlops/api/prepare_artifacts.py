"""Copy the minimal frozen files the API needs into ./artifacts (read-only
copies; this never modifies beataml_modeling). Run once before `docker build`,
and again if the frozen pipelines are ever regenerated upstream."""
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BEAT = ROOT.parent.parent / 'beataml_modeling'
ARTIFACTS = ROOT / 'artifacts'

FILES = [
    'workflow.py',                              # exposes LogCPM/FoldGeneFilter for joblib to unpickle
    'xlsx_read.py',                             # workflow.py imports this at module level
    'results/gene_annotation.csv',               # defines the exact 51,016-gene input order
    'results/FLT3-ITD_elastic_net.joblib',
    'results/NPM1_elastic_net.joblib',
]


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ARTIFACTS.mkdir(exist_ok=True)
    for rel in FILES:
        src = BEAT / rel
        if not src.exists():
            raise FileNotFoundError(f'Missing frozen artifact: {src}')
        dst = ARTIFACTS / Path(rel).name
        shutil.copy2(src, dst)
        print(f'{sha256(dst)}  {dst.relative_to(ROOT)}')
    print(f'\nCopied {len(FILES)} files into {ARTIFACTS}. Now run: docker build -t beataml-api .')


if __name__ == '__main__':
    main()
