"""Reproduce public downloads with curl; cached files are integrity-checked.

Network use may require permission in Codex. Never accesses BEAT-AML files.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'
SOURCES={
 'cbio_clinical_patient.txt':'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/laml_tcga/data_clinical_patient.txt',
 'cbio_clinical_sample.txt':'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/laml_tcga/data_clinical_sample.txt',
 'cbio_mutations.txt':'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/laml_tcga/data_mutations.txt',
 'cbio_cases_sequenced.txt':'https://raw.githubusercontent.com/cBioPortal/datahub/master/public/laml_tcga/case_lists/cases_sequenced.txt',
 'TCGA_AML_Table_S1.xlsx':'https://api.gdc.cancer.gov/data/b9196563-a05d-40b8-80dc-640ec712eb06',
 'TCGA_full_sample_list.tsv':'https://api.gdc.cancer.gov/data/24ffea56-2e7e-40d1-b7a2-b14ee4efbe6e',
 'TCGA_Genetic_crosscheck.txt':'https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/TCGA_Genetic.txt',
}

def digest(p,kind='sha256'):
 h=hashlib.new(kind)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def curl(url,path):
 subprocess.run(['curl','-sS','-L','--fail','--retry','3','--max-time','180',url,'-o',str(path)],check=True)

def main():
 parser=argparse.ArgumentParser(); parser.add_argument('--verify-only',action='store_true');args=parser.parse_args()
 DATA.mkdir(exist_ok=True);(DATA/'counts').mkdir(exist_ok=True)
 if not (DATA/'gdc_files.json').exists():
  if args.verify_only:raise FileNotFoundError('gdc_files.json')
  subprocess.run(['curl','-sS','-L','--fail','--get','--data-urlencode','filters@'+str(DATA/'gdc_filter.json'),
   '--data-urlencode','size=1000','--data-urlencode',
   'fields=file_id,file_name,file_size,md5sum,access,data_type,analysis.workflow_type,cases.submitter_id,cases.case_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.sample_id',
   'https://api.gdc.cancer.gov/files','-o',str(DATA/'gdc_files.json')],check=True)
 hits=json.loads((DATA/'gdc_files.json').read_text())['data']['hits']
 config=['parallel','parallel-max = 6','location','fail','silent','show-error','retry = 3','max-time = 180']
 needed=0
 for h in hits:
  path=DATA/'counts'/(h['file_id']+'.tsv')
  if path.exists() and path.stat().st_size==h['file_size'] and digest(path,'md5')==h['md5sum']:continue
  if args.verify_only:raise ValueError('Missing/invalid STAR file: '+str(path))
  config += ['url = "https://api.gdc.cancer.gov/data/'+h['file_id']+'"','output = '+json.dumps(str(path))]
  needed+=1
 if needed:
  cfg=DATA/'download_counts.curl';cfg.write_text('\n'.join(config)+'\n')
  subprocess.run(['curl','--config',str(cfg)],check=True)
 for name,url in SOURCES.items():
  if not (DATA/name).exists():
   if args.verify_only:raise FileNotFoundError(name)
   curl(url,DATA/name)
  with (DATA/name).open('rb') as f: head=f.read(200)
  if b'version https://git-lfs.github.com/spec/v1' in head:raise ValueError('LFS pointer, not data: '+name)
 for h in hits:
  path=DATA/'counts'/(h['file_id']+'.tsv')
  assert path.stat().st_size==h['file_size'] and digest(path,'md5')==h['md5sum']
 # Freeze-dated external metadata can change upstream. Fail rather than silently change labels.
 saved=ROOT/'results/source_hashes.json'
 if saved.exists():
  hashes=json.loads(saved.read_text())
  for name,h in hashes['metadata_and_labels'].items():
   if (ROOT/name).exists():assert digest(ROOT/name)==h,'Source revision changed: '+name
 print(f'Verified {len(hits)} STAR files and curated label sources.')

if __name__=='__main__':main()
