"""Check frozen identities and regenerated results, without protein-model imports."""
from pathlib import Path
import csv,hashlib,json,math,sys
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(n):
 with (ROOT/'data'/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
provenance=json.loads((ROOT/'provenance/checkpoint-b-data-provenance.json').read_text())
for s in provenance['sources']:assert sha(ROOT/'repro/inputs'/s['path'])==s['sha256'],s['path']
for s in json.loads((ROOT/'provenance/source-files.json').read_text()):assert sha(ROOT/s['destination'])==s['sha256'],s['destination']
checkpoint=json.loads((ROOT/'provenance/checkpoint-b-verification.json').read_text())
for s in checkpoint['artifacts']:
 if s['path'].startswith('data/') or (s['path'].startswith('figures/') and s['path'].endswith('.pdf')):
  assert sha(ROOT/s['path'])==s['sha256'],s['path']
files=list((ROOT/'data').glob('*.tsv'));assert len(files)==11
for f in files:assert f.read_bytes()==(ROOT/'repro/data'/f.name).read_bytes(),f.name
models={}
for s in provenance['sources']:
 if s['path'].endswith('.json') and 'analysis' in s['path']:
  data=json.loads((ROOT/'repro/inputs'/s['path']).read_text())
  for k,v in data['models'].items():
   if k.startswith('esm2_'):models['650M' if '650m' in k else '3B']=v
for r in rows('esm-all-fixed-blocks.tsv'):
 b=models[r['model']][r['block']];v=b['genotype_native'][r['genotype']];eta=b['eta'];m=float(r['M']);d=float(r['D'])
 for k in ['M','D','B','T','M_over_D','L','U']:assert float(r[k])==v[k]
 assert d>4*eta
 assert math.isclose(float(r['M_over_D']),m/d,abs_tol=1e-12)
 assert math.isclose(float(r['L']),max(0,m-8*eta)/(d+4*eta),abs_tol=1e-12)
 assert math.isclose(float(r['U']),(m+8*eta)/(d-4*eta),abs_tol=1e-12)
counts={}
for b in ['8','11','16']:
 pairs=[r for r in rows('profam-matched-norm-differences.tsv') if r['block']==b];assert len(pairs)==10
 for r in pairs:assert math.isclose(float(r['difference_of_norms']),float(r['hom_norm'])-float(r['nonhom_norm']),abs_tol=1e-12)
 counts[b]=sum(float(r['difference_of_norms'])>0 for r in pairs)
assert counts=={'8':0,'11':0,'16':8}
run=json.loads((ROOT/'historical-execution/original/design/run-plan.json').read_text())
assert run['masks']['genotype_primary']==[i for i in range(337,431) if i not in [401,403,404,423,425,428]]
assert run['masks']['WT_primary']==list(range(337,431))
for model,blocks in models.items():
 for key,b in blocks.items():
  for genotype,g in b['genotype_native'].items():
   count=sum(x['M']>=g['M'] for x in b['genotype_shuffled'][genotype])
   assert count==b['interpretation']['shuffle_counts_at_or_above_native'][genotype]
report={'status':'passed','scope':'saved-summary reconstruction and arithmetic; no new inference','frozen_inputs':7,'historical_source_files':len(json.loads((ROOT/'provenance/source-files.json').read_text())),'tables_byte_identical':11,'fixed_ESM_metric_rows':30,'profam_positive_matched_differences':counts,'mask_sizes':[88,94],'new_model_forwards':0,'python':sys.version.split()[0]}
(ROOT/'provenance/reproduction-check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
