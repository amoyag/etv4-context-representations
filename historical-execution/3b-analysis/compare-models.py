import json,csv,hashlib
from pathlib import Path
p=Path('/Users/aurelio/KnowledgeSystem/areas/research/ministerio-PID2026/pilot/ets-regulatory-context')
r=p/'followups/esm2-3b-20261111/job-3647838'
a=p/'results/E2-analysis-v1.json';b=r/'e2/results/E2-3B-analysis.json'
old=json.loads(a.read_text());new=json.loads(b.read_text());rows=[];details={}
for label,data,model,primary in [('650M',old,'esm2_650m','22'),('3B',new,'esm2_3b','24')]:
 for layer,x in data['models'][model].items():
  y=x['genotype_native']['Y401N_Y403G']; controls={k:v for k,v in x['genotype_native'].items() if k!='Y401N_Y403G'};sh=x['genotype_shuffled']['Y401N_Y403G']
  row=dict(model=label,layer=layer,role='primary' if layer==primary else 'sensitivity',classification=x['interpretation']['classification'],**x['interpretation']['conditions'],eta=x['eta'],Y_M=y['M'],Y_D=y['D'],Y_ratio=y['M_over_D'],Y_ratio_lower=y['L'],Y_ratio_upper=y['U'],max_control_ratio_upper=max(v['U'] for v in controls.values()),raw_M_margin=y['M']-max(v['M'] for v in controls.values()),native_minus_max_shuffle=y['M']-max(v['M'] for v in sh),required_margin=8*x['eta'],shuffle_count=x['interpretation']['shuffle_counts_at_or_above_native']['Y401N_Y403G'])
  rows.append(row)
  details[f'{label}/{layer}']={'genotype_native':x['genotype_native'],'interpretation':x['interpretation'],'factorial_native_rms':x['factorial_native']['rms'],'factorial_shuffle_range':[min(v['rms'] for v in x['factorial_shuffled']),max(v['rms'] for v in x['factorial_shuffled'])],'factorial_shuffle_count':sum(v['rms']>=x['factorial_native']['rms'] for v in x['factorial_shuffled']),'native_N_response':x['native_N_response_with_C']['rms'],'N_shuffle_range':[min(v['rms'] for v in x['N_only_shuffled_response']),max(v['rms'] for v in x['N_only_shuffled_response'])]}
out=r/'local-analysis'
with (out/'model-comparison.tsv').open('x') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(out/'model-comparison.json').write_text(json.dumps({'sources':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [a,b]},'rows':rows,'details':details},indent=2)+'\n')
print(json.dumps(rows,indent=2))
