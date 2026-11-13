#!/usr/bin/env python3
"""Verify E0 data/design invariants without loading models or executing inference."""
import csv, hashlib, json, math
from pathlib import Path
from collections import Counter
R=Path(__file__).resolve().parents[1]
def load(p):return json.loads((R/p).read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sha(s):return hashlib.sha256(s.encode()).hexdigest()
def rows(p):return list(csv.DictReader((R/p).open(),delimiter='\t'))
checks=[]
def check(name,condition):
    if not condition: raise AssertionError(name)
    checks.append(name)
reg=rows('data/derived/construct-registry.tsv');obs=rows('data/derived/measurement-registry.tsv');inputs=load('data/derived/model-inputs.json');prompts=load('data/derived/prompt-panel.json');plan=load('design/run-plan.json')
ref=load('data/raw/E0-20261110/P43268.json')['sequence']['value']
check('17 distinct native sequences',len(reg)==len({r['sequence'] for r in reg})==17)
for r in reg:
    a,b=int(r['start']),int(r['end']);s=r['sequence']; expected=list(ref[a-1:b])
    if r['genotype']!='WT':
        for sub in r['genotype'].split('_'):
            pos=int(sub[1:-1]);check('WT amino acid '+sub,ref[pos-1]==sub[0]);expected[pos-a]=sub[-1]
    check('exact mutation and bounds '+r['construct_id'],s==''.join(expected) and len(s)==b-a+1 and sha(s)==r['sequence_sha256'])
    check('input identity '+r['construct_id'],inputs[r['construct_id']]['sequence']==s)
expected={'337_430':(6.1,.6,25),'337_436':(12,2,23),'337_484':(11,2,4),'165_430':(12,1,3),'165_436':(60,10,11),'165_484':(66,9,18),'1_484':(83,8,35)}
kd=[o for o in obs if o['endpoint']=='equilibrium_KD'];check('seven numeric KD only',len(kd)==7)
for o in kd:
    key='_'.join(o['construct_id'].split('_')[1:3]);check('S3 value units SEM n '+key,(float(o['value']),float(o['uncertainty']),int(o['n']))==expected[key] and o['units']=='1e-11 M' and o['uncertainty_type']=='SEM' and o['evidence_type']=='numeric')
check('no numeric mutant labels',all(o['evidence_type']!='numeric' for o in obs if o['genotype']!='WT'))
check('all observations resolve to measured constructs',all(o['construct_id'] in {r['construct_id'] for r in reg} for o in obs))
check('genotype mask 88 excludes mutations',len(plan['masks']['genotype_primary'])==88 and not set(plan['masks']['genotype_primary']) & {401,403,404,423,425,428,432})
for k in range(1,21):
    n=inputs[f'Nshuffle_{k}_430']['sequence'][:172]; c=inputs[f'Cshuffle_{k}']['sequence'][94:]
    check(f'N composition endpoints {k}',Counter(n)==Counter(ref[164:336]) and n[0]==ref[164] and n[-1]==ref[335])
    check(f'C composition endpoints {k}',Counter(c)==Counter(ref[430:436]) and c[0]==ref[430] and c[-1]==ref[435])
    check(f'paired factorial construction {k}',inputs[f'NCshuffle_{k}']['sequence']==n+ref[336:430]+c)
    check(f'native C Nshuffle {k}',inputs[f'Nshuffle_{k}_436']['sequence']==n+ref[336:436])
    prefix=inputs[f'FLshuffle_{k}_WT']['sequence'][:336]
    check(f'FL composition {k}',Counter(prefix)==Counter(ref[:336]) and prefix[0]==ref[0] and prefix[-1]==ref[335])
    for g in ['WT','Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']:
        x=inputs[f'FLshuffle_{k}_{g}'];check(f'shared FL prefix intact genotype suffix {k} {g}',x['sequence']==prefix+inputs[f'ETV4_1_484_{g}']['sequence'][336:] and x['evidence_type']=='model-only')
check('31 prompts',len(prompts)==31)
for j in range(1,6):
    for order in ['ab','ba']:
        h,c,s=[prompts[f'{group}_{j}_{order}'] for group in ['hom','nonhom','homshuffle']]
        check(f'exact matched ordered lengths {j} {order}',[len(x) for x in h['sequences']]==[len(x) for x in c['sequences']]==[len(x) for x in s['sequences']])
        check(f'prompt shuffle compositions {j} {order}',all(Counter(a)==Counter(b) for a,b in zip(h['sequences'],s['sequences'])))
        check(f'no focal duplication {j} {order}',ref not in h['sequences'] and len(set(h['sequences']))==2)
    for group in ['hom','nonhom','homshuffle']:
        check(f'reversed exact order {group} {j}',prompts[f'{group}_{j}_ab']['sequences']==list(reversed(prompts[f'{group}_{j}_ba']['sequences'])))
for id,x in inputs.items():check('stored sequence length hash '+id,len(x['sequence'])==x['length'] and sha(x['sequence'])==x['sha256'])
runs=plan['runs']; check('494 unique named runs',len(runs)==len({x['run_id'] for x in runs})==494)
count=Counter(x['model'] for x in runs);check('planned counts',dict(count)==plan['run_counts']=={'esm2_650m':217,'profam1':277})
for x in runs:
    s=inputs[x['input_id']]['sequence'];p=prompts[x['prompt_id']];is_esm=x['model']=='esm2_650m'
    nt=len(s)+2 if is_esm else 2+sum(len(v)+1 for v in p['sequences'])+len(s)+1
    check('token budget and input hash '+x['run_id'],nt==x['expected_tokens'] and nt<=(1024 if is_esm else 2048) and sha(s)==x['input_sha256'])
distinct=Counter()
for model in count:
    unique={(inputs[x['input_id']]['sequence'],tuple(prompts[x['prompt_id']]['sequences'])) for x in runs if x['model']==model}
    distinct[model]=len(unique)
check('scheduled versus distinct combinations disclosed',dict(distinct)=={'esm2_650m':211,'profam1':253})
q=load('design/qualification-inputs.json')
check('five fixed artificial ProFam cases',len(q['cases'])==5 and q['repeats']==3)
qmap={x['id']:x for x in q['cases']}
for k in ['Q1','Q2']:
    check('artificial causal extension '+k,qmap[k+'_C']['target']==qmap[k]['target']+'ACDEFG' and qmap[k+'_C']['prefixes']==qmap[k]['prefixes'])
check('budget includes qualification reserve',count['esm2_650m']+plan['E1_qualification_reserve']['esm_embedding_passes']<=250 and count['profam1']+plan['E1_qualification_reserve']['profam_forward_passes']<=1200)
for model in ['esm2_650m','profam1']:
    L=plan['models'][model]['layers'];check('fixed layers '+model,plan['models'][model]['primary_block']==math.ceil(2*L/3) and plan['models'][model]['sensitivity_blocks']==[math.ceil(L/2),L])
for source in load('sources/manifest.json')['sources']:check('source hash '+source['id'],digest(Path(source['path']))==source['sha256'])
f=R/'design/frozen-protocol.json'
if f.exists():
    for item in load('design/frozen-protocol.json')['files']:check('frozen hash '+item['path'],digest(R/item['path'])==item['sha256'])
print(json.dumps({'passed':True,'checks_executed':len(checks),'native_constructs':len(reg),'observation_records':len(obs),'E2_planned_passes':dict(count),'E2_executed_passes':0,'weights_loaded':False,'max_expected_ProFam_tokens':max(x['expected_tokens'] for x in runs if x['model']=='profam1'),'WT_factorial_log_ratio':math.log(60*6.1/144),'checks':checks},indent=2))
