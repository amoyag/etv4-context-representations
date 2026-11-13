#!/usr/bin/env python3
"""E0 source reconciliation and deterministic design generation. No model imports/runs."""
import csv, hashlib, itertools, json, math, random, re, subprocess
from pathlib import Path
from collections import Counter
from fractions import Fraction
from Bio.Align import PairwiseAligner
ROOT = Path(__file__).resolve().parents[1]
if (ROOT / 'design/frozen-protocol.json').exists():
    raise SystemExit('E0 is frozen. Use scripts/verify_e0.py; preserve the old freeze before an explicitly versioned design amendment.')
RAW = ROOT / 'data/raw/E0-20261110'
AA = 'ACDEFGHIKLMNPQRSTVWY'
def sha(s): return hashlib.sha256(s.encode()).hexdigest()
def read(p): return json.loads(p.read_text())
def writej(p,x): p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def tsv(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def seq(r): return r['sequence']['value']
def freqs(s): return [s.count(a)/len(s) for a in AA]
def dist(a,b): return sum(abs(x-y) for x,y in zip(freqs(a),freqs(b)))
def exact_composition_cost(a,b):
    assert len(a)==len(b)
    return Fraction(sum(abs(a.count(x)-b.count(x)) for x in AA),len(a))
def charge(s): return (sum(s.count(a) for a in 'KR')-sum(s.count(a) for a in 'DE'))/len(s)
def shuf(s,seed,endpoints=True):
    a=list(s[1:-1] if endpoints else s);random.Random(seed).shuffle(a)
    return (s[0]+''.join(a)+s[-1]) if endpoints else ''.join(a)
ref=read(RAW/'P43268.json'); wt=seq(ref); assert len(wt)==484
s1=(RAW/'supplement-S1-pdftotext-raw.txt').read_text()
parts=list(dict.fromkeys(re.findall(r'^ETV4 ([A-Z-]+)$',s1,re.M)))
assert ''.join(x.replace('-','') for x in parts)==wt
for p,a in {401:'Y',403:'Y',404:'E',423:'E',425:'E',428:'F',432:'F'}.items(): assert wt[p-1]==a
source_checks=[]
for s in read(ROOT/'sources/manifest.json')['sources']:
    p=Path(s['path']); h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==s['sha256'];source_checks.append({'id':s['id'],'sha256':h,'passed':True})
# Source tables: reported values retain their original scale and SEM.
values=[(337,430,6.1,.6,1,.1,25),(337,436,12,2,2,.4,23),(337,484,11,2,1.8,.4,4),(165,430,12,1,2.1,.2,3),(165,436,60,10,10,3,11),(165,484,66,9,11,2,18),(1,484,83,8,14,2,35)]
mutants={'WT':[], 'Y401N_Y403G':[(401,'Y','N'),(403,'Y','G')], 'E404K':[(404,'E','K')], 'E423K':[(423,'E','K')], 'E425K':[(425,'E','K')], 'F428A_F432A':[(428,'F','A'),(432,'F','A')]}
constructs=[]; measurements=[]; inputs={}
def cid(a,b,g='WT'): return f'ETV4_{a}_{b}_{g}'
def add_input(id,s,kind,origin,seed=None,native_start=None):
    assert set(s)<=set(AA)
    inputs[id]={'sequence':s,'length':len(s),'sha256':sha(s),'evidence_type':kind,'origin':origin,'seed':seed,'native_start':native_start}
for g,ms in mutants.items():
    full=list(wt)
    for p,a,b in ms: assert full[p-1]==a;full[p-1]=b
    full=''.join(full)
    bounds=[(a,b) for a,b,*_ in values] if g=='WT' else [(1,484),(337,436)]
    for a,b in bounds:
        id=cid(a,b,g);s=full[a-1:b];add_input(id,s,'published_native_construct','Currie_S3' if g=='WT' else 'Currie_S10',native_start=a)
        constructs.append(dict(construct_id=id,species='Homo sapiens',taxon_id=9606,accession='P43268',isoform='canonical; exact S1 match',sequence_version=3,start=a,end=b,genotype=g,context='FL' if a==1 else f'{a}-{b}',length=len(s),sequence=s,sequence_sha256=sha(s),source_id='Currie_NAR_2017_supplement',locator='S1 p5; S3 p4' if g=='WT' else 'S1 p5; S10 p14',input_convention='native amino acids only',expression='E. coli; pET28; His tag cleaved',scar_status='residual EMSA cloning scar unspecified',measurement_support='numeric WT KD' if g=='WT' else 'qualitative mutant comparison; no tabulated absolute KD',limitations='Not an exact reconstruction of any unreported vector scar; NMR ligation scars excluded'))
def measurement(id,construct,g,context,endpoint,units,value,uncertainty,n,etype,locator,limits='',utype='SEM',comparison=''):
    measurements.append(dict(observation_id=id,source_id='Currie_NAR_2017_supplement',locator=locator,construct_id=construct,genotype=g,context=context,endpoint=endpoint,units=units,value=value,uncertainty_type=utype,uncertainty=uncertainty,n=n,evidence_type=etype,comparison_construct_id=comparison,assay_id='Currie_EMSA_27bp',limitations=limits))
for a,b,k,e,fold,fe,n in values:
    id=cid(a,b)
    measurement(id+'_KD',id,'WT',f'{a}-{b}','equilibrium_KD','1e-11 M',k,e,n,'numeric','Table S3 p4','Reported aggregated replicates; raw titrations/covariances unavailable')
    measurement(id+'_fold',id,'WT',f'{a}-{b}','fold_inhibition_vs_337_430','dimensionless',fold,fe,n,'numeric','Table S3 p4','Do not treat shared reference ratios as independent',comparison=cid(337,430))
categories={'WT':('autoinhibited','reference'), 'Y401N_Y403G':('NID_associated_inhibition_relieved','DNA_affinity_weakened'), 'E404K':('autoinhibition_reported_preserved','DNA_affinity_strengthened_in_both_contexts'), 'E423K':('autoinhibition_reported_preserved','no_DNA_binding_effect_reported'), 'E425K':('autoinhibition_reported_preserved','no_DNA_binding_effect_reported'), 'F428A_F432A':('autoinhibition_reported_preserved','no_DNA_binding_effect_reported')}
for g,n in zip(mutants,[23,5,4,4,4,6]):
    for endpoint,cat in zip(['FL_to_CID_regulation_category','absolute_affinity_change_category'],categories[g]):
        measurement(g+'_'+endpoint,cid(1,484,g),g,'FL vs 337-436',endpoint,'category',cat,'',n if endpoint.startswith('FL') else '', 'qualitative','Figure S10B and caption p14','No absolute mutant KD transcribed; no equivalence test. Affinity assertions for E404K and remaining controls are caption statements with data not shown. n belongs to S10B ratio, not each construct.',utype='not_extracted',comparison=cid(337,436,g))
# Preserve distinct sample sizes without turning a plotted curve into a new numeric affinity.
for a,b,g,n in [(337,436,'WT',7),(1,484,'WT',14),(337,436,'Y401N_Y403G',3),(1,484,'Y401N_Y403G',4)]:
    measurement(cid(a,b,g)+'_S10A',cid(a,b,g),g,f'{a}-{b}','binding_isotherm_displayed','fraction_bound','not_digitised','',n,'qualitative','Figure S10A p14','Curve-level replicate count, not the S3 or S10B n; no fitted absolute KD extracted',utype='SEM_plotted_not_extracted')
assert len(constructs)==17 and len({c['sequence_sha256'] for c in constructs})==17
# Conditioning panel, sequence-only criteria.
domain=wt[338:420]
def bestmatch(s):
    k=max(range(len(s)-81),key=lambda k:sum(a==b for a,b in zip(domain,s[k:k+82])))
    return k+1,sum(a==b for a,b in zip(domain,s[k:k+82]))/82
align=PairwiseAligner(mode='global',match_score=1,mismatch_score=0,open_gap_score=-2,extend_gap_score=-.5)
def global_id(s):
    a=align.align(wt,s)[0];return sum(x==y and x!='-' for x,y in zip(a[0],a[1]))/max(len(wt),len(s))
hom=[];hqc=[];seen=set()
for r in sorted(read(RAW/'homolog-candidates.json')['results'],key=lambda r:r['primaryAccession']):
    s=seq(r);pos,identity=bestmatch(s);fid=global_id(s)
    keep=400<=len(s)<=550 and set(s)<=set(AA) and not r.get('proteinDescription',{}).get('flag')=='Fragment' and .6<=identity<1 and fid<.9 and s not in seen
    hqc.append({'accession':r['primaryAccession'],'domain_start':pos,'domain_identity':identity,'full_identity':fid,'length':len(s),'retained':keep,'species':r['organism']['scientificName'],'sequence_version':r['entryAudit']['sequenceVersion']})
    if keep: hom.append(r);seen.add(s)
controls=[];cqc=[]
for r in read(RAW/'control-candidates.json')['results']:
    s=seq(r);pos,identity=bestmatch(s)
    annotation=json.dumps([r.get('proteinDescription'),r.get('comments'),r.get('uniProtKBCrossReferences')]).lower()
    keep=set(s)<=set(AA) and identity<.3 and 'ets' not in re.findall(r'\b[a-z]+\b',annotation) and r.get('proteinDescription',{}).get('flag')!='Fragment'
    cqc.append({'accession':r['primaryAccession'],'length':len(s),'max_ungapped_domain_identity':identity,'retained':keep})
    if keep: controls.append(r)
pairs=random.Random(20261110).sample(list(itertools.combinations(hom,2)),5)
prompts={'empty':{'accessions':[],'sequences':[],'kind':'empty','set':0,'order':'none','token_prefix_length':2}}
matching=[]
for j,pair in enumerate(pairs,1):
    eligible=[[r for r in controls if len(seq(r))==len(seq(h))] for h in pair]
    candidates=[(exact_composition_cost(seq(pair[0]),seq(a))+exact_composition_cost(seq(pair[1]),seq(b)),a['primaryAccession'],b['primaryAccession'],a,b) for a in eligible[0] for b in eligible[1] if a['primaryAccession']!=b['primaryAccession']]
    assert candidates,'No exactly length-matched distinct controls'
    best=min(candidates,key=lambda x:x[:3]); cp=best[3:]
    matching.append({'set':j,'homologues':[r['primaryAccession'] for r in pair],'controls':[r['primaryAccession'] for r in cp],'lengths':[len(seq(r)) for r in pair],'composition_L1':[dist(seq(h),seq(c)) for h,c in zip(pair,cp)],'charge_proxy_difference':[charge(seq(h))-charge(seq(c)) for h,c in zip(pair,cp)],'taxonomy_homologues':[r['organism']['scientificName'] for r in pair],'taxonomy_controls':[r['organism']['scientificName'] for r in cp]})
    for kind,rs in [('hom',pair),('nonhom',cp),('homshuffle',pair)]:
        ss=[shuf(seq(r),14000+10*j+k,False) if kind=='homshuffle' else seq(r) for k,r in enumerate(rs)]
        for order,ix in [('ab',[0,1]),('ba',[1,0])]:
            prompts[f'{kind}_{j}_{order}']={'accessions':[rs[k]['primaryAccession'] for k in ix],'sequences':[ss[k] for k in ix],'sequence_sha256':[sha(ss[k]) for k in ix],'kind':kind,'set':j,'order':order,'token_prefix_length':2+sum(len(s)+1 for s in ss),'seeds':[14000+10*j+k for k in ix] if kind=='homshuffle' else []}
# Model-only flank controls; retain coordinate map by origin and unchanged lengths.
for k in range(20):
    n=shuf(wt[164:336],11001+k);c=shuf(wt[430:436],12001+k); d=wt[336:430]
    add_input(f'Nshuffle_{k+1}_430',n+d,'model-only','N=165:336; domain=337:430',11001+k,165)
    add_input(f'Cshuffle_{k+1}',d+c,'model-only','domain=337:430; C=431:436',12001+k,337)
    add_input(f'NCshuffle_{k+1}',n+d+c,'model-only','paired N and C controls',[11001+k,12001+k],165)
    add_input(f'Nshuffle_{k+1}_436',n+wt[336:436],'model-only','N shuffled; 337:436 unchanged',11001+k,165)
    nfull=shuf(wt[:336],13001+k)
    for g in mutants:
        suffix=inputs[cid(1,484,g)]['sequence'][336:]
        add_input(f'FLshuffle_{k+1}_{g}',nfull+suffix,'model-only','1:336 shuffled; 337:484 native genotype',13001+k,1)
# Run list stores each named input/prompt task once, with all consuming blocks.
runs={}
def run(model,input_id,prompt='empty',block='primary',required=True):
    key=(model,input_id,prompt)
    if key in runs:
        if block not in runs[key]['blocks']:runs[key]['blocks'].append(block)
        return
    x=inputs[input_id];pr=prompts[prompt]; nt=x['length']+2 if model=='esm2_650m' else pr['token_prefix_length']+x['length']+1
    assert nt<= (1024 if model=='esm2_650m' else 2048)
    runs[key]={'run_id':f'{model}__{input_id}__{prompt}','model':model,'input_id':input_id,'prompt_id':prompt,'blocks':[block],'required':required,'expected_tokens':nt,'input_sha256':x['sha256'],'layers':[17,22,33] if model=='esm2_650m' else [8,11,16],'weight_execution_status':'not_run'}
for c in constructs:run('esm2_650m',c['construct_id'],block='measured_panel')
for k in range(1,21):
    for x in [f'Nshuffle_{k}_430',f'Cshuffle_{k}',f'NCshuffle_{k}',f'Nshuffle_{k}_436']:run('esm2_650m',x,block='WT_flank_controls')
    for g in mutants:run('esm2_650m',f'FLshuffle_{k}_{g}',block='genotype_context_controls')
# 217 ESM passes: 17 + 80 + 120; extra E1 qualification budget is separate and bounded below.
wtfour=[cid(337,430),cid(165,430),cid(337,436),cid(165,436)]
for pr in prompts:
    for x in wtfour:run('profam1',x,pr,'WT_N_context_and_visibility')
for c in constructs:run('profam1',c['construct_id'],block='secondary_partial_context_panel')
for pr in ['empty','hom_1_ab','hom_1_ba','nonhom_1_ab','nonhom_1_ba','homshuffle_1_ab','homshuffle_1_ba']:
    for k in range(1,21):run('profam1',f'Nshuffle_{k}_436',pr,'N_flank_specificity')
counts=Counter(r['model'] for r in runs.values())
assert counts['esm2_650m']<=250 and counts['profam1']<=1200
masks={'genotype_primary':[i for i in range(337,431) if i not in [401,403,404,423,425,428,432]],'WT_primary':list(range(337,431)),'direct_substitutions_in_shared_segment':[401,403,404,423,425,428],'substitutions_outside_primary_mask':[432],'structural_ETS_domain':list(range(339,421)),'source_surface_markers':[401,403,407,420,430]}
plan={'schema_version':1,'stage':'E0','status':'specified_preinference','date':'2026-11-10','no_model_execution':True,'models':{'esm2_650m':{'layers':33,'primary_block':22,'sensitivity_blocks':[17,33]},'profam1':{'layers':16,'primary_block':11,'sensitivity_blocks':[8,16],'config':'models/profam-1/.hydra/config.yaml'}},'masks':masks,'random_generator':'Python random.Random MT19937; input sequences materialized and hashed','seeds':{'prompt_pairs':20261110,'N':list(range(11001,11021)),'C':list(range(12001,12021)),'FL_N':list(range(13001,13021))},'inputs_file':'data/derived/model-inputs.json','prompts_file':'data/derived/prompt-panel.json','run_counts':dict(counts),'caps':{'profam_forward_passes':1200,'esm_embedding_passes':250,'masked_mutation_passes':100},'E1_qualification_reserve':{'profam_forward_passes':40,'esm_embedding_passes':30,'uses_only_unlabelled_artificial_sequences':True},'secondary_mutation_scores':{'enabled':False,'reason':'Prespecified omission; use representation and CID-construct controls; no score fishing'},'optional_models_and_ETV1_ETV5_biological_tests':{'enabled':False},'epsilon_activation_normalization':1e-12,'zero_vector_policy':'Fail readout if activation norm <= epsilon; no inflated denominator epsilon','numerical_floor_rule':'E1 establishes eta per model/layer as max(1e-7,10 times maximal per-residue normalized-vector L2 discrepancy among three repeats of each artificial input). Undefined T when WT RMS context denominator <= 4*eta; effect below 8*eta cannot support specificity. Freeze backend and eta in E1 before any biological input.','metrics_document':'design/protocol-v1.md','runs':list(runs.values())}
tsv(ROOT/'data/derived/construct-registry.tsv',constructs);tsv(ROOT/'data/derived/measurement-registry.tsv',measurements)
writej(ROOT/'data/derived/model-inputs.json',inputs);writej(ROOT/'data/derived/prompt-panel.json',prompts)
writej(ROOT/'data/derived/prompt-selection.json',{'rules':'design/E0-selection-rules.md','homologue_candidates':hqc,'control_candidates':cqc,'matching':matching,'pretraining_overlap':'Not excluded. Focal and prompt accessions may occur in UniRef/AFDB/ProFam Atlas; no held-out generalization claim. Atlas membership not exhaustively queried.'})
writej(ROOT/'design/run-plan.json',plan)
writej(ROOT/'data/derived/source-reconciliation.json',{'source_hash_checks':source_checks,'ETV4_S1_unique_chunks':parts,'S1_matches_UniProt':True,'accession':'P43268','sequence_version':3,'sequence_sha256':sha(wt),'number_native_constructs':len(constructs),'number_measurement_records':len(measurements),'WT_factorial_log_interaction':math.log(60*6.1/(12*12)),'WT_factorial_ratio':60*6.1/(12*12),'WT_factorial_error':'No SE computed: covariance/raw data unavailable','EMSA_scars':'Not specified; no exact scar claim','ETV1_optional_excluded':'Source 479 residues; current canonical P50549 477; no substitution of coordinates'})
print(json.dumps({'constructs':len(constructs),'observations':len(measurements),'inputs':len(inputs),'prompts':len(prompts),'run_counts':dict(counts),'homologues':len(hom),'max_tokens':max(r['expected_tokens'] for r in runs.values())},indent=2))
