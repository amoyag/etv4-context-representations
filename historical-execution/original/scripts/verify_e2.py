"""Separate float64 reconstruction from E2 NPZ files; no model loads.
Does not import production analysis metrics. Reports numerical comparisons,
full finite/coordinate checks and plan-derived coverage.
"""
import argparse, collections, csv, datetime, json, math
from pathlib import Path
import numpy as np
from run_experiments import ROOT, read, digest, write_json

def unit(a):
    a=np.asarray(a,dtype=np.float64)
    lengths=np.sqrt(np.einsum('ij,ij->i',a,a))
    assert np.isfinite(a).all() and np.all(lengths>1e-12)
    return a/lengths[:,None]
def magnitude(a): return float(np.sqrt(np.einsum('ij,ij->',a,a)/len(a)))
def profile(a): return np.sqrt(np.einsum('ij,ij->i',a,a)).tolist()
def ctx(f,c,wf,wc,eta):
    q=unit(f)-unit(c); qw=unit(wf)-unit(wc)
    M=magnitude(q-qw);B=magnitude(qw);D=magnitude(unit(c)-unit(wc))
    return dict(M=M,B=B,D=D,T=M/B if B>4*eta else None,M_over_D=M/D if D>4*eta else None,
        L=max(0,M-8*eta)/(D+4*eta) if D>4*eta else None,U=(M+8*eta)/(D-4*eta) if D>4*eta else None)

def verify(manifests,analysis):
    plan=read('design/run-plan.json');config=read('design/execution-config.json')
    inputs=read('data/derived/model-inputs.json');prompts=read('data/derived/prompt-panel.json')
    actual=json.loads(Path(analysis).read_text());cfg_hash=digest(ROOT/'design/execution-config.json')
    assert actual['execution_config_sha256']==cfg_hash
    records={};coverage=[];arrays={};token_alphabets={};specials={}
    for model in config['models']:
        qual=read(config['models'][model].get('selected_report',read('environment/model-manifest.json')['models'][model]['selected_report']))
        qcases={r['id']:r for r in read('design/qualification-inputs.json')['cases']}
        token_alphabets[model]={}
        for row in qual['forwards']:
            q=qcases[row['case']];mp=row['mapping'];seq=q['target']
            for aa,tok in zip(seq,mp['token_ids'][mp['token_start']:]):
                assert token_alphabets[model].setdefault(aa,tok)==tok
        assert len(token_alphabets[model])==20
        specials[model]=(mp['token_ids'][:2] if model=='profam1' else mp['token_ids'][:1],mp['token_ids'][-1])
    for path in manifests:
        for r in json.loads(Path(path).read_text()):
            assert r['run_id'] not in records
            records[r['run_id']]=r
    expected_ids={r['run_id'] for r in plan['runs']};assert set(records)==expected_ids
    assert actual['fixed_run_count']==len(records)==494
    for p in plan['runs']:
        r=records[p['run_id']];model=p['model'];x=inputs[p['input_id']];spec=config['models'][model]
        for key,value in p.items():assert r[key]==value,(key,p['run_id'])
        assert r['execution_config_sha256']==cfg_hash and r['native_start']==x['native_start']
        assert digest(ROOT/r['activation_file'])==r['sha256']
        assert __import__('hashlib').sha256(x['sequence'].encode()).hexdigest()==p['input_sha256']
        prefixes=prompts[p['prompt_id']]['sequences']
        special,last=specials[model];toks=list(special)
        for seq in prefixes+[x['sequence']]:
            toks.extend(token_alphabets[model][a] for a in seq);toks.append(last)
        start=(2+sum(len(s)+1 for s in prefixes)) if model=='profam1' else 1
        m=r['mapping'];assert m['token_ids']==toks and m['tokens']==len(toks)==p['expected_tokens']
        assert m['token_start']==start and m['target_length']==x['length']
        assert m['position_ids']==list(range(len(toks))) and m['target_token_indices']==list(range(start,start+x['length']))
        with np.load(ROOT/r['activation_file'],allow_pickle=False) as f:
            assert set(f.files)==set(map(str,spec['layers']))
            aa={}
            for layer in f.files:
                a=f[layer];assert a.dtype==np.float32 and a.shape==(x['length'],spec['width'])
                unit(a);aa[layer]=a
        arrays[p['run_id']]=aa
        coverage.append({'run_id':p['run_id'],'model':model,'status':'complete','activation_file':r['activation_file'],'sha256':r['sha256'],'finite':True,'coordinates':True})
    checks=0;max_abs=0.;checked_paths=set()
    def cmp(value,expected,path):
        nonlocal checks,max_abs
        if isinstance(expected,dict):
            for k,v in expected.items():cmp(value[k],v,path+'/'+str(k))
        elif isinstance(expected,list):
            assert len(value)==len(expected),path
            for i,v in enumerate(expected):cmp(value[i],v,path+'/'+str(i))
        elif isinstance(expected,(float,int)) and not isinstance(expected,bool):
            assert value is not None and math.isfinite(value),path
            err=abs(value-expected);max_abs=max(max_abs,err)
            assert math.isclose(value,expected,rel_tol=2e-12,abs_tol=2e-14),(path,value,expected)
            checks+=1;checked_paths.add(path)
        else: assert value==expected,(path,value,expected);checked_paths.add(path)
    gen=['Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']
    for model,spec in config['models'].items():
      for layer in spec['layers']:
        layer=str(layer);eta=spec['eta'][layer];block=actual['models'][model][layer];base=f'/models/{model}/{layer}'
        def compare(key,expected):cmp(block[key],expected,base+'/'+key)
        def h(ident,prompt='empty',mask=None):
            if mask is None:mask=plan['masks']['WT_primary']
            x=inputs[ident];ix=np.array(mask)-x['native_start']
            assert np.all((ix>=0)&(ix<x['length']))
            return arrays[f'{model}__{ident}__{prompt}'][layer][ix].astype(np.float64)
        mask=plan['masks']['genotype_primary']
        def c(g):return h(f'ETV4_337_436_{g}',mask=mask)
        def f(g):return h(f'ETV4_1_484_{g}',mask=mask)
        native={g:ctx(f(g),c(g),f('WT'),c('WT'),eta) for g in gen}
        compare('eta',eta);compare('genotype_native',native)
        direct=[401,403,404,423,425,428,432]
        dqw=unit(h('ETV4_1_484_WT',mask=direct))-unit(h('ETV4_337_436_WT',mask=direct))
        for g in gen:
            df=unit(h(f'ETV4_1_484_{g}',mask=direct));dc=unit(h(f'ETV4_337_436_{g}',mask=direct))
            cmp(block['substitution_site_diagnostics'][g],{'native_positions':direct,'context_change_norm':profile(df-dc-dqw),'CID_change_norm':profile(dc-unit(h('ETV4_337_436_WT',mask=direct)))},base+'/substitution_site_diagnostics/'+g)
        if model=='esm2_650m':
            sh={g:[ctx(h(f'FLshuffle_{s}_{g}',mask=mask),c(g),h(f'FLshuffle_{s}_WT',mask=mask),c('WT'),eta) for s in range(1,21)] for g in gen}
            compare('genotype_shuffled',sh)
            Y=gen[0];K=gen[1:];resolved=all(m['B']>4*eta and m['D']>4*eta for m in native.values())
            cond={'denominators_resolved':resolved,'M_separation':native[Y]['M']>max(native[g]['M'] for g in K)+8*eta,
                  'ratio_bounds_separation':resolved and native[Y]['L']>max(native[g]['U'] for g in K),
                  'native_order_specificity':native[Y]['M']>max(x['M'] for x in sh[Y])+8*eta}
            label='numerically_unresolved' if not resolved else ('source_consistent_native_order_specificity' if all(cond.values()) else 'specificity_not_established')
            compare('interpretation',{'classification':label,'conditions':cond,'shuffle_counts_at_or_above_native':{g:sum(v['M']>=native[g]['M'] for v in sh[g]) for g in gen}})
        zero,n,c0,nc=[h(f'ETV4_{a}_{b}_WT') for a,b in [(337,430),(165,430),(337,436),(165,436)]]
        def fac(z,n,c,nc):
            j=unit(nc)-unit(n)-unit(c)+unit(z)
            return {'rms':magnitude(j),'profile':profile(j)}
        compare('factorial_native',fac(zero,n,c0,nc));q=unit(nc)-unit(c0)
        compare('native_N_response_with_C',{'rms':magnitude(q),'profile':profile(q)})
        if model=='esm2_650m':
            compare('factorial_shuffled',[{'seed':s,**fac(zero,h(f'Nshuffle_{s}_430'),h(f'Cshuffle_{s}'),h(f'NCshuffle_{s}'))} for s in range(1,21)])
            compare('N_only_shuffled_response',[{'seed':s,'rms':magnitude(unit(h(f'Nshuffle_{s}_436'))-unit(c0))} for s in range(1,21)])
        else:
            qs={};pr={}
            for p in prompts:
                pn,pc=[h(ident,p) for ident in ['ETV4_165_436_WT','ETV4_337_436_WT']]
                zn,zc=[h(ident,p) for ident in ['ETV4_165_430_WT','ETV4_337_430_WT']]
                q=unit(pn)-unit(pc);qs[p]=q
                vis=max(max(profile(unit(pn)-unit(zn))),max(profile(unit(pc)-unit(zc))))
                assert vis<=eta,(p,layer,vis,eta)
                pr[p]={'rms':magnitude(q),'profile':profile(q),'causal_extension_error':vis}
            compare('prompt_response',pr)
            compare('prompt_contrasts',[{'set':j,'order':o,'hom_minus_nonhom':magnitude(qs[f'hom_{j}_{o}']-qs[f'nonhom_{j}_{o}']),'hom_minus_homshuffle':magnitude(qs[f'hom_{j}_{o}']-qs[f'homshuffle_{j}_{o}'])} for j in range(1,6) for o in ['ab','ba']])
            compare('against_empty',{p:magnitude(q-qs['empty']) for p,q in qs.items()})
            compare('order_contrasts',{f'{g}_{j}':magnitude(qs[f'{g}_{j}_ab']-qs[f'{g}_{j}_ba']) for g in ['hom','nonhom','homshuffle'] for j in range(1,6)})
            ns={}
            for p in ['empty']+[f'{g}_1_{o}' for g in ['hom','nonhom','homshuffle'] for o in ['ab','ba']]:
                ns[p]=[]
                for s in range(1,21):
                    q=unit(h(f'Nshuffle_{s}_436',p))-unit(h('ETV4_337_436_WT',p))
                    ns[p].append({'seed':s,'rms':magnitude(q),'profile':profile(q)})
            compare('N_shuffles',ns)
        compare('boundary_response_vs_CID',{i:{'rms':magnitude(unit(h(i))-unit(c0))} for i in ['ETV4_337_484_WT','ETV4_165_484_WT','ETV4_1_484_WT']})
    def baseline(seq):
        counts=collections.Counter(seq);return {'length':len(seq),'frequencies':{a:counts[a]/len(seq) for a in 'ACDEFGHIKLMNPQRSTVWY'},'charge_proxy':(counts['K']+counts['R']-counts['D']-counts['E'])/len(seq)}
    ref=inputs['ETV4_1_484_WT']['sequence']
    for ident,x in inputs.items():
        seq=x['sequence'];start=x['native_start'];a=339-start;b=421-start
        cmp(actual['baselines'][ident],{'whole':baseline(seq),'ETS_339_420':baseline(seq[a:b]),'N_flank':baseline(seq[:a]) if a else None,'C_flank':baseline(seq[b:]) if b<len(seq) else None,
            'native_position_identity':sum(seq[i]==ref[start-1+i] for i in range(len(seq)))/len(seq)},'/baselines/'+ident)
    rows=list(csv.DictReader((ROOT/'data/derived/measurement-registry.tsv').open(),delimiter='\t'))
    # Confirm the frozen source means independently of the hard-coded production expression.
    kd={r['construct_id']:float(r['value']) for r in rows if r['endpoint']=='equilibrium_KD'}
    if not kd:
        raise ValueError('Measurement field names require explicit reconciliation')
    expected_log=math.log(kd['ETV4_165_436_WT']*kd['ETV4_337_430_WT']/kd['ETV4_165_430_WT']/kd['ETV4_337_436_WT'])
    cmp(actual['biological_factorial_log_ratio'],expected_log,'/biological_factorial_log_ratio')
    # Every reported numerical leaf must have undergone separate arithmetic comparison.
    def numeric_paths(x,p=''):
        if isinstance(x,dict):
            for k,v in x.items():yield from numeric_paths(v,p+'/'+k)
        elif isinstance(x,list):
            for i,v in enumerate(x):yield from numeric_paths(v,p+'/'+str(i))
        elif isinstance(x,(int,float)) and not isinstance(x,bool):yield p
    missing=set(numeric_paths(actual))-checked_paths-{'/fixed_run_count'}
    assert not missing,sorted(missing)
    out={'passed':True,'actor':'Codex /root; separately implemented arithmetic, same agent','checked_at':datetime.datetime.now().astimezone().isoformat(),
      'numeric_comparisons':checks,'max_absolute_difference':max_abs,'comparison_tolerance':'2e-12 relative, 2e-14 absolute for float64 reconstruction only; scientific eta unchanged',
      'all_numeric_analysis_fields_checked':True,'source_analysis':str(Path(analysis).relative_to(ROOT)),'source_analysis_sha256':digest(analysis),'raw_arrays':len(records),
      'manifest_hashes':[{'path':str(Path(p).relative_to(ROOT)),'sha256':digest(p)} for p in manifests]}
    cov={'core_complete':True,'missing_core':sorted(expected_ids-set(records)),'completed_runs':len(records),'counts_by_model':dict(collections.Counter(r['model'] for r in records.values())),
      'run_plan_sha256':digest(ROOT/'design/run-plan.json'),'execution_config_sha256':cfg_hash,'runs':coverage}
    return out,cov
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifests',nargs='+',required=True);p.add_argument('--analysis',required=True);p.add_argument('--output',required=True);p.add_argument('--coverage',required=True);a=p.parse_args()
    out,cov=verify([str(Path(x).resolve()) for x in a.manifests],str(Path(a.analysis).resolve()))
    for path,obj in [(a.output,out),(a.coverage,cov)]:
        assert not Path(path).exists(),path
        write_json(path,obj)
    print(json.dumps(out,indent=2))
