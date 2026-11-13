#!/usr/bin/env python3
"""Recompute E1 qualification from retained artificial arrays; no model is loaded."""
import argparse
from collections import Counter
from datetime import datetime
import itertools
from pathlib import Path
import shutil
from zoneinfo import ZoneInfo
import numpy as np
from run_experiments import ROOT, read, digest, write_json, LAYERS


def unit(x):
    # Deliberately independent arithmetic from analyse_results.normalized.
    x=np.asarray(x,np.float64)
    lengths=np.sqrt(np.einsum('ij,ij->i',x,x))
    if not np.isfinite(x).all() or np.any(lengths<=1e-12): raise ValueError('Invalid activation')
    return x/lengths[:,None]


def discrepancy(a,b):
    d=unit(a)-unit(b)
    return float(np.sqrt(np.einsum('ij,ij->i',d,d)).max())


def verify():
    reports=[read(str(p.relative_to(ROOT))) for p in sorted((ROOT/'results/qualification').glob('*/report.json'))]
    cases={q['id']:q for q in read('design/qualification-inputs.json')['cases']}
    counts=Counter(); group={}; diagnostics={}; checked_files=[]
    snapshots=[ROOT/'scripts/run_experiments.py',ROOT/'scripts/qualify_models.py']+list((ROOT/'handoff/E1-history').glob('source-*/*.py'))
    known_hashes={digest(p) for p in snapshots}
    for report in reports:
        assert report['biological_forwards']==0
        assert report['inputs_sha256']==digest(ROOT/'design/qualification-inputs.json')
        assert report['adapter_sha256'] in known_hashes, 'Missing exact adapter snapshot: '+report['run_id']
        assert report['qualifier_sha256'] in known_hashes, 'Missing exact qualifier snapshot: '+report['run_id']
        counts[report['model']]+=len(report['forwards'])
        if report['status']=='failed':
            assert not report['forwards'], 'Unreconciled partial failed qualification'
            continue
        assert report['status']=='measured'
        arrays={}; tokens_by_aa={}
        for row in report['forwards']:
            path=ROOT/row['path']; assert digest(path)==row['sha256']
            checked_files.append({'path':row['path'],'sha256':row['sha256']})
            q=cases[row['case']];m=row['mapping'];target=q['target']
            prefixes=q['prefixes'] if report['model']=='profam1' else []
            start=2+sum(len(s)+1 for s in prefixes) if report['model']=='profam1' else 1
            assert m['token_start']==start and m['tokens']==start+len(target)+1
            assert m['target_token_indices']==list(range(start,start+len(target)))
            assert m['position_ids']==list(range(m['tokens']))
            assert row['post_block_pre_norm_check'] and row['all_residue_token_identities_checked']
            # Check every residue in every prefix and target against one consistent token alphabet.
            offset=2 if report['model']=='profam1' else 1
            for seq in prefixes+[target]:
                for aa,token in zip(seq,m['token_ids'][offset:offset+len(seq)]):
                    assert tokens_by_aa.setdefault(aa,token)==token
                offset+=len(seq)+1
            assert len(set(tokens_by_aa.values()))==len(tokens_by_aa)
            with np.load(path,allow_pickle=False) as f:
                a={k:f[k] for k in f.files}
            assert set(a)==set(map(str,LAYERS[report['model']]))
            for x in a.values():
                assert x.dtype==np.float32 and x.shape==(len(target),report['feature_width'])
                unit(x)
            arrays.setdefault(row['case'],[]).append(a)
        record=(report,arrays)
        if report.get('diagnostics_only'):
            diagnostics[report['model']]=record
        else:
            assert (report['model'],report['device']) not in group
            group[report['model'],report['device']]=record
    assert counts['profam1']<=40 and counts['esm2_650m']<=30
    output={'passed':True,'scope':'Artificial E1 only','forwards':dict(counts),'biological_forwards':0,
            'models':{},'raw_array_count':len(checked_files),'raw_arrays':checked_files,
            'frozen_E0_sha256':digest(ROOT/'design/frozen-protocol.json'),
            'failed_attempts':[r['run_id'] for r in reports if r['status']=='failed']}
    plan=read('design/run-plan.json')
    for model,layers in LAYERS.items():
        cpu,ca=group[model,'cpu'];mps,ma=group[model,'mps']
        assert set(ca)==set(ma)==({'Q1','Q2','Q3','Q1_C','Q2_C'} if model=='profam1' else {'Q1','Q2','Q3'})
        assert all(len(v)==3 for v in list(ca.values())+list(ma.values()))
        eta={}; errors={}; causal={}
        for layer in map(str,layers):
            repeat=max(discrepancy(a[layer],b[layer]) for arrays in [ca,ma] for reps in arrays.values() for a,b in itertools.combinations(reps,2))
            cross=max(discrepancy(a[layer],b[layer]) for q in ca for a in ca[q] for b in ma[q])
            eta[layer]=max(1e-7,10*max(repeat,cross))
            errors[layer]={'repeat_max':repeat,'cpu_mps_max':cross,'eta':eta[layer]}
            if model=='profam1':
                causal[layer]={device:max(discrepancy(arrays[q][i][layer],arrays[q+'_C'][j][layer][:len(arrays[q][i][layer])]) for q in ['Q1','Q2'] for i in range(3) for j in range(3)) for device,arrays in [('cpu',ca),('mps',ma)]}
                assert all(v<=eta[layer] for v in causal[layer].values()), 'Causal failure under independently calibrated floor'
        # Select only after comparing measured reference and accelerated run cost.
        cpu_long=max(r['end_to_end_seconds'] for r in cpu['forwards'] if r['case'] in ['Q2','Q3'])
        mps_long=max(r['end_to_end_seconds'] for r in mps['forwards'] if r['case'] in ['Q2','Q3'])
        device='mps' if mps_long<cpu_long else 'cpu'
        selected=mps if device=='mps' else cpu
        # Bound every planned pass by the slowest observed artificial pass, including cold/IO costs.
        per_pass=max(r['end_to_end_seconds'] for r in selected['forwards'])
        seconds=2*(selected['load_seconds']+per_pass*plan['run_counts'][model])
        diagnostic,da=diagnostics[model]
        assert diagnostic['device']==device
        if model=='esm2_650m':
            visibility={k:discrepancy(da['Q1'][0][k],da['Q1_C'][0][k][:100]) for k in eta}
            assert all(v>8*eta[k] for k,v in visibility.items()), 'ESM did not see C extension above error floor'
        else: visibility=causal
        rss=max(r['peak_process_rss_bytes'] for r in [cpu,mps,diagnostic])
        driver=max(r.get('mps_driver_allocated_bytes',0) for report in [mps,diagnostic] for r in report['forwards'])
        assert driver>0, 'MPS allocations not measured'
        memory_bound=rss+driver
        assert memory_bound<12*2**30, 'Insufficient 16-GiB host headroom'
        storage=sum((r['expected_tokens']-2 if model=='esm2_650m' else
                     read('data/derived/model-inputs.json')[r['input_id']]['length'])*selected['feature_width']*3*4
                    for r in plan['runs'] if r['model']==model)
        output['models'][model]={'device':device,'dtype':'float32','width':selected['feature_width'],'layers':layers,
            'eta':eta,'errors':errors,'visibility':visibility,'parameter_count':selected['parameter_count'],
            'max_measured_end_to_end_seconds':per_pass,'projected_seconds_2x':seconds,
            'RSS_peak_bytes':rss,'MPS_driver_max_sample_bytes':driver,'conservative_RSS_plus_driver_bytes':memory_bound,
            'activation_bytes_uncompressed':storage,'cpu_long_max_seconds':cpu_long,'mps_long_max_seconds':mps_long,
            'selected_report':str((ROOT/'results/qualification'/selected['run_id']/'report.json').relative_to(ROOT))}
    now=datetime.now(ZoneInfo('Europe/Madrid'));deadline=datetime(2026,11,12,16,0,tzinfo=ZoneInfo('Europe/Madrid'))
    cost=sum(m['projected_seconds_2x'] for m in output['models'].values())
    # Reserve one hour for E2 arithmetic/figures/integrity work; this is a planning allowance, not a measured runtime.
    output['resources']={'estimated_at':now.isoformat(),'export_freeze':deadline.isoformat(),'model_seconds_2x':cost,
        'analysis_planning_allowance_seconds':3600,'projected_core_seconds_including_allowance':cost+3600,
        'available_seconds':(deadline-now).total_seconds(),'fits_before_freeze':cost+3600<(deadline-now).total_seconds(),
        'disk_free_bytes':shutil.disk_usage(ROOT).free,'activation_storage_bytes_2x':2*sum(m['activation_bytes_uncompressed'] for m in output['models'].values()),
        'assumption':'Prompt E2 activation after author review and otherwise available local host; timing does not schedule or authorize E2.'}
    assert output['resources']['fits_before_freeze']
    assert output['resources']['activation_storage_bytes_2x']<output['resources']['disk_free_bytes']
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output');args=parser.parse_args()
    result=verify()
    if args.output: write_json(ROOT/args.output,result)
    print(__import__('json').dumps({k:v for k,v in result.items() if k!='raw_arrays'},indent=2))
