#!/usr/bin/env python3
"""Recompute the fixed ETV4 summaries from archived NPZ arrays (CPU, NumPy only).

Adapted from historical-execution/original/scripts/analyse_results.py and the
3B variant. Scientific calculations are unchanged; file access and command-line
handling are portable. Original execution guards remain intact in the archive.
"""
import argparse
import json
import numpy as np
from pathlib import Path
import hashlib


def normalized(h):
    h = np.asarray(h, dtype=np.float64)
    if h.ndim != 2 or not h.size or not np.isfinite(h).all():
        raise ValueError('Expected finite nonempty residue-by-feature array')
    norm = np.linalg.norm(h, axis=-1, keepdims=True)
    if np.any(norm <= 1e-12):
        raise ValueError('Uninterpretable zero-norm activation')
    return h / norm


def rms(v):
    v = np.asarray(v, dtype=np.float64)
    if v.ndim != 2 or not v.size or not np.isfinite(v).all():
        raise ValueError('Expected finite nonempty residue-by-feature array')
    return float(np.sqrt(np.mean(np.sum(v * v, axis=-1))))


def max_discrepancy(a, b):
    if a.shape != b.shape:
        raise ValueError('Unaligned arrays')
    return float(np.max(np.linalg.norm(normalized(a) - normalized(b), axis=1)))


def context_metrics(fl, cid, wt_fl, wt_cid, eta):
    if eta < 1e-7 or len({a.shape for a in (fl, cid, wt_fl, wt_cid)}) != 1:
        raise ValueError('Invalid numerical floor or unmatched coordinates')
    q, qw = normalized(fl) - normalized(cid), normalized(wt_fl) - normalized(wt_cid)
    M, B, D = rms(q - qw), rms(qw), rms(normalized(cid) - normalized(wt_cid))
    return {'M': M, 'B': B, 'D': D, 'T': M/B if B > 4*eta else None,
            'M_over_D': M/D if D > 4*eta else None,
            'L': max(0., M-8*eta)/(D+4*eta) if D > 4*eta else None,
            'U': (M+8*eta)/(D-4*eta) if D > 4*eta else None}


def specificity(metrics, shuffled_M, eta):
    Y = 'Y401N_Y403G'
    controls = ['E404K', 'E423K', 'E425K', 'F428A_F432A']
    expected = {Y, *controls}
    if set(metrics) != expected or set(shuffled_M) != expected or any(len(v) != 20 for v in shuffled_M.values()):
        raise ValueError('Complete fixed genotype and 20-shuffle panel required')
    resolved = all(m['B'] > 4*eta and m['D'] > 4*eta for m in metrics.values())
    conditions = {'denominators_resolved': resolved,
                  'M_separation': metrics[Y]['M'] > max(metrics[g]['M'] for g in controls) + 8*eta,
                  'ratio_bounds_separation': resolved and metrics[Y]['L'] > max(metrics[g]['U'] for g in controls),
                  'native_order_specificity': metrics[Y]['M'] > max(shuffled_M[Y]) + 8*eta}
    label = ('numerically_unresolved' if not resolved else
             'source_consistent_native_order_specificity' if all(conditions.values()) else 'specificity_not_established')
    return {'classification': label, 'conditions': conditions,
            'shuffle_counts_at_or_above_native': {g: sum(v >= metrics[g]['M'] for v in shuffled_M[g]) for g in metrics}}


def factorial(zero, n, c, nc):
    if len({a.shape for a in (zero, n, c, nc)}) != 1:
        raise ValueError('Unaligned factorial states')
    j = normalized(nc) - normalized(n) - normalized(c) + normalized(zero)
    return {'rms': rms(j), 'profile': np.linalg.norm(j, axis=1).tolist()}


def prompt_response(n, c):
    if n.shape != c.shape:
        raise ValueError('Unaligned prompt states')
    return normalized(n) - normalized(c)


def sequence_baseline(seq):
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    if not seq or set(seq) - set(aa):
        raise ValueError('Invalid sequence')
    return {'length': len(seq), 'frequencies': {a: seq.count(a)/len(seq) for a in aa},
            'charge_proxy': sum(seq.count(a)*s for a, s in [('K',1),('R',1),('D',-1),('E',-1)])/len(seq)}


def analyse_export(group, arrays_root, destination):
    """Reanalyse archived arrays without importing an inference adapter."""
    root = Path(__file__).resolve().parents[1]
    source = root / 'historical-execution' / group
    read = lambda name: json.loads((source / name).read_text())
    digest = lambda path: hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError('Output already exists; choose a new output path')
    plan, config = read('design/run-plan.json'), read('design/execution-config.json')
    names = ['650m', 'profam'] if group == 'original' else ['3b']
    manifest_paths = [root / 'provenance/activation-manifests' / (n + '.json') for n in names]
    model_folder = {'esm2_650m':'650m', 'profam1':'profam', 'esm2_3b':'3b'}
    def array_path(record):
        return Path(arrays_root) / model_folder[record['model']] / Path(record['activation_file']).name
    records = {}
    for manifest in manifest_paths:
        for row in json.loads(Path(manifest).read_text()):
            if row['run_id'] in records: raise ValueError('Duplicate run ID')
            records[row['run_id']] = row
    if set(records) != {r['run_id'] for r in plan['runs']}:
        raise ValueError('Complete fixed run matrix required')
    cfg_hash = digest(source / 'design/execution-config.json')
    for row in plan['runs']:
        record = records[row['run_id']]
        for key in ['model','input_id','prompt_id','input_sha256','expected_tokens','layers']:
            if record[key] != row[key]: raise ValueError('Frozen plan mismatch: '+key)
        path = array_path(record)
        if digest(path) != record['sha256']:
            raise ValueError('Activation identity mismatch')
        if record['execution_config_sha256'] != cfg_hash:
            raise ValueError('Execution configuration mismatch')
    inputs = read('data/derived/model-inputs.json')
    prompts = read('data/derived/prompt-panel.json')
    def h(model, ident, prompt, layer, mask):
        row = records[f'{model}__{ident}__{prompt}']
        source = inputs[ident]
        if row['native_start'] != source['native_start']: raise ValueError('Native origin mismatch')
        ix = native_indices(mask, source['native_start'], source['length'])
        with np.load(array_path(row), allow_pickle=False) as arrays:
            a = arrays[str(layer)]
            if a.shape != (source['length'],config['models'][model]['width']): raise ValueError('Activation shape mismatch')
            return a[ix].astype(np.float64)
    genotypes=['Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']
    result={'execution_config_sha256':cfg_hash,'fixed_run_count':len(records),'models':{},'baselines':{}}
    for model, spec in config['models'].items():
        result['models'][model]={}
        for layer in spec['layers']:
            eta=spec['eta'][str(layer)]; mask=plan['masks']['genotype_primary']
            def cid(g): return h(model,f'ETV4_337_436_{g}','empty',layer,mask)
            def fl(g): return h(model,f'ETV4_1_484_{g}','empty',layer,mask)
            wt_fl,wt_cid=fl('WT'),cid('WT')
            native={g:context_metrics(fl(g),cid(g),wt_fl,wt_cid,eta) for g in genotypes}
            block={'genotype_native':native,'eta':eta}
            direct_positions=[401,403,404,423,425,428,432]
            block['substitution_site_diagnostics']={}
            dwf=h(model,'ETV4_1_484_WT','empty',layer,direct_positions)
            dwc=h(model,'ETV4_337_436_WT','empty',layer,direct_positions)
            for g in genotypes:
                df=h(model,f'ETV4_1_484_{g}','empty',layer,direct_positions)
                dc=h(model,f'ETV4_337_436_{g}','empty',layer,direct_positions)
                delta=prompt_response(df,dc)-prompt_response(dwf,dwc)
                block['substitution_site_diagnostics'][g]={'native_positions':direct_positions,
                    'context_change_norm':np.linalg.norm(delta,axis=1).tolist(),
                    'CID_change_norm':np.linalg.norm(normalized(dc)-normalized(dwc),axis=1).tolist(),
                    'scope':'Separate descriptive sites; not independent biological replicates.'}
            if model in ('esm2_650m', 'esm2_3b'):
                shuffled={g:[] for g in genotypes}
                for seed in range(1,21):
                    wf=h(model,f'FLshuffle_{seed}_WT','empty',layer,mask)
                    for g in genotypes:
                        gf=h(model,f'FLshuffle_{seed}_{g}','empty',layer,mask)
                        shuffled[g].append(context_metrics(gf,cid(g),wf,wt_cid,eta))
                block['genotype_shuffled']=shuffled
                block['interpretation']=specificity(native,{g:[x['M'] for x in v] for g,v in shuffled.items()},eta)
            else:
                block['genotype_scope']='Secondary causal partial-context diagnostic; cannot overturn ESM primary comparison.'
            mask=plan['masks']['WT_primary']
            def wt(ident,prompt='empty'): return h(model,ident,prompt,layer,mask)
            zero,n,c,nc=[wt(f'ETV4_{a}_{b}_WT') for a,b in [(337,430),(165,430),(337,436),(165,436)]]
            block['factorial_native']=factorial(zero,n,c,nc)
            block['native_N_response_with_C']={'rms':rms(prompt_response(nc,c)),
                'profile':np.linalg.norm(prompt_response(nc,c),axis=1).tolist()}
            if model in ('esm2_650m', 'esm2_3b'):
                block['factorial_shuffled']=[{'seed':s,**factorial(zero,wt(f'Nshuffle_{s}_430'),wt(f'Cshuffle_{s}'),wt(f'NCshuffle_{s}'))} for s in range(1,21)]
                block['N_only_shuffled_response']=[{'seed':s,'rms':rms(prompt_response(wt(f'Nshuffle_{s}_436'),c))} for s in range(1,21)]
            else:
                responses={}
                block['prompt_response']={}
                for prompt in prompts:
                    pn,pc=wt('ETV4_165_436_WT',prompt),wt('ETV4_337_436_WT',prompt)
                    q=prompt_response(pn,pc);responses[prompt]=q
                    zn,zc=wt('ETV4_165_430_WT',prompt),wt('ETV4_337_430_WT',prompt)
                    visibility=max(max_discrepancy(pn,zn),max_discrepancy(pc,zc))
                    if visibility>eta: raise ValueError('E2 ProFam causal-invariance failure: '+prompt)
                    block['prompt_response'][prompt]={'rms':rms(q),'profile':np.linalg.norm(q,axis=1).tolist(),'causal_extension_error':visibility}
                block['prompt_contrasts']=[]
                for j in range(1,6):
                    for order in ['ab','ba']:
                        hom=f'hom_{j}_{order}';non=f'nonhom_{j}_{order}';sh=f'homshuffle_{j}_{order}'
                        block['prompt_contrasts'].append({'set':j,'order':order,'hom_minus_nonhom':rms(responses[hom]-responses[non]),'hom_minus_homshuffle':rms(responses[hom]-responses[sh])})
                block['against_empty']={p:rms(q-responses['empty']) for p,q in responses.items()}
                block['order_contrasts']={f'{group}_{j}':rms(responses[f'{group}_{j}_ab']-responses[f'{group}_{j}_ba']) for group in ['hom','nonhom','homshuffle'] for j in range(1,6)}
                block['N_shuffles']={}
                for p in ['empty']+[f'{group}_1_{order}' for group in ['hom','nonhom','homshuffle'] for order in ['ab','ba']]:
                    pc=wt('ETV4_337_436_WT',p); values=[]
                    for seed in range(1,21):
                        q=prompt_response(wt(f'Nshuffle_{seed}_436',p),pc)
                        values.append({'seed':seed,'rms':rms(q),'profile':np.linalg.norm(q,axis=1).tolist()})
                    block['N_shuffles'][p]=values
            # Prespecified boundary effects relative to CID, shared 337–430 positions.
            block['boundary_response_vs_CID']={ident:{'rms':rms(prompt_response(wt(ident),c))} for ident in ['ETV4_337_484_WT','ETV4_165_484_WT','ETV4_1_484_WT']}
            result['models'][model][str(layer)]=block
    reference=inputs['ETV4_1_484_WT']['sequence']
    for ident,x in inputs.items():
        seq,start=x['sequence'],x['native_start']; domain=native_indices(range(339,421),start,len(seq))
        a,b=domain[0],domain[-1]+1
        result['baselines'][ident]={'whole':sequence_baseline(seq),'ETS_339_420':sequence_baseline(seq[a:b]),
            'N_flank':sequence_baseline(seq[:a]) if a else None,'C_flank':sequence_baseline(seq[b:]) if b<len(seq) else None,
            'native_position_identity':sum(u==v for u,v in zip(seq,reference[start-1:start-1+len(seq)]))/len(seq)}
    result['biological_factorial_log_ratio']=float(np.log(60*6.1/(12*12)))
    result['biological_factorial_source']='Frozen S3 means; no covariance/error bars or calibration to representation norms.'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    return result


def native_indices(native_positions, native_start, target_length):
    positions = list(native_positions)
    if len(set(positions)) != len(positions):
        raise ValueError('Duplicate native coordinate')
    offsets = [i - native_start for i in positions]
    if any(i < 0 or i >= target_length for i in offsets):
        raise ValueError('Native position outside target')
    return offsets


def compare(actual, expected, path='root'):
    """Check every saved field, including profiles, labels and negative findings.

    The tolerance covers float64 reduction/library differences only. It is not
    a changed scientific threshold and does not replace any protocol eta.
    """
    if isinstance(expected, dict):
        if actual.keys() != expected.keys(): raise ValueError('Keys differ: '+path)
        return sum(compare(actual[k], v, path+'/'+k) for k,v in expected.items())
    if isinstance(expected, list):
        if len(actual) != len(expected): raise ValueError('Length differs: '+path)
        return sum(compare(a,b,path+'/'+str(i)) for i,(a,b) in enumerate(zip(actual,expected)))
    if isinstance(expected, (float,int)) and not isinstance(expected,bool):
        if not np.isclose(actual, expected, rtol=1e-12, atol=1e-14):
            raise ValueError(f'Numeric mismatch at {path}: {actual} != {expected}')
        return 1
    if actual != expected: raise ValueError('Value differs: '+path)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arrays', required=True, type=Path,
                        help='Directory containing extracted 650m/, profam/ and 3b/ folders')
    parser.add_argument('--output', required=True, type=Path, help='New output directory')
    args = parser.parse_args()
    if args.output.exists(): parser.error('Choose a new output directory')
    root = Path(__file__).resolve().parents[1]
    refs = {
        'original': 'repro/inputs/exports/results-v1/results/E2-analysis-v1.json',
        '3b': 'repro/inputs/followups/esm2-3b-20261111/job-3647838/exports/results-v1/e2/results/E2-3B-analysis.json'}
    report = {'model_inference':False, 'comparison_rtol':1e-12, 'comparison_atol':1e-14, 'groups':{}}
    for group, reference in refs.items():
        result = analyse_export(group,args.arrays,args.output/(group+'-analysis.json'))
        count = compare(result,json.loads((root/reference).read_text()))
        report['groups'][group] = {'numeric_values_checked':count,'all_fields_match':True,'activation_hashes_verified':True}
        print(group, ': all fields match;', count, 'numeric values checked',flush=True)
    (args.output/'verification.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__ == '__main__':
    main()
