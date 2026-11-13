#!/usr/bin/env python3
"""E1 only: frozen artificial input tests, immutable raw arrays and timing logs."""
import argparse
from datetime import datetime, timezone
import fcntl
import importlib.metadata
import itertools
import json
from pathlib import Path
import platform
import resource
import time
import traceback
import numpy as np
from run_experiments import Adapter, ROOT, LAYERS, read, digest, write_json, native_indices
from analyse_results import max_discrepancy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=LAYERS, required=True)
    parser.add_argument('--device', choices=['cpu', 'mps'], default='cpu')
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--diagnostics', action='store_true', help='One Q3 memory pass for ProFam; Q1/Q1_C visibility pair for ESM')
    args = parser.parse_args()
    state = read('WORK_STATE.json')
    if state.get('activation_scope') != ['E1'] or not state.get('active_execution'):
        raise RuntimeError('Standalone E1 is not active')
    if Path(args.run_id).name != args.run_id or args.run_id in ('.', '..'):
        raise ValueError('Invalid run ID')
    lock = (ROOT / '.model-execution.lock.tmp').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    dest = ROOT / 'results/qualification' / args.run_id
    dest.mkdir(parents=True, exist_ok=False)
    report = {'run_id': args.run_id, 'model': args.model, 'device': args.device,
              'started_at': datetime.now(timezone.utc).isoformat(), 'status': 'running',
              'inputs_sha256': digest(ROOT / 'design/qualification-inputs.json'),
              'adapter_sha256': digest(ROOT / 'scripts/run_experiments.py'),
              'qualifier_sha256': digest(Path(__file__)), 'forwards': [], 'biological_forwards': 0,
              'host': platform.node(), 'platform': platform.platform(),
              'versions': {p: importlib.metadata.version(p) for p in ['torch','profam','fair-esm','transformers','numpy']},
              'python': platform.python_version(), 'dtype': 'float32', 'threads': 4, 'interop_threads': 1,
              'warmups': 0, 'batch_size': 1, 'storage_dtype': 'float32', 'summaries_dtype': 'float64'}
    report['diagnostics_only'] = args.diagnostics
    def save():
        report['peak_process_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if platform.system() == 'Darwin' else 1024)
        write_json(dest / 'report.json', report)
    save()
    try:
        adapter = Adapter(args.model, args.device)
        report.update(load_seconds=adapter.load_seconds, parameter_count=adapter.parameter_count,
                      feature_width=adapter.width, asset_checks=adapter.asset_checks,
                      state_dict_audit=getattr(adapter, 'state_dict_audit', 'fair-esm strict core load'),
                      mps_available=adapter.torch.backends.mps.is_available())
        if args.model == 'profam1':
            report['rope_type'] = adapter.model.model.model.rotary_emb.rope_type
            report['rope_scaling'] = adapter.model.model.config.rope_scaling
        cases = read('design/qualification-inputs.json')['cases']
        if args.diagnostics:
            wanted = ['Q3'] if args.model == 'profam1' else ['Q1','Q1_C']
            cases = [q for q in cases if q['id'] in wanted]
        elif args.model == 'esm2_650m':
            cases = [q for q in cases if not q['id'].endswith('_C')]
        arrays = {}
        for q in cases:
            prefixes = q['prefixes'] if args.model == 'profam1' else []
            arrays[q['id']] = []
            for repeat in range(1 if args.diagnostics else 3):
                print(f"START {args.model} {args.device} {q['id']} repeat={repeat}", flush=True)
                began = time.perf_counter()
                h, mapping, seconds = adapter.forward(q['target'], prefixes)
                path = dest / f"{q['id']}-{repeat}.npz"
                np.savez(path, **h)
                elapsed = time.perf_counter() - began
                assert native_indices(range(337,337+len(q['target'])),337,len(q['target']),mapping['token_start']) == mapping['target_token_indices']
                row = {'case': q['id'], 'repeat': repeat, 'forward_seconds': seconds, 'end_to_end_seconds': elapsed,
                       'path': str(path.relative_to(ROOT)), 'sha256': digest(path), 'mapping': mapping,
                       'post_block_pre_norm_check': True, 'all_residue_token_identities_checked': True}
                if args.device == 'mps':
                    row['mps_current_allocated_bytes'] = adapter.torch.mps.current_allocated_memory()
                    row['mps_driver_allocated_bytes'] = adapter.torch.mps.driver_allocated_memory()
                report['forwards'].append(row)
                arrays[q['id']].append(h)
                save()
                print(f"DONE {q['id']} repeat={repeat} {seconds:.3f}s peak_RSS={report['peak_process_rss_bytes']/2**30:.3f}GiB", flush=True)
        if args.diagnostics:
            if args.model == 'esm2_650m':
                report['bidirectional_extension_effect'] = {layer:max_discrepancy(arrays['Q1'][0][layer],arrays['Q1_C'][0][layer][:100]) for layer in map(str,adapter.layers)}
            report['status'] = 'measured'
            report['finished_at'] = datetime.now(timezone.utc).isoformat()
            save()
            print(json.dumps({k:v for k,v in report.items() if k in ['status','bidirectional_extension_effect','peak_process_rss_bytes']},indent=2),flush=True)
            return
        repeatability = {}
        for layer in map(str, adapter.layers):
            d = max(max_discrepancy(a[layer], b[layer]) for reps in arrays.values() for a,b in itertools.combinations(reps,2))
            repeatability[layer] = {'max_pairwise_unit_vector_l2': d, 'eta_repeat_only': max(1e-7,10*d)}
        report['repeatability'] = repeatability
        if args.model == 'profam1':
            report['causal_extension'] = {q: {layer: max_discrepancy(arrays[q][0][layer],arrays[q+'_C'][0][layer][:len(arrays[q][0][layer])]) for layer in repeatability} for q in ['Q1','Q2']}
            report['causal_pass_repeat_only'] = all(v <= repeatability[k]['eta_repeat_only'] for row in report['causal_extension'].values() for k,v in row.items())
        report['status'] = 'measured'
        report['finished_at'] = datetime.now(timezone.utc).isoformat()
        save()
        print(json.dumps({k:v for k,v in report.items() if k in ['status','repeatability','causal_extension','causal_pass_repeat_only','load_seconds','parameter_count','peak_process_rss_bytes']}, indent=2), flush=True)
    except BaseException as exc:
        report['status'] = 'failed'
        report['error'] = repr(exc)
        report['traceback'] = traceback.format_exc()
        save()
        raise


if __name__ == '__main__':
    main()
