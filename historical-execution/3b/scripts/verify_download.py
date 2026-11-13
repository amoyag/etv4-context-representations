"""Model-free verification after downloading an E2 workspace to the Mac."""
import argparse
import json
import os
from pathlib import Path
import hashlib
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

p=argparse.ArgumentParser();p.add_argument('workspace');args=p.parse_args();root=Path(args.workspace).resolve()
plan=json.loads((root/'design/run-plan.json').read_text())
cfg=json.loads((root/'design/execution-config.json').read_text())
receipt=json.loads((root/'qualification-receipt.json').read_text())
assert sha(root/'design/execution-config.json')==receipt['config_sha256'] and receipt['passed']
for rel,h in receipt['input_hashes'].items():assert sha(root/rel)==h
rows=json.loads((root/'results/activations/esm2_3b/manifest.json').read_text())
assert len(rows)==217 and len({r['run_id'] for r in rows})==217
assert {r['run_id'] for r in rows}=={r['run_id'] for r in plan['runs']}
inputs=json.loads((root/'data/derived/model-inputs.json').read_text())
for r in rows:
    x=inputs[r['input_id']];path=(root/r['activation_file']).resolve()
    assert path.is_relative_to(root/'results/activations') and sha(path)==r['sha256']
    assert r['mapping']['tokens']==x['length']+2
    assert r['mapping']['target_token_indices']==list(range(1,x['length']+1))
    assert r['mapping']['position_ids']==list(range(x['length']+2))
    assert r['execution_config_sha256']==receipt['config_sha256']
    with np.load(path,allow_pickle=False) as data:
        assert set(data.files)=={'19','24','36'}
        for a in data.values():
            assert a.dtype==np.float32 and a.shape==(x['length'],2560) and np.isfinite(a).all()
            assert np.all(np.sqrt(np.einsum('ij,ij->i',a.astype('float64'),a.astype('float64')))>1e-12)
print(json.dumps({'passed':True,'rows':len(rows),'scope':'download integrity, dtype, shape, coordinates; not independent scientific review'}))
