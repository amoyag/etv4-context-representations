"""Versioned ESM-2 3B follow-up. No model executes without a Slurm allocation."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('E2_WORKSPACE', PACKAGE))
LAYERS = [19, 24, 36]

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def read(name):
    return json.loads((ROOT/name).read_text())

def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')

def native_indices(positions, start, length, token_start=0):
    positions = list(positions)
    assert len(set(positions)) == len(positions)
    offsets = [x-start for x in positions]
    assert offsets and min(offsets) >= 0 and max(offsets) < length
    return [x+token_start for x in offsets]

def verify_package():
    for rel, sha in json.loads((PACKAGE/'package-manifest.json').read_text()).items():
        assert digest(PACKAGE/rel) == sha, 'Package changed: '+rel

def require_e2():
    cfg = read('design/execution-config.json')
    assert cfg['status'] == 'frozen' and set(cfg['models']) == {'esm2_3b'}
    receipt = read('qualification-receipt.json')
    assert receipt['passed'] and digest(ROOT/'design/execution-config.json') == receipt['config_sha256']
    for rel, sha in receipt['input_hashes'].items():
        assert digest(ROOT/rel) == sha, 'Qualified input changed: '+rel

def initialize(workspace):
    global ROOT
    ROOT = Path(workspace).resolve()
    ROOT.mkdir(parents=True, exist_ok=False)
    for part in ['design', 'data']:
        shutil.copytree(PACKAGE/part, ROOT/part)
    (ROOT/'results/activations').mkdir(parents=True)
    shutil.copyfile(PACKAGE/'package-manifest.json', ROOT/'package-manifest.json')

def environment_check():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('GPU/model execution requires an explicitly submitted Slurm job')
    import importlib.metadata as md
    import torch
    assert torch.cuda.is_available(), 'CUDA runtime/driver unavailable'
    assert torch.cuda.device_count() == 1 and 'A100' in torch.cuda.get_device_name(0)
    installed={'torch': str(torch.__version__), 'fair-esm': md.version('fair-esm'), 'numpy': md.version('numpy')}
    expected=json.loads((PACKAGE/'design/preparation.json').read_text())['versions']
    assert installed == expected, (installed,expected)
    assert torch.version.cuda == '13.2', ('Unexpected CUDA build', torch.version.cuda)
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    torch.manual_seed(19001)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.use_deterministic_algorithms(True)
    torch.cuda.reset_peak_memory_stats()
    return dict(torch=torch.__version__, cuda=torch.version.cuda,
                gpu=torch.cuda.get_device_name(0), vram=torch.cuda.get_device_properties(0).total_memory,
                job_id=os.environ['SLURM_JOB_ID'])

class Adapter:
    def __init__(self):
        import torch, esm
        self.torch = torch
        assets=json.loads((PACKAGE/'checkpoint-manifest.json').read_text())
        for row in assets['files']:
            assert digest(PACKAGE/row['path']) == row['sha256'], row['path']
        with torch.serialization.safe_globals([argparse.Namespace]):
            self.model, self.alphabet=esm.pretrained.load_model_and_alphabet_local(str(PACKAGE/'models/esm2_t36_3B_UR50D.pt'))
        self.model.eval().requires_grad_(False).to(dtype=torch.float32)
        assert len(self.model.layers)==36 and self.model.embed_dim==2560
        assert not self.model.training and not any(p.requires_grad for p in self.model.parameters())
        self.device='cpu'
    def to(self, device):
        self.model.to(device=device, dtype=self.torch.float32)
        self.device=device
    def forward(self, sequence):
        import numpy as np
        t=self.torch
        assert sequence and not set(sequence)-set('ACDEFGHIKLMNPQRSTVWY')
        _,_,tokens=self.alphabet.get_batch_converter()([('target',sequence)])
        expected=[self.alphabet.cls_idx]+[self.alphabet.get_idx(a) for a in sequence]+[self.alphabet.eos_idx]
        assert tokens.tolist()[0]==expected
        tokens=tokens.to(self.device)
        captured={}; final={}; handles=[]
        def hook(layer):
            def capture(module, args, out):
                value=out[0] if isinstance(out,tuple) else out
                captured[layer]=value.transpose(0,1).detach().clone()
            return capture
        for layer in LAYERS:
            handles.append(self.model.layers[layer-1].register_forward_hook(hook(layer)))
        def capture_final(module,args):
            final['x']=args[0].detach().clone()
        handles.append(self.model.emb_layer_norm_after.register_forward_pre_hook(capture_final))
        if self.device=='cuda': t.cuda.synchronize()
        start=time.perf_counter()
        try:
            with t.inference_mode():
                out=self.model(tokens,repr_layers=LAYERS[:-1],return_contacts=False)
            # Both tensors must use (batch, tokens, features).
            assert t.equal(captured[36],final['x'].transpose(0,1))
            for layer in LAYERS[:-1]: assert t.equal(captured[layer],out['representations'][layer])
            if self.device=='cuda': t.cuda.synchronize()
            elapsed=time.perf_counter()-start
            arrays={str(k):v[0,1:len(sequence)+1].cpu().numpy().copy() for k,v in captured.items()}
            for x in arrays.values():
                assert x.shape==(len(sequence),2560) and x.dtype==np.float32 and np.isfinite(x).all()
                assert np.all(np.linalg.norm(x.astype(np.float64),axis=1)>1e-12)
        finally:
            for h in handles:h.remove()
        mapping=dict(token_start=1,target_length=len(sequence),tokens=len(expected),token_ids=expected,
                     position_ids=list(range(len(expected))),target_token_indices=list(range(1,len(sequence)+1)))
        return arrays,mapping,elapsed

def qualification(workspace):
    import numpy as np
    from analyse_results import max_discrepancy
    env=environment_check();verify_package();initialize(workspace)
    adapter=Adapter();cases={q['id']:q for q in read('design/qualification-inputs.json')['cases']}
    arrays={};records=[]
    for device in ['cpu','cuda']:
        adapter.to(device)
        for case in ['Q1','Q2','Q3']:
            arrays[device,case]=[]
            for repeat in range(3):
                a,m,seconds=adapter.forward(cases[case]['target'])
                path=ROOT/f'{device}-{case}-{repeat}.npz';np.savez(path,**a)
                records.append(dict(path=path.name,sha256=digest(path),seconds=seconds,mapping=m,device=device))
                arrays[device,case].append(a)
                write_json(ROOT/'progress.json',records)
                print(device,case,repeat,seconds,flush=True)
    extra={}
    for case in ['Q1','Q1_C']:
        a,m,seconds=adapter.forward(cases[case]['target']);extra[case]=a
        path=ROOT/f'cuda-{case}-extension.npz';np.savez(path,**a)
        records.append(dict(path=path.name,sha256=digest(path),seconds=seconds,mapping=m))
    eta={};discrepancies={}
    for k in map(str,LAYERS):
        repeat=max(max_discrepancy(a[k],b[k]) for reps in arrays.values() for a,b in itertools.combinations(reps,2))
        cross=max(max_discrepancy(a[k],b[k]) for q in ['Q1','Q2','Q3'] for a in arrays['cpu',q] for b in arrays['cuda',q])
        eta[k]=max(1e-7,10*max(repeat,cross))
        extension=max_discrepancy(extra['Q1'][k],extra['Q1_C'][k][:len(extra['Q1'][k])])
        assert extension>eta[k], 'Bidirectional extension below numerical floor'
        discrepancies[k]=dict(repeat=repeat,cpu_cuda=cross,extension=extension)
    config=dict(status='frozen',model='esm2_3b',models={'esm2_3b':dict(device='cuda',dtype='float32',width=2560,layers=LAYERS,eta=eta)},primary_block=24,
                batch_size=1,seed=19001,tf32=False,gradients=False,hardware=env,
                package_sha256=digest(PACKAGE/'package-manifest.json'),checkpoint_manifest_sha256=digest(PACKAGE/'checkpoint-manifest.json'))
    write_json(ROOT/'design/execution-config.json',config)
    import torch
    projected_e2_seconds_2x=217*2*max(x['seconds'] for x in records if x.get('device')=='cuda')
    assert projected_e2_seconds_2x < 3000, 'Measured envelope exceeds one-hour E2 request; revise resources explicitly'
    assert torch.cuda.max_memory_reserved() < 0.9*env['vram'], 'Insufficient GPU headroom'
    assert shutil.disk_usage(ROOT).free > 2*2484756480, 'Insufficient scratch headroom'
    input_hashes={str(p.relative_to(ROOT)):digest(p) for part in ['data','design'] for p in (ROOT/part).rglob('*.json') if p.name!='execution-config.json'}
    receipt=dict(passed=True,scope='artificial_only',forwards=len(records),eta=eta,discrepancies=discrepancies,
                 config_sha256=digest(ROOT/'design/execution-config.json'),input_hashes=input_hashes,
                 package_sha256=digest(PACKAGE/'package-manifest.json'),files=records,
                 peak_gpu_reserved_bytes=torch.cuda.max_memory_reserved(),projected_e2_seconds_2x=projected_e2_seconds_2x,hardware=env)
    write_json(ROOT/'qualification-receipt.json',receipt)
    print(json.dumps(receipt,indent=2),flush=True)

def experiment(workspace, qualified):
    import numpy as np
    env=environment_check();verify_package()
    q=Path(qualified).resolve();receipt=json.loads((q/'qualification-receipt.json').read_text())
    assert receipt['passed'] and receipt['package_sha256']==digest(PACKAGE/'package-manifest.json')
    assert digest(q/'design/execution-config.json')==receipt['config_sha256']
    for row in receipt['files']:assert digest(q/row['path'])==row['sha256']
    for rel,sha in receipt['input_hashes'].items():assert digest(PACKAGE/rel)==sha
    initialize(workspace)
    shutil.copyfile(q/'design/execution-config.json',ROOT/'design/execution-config.json')
    shutil.copyfile(q/'qualification-receipt.json',ROOT/'qualification-receipt.json')
    require_e2()
    config=read('design/execution-config.json')
    for key in ['torch','cuda','gpu','vram']:assert env[key]==config['hardware'][key], 'Hardware/runtime changed'
    adapter=Adapter();adapter.to('cuda');inputs=read('data/derived/model-inputs.json');rows=read('design/run-plan.json')['runs']
    records=[];dest=ROOT/'results/activations/esm2_3b';dest.mkdir()
    assert len(rows)==217
    for row in rows:
        x=inputs[row['input_id']]
        assert hashlib.sha256(x['sequence'].encode()).hexdigest()==row['input_sha256']
        assert row['layers']==LAYERS and row['model']=='esm2_3b' and row['prompt_id']=='empty'
        a,m,seconds=adapter.forward(x['sequence']);assert m['tokens']==row['expected_tokens']
        path=dest/(row['run_id']+'.npz');np.savez(path,**a)
        records.append({**row,'activation_file':str(path.relative_to(ROOT)),'sha256':digest(path),'seconds':seconds,
                        'mapping':m,'native_start':x['native_start'],'execution_config_sha256':digest(ROOT/'design/execution-config.json')})
        write_json(dest/'manifest.json',records)
        print(len(records),row['run_id'],seconds,flush=True)
    import torch
    write_json(ROOT/'results/inference-complete.json',dict(rows=len(records),hardware=env,peak_gpu_reserved_bytes=torch.cuda.max_memory_reserved(),
               status='inference_complete_analysis_pending',manifest_sha256=digest(dest/'manifest.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['qualify','e2']);p.add_argument('--workspace',required=True);p.add_argument('--qualified')
    args=p.parse_args()
    if args.stage=='qualify':qualification(args.workspace)
    else:
        if not args.qualified:p.error('--qualified required for e2')
        experiment(args.workspace,args.qualified)
