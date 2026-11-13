"""Read-only E0/E1 integrity and host inspection; no model load."""
import datetime, importlib.metadata, json, platform, shutil, subprocess
from pathlib import Path
import torch
from run_experiments import ROOT, read, digest, write_json
from check_workspace import check_stage

def main():
    check_stage(ROOT,'E1',read('design/stages.json')['stages'])
    config=read('design/execution-config.json'); manifest=read('environment/model-manifest.json')
    checked=[]
    for x in config['frozen_files']:
        assert digest(ROOT/x['path'])==x['sha256'],x['path']
        checked.append(x)
    versions={k:importlib.metadata.version(k) for k in manifest['versions']}
    assert versions==manifest['versions']
    assert platform.python_version()==manifest['python']
    hardware=subprocess.check_output(['sysctl','-n','hw.model','hw.memsize'],text=True).splitlines()
    assert hardware==['Macmini9,1','17179869184']
    assert torch.backends.mps.is_available()
    out={'checked_at':datetime.datetime.now().astimezone().isoformat(),'host':platform.node(),'hardware':hardware,
      'python':platform.python_version(),'versions':versions,'mps_available':True,'mps_recommended_bytes':torch.mps.recommended_max_memory(),
      'disk_free_bytes':shutil.disk_usage(ROOT).free,'checked_frozen_files':checked,'model_forwards':0,
      'receipts':[{'path':f'handoff/{s}.json','sha256':digest(ROOT/f'handoff/{s}.json')} for s in ['E0','E1']]}
    path=ROOT/'logs/E2-host-preflight.json'
    assert not path.exists()
    write_json(path,out)
    print('PASS: E0/E1 receipts, frozen code/weights, exact versions and Macmini MPS; zero model forwards.')
if __name__=='__main__': main()
