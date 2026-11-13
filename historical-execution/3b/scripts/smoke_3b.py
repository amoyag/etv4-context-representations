"""One artificial GPU forward; no biological input or qualification receipt."""
import json, os, time
from pathlib import Path
import run_experiments as r
import torch
env=r.environment_check();r.verify_package()
a=r.Adapter();a.to('cuda')
case=json.loads((r.PACKAGE/'design/qualification-inputs.json').read_text())['cases'][0]
arrays,mapping,seconds=a.forward(case['target'])
torch.cuda.synchronize()
result=dict(passed=True,scope='one artificial Q1 forward only',hardware=env,layers=sorted(arrays),shapes={k:list(v.shape) for k,v in arrays.items()},seconds=seconds,peak_allocated=torch.cuda.max_memory_allocated(),peak_reserved=torch.cuda.max_memory_reserved(),biological_forwards=0,qualification_complete=False)
print(json.dumps(result,indent=2),flush=True)
