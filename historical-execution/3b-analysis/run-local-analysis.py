import os,subprocess,json,sys
from pathlib import Path
project=Path('/Users/aurelio/KnowledgeSystem/areas/research/ministerio-PID2026')
base=project/'pilot/ets-regulatory-context/followups/esm2-3b-20261111'
root=base/'job-3647838/e2'
pkg=project/'infrastructure/picasso/esm2-3b-20261111-v3/scripts'
post=base/'postprocess'
env={**os.environ,'E2_WORKSPACE':str(root),'PYTHONPATH':str(pkg),'MPLCONFIGDIR':'/private/tmp/pid2026-mpl-3647838'}
py='/private/tmp/pid2026-3b-analysis-3647838/bin/python'
man=str(root/'results/activations/esm2_3b/manifest.json');analysis=str(root/'results/E2-3B-analysis.json')
logs=base/'job-3647838/local-analysis';logs.mkdir(exist_ok=True)
commands=[('download-verification',[str(pkg/'verify_download.py'),str(root)]),('analysis',[str(pkg/'analyse_results.py'),'--manifests',man,'--output',analysis]),('arithmetic-verification',[str(post/'verify_e2.py'),'--manifests',man,'--analysis',analysis,'--output',str(root/'results/arithmetic-verification.json'),'--coverage',str(root/'results/coverage.json')]),('presentation',[str(post/'package_e2_results.py'),'--analysis',analysis])]
for name,args in commands:
    print('Starting',name,flush=True)
    with (logs/(name+'.log')).open('x') as out:
        result=subprocess.run([py,*args],env=env,stdout=out,stderr=subprocess.STDOUT)
    if result.returncode:
        print((logs/(name+'.log')).read_text());sys.exit(result.returncode)
    print('Passed',name,flush=True)
print(json.dumps({k:v['interpretation'] for k,v in json.loads(Path(analysis).read_text())['models']['esm2_3b'].items()},indent=2))
