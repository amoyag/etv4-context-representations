"""Execute the frozen run_experiments entry point with resource-only instrumentation.
No change to model inputs, inference, layers, precision or analysis thresholds.
"""
import datetime, json, resource, runpy, sys, time, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from check_stage_transition import check_transition
from run_experiments import digest, write_json

def main():
    check_transition(ROOT,'E2')
    run_id=sys.argv[sys.argv.index('--run-id')+1]
    report_path=ROOT/'logs'/f'{run_id}-resources.json'
    if report_path.exists(): raise ValueError('Resource report already exists')
    code=ROOT/'scripts/run_experiments.py'
    ns=runpy.run_path(str(code),run_name='frozen_e2_execution')
    cls=ns['Adapter']; original=cls.forward
    report={'run_id':run_id,'started_at':datetime.datetime.now().astimezone().isoformat(),
      'runner_sha256':digest(code),'monitor_sha256':digest(__file__),'samples':[],
      'memory_scope':'Process peak RSS plus synchronized post-forward MPS allocation samples; sum may double-count unified memory and is not a system-wide peak.'}
    start=time.perf_counter()
    def observed(self,target,prefixes=()):
        t=time.perf_counter()
        arrays,mapping,seconds=original(self,target,prefixes)
        torch=self.torch
        report['load_seconds']=self.load_seconds
        report['model']=self.model_id;report['device']=self.device
        report['samples'].append({'index':len(report['samples'])+1,'tokens':mapping['tokens'],
          'forward_seconds':seconds,'adapter_call_seconds':time.perf_counter()-t,
          'rss_peak_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          'mps_current_bytes':torch.mps.current_allocated_memory(),'mps_driver_bytes':torch.mps.driver_allocated_memory(),
          'finite_shapes_and_coordinates_checked':True})
        write_json(report_path,report)
        return arrays,mapping,seconds
    cls.forward=observed
    try:
        ns['main']()
        report['status']='complete'
    except BaseException:
        report['status']='failed';report['error']=traceback.format_exc()
        raise
    finally:
        report['ended_at']=datetime.datetime.now().astimezone().isoformat()
        report['wall_seconds']=time.perf_counter()-start
        report['rss_peak_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        report['completed_forwards']=len(report['samples'])
        write_json(report_path,report)
if __name__=='__main__': main()
