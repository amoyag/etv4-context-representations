"""Supplement full-target causal and repeated-input checks without new forwards."""
import collections,json
from pathlib import Path
import numpy as np
from run_experiments import ROOT,read,digest,write_json
from verify_e2 import unit,profile

def main():
    config=read('design/execution-config.json');inputs=read('data/derived/model-inputs.json');prompts=read('data/derived/prompt-panel.json')
    records={r['run_id']:r for p in (ROOT/'results/activations').glob('*/manifest.json') for r in json.loads(p.read_text())}
    checks=[]
    def arrays(ident):return np.load(ROOT/records[ident]['activation_file'],allow_pickle=False)
    def compare(a,b,kind,shared_prefix=False):
        model=records[a]['model'];assert records[b]['model']==model
        assert records[a]['mapping']['token_start']==records[b]['mapping']['token_start']
        with arrays(a) as aa,arrays(b) as bb:
            for layer in aa.files:
                x=aa[layer];y=bb[layer][:len(x)] if shared_prefix else bb[layer]
                assert x.shape==y.shape
                err=max(profile(unit(x)-unit(y)));eta=config['models'][model]['eta'][layer]
                assert err<=eta,(kind,a,b,layer,err,eta)
                checks.append({'kind':kind,'a':a,'b':b,'layer':int(layer),'residues':len(x),'max_normalized_error':err,'eta':eta,'bitwise_identical':bool(np.array_equal(x,y))})
    for prompt in prompts:
        for start in [165,337]:compare(f'profam1__ETV4_{start}_430_WT__{prompt}',f'profam1__ETV4_{start}_436_WT__{prompt}','full_target_C_extension',True)
    for start in [165,337]:compare(f'profam1__ETV4_{start}_436_WT__empty',f'profam1__ETV4_{start}_484_WT__empty','empty_prompt_distal_tail_extension',True)
    groups=collections.defaultdict(list)
    for k,r in records.items():groups[(r['model'],inputs[r['input_id']]['sequence'],tuple(prompts[r['prompt_id']]['sequences']))].append(k)
    for keys in groups.values():
        for key in keys[1:]:compare(keys[0],key,'scheduled_identical_input_repeat')
    out={'passed':True,'checks':checks,'comparison_count':len(checks),'counts_by_kind':dict(collections.Counter(x['kind'] for x in checks)),
      'max_normalized_error':max(x['max_normalized_error'] for x in checks),'all_bitwise_identical':all(x['bitwise_identical'] for x in checks),'new_model_forwards':0,
      'scope':'Checks all earlier target residues in the 430/436 causal pairs, empty-prompt distal-tail extensions, and every repeated scheduled input. Error thresholds are frozen eta.'}
    dest=ROOT/'results/E2-full-control-verification.json';assert not dest.exists();write_json(dest,out)
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
if __name__=='__main__':main()
