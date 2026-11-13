import unittest
import numpy as np
from analyse_results import normalized,rms,context_metrics,specificity,factorial,prompt_response,sequence_baseline,max_discrepancy
class Arithmetic(unittest.TestCase):
    def test_normalization_float64_and_zero_guard(self):
        np.testing.assert_allclose(normalized(np.array([[3,4]],np.float32)),[[.6,.8]])
        self.assertEqual(normalized(np.array([[3,4]],np.float32)).dtype,np.float64)
        for arr in (np.zeros((2,2)),np.array([[np.nan,1.]]),np.zeros((0,2))):
            with self.assertRaises(ValueError): normalized(arr)
    def test_rms_not_feature_mean(self):
        self.assertEqual(rms([[3,4],[3,4]]),5.)
    def test_context_hand_calculation(self):
        # q_WT=(-1,1), q_g=(0,0); M=sqrt(2), B=sqrt(2), D=sqrt(2).
        e1=np.array([[1.,0.]]);e2=np.array([[0.,1.]])
        m=context_metrics(e2,e2,e2,e1,1e-7)
        for key in ['M','B','D']: self.assertAlmostEqual(m[key],2**.5)
        self.assertAlmostEqual(m['T'],1.);self.assertAlmostEqual(m['M_over_D'],1.)
        self.assertLess(m['L'],1.);self.assertGreater(m['U'],1.)
    def test_undefined_denominators_not_regularized(self):
        e1=np.array([[1.,0.]]);e2=np.array([[0.,1.]])
        m=context_metrics(e2,e1,e1,e1,1e-7)
        self.assertIsNone(m['T']);self.assertIsNone(m['M_over_D']);self.assertIsNone(m['U'])
    def test_factorial_cancellation_and_nonzero(self):
        a=np.array([[1.,0.]]);b=np.array([[0.,1.]])
        self.assertEqual(factorial(a,b,a,b)['rms'],0.)
        self.assertAlmostEqual(factorial(a,a,a,b)['rms'],2**.5)
        self.assertAlmostEqual(rms(prompt_response(a,b)),2**.5)
    def test_fixed_rule_ties_controls_and_unresolved(self):
        names=['Y401N_Y403G','E404K','E423K','E425K','F428A_F432A']
        m={n:{'M':.1,'B':1.,'D':1.,'L':.09,'U':.11} for n in names}
        m[names[0]]={'M':.5,'B':1.,'D':1.,'L':.49,'U':.51}
        sh={n:[.2]*20 for n in names}
        self.assertEqual(specificity(m,sh,1e-7)['classification'],'source_consistent_native_order_specificity')
        sh[names[0]][0]=.5
        self.assertEqual(specificity(m,sh,1e-7)['classification'],'specificity_not_established')
        m[names[1]]['D']=4e-7;m[names[1]]['U']=None
        self.assertEqual(specificity(m,sh,1e-7)['classification'],'numerically_unresolved')
        sh[names[0]].pop()
        with self.assertRaises(ValueError): specificity(m,sh,1e-7)
    def test_baseline_and_error_units(self):
        b=sequence_baseline('AKKRED');self.assertEqual(b['length'],6)
        self.assertAlmostEqual(b['charge_proxy'],1/6)
        self.assertAlmostEqual(sum(b['frequencies'].values()),1.)
        self.assertAlmostEqual(max_discrepancy(np.array([[1.,0.]]),np.array([[0.,1.]])),2**.5)





class Preparation(unittest.TestCase):
    def test_plan_matches_original_sequences_and_masks(self):
        import json, hashlib
        from pathlib import Path
        package=Path(__file__).resolve().parents[1]
        plan=json.loads((package/'design/run-plan.json').read_text())
        inputs=json.loads((package/'data/derived/model-inputs.json').read_text())
        self.assertEqual(len(plan['runs']),217)
        self.assertEqual(len({r['run_id'] for r in plan['runs']}),217)
        for r in plan['runs']:
            x=inputs[r['input_id']]
            self.assertEqual(hashlib.sha256(x['sequence'].encode()).hexdigest(),r['input_sha256'])
            self.assertEqual(len(x['sequence'])+2,r['expected_tokens'])
            self.assertEqual(r['layers'],[19,24,36])
            self.assertEqual(r['prompt_id'],'empty')
    def test_login_guard(self):
        import os
        from unittest.mock import patch
        import run_experiments as r
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaisesRegex(RuntimeError,'Slurm'):r.environment_check()
    def test_complete_synthetic_analysis_and_tamper_rejection(self):
        import json, tempfile, shutil, hashlib
        from pathlib import Path
        import run_experiments as r
        from analyse_results import analyse_export
        package=Path(__file__).resolve().parents[1]
        old=r.ROOT
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp).resolve();r.ROOT=root
                for part in ['data','design']:shutil.copytree(package/part,root/part)
                dest=root/'results/activations/test';dest.mkdir(parents=True)
                cfg={'status':'frozen','models':{'esm2_3b':{'width':2,'layers':[19,24,36],'eta':{str(k):1e-7 for k in [19,24,36]}}}}
                r.write_json(root/'design/execution-config.json',cfg)
                cfg_hash=r.digest(root/'design/execution-config.json')
                r.write_json(root/'qualification-receipt.json',{'passed':True,'config_sha256':cfg_hash,'input_hashes':{}})
                inputs=r.read('data/derived/model-inputs.json');records=[]
                for row in r.read('design/run-plan.json')['runs']:
                    x=inputs[row['input_id']]
                    rng=np.random.default_rng(int(hashlib.sha256(x['sequence'].encode()).hexdigest()[:8],16))
                    arrays={str(k):rng.normal(size=(x['length'],2)).astype(np.float32) for k in [19,24,36]}
                    path=dest/(row['run_id']+'.npz');np.savez(path,**arrays)
                    records.append({**row,'activation_file':str(path.relative_to(root)),'sha256':r.digest(path),'native_start':x['native_start'],'execution_config_sha256':cfg_hash})
                manifest=dest/'manifest.json';r.write_json(manifest,records)
                output=root/'results/synthetic.json';analyse_export([manifest],output)
                result=json.loads(output.read_text())
                self.assertEqual(result['fixed_run_count'],217)
                self.assertEqual(set(result['models']['esm2_3b']),{'19','24','36'})
                self.assertEqual(len(result['models']['esm2_3b']['24']['factorial_shuffled']),20)
                path.write_bytes(b'tampered')
                with self.assertRaisesRegex(ValueError,'identity'):analyse_export([manifest],root/'results/bad.json')
        finally:r.ROOT=old

if __name__ == '__main__': unittest.main()
