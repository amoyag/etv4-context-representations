"""Model-free tests with analytical expectations; actual-model tests are separate."""
import importlib.util
import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_experiments import native_indices, validate_sequences, require_e2
from analyse_results import normalized, rms, context_metrics, specificity, factorial, prompt_response, sequence_baseline, max_discrepancy


class Coordinates(unittest.TestCase):
    def test_native_to_tokens_with_prompt_and_n_prefix(self):
        # Two 477-residue family sequences each followed by SEP; target begins at 165.
        self.assertEqual(native_indices([337,401,430],165,272,958),[1130,1194,1223])
        self.assertEqual(native_indices([337,430],337,100,1),[1,94])
    def test_invalid_coordinates(self):
        for p in ([336],[437],[337,337]):
            with self.assertRaises(ValueError): native_indices(p,337,100)
    def test_no_implicit_unknown_truncation_or_empty(self):
        for seq in ('','AX','ac','A[SEP]'):
            with self.assertRaises(ValueError): validate_sequences(seq,[])
    def test_e2_is_not_authorized(self):
        with patch('run_experiments.read', return_value={'active_execution':True,'activation_scope':['E1']}):
            with self.assertRaisesRegex(RuntimeError,'E2 is not'): require_e2()


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


if __name__ == '__main__': unittest.main()
