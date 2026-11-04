"""Independent numerical and biological-contract checks for calibration."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from scipy.stats import permutation_test

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('calibration',ROOT/'scripts/10_functional_network/10M_CALIBRATE_PROCESS_COHERENCE.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class CalibrationTests(unittest.TestCase):
    def test_exact_flips_agree_with_scipy_including_ties(self):
        for x in [np.array([.3,-.1,.2,.1]),np.array([1.,1.,0.,-1.])]:
            _,p=m.signflip_test(x)
            ref=permutation_test((x,),np.mean,permutation_type='samples',
                                 n_resamples=np.inf,alternative='greater')
            self.assertAlmostEqual(p[0],ref.pvalue,places=13)

    def test_no_pseudoreplication_of_coordinate_groups(self):
        # Four coordinate observations form two blocks, hence four flips.
        y=np.array([[0.,0.],[1.,1.],[0.,0.],[2.,2.]])
        c,s,_,_=m.block_contributions(y,['a','a','b','b'],np.array([[0,1]]),np.ones((1,1)))
        _,p=m.signflip_test(c)
        self.assertEqual(len(s),2)
        self.assertEqual(m.exact_signs(len(s)).shape,(4,2))
        self.assertEqual(p[0],.25)

    def test_holm_keeps_complete_declared_family(self):
        np.testing.assert_allclose(m.holm([.01,.04,.03,1]),[.04,.09,.09,1])
        self.assertEqual(m.holm([1/1024]+[1]*7)[0],8/1024)

    def test_coherence_equals_pooled_pair_proportionality(self):
        y=np.array([[0.,1.,3.],[2.,3.,2.],[5.,6.,1.],[0.,2.,4.],[1.,3.,2.]])
        pairs=np.array([[0,1],[1,2]])
        c,_,_,_=m.block_contributions(y,['a','a','a','b','b'],pairs,np.ones((2,1))*.5)
        expected,_,_=m.pilot.pair_rho(y,pairs[:,0],pairs[:,1],['a','a','a','b','b'],True)
        self.assertAlmostEqual(c.mean(),expected.mean(),places=13)

    def test_invariant_markers_do_not_create_perfect_coherence(self):
        y=np.zeros((6,2))
        c,_,rho,ok=m.block_contributions(y,['a']*3+['b']*3,np.array([[0,1]]),np.ones((1,1)))
        self.assertFalse(ok[0]);self.assertEqual(rho[0],0)
        self.assertEqual(m.signflip_test(c)[1][0],1)

    def test_shared_marker_across_steps_fails_instead_of_self_edge(self):
        rules=pd.DataFrame({'expression':['K00001','K00001']},index=['r1','r2'])
        chains=[{'id':'c','steps':[{'rule':'r1','label':'one'},{'rule':'r2','label':'two'}]}]
        with self.assertRaisesRegex(ValueError,'Shared KO'):
            m.chain_pairs(chains,rules,['K00001'])

    def test_alternatives_remain_in_same_step(self):
        rules=pd.DataFrame({'expression':['(K00001,K00002)','K00003']},index=['r1','r2'])
        chains=[{'id':'c','steps':[{'rule':'r1','label':'one'},{'rule':'r2','label':'two'}]}]
        pairs,w,ok,_=m.chain_pairs(chains,rules,['K00001','K00002','K00003'])
        self.assertEqual(pairs.tolist(),[[0,2],[1,2]])
        np.testing.assert_allclose(w[:,0],[.5,.5]);self.assertTrue(ok[0])


if __name__=='__main__':unittest.main()
