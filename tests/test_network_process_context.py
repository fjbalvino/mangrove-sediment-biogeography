"""Independent biological-logic and sampling checks for the network pilot."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('network_pilot', ROOT / 'scripts/10_functional_network/10L_AUDIT_PROCESS_CONTEXT.py')
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)
import biogeochemical_context as context


class ProcessContextTests(unittest.TestCase):
    def test_missing_from_network_is_not_missing_from_matrix(self):
        expr = 'K00001+K00002'
        self.assertEqual(context.module_block_state(expr, {'K00001'}, {'K00001','K00002'}), 'partial_representation')
        self.assertEqual(context.module_block_state(expr, set(), {'K00001','K00002'}), 'not_represented_in_module')
        self.assertEqual(context.module_block_state(expr, set(), set()), 'not_resolvable_from_matrix')

    def test_alternatives_and_optional_components(self):
        self.assertEqual(context.module_block_state('(K00001,K00002)-K00003', {'K00002'}, {'K00002'}), 'logical_block_represented')
        self.assertEqual(context.get_blocks('(K00001+K00002) (K00003,K00004)'), ['(K00001+K00002)','(K00003,K00004)'])

    def test_adjacent_steps_must_cooccur_in_the_same_sample(self):
        values = {'K00001':np.array([1,0]), 'K00002':np.array([0,1])}
        np.testing.assert_array_equal(context.logic.evaluate(context.logic.parse('K00001 K00002'), values, 2), [0,0])

    def test_proportionality_against_log_ratio_variance_identity(self):
        x = np.array([[1.,3.],[2.,2.],[4.,7.],[8.,9.]])
        observed,_,_ = pilot.pair_rho(x, np.array([0]), np.array([1]))
        expected = 1 - np.var(x[:,0]-x[:,1],ddof=1)/(np.var(x[:,0],ddof=1)+np.var(x[:,1],ddof=1))
        self.assertAlmostEqual(observed[0], expected, places=13)

    def test_between_study_shift_does_not_prove_within_study_covariation(self):
        x = np.array([[0.,2.],[2.,0.],[10.,12.],[12.,10.]])
        pooled,_,_ = pilot.pair_rho(x, np.array([0]), np.array([1]))
        within,n,ns = pilot.pair_rho(x, np.array([0]), np.array([1]), ['a','a','b','b'], True)
        self.assertGreater(pooled[0], .9)
        self.assertAlmostEqual(within[0], -1.)
        self.assertEqual((n,ns),(4,2))

    def test_single_coordinate_study_adds_no_within_study_replication(self):
        x = np.array([[0.,0.],[1.,1.],[100.,-100.]])
        result,n,ns = pilot.pair_rho(x, np.array([0]), np.array([1]), ['a','a','b'], True)
        self.assertEqual((n,ns),(2,1))
        self.assertAlmostEqual(result[0],1.)

    def test_repeated_runs_do_not_upweight_a_biosample(self):
        y = pd.DataFrame([[0.,2.,10.]],index=['K00001'],columns=['r1','r2','r3'])
        meta = pd.DataFrame({'analysis_group':['g']*3,'biosample_group':['a','a','b'],
            'study_id':['s']*3,'coordinate_group':['c']*3,'lat':[0]*3,'lon':[0]*3},index=y.columns)
        grouped,_ = pilot.grouped_profiles(y,meta)
        self.assertAlmostEqual(grouped.iloc[0,0],5.5)


if __name__ == '__main__':
    unittest.main()
