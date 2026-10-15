"""Independent checks for the functional atlas's statistical contracts."""
from pathlib import Path
import importlib.util
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'functional_atlas', ROOT/'scripts/08_functional_atlas/08A_ANALYSE_FUNCTIONAL_ATLAS.py')
atlas = importlib.util.module_from_spec(spec)
spec.loader.exec_module(atlas)


class FunctionalStatistics(unittest.TestCase):
    def test_unbalanced_euclidean_permanova_against_group_sums(self):
        y = np.array([[0,1],[2,0],[1,3],[4,2],[7,3],[6,5],[8,4]],dtype=float)
        groups = np.array(['A']*3+['B']*4)
        total = ((y-y.mean(axis=0))**2).sum()
        residual = sum(((y[groups==g]-y[groups==g].mean(axis=0))**2).sum() for g in ['A','B'])
        fit = atlas.permutation_model(y,groups,17,nperm=19)
        self.assertAlmostEqual(fit['R2'],(total-residual)/total,places=12)
        self.assertAlmostEqual(fit['pseudo_F'],(total-residual)/(residual/5),places=12)
        self.assertEqual((fit['df'],fit['residual_df']),(1,5))

    def test_conditional_permanova_against_direct_least_squares(self):
        y = np.array([[0,1],[2,0],[1,3],[4,2],[7,3],[6,5],[8,4]],dtype=float)
        groups = np.array(['A']*3+['B']*4)
        x = np.array([0,1,3,2,5,4,7],dtype=float)
        reduced = np.column_stack([np.ones(7),x])
        full = np.column_stack([reduced,groups=='B'])
        def sse(design):
            residual = y-design@np.linalg.lstsq(design,y,rcond=None)[0]
            return (residual**2).sum()
        term = sse(reduced)-sse(full)
        total = ((y-y.mean(axis=0))**2).sum()
        fit = atlas.permutation_model(y,groups,17,nperm=19,covariates=x)
        self.assertAlmostEqual(fit['sum_of_squares'],term,places=10)
        self.assertAlmostEqual(fit['R2'],term/total,places=12)
        self.assertAlmostEqual(fit['partial_R2'],term/(term+sse(full)),places=12)
        self.assertAlmostEqual(fit['pseudo_F'],term/(sse(full)/4),places=12)
        self.assertEqual((fit['df'],fit['residual_df']),(1,4))

    def test_persistence_priority_and_fraction_threshold(self):
        ids = [f's{i}' for i in range(20)]
        counts = pd.DataFrame([
            [100000]*20, [1]*19+[0], [50]*10+[0]*10,
            [1]*10+[0]*10, [50]+[0]*19,
        ], index=['absolute','near','realm','too_small','one_locality'],columns=ids)
        meta = pd.DataFrame({'REALM':['R1']*10+['R2']*10,
            'core_locality':['L1']*8+['L2','L3']+['L4']*8+['L5','L6']},index=ids)
        annotations,summary,support,_ = atlas.classification(counts,meta)
        classes = annotations.set_index('KO')['class'].to_dict()
        self.assertEqual(classes,{'absolute':'Absolute core','near':'Near-global',
            'realm':'Realm-core','too_small':'Periphery','one_locality':'Periphery'})
        self.assertEqual(int(summary.KOs.sum()),5)
        self.assertAlmostEqual(summary.KO_assignment_percent.sum(),100)
        row = support[(support.KO=='too_small') & (support.realm=='R1')].iloc[0]
        self.assertEqual(row.localities_present,3)
        self.assertLess(row.mean_relative_abundance,0.0001)
        self.assertGreater(row.mean_relative_abundance,0.000001)

    def test_site_means_give_equal_weight_to_profiles(self):
        observed = atlas.group_mean(np.array([[0.,2.],[2.,4.],[10.,8.]]),
                                    np.array(['A','A','B']),['A','B'])
        np.testing.assert_allclose(observed,[[1.,3.],[10.,8.]])


if __name__ == '__main__':
    unittest.main()
