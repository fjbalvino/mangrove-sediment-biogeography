"""Scientific invariants introduced by the upstream implementation audit."""
from pathlib import Path
import importlib.util
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT/'scripts/12_gene_environment/template'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


rules = module('rules_audit', TEMPLATE/'evaluate.py')
models = module('models_audit', TEMPLATE/'analysis_nitrogen/scripts/12_gene_environment/12B_BUILD_NITROGEN_PILOT.py')


class AuditContracts(unittest.TestCase):
    def test_partial_module_does_not_relax_complex(self):
        values = {'K00001': np.array([1]), 'K00002': np.array([0]), 'K00003': np.array([1])}
        x = rules.block_evidence('K00001+K00002 K00003 K00004', values, 1)
        self.assertEqual(x['n_blocks'], 3)
        self.assertEqual(x['states'].ravel().tolist(), [0,1,2])
        self.assertAlmostEqual(x['coverage_lower'][0], 1/3)
        self.assertAlmostEqual(x['coverage_upper'][0], 2/3)
        self.assertEqual(rules.evaluate(rules.parse('K00001+K00002 K00003 K00004'), values, 1)[0], 0)

    def test_unknown_optional_and_alternative(self):
        x = rules.block_evidence('(K00001,K00002) K00003-K00004', {'K00001':np.array([1]), 'K00003':np.array([1])}, 1)
        np.testing.assert_array_equal(x['coverage_lower'], [1])
        np.testing.assert_array_equal(x['coverage_upper'], [1])
        x = rules.block_evidence('K00001+K00002', {'K00001':np.array([1])}, 1)
        self.assertEqual(x['n_blocks'], 1)
        self.assertEqual(x['n_unknown'][0], 1)

    def test_partial_information_matches_fixed_effect_geometry(self):
        rng = np.random.default_rng(42)
        study = np.repeat(np.arange(4), 10)
        z = rng.normal(size=(40, 2))
        z[:,0] += 0.95*z[:,1]
        x = z + np.repeat(np.arange(4)*9,10)[:,None]
        y = x @ [2,-1] + study*5 + rng.normal(size=40)
        fit = models.within_fit(x,y,study)
        design = np.column_stack([np.eye(4)[study],x])
        expected = np.linalg.lstsq(design,y,rcond=None)[0][-2:]
        np.testing.assert_allclose(fit['beta'], expected, atol=1e-11)
        for j in range(2):
            residual = fit['partial_predictors'][:,j]
            self.assertAlmostEqual(residual@y/(residual@residual),fit['beta'][j])
            self.assertLessEqual(fit['conditional_variation_fraction'][j],1+1e-12)
            self.assertGreaterEqual(fit['effective_information_studies'][j],1-1e-12)
            self.assertLessEqual(fit['effective_information_studies'][j],4+1e-12)
        self.assertFalse(np.allclose(fit['max_study_information_share'],fit['max_conditional_study_information_share']))

    def test_equal_study_information_and_scale_invariance(self):
        x = np.tile(np.array([[-2,-1],[-1,2],[1,-2],[2,1.]]),(4,1))
        study = np.repeat(np.arange(4),4)
        y = x @ [3,2] + study
        a = models.within_fit(x,y,study)
        b = models.within_fit(x*[100,0.01],y,study)
        np.testing.assert_allclose(a['effective_information_studies'],4)
        np.testing.assert_allclose(a['max_conditional_study_information_share'],0.25)
        np.testing.assert_allclose(a['conditional_variation_fraction'],b['conditional_variation_fraction'])
        self.assertIsNone(models.within_fit(np.column_stack([x[:,0],x[:,0]]),y,study))


if __name__ == '__main__':
    unittest.main()
