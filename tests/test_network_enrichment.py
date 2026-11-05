"""Check exact conditional enrichment against exhaustive allocations."""
import importlib.util
import itertools
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, false_discovery_control

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('enrichment',ROOT/'scripts/10_functional_network/10P_ENRICH_NETWORK_FUNCTIONS.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class EnrichmentTests(unittest.TestCase):
    def test_exact_stratified_distribution_matches_exhaustive_draws(self):
        counts=[]
        for a in itertools.combinations(range(5),2):
            for b in itertools.combinations(range(4),2):
                counts.append(sum(i<3 for i in a)+sum(i<1 for i in b))
        expected=np.bincount(counts,minlength=4)/len(counts)
        distribution,mean,var=m.exact_null([(5,3,2),(4,1,2)])
        np.testing.assert_allclose(distribution,expected,atol=1e-14)
        self.assertAlmostEqual(mean,np.mean(counts))
        self.assertAlmostEqual(var,np.var(counts))

    def test_unstratified_tail_matches_scipy(self):
        d,_,_=m.exact_null([(100,20,15)])
        for observed in range(16):
            self.assertAlmostEqual(d[observed:].sum(),hypergeom.sf(observed-1,100,20,15),places=13)

    def test_constant_blocks_do_not_create_evidence(self):
        d,mean,var=m.exact_null([(5,5,2),(4,0,1)])
        self.assertAlmostEqual(d[2],1.)
        self.assertEqual(mean,2.)
        self.assertAlmostEqual(var,0.)

    def test_permutations_preserve_blockwise_module_counts(self):
        labels=np.array([0,0,1,1,2,2]); blocks=np.array([0,0,0,1,1,1])
        draws=m.permutation_labels(labels,blocks,np.random.default_rng(2),100)
        for draw in draws:
            for block in [0,1]:
                np.testing.assert_array_equal(np.sort(draw[blocks==block]),np.sort(labels[blocks==block]))

    def test_BY_uses_complete_family_with_zero_count_tests(self):
        p=np.array([.001,.01,1.,1.])
        order=np.argsort(p);ranks=np.arange(1,len(p)+1)
        expected=np.minimum.accumulate((p[order]*len(p)*sum(1/ranks)/ranks)[::-1])[::-1].clip(0,1)
        np.testing.assert_allclose(false_discovery_control(p,method='by')[order],expected)
        self.assertGreater(false_discovery_control(p,method='by')[0],false_discovery_control(p[:2],method='by')[0])

    def test_covariate_blocks_do_not_use_module_or_function_labels(self):
        frame=pd.DataFrame({'x':np.linspace(0,10,60),'y':np.sin(np.arange(60)),
            'module_label':['a']*30+['b']*30},index=['K%05d'%i for i in range(60)])
        config={'covariates':['x','y'],'block_target_size':15,'block_min_size':8,
                'kmeans_n_init':50,'kmeans_max_iter':200,'seed':1}
        a=m.match_blocks(frame,config)
        frame.module_label=frame.module_label.iloc[::-1].to_numpy()
        np.testing.assert_array_equal(a,m.match_blocks(frame,config))


if __name__=='__main__':unittest.main()
