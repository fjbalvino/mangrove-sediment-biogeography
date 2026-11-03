import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('nitrogen_pilot', ROOT/'scripts/12_gene_environment/12B_BUILD_NITROGEN_PILOT.py')
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


class NitrogenContracts(unittest.TestCase):
    def test_id_errors(self):
        for values in [['A','A'],['A',''],[' A','B'],['A',None]]:
            with self.assertRaises(ValueError):mod.validate_ids(values,'test')

    def test_grouping_by_id_not_position(self):
        f=pd.DataFrame([[2.,4.,10.]],index=['K00001'],columns=['a','b','c'])
        meta=pd.DataFrame({'analysis_group':['s','s','s'],'biosample_group':['B','A','A']},index=['c','b','a'])
        self.assertEqual(mod.hierarchical_summary(f,meta).loc['s','K00001'],6.5)

    def test_unavailable_not_zero(self):
        f=pd.DataFrame([[np.nan,np.nan],[0.,0.]],index=['missing','zero'],columns=['a','b'])
        meta=pd.DataFrame({'analysis_group':['s','s'],'biosample_group':['A','B']},index=['a','b'])
        r=mod.hierarchical_summary(f,meta)
        self.assertTrue(np.isnan(r.loc['s','missing']))
        self.assertEqual(r.loc['s','zero'],0.)

    def test_within_matches_full_dummy_OLS(self):
        rng=np.random.default_rng(726);s=np.repeat(np.arange(10),5);x=rng.normal(size=(50,3))
        y=x@np.array([.2,-.4,.7])+s*10+rng.normal(size=50)*.1
        f=mod.within_fit(x,y,s)
        design=np.column_stack([np.eye(10)[s],x]);b=np.linalg.lstsq(design,y,rcond=None)[0][-3:]
        np.testing.assert_allclose(f['beta'],b,rtol=1e-10,atol=1e-10)

    def test_confounded_predictor_rejected(self):
        s=np.repeat(np.arange(5),3);x=s[:,None].astype(float)
        self.assertIsNone(mod.within_fit(x,s.astype(float),s))

    def test_study_intercept_invariance(self):
        rng=np.random.default_rng(826);s=np.repeat(np.arange(8),8);x=rng.normal(size=(64,3));y=x@np.array([.5,-1.,2.])
        np.testing.assert_allclose(mod.within_fit(x,y+s*100,s)['beta'],[.5,-1,2],atol=1e-12)

    def test_load_reorders_by_run_and_rejects_duplicate_headers(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'ids.txt').write_text('B\nA\n\n')
            meta=pd.DataFrame({'run_accession':['A','B'],'study_id':['s','s'],'sample_accession':['x','y'],'lat':[1.,1.],'lon':[2.,2.]})
            meta.to_csv(d/'meta.csv',index=False)
            (d/'matrix.tsv').write_text('KO\tA\tB\nK00001\t1\t2\n')
            m,v=mod.load_inputs(d/'meta.csv',d/'matrix.tsv',d/'ids.txt',2)
            self.assertEqual(list(m.index),['B','A']);self.assertEqual(v.loc['K00001','B'],2)
            (d/'matrix.tsv').write_text('KO\tA\tA\nK00001\t1\t2\n')
            with self.assertRaises(ValueError):mod.load_inputs(d/'meta.csv',d/'matrix.tsv',d/'ids.txt',2)


if __name__=='__main__':unittest.main()
