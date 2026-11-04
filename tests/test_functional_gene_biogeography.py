"""Portability and measurement contracts using synthetic, temporary inputs."""
from pathlib import Path
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DRIVER=ROOT/'workflow/08_RUN_GENE_BIOGEOGRAPHY.py'


class GeneBiogeography(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='gene workflow ')
        self.addCleanup(self.tmp.cleanup);self.base=Path(self.tmp.name)
        self.ids=['run_'+str(i) for i in range(16)]
        rng=np.random.default_rng(142)
        self.meta=pd.DataFrame({'run_accession':self.ids,'study_id':np.repeat(['A','B','C','D'],4),
            'sample_accession':self.ids,'lat':np.repeat([10.,11.,12.,13.],4),'lon':np.tile([20.,21.,22.,23.],4),
            'Country':'Synthetic','Locality':'Synthetic','REALM':'Synthetic','bio1':rng.normal(22,2,16),
            'soc_mean':rng.uniform(50,400,16),'phh2o_mean':rng.uniform(50,80,16)})
        # An extra metadata row must never be appended to the selected cohort.
        extra=self.meta.iloc[[0]].copy();extra.run_accession='unselected'
        pd.concat([self.meta,extra]).sample(frac=1,random_state=2).to_csv(self.base/'meta.csv',index=False)
        self.mat=pd.DataFrame(rng.integers(1,50,(8,16)),index=['K02588','K02586','K02591','K00855','K00399','K02040','K04759','K00001'],columns=self.ids)
        self.mat.loc['K00855',self.ids[0]]=0
        self.mat.to_csv(self.base/'counts.tsv',sep='\t',index_label='KO')
        self.env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')

    def run_command(self,*args,ok=True):
        result=subprocess.run([sys.executable,str(DRIVER),*map(str,args)],capture_output=True,text=True,env=self.env,timeout=90)
        if ok:self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        else:self.assertNotEqual(result.returncode,0)
        return result

    def args(self):
        return ['--matrix',self.base/'counts.tsv','--metadata',self.base/'meta.csv','--all-matrix-samples',
                '--unit','predicted_protein_records','--run-dir',self.base/'run']

    def test_dynamic_cohort_id_matching_and_unknowns(self):
        self.run_command(*self.args(),'--execute')
        run=self.base/'run'
        manifest=json.loads((run/'run_manifest.json').read_text())
        self.assertEqual(manifest['profiles'],16)
        rel=pd.read_csv(run/'atlas/results/relative_representation.tsv.gz',sep='\t',index_col=0)
        self.assertEqual(list(rel.columns),self.ids)
        self.assertEqual(rel.loc['K00855',self.ids[0]],0)
        self.assertTrue(rel.loc['K01601'].isna().all())
        expected=self.mat.loc['K02588']/self.mat.sum()*1e6
        np.testing.assert_allclose(rel.loc['K02588'],expected)
        self.assertNotIn('unselected',rel.columns)
        states=pd.read_csv(run/'outputs/rule_block_states_ge1.tsv.gz',sep='\t').set_index(['rule_id','block'])
        self.assertEqual(list(states.columns),self.ids)
        self.assertTrue((states.loc['N01_nifHDK'].to_numpy()==1).all())
        self.assertEqual(states.loc['N01_nifHDKENB',self.ids[0]].tolist(),[1,2,2,2])
        preview=self.run_command('--stage','atlas','--atlas-cycle','N','--run-dir',run)
        self.assertIn('--cycle N',preview.stdout)
        self.assertNotIn('plot_supplementary_panels.py',preview.stdout)
        self.assertFalse((run/'atlas/figures/nitrogen/viewer.html').exists())
        self.run_command(*self.args(),'--execute',ok=False)
        with (run/'inputs/metadata_functional.csv').open('a') as handle:handle.write('\n')
        bad=self.run_command('--stage','atlas','--run-dir',run,'--execute',ok=False)
        self.assertIn('Prepared input/source changed',bad.stderr)

    def test_duplicate_ids_and_normalized_values_rejected(self):
        self.meta.loc[1,'run_accession']=self.meta.loc[0,'run_accession']
        self.meta.to_csv(self.base/'meta.csv',index=False)
        bad=self.run_command(*self.args(),'--execute',ok=False)
        self.assertIn('duplicate identifier',bad.stderr)
        self.meta.loc[1,'run_accession']=self.ids[1];self.meta.to_csv(self.base/'meta.csv',index=False)
        self.mat=self.mat.astype(float);self.mat.iloc[0,0]=1.5;self.mat.to_csv(self.base/'counts.tsv',sep='\t')
        bad=self.run_command(*self.args(),'--execute',ok=False)
        self.assertIn('never round normalized data',bad.stderr)

    def test_preview_is_read_only(self):
        self.run_command(*self.args())
        self.assertFalse((self.base/'run').exists())

    def test_constant_predictors_are_not_estimable(self):
        p=ROOT/'scripts/12_gene_environment/template/analysis_nitrogen/scripts/12_gene_environment/12B_BUILD_NITROGEN_PILOT.py'
        spec=importlib.util.spec_from_file_location('model',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        self.assertIsNone(m.within_fit(np.ones((8,3)),np.arange(8),['a']*4+['b']*4))


if __name__=='__main__':unittest.main()
