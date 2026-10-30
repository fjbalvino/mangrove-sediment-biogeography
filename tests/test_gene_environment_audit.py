"""Small synthetic checks of provenance-sensitive audit behavior."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/12_gene_environment/12A_AUDIT_GENE_ENVIRONMENT_DESIGN.py'
SPEC = importlib.util.spec_from_file_location('gene_environment_audit', SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.meta = self.root/'meta.csv'; self.matrix=self.root/'matrix.tsv'
        self.ids=self.root/'ids.txt'; self.out=self.root/'out'
        pd.DataFrame(dict(run_accession=['s1','s2','s3','s4'], study_id=['A','A','B','B'],
            sampling_site_id=['p1','p1','p2','p3'], sample_accession=['b1','b1','b2','b3'],
            lat=[1,1,2,3],lon=[4,4,5,6], REALM=['R1','R1','R2','R2'],
            bio1=[24,24,25,26],soc_mean=[50,50,70,90],phh2o_mean=[60,60,65,70],
            human_footprint_2020=[1,2,3,4])).to_csv(self.meta,index=False)
        pd.DataFrame({'KO':['K00001','K00002'],'s1':[2,0],'s2':[3,1],'s3':[2,1],'s4':[1,4]}).to_csv(self.matrix,sep='\t',index=False)
        self.ids.write_text('s1\ns2\ns3\ns4\n')

    def tearDown(self):
        self.tmp.cleanup()

    def run_audit(self):
        argv=[str(SCRIPT),'--metadata',str(self.meta),'--matrix',str(self.matrix),
              '--samples',str(self.ids),'--unit','predicted_protein_records','--outdir',str(self.out)]
        with patch('sys.argv',argv):
            MODULE.main()

    def test_preserves_inputs_and_reports_dependence(self):
        before=self.meta.read_bytes(); self.run_audit()
        d=json.loads((self.out/'audit_summary.json').read_text())
        self.assertEqual(d['repeated_BioSamples'],1)
        self.assertEqual(d['rank_increment_realm_given_study'],0)
        self.assertEqual(self.meta.read_bytes(),before)
        self.assertTrue(any('exceeds 14' in x for x in d['unit_flags']))
        self.assertEqual(len(pd.read_csv(self.out/'coordinate_value_conflicts.tsv',sep='\t')),1)

    def test_negative_values_rejected(self):
        x=pd.read_csv(self.matrix,sep='\t');x.loc[0,'s1']=-1;x.to_csv(self.matrix,sep='\t',index=False)
        with self.assertRaises(ValueError): self.run_audit()
        self.assertFalse(self.out.exists())

    def test_missing_selected_accession_rejected(self):
        self.ids.write_text('s1\nmissing\n')
        with self.assertRaises(ValueError): self.run_audit()

    def test_no_overwrite(self):
        self.out.mkdir();(self.out/'keep.txt').write_text('keep')
        with self.assertRaises(ValueError):self.run_audit()
        self.assertEqual((self.out/'keep.txt').read_text(),'keep')


if __name__ == '__main__':
    unittest.main()
