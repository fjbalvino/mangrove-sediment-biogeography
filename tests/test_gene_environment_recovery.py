"""Protect scientific identity, pre-existing data and georeference provenance."""
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1] / 'scripts/12_gene_environment/12B_APPLY_METADATA_RECOVERY.py'
SPEC = importlib.util.spec_from_file_location('recovery', PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.metadata = [{'run_accession':'run_a','experiment_accession':'experiment_a',
                          'lat':'','lon':'','soc_mean':'100','Country':'A'}]
        self.evidence = {key:'' for key in MODULE.REQUIRED}
        self.evidence.update(run_accession='run_a',experiment_accession='experiment_a',
                             recovery_status='published_site_coordinate',
                             recovery_lat='22.5',recovery_lon='39.1',
                             recovery_source_url='https://example.org/paper',
                             recovery_identity_basis='explicit sample label',
                             recovery_locality='sampling station',recovery_datum='unknown')

    def call(self):
        return MODULE.merge(list(self.metadata[0]),self.metadata,list(self.evidence),[self.evidence])

    def test_gazetteer_point_never_becomes_sampling_coordinate(self):
        self.evidence.update(recovery_status='locality_resolved_position_pending',
                             recovery_lat='',recovery_lon='',
                             recovery_gazetteer_lat='22.5',recovery_gazetteer_lon='39.1')
        _, rows, changed = self.call()
        self.assertEqual((rows[0]['lat'], rows[0]['lon']), ('',''))
        self.assertEqual(changed, [])

    def test_identity_mismatch_rejected(self):
        self.evidence['experiment_accession']='another_experiment'
        with self.assertRaisesRegex(ValueError,'identity mismatch'):
            self.call()

    def test_preexisting_coordinate_cannot_be_replaced(self):
        self.metadata[0].update(lat='20',lon='40')
        with self.assertRaisesRegex(ValueError,'cannot be overwritten'):
            self.call()

    def test_original_schema_preserved_without_importing_other_measurements(self):
        self.evidence['recovery_observed_TOC_g_kg']='9.83'
        fields,rows,_=self.call()
        self.assertEqual(fields,list(self.metadata[0]))
        self.assertEqual(list(rows[0]),list(self.metadata[0]))
        self.assertEqual(rows[0]['soc_mean'],'100')
        self.assertNotIn('recovery_observed_TOC_g_kg',rows[0])
        self.assertEqual(self.metadata[0]['lat'],'')


if __name__=='__main__':
    unittest.main()
