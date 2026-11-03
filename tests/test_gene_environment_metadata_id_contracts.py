"""Sample identity must dominate row order and shared geographic coordinates."""
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1] / 'scripts/02_environmental_covariates/recovery/id_contracts.py'
SPEC = importlib.util.spec_from_file_location('metadata_id_contracts', PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(run, value='', lat='20', lon='100'):
    return {'run_accession': run, 'lat': lat, 'lon': lon, 'bio1': value, 'original_note': run}


class SampleIdentityTests(unittest.TestCase):
    def test_shared_coordinates_and_reordered_rows_join_only_by_id(self):
        master = [row('A', '10'), row('B', '20')]
        target = [row('B'), row('A')]
        result = MODULE.transfer_by_id(master, target, ['bio1'])
        self.assertEqual([r['run_accession'] for r in result], ['B', 'A'])
        self.assertEqual([r['bio1'] for r in result], ['20', '10'])
        self.assertEqual([r['original_note'] for r in result], ['B', 'A'])
        self.assertEqual(target, [row('B'), row('A')])
        self.assertEqual([list(r) for r in result], [list(r) for r in target])

    def test_duplicate_ids_in_either_database_are_rejected(self):
        for master, target in [([row('A'), row('A')], [row('A')]), ([row('A')], [row('A'), row('A')])]:
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                MODULE.transfer_by_id(master, target, ['bio1'])

    def test_unknown_and_invalid_ids_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'absent from master'):
            MODULE.transfer_by_id([row('A')], [row('B')], ['bio1'])
        for invalid in ['', ' A ', None, 'NA']:
            with self.assertRaises(ValueError):
                MODULE.index_by_id([row(invalid)], 'test')

    def test_coordinate_conflict_cannot_be_fixed_by_another_id(self):
        with self.assertRaisesRegex(ValueError, 'coordinate conflict'):
            MODULE.transfer_by_id([row('A', '10')], [row('A', lat='21')], ['bio1'])

    def test_existing_values_remain_unchanged_and_conflicts_fail(self):
        self.assertEqual(MODULE.transfer_by_id([row('A', '10')], [row('A', '10.0')], ['bio1'])[0]['bio1'], '10.0')
        with self.assertRaisesRegex(ValueError, 'conflicting existing value'):
            MODULE.transfer_by_id([row('A', '10')], [row('A', '11')], ['bio1'])

    def test_count_and_exclusion_contracts(self):
        with self.assertRaisesRegex(ValueError, 'expected 401 samples'):
            MODULE.index_by_id([row('A')], 'master', 401)
        target = [row('A')]
        self.assertEqual(MODULE.transfer_by_id([row('A', '10')], target, ['bio1'], excluded={'A'}), target)


if __name__ == '__main__':
    unittest.main()
