"""Direction, distance, completeness and preservation contracts for NA recovery."""
import csv
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts/02_environmental_covariates"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("cardinal", SCRIPTS / "02J_RECOVER_BIOORACLE_SALINITY_NA.py")
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def periods(values):
    return dict(zip(m.bio.TIMES, values))


class FakeSource:
    def __init__(self, valid):
        self.cells = {index: periods(values) for index, values in valid.items()}
        self.records = {}

    def ensure(self, indices):
        for index in indices:
            self.cells.setdefault(index, periods([(None, None), (None, None)]))


class CardinalContracts(unittest.TestCase):
    def search(self, point, source, steps=(5000, 10000), bound=20000):
        with redirect_stdout(io.StringIO()):
            return m.cardinal_search({point}, source, steps, bound)

    def test_original_point_precedes_offsets_and_zero_is_data(self):
        p = (.025, .025)
        source = FakeSource({m.bio.cell_index(*p): [(0, 0), (0, 0)]})
        selected, attempts = self.search(p, source)
        self.assertEqual(selected[p]['direction'], 'origin')
        self.assertEqual(selected[p]['values'], (0, 0))
        self.assertEqual(len(attempts), 1)

    def test_north_precedes_east_at_same_radius(self):
        p = (.025, .025)
        north = m.bio.cell_index(*m.offset_point(*p, 5000, 0))
        east = m.bio.cell_index(*m.offset_point(*p, 5000, 90))
        source = FakeSource({north: [(10, 2), (20, 4)], east: [(30, 6), (40, 8)]})
        selected, _ = self.search(p, source)
        self.assertEqual(selected[p]['direction'], 'N')
        self.assertEqual(selected[p]['values'], (15, 3))

    def test_increasing_radii_and_actual_center_distance_are_distinct(self):
        p = (.025, .025)
        east = m.bio.cell_index(*m.offset_point(*p, 10000, 90))
        selected, attempts = self.search(p, FakeSource({east: [(1, 2), (3, 4)]}), steps=(500, 1000, 5000, 10000))
        self.assertEqual(selected[p]['direction'], 'E')
        self.assertEqual(selected[p]['radius_m'], 10000)
        self.assertGreater(selected[p]['cell_distance_km'], 10)
        self.assertTrue(any(a['result'] == 'repeated_cell' for a in attempts))

    def test_actual_center_bound_rejects_nominal_radius_hit(self):
        p = (.025, .025)
        east = m.bio.cell_index(*m.offset_point(*p, 10000, 90))
        selected, attempts = self.search(p, FakeSource({east: [(1, 2), (3, 4)]}), bound=10000)
        self.assertNotIn(p, selected)
        self.assertTrue(any(a['result'] == 'outside_cell_center_distance_bound' for a in attempts))

    def test_both_decades_and_both_variables_required(self):
        p = (.025, .025)
        north = m.bio.cell_index(*m.offset_point(*p, 5000, 0))
        east = m.bio.cell_index(*m.offset_point(*p, 5000, 90))
        selected, _ = self.search(p, FakeSource({north: [(1, 2), (3, None)], east: [(10, 2), (20, 4)]}))
        self.assertEqual(selected[p]['direction'], 'E')

    def test_offsets_wrap_dateline_and_keep_geodesic_distance(self):
        p = (20, 179.99)
        target = m.offset_point(*p, 10000, 90)
        self.assertLess(target[1], -179)
        self.assertAlmostEqual(m.bio.distance_km(p, target), 10, places=7)

    def test_run_preserves_old_values_and_warnings_while_filling_na(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fields = ['run_accession', 'lat', 'lon', 'original', *m.bio.COLUMNS]
            rows = [dict(zip(fields, v)) for v in [('A', '0.025', '0.025', '0001', '1.2300', '2.000'), ('B', '0.025', '0.025', 'NA', 'NA', 'NA'), ('C', '', '', 'quoted,"', 'NA', 'NA')]]
            m.bio.write_csv(root/'input.csv', fields, rows)
            prior = []
            for r in rows:
                d = dict(r, status='assigned_native' if r['run_accession']=='A' else 'held_existing_coordinate_warning', reason='original warning', cell_lat='', cell_lon='', distance_km='')
                for v in m.bio.VARIABLES:
                    for decade in ['2000_2009', '2010_2019']:d[v+'_'+decade] = ''
                prior.append(d)
            m.bio.write_csv(root/'prior.csv', list(prior[0]), prior)
            args = SimpleNamespace(base=root/'input.csv', previous_audit=root/'prior.csv', out=root/'output.csv', audit_dir=root/'audit', cache_dir=root/'cache', offline=True, workers=1, radial_steps_m='5000,10000', max_anchor_distance_m=20000)
            source = FakeSource({m.bio.cell_index(.025,.025): [(10,2),(20,4)]})
            with patch.object(m.bio,'fetch',return_value=(b'{}',{'url':'fixture'})), patch.object(m.bio,'validate_source'), patch.object(m,'CellSource',return_value=source), redirect_stdout(io.StringIO()):
                m.run(args)
            new_fields, new = m.bio.read_metadata(args.out)
            self.assertEqual(new_fields, fields)
            self.assertEqual(new[0], rows[0])
            self.assertEqual(new[2], rows[2])
            self.assertEqual(new[1][m.bio.COLUMNS[0]], '15')
            self.assertEqual(new[1]['original'], 'NA')
            with (root/'audit/sample_salinity_audit.csv').open() as h: audit=list(csv.DictReader(h))
            self.assertEqual(audit[1]['previous_reason'],'original warning')
            self.assertEqual(audit[1]['status'],'assigned_recovered_native')
            bad = list(prior);bad[1]=dict(bad[1],lat='.026')
            m.bio.write_csv(root/'prior.csv',list(prior[0]),bad)
            with self.assertRaisesRegex(ValueError,'Stale'):
                m.run(args)

    def test_partial_pairs_fail_instead_of_mixing_cells(self):
        with self.assertRaises(ValueError):
            m.missing_pair(dict(zip(m.bio.COLUMNS, ['1', 'NA'])))


if __name__ == '__main__':
    unittest.main()
