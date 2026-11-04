"""Biological and bookkeeping invariants for repertoire-wide reference mapping."""
import importlib.util
import gzip
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('broad_annotation',
    ROOT / 'scripts/10_functional_network/10O_EXTEND_FUNCTIONAL_ANNOTATION.py')
broad = importlib.util.module_from_spec(spec)
spec.loader.exec_module(broad)


class BroadAnnotationTests(unittest.TestCase):
    def test_shared_KO_retains_multiple_paths_and_disease_context(self):
        tree = {'name': 'ko00001', 'children': [
            {'name': root, 'children': [{'name': cat, 'children': [
                {'name': term, 'children': [{'name': 'K00001 example'}]}]}]}
            for root, cat, term in [
                ('09100 Metabolism', '09101 Carbohydrate metabolism', '00010 Glycolysis [PATH:ko00010]'),
                ('09100 Metabolism', '09102 Energy metabolism', '00680 Methane metabolism [PATH:ko00680]'),
                ('09160 Human Diseases', '09161 Cancer', '05200 Cancer [PATH:ko05200]')]]}
        frame = broad.hierarchy_rows(tree)
        self.assertEqual(len(frame), 3)
        self.assertEqual(broad.biological_pathways(frame), {'ko00010', 'ko00680'})

    def test_duplicate_memberships_do_not_inflate_module_counts(self):
        frame = pd.DataFrame([dict(KO='K00001', scope='protein_family', category='families', term='transport')]*3)
        out = broad.summarize_membership(frame, {'M01': {'K00001', 'K00002'}})
        self.assertEqual(out.iloc[0].n_KOs, 1)
        self.assertEqual(out.iloc[0].fraction_module_KOs, .5)

    def test_pathway_map_and_ko_links_are_one_membership(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'links.tsv'
            p.write_text('ko:K00001\tpath:map00010\nko:K00001\tpath:ko00010\n')
            out = broad.link_rows(p, 'pathway_id')
            self.assertEqual(len(out), 1)
            self.assertEqual(out.iloc[0].pathway_id, 'ko00010')

    def test_logical_blocks_and_same_sample_detection(self):
        entries = {'M00001': {'NAME': 'test', 'DEFINITION': '(K00001,K00002)+K00003-K00004 K00005'}}
        compiled, audit = broad.compile_modules(entries)
        values = {'K00001': np.array([1, 0]), 'K00002': np.array([0, 1]),
                  'K00003': np.array([1, 1]), 'K00005': np.array([0, 1])}
        state, found, unknown, n = broad.evidence_states(compiled['M00001'], values, 2)
        np.testing.assert_array_equal(state, [0, 1])
        np.testing.assert_array_equal(found, [1, 2])
        self.assertEqual(n, 2)
        np.testing.assert_array_equal(unknown, [0, 0])

    def test_missing_rows_are_unknown_not_zero(self):
        compiled, _ = broad.compile_modules({'M00001': {'DEFINITION': 'K00001 K00002'}})
        state, found, unknown, n = broad.evidence_states(compiled['M00001'], {'K00001': np.array([1])}, 1)
        self.assertEqual(state[0], 2)
        self.assertEqual((found[0], unknown[0], n), (1, 1, 2))

    def test_standalone_optional_block_does_not_inflate_coverage(self):
        compiled, _ = broad.compile_modules({'M00001': {'DEFINITION': 'K00001 -K00002'}})
        state, found, unknown, n = broad.evidence_states(compiled['M00001'], {'K00001': np.array([0])}, 1)
        self.assertEqual(state[0], 0)
        self.assertEqual((found[0], unknown[0], n), (0, 0, 1))

    def test_unassigned_reference_steps_and_cycles_are_not_invented(self):
        entries = {'M00001': {'DEFINITION': '-- K00001'},
                   'M00002': {'DEFINITION': 'M00003'}, 'M00003': {'DEFINITION': 'M00002'}}
        compiled, audit = broad.compile_modules(entries)
        self.assertFalse(compiled)
        self.assertTrue(audit.evaluation_status.eq('not_evaluated_definition').all())

    def test_nested_module_expansion_preserves_AND_OR(self):
        entries = {'M00001': {'DEFINITION': 'K00001+K00002'},
                   'M00002': {'DEFINITION': 'M00001,K00003'}}
        compiled, _ = broad.compile_modules(entries)
        state, *_ = broad.evidence_states(compiled['M00002'], {'K00003': np.array([1])}, 1)
        self.assertEqual(state[0], 1)

    def test_reference_tampering_stops_annotation(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'data.tsv').write_text('original')
            broad.write_json(p/'manifest.json', {'files': {'data.tsv': {'sha256': broad.digest(p/'data.tsv')}}})
            (p/'data.tsv').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                broad.verify_reference(p)

    def test_historical_function_is_not_silently_made_unknown(self):
        self.assertEqual(broad.description_status('dnaA; replication protein', False), 'historical_description_only')
        self.assertEqual(broad.description_status('uncharacterized protein', True), 'explicitly_uncharacterized_description')

    def test_compressed_exports_are_complete_and_round_trip(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'evidence.tsv.gz'
            expected = pd.DataFrame({'KO':['K00001']*10000, 'state':np.arange(10000)%3})
            broad.write_tsv(expected, path)
            self.assertTrue(gzip.decompress(path.read_bytes()).endswith(b'\n'))
            pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'), expected)
            self.assertEqual(list(Path(d).glob('.writing-*')), [])


if __name__ == '__main__':
    unittest.main()
