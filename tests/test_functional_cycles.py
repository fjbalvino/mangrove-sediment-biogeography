"""Synthetic checks of ternary logic, missing KO rows and portable cycle runs."""
import importlib.util
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BLOCK = ROOT / 'scripts/11_biogeochemical_cycles'
spec = importlib.util.spec_from_file_location('cycle_rules', BLOCK / '11A_EVALUATE_BIOGEOCHEMICAL_RULES.py')
rules = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rules)
spec_verify = importlib.util.spec_from_file_location('cycle_verify', BLOCK / '11C_VERIFY_BIOGEOCHEMICAL_OUTPUTS.py')
verifier = importlib.util.module_from_spec(spec_verify)
spec_verify.loader.exec_module(verifier)


class TestTernaryLogic(unittest.TestCase):
    def test_truth_tables(self):
        for a, b in itertools.product([0, 1, 2], repeat=2):
            values = {'K00001': np.array([a]), 'K00002': np.array([b])}
            conjunction = 0 if 0 in (a, b) else 2 if 2 in (a, b) else 1
            disjunction = 1 if 1 in (a, b) else 2 if 2 in (a, b) else 0
            self.assertEqual(rules.evaluate(rules.parse('K00001+K00002'), values, 1)[0], conjunction)
            self.assertEqual(rules.evaluate(rules.parse('K00001,K00002'), values, 1)[0], disjunction)

    def test_missing_optional_dependency_and_invalid_syntax(self):
        self.assertGreaterEqual(rules.self_tests(), 44)
        self.assertEqual(rules.evaluate(rules.parse('K00001'), {}, 1)[0], 2)
        with self.assertRaisesRegex(ValueError, 'Recursive'):
            rules.parse('M00001', {'M00001': 'M00002', 'M00002': 'M00001'})
        with self.assertRaisesRegex(ValueError, 'Unresolved'):
            rules.parse('M00001', {})
        for expression in ['K00001,,K00002', ',K00001', 'K00001,']:
            with self.assertRaisesRegex(ValueError, 'Empty alternative'):
                rules.parse(expression)

    def test_matrix_rejects_duplicate_headers_and_noncounts(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'counts.tsv'
            for content in ['KO\ta\ta\nK00001\t1\t2\n', 'KO\ta\nK00001\t-1\n',
                            'KO\ta\nK00001\t0.5\n', 'KO\ta\nK00001\tNA\n',
                            'KO\ta\nK00001\t1\nK00001\t2\n']:
                path.write_text(content)
                with self.assertRaises(ValueError):
                    rules.read_counts(path, ['a'])


class TestCycleWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.folder = Path(cls.temp.name)
        cls.matrix = cls.folder / 'counts.tsv'
        # Deliberately omit nifK and all MCR/Nrf/ACL rows to exercise UNKNOWN handling.
        cls.matrix.write_text('KO\ts1\ts2\ts3\nK02588\t2\t1\t0\nK02586\t2\t1\t0\nK00001\t0\t0\t0\n')
        cls.meta = cls.folder / 'metadata.csv'
        cls.meta.write_text('run_accession,study_id,sampling_site_id,REALM,lat,lon\ns1,A,site1,R1,1,2\ns2,B,site2,R2,,\ns3,B,site2,R2,,\n')
        cls.samples = cls.folder / 'samples.txt'
        cls.samples.write_text('s3\ns1\ns2\n')
        cls.output = cls.folder / 'output'
        cls.command = [sys.executable, str(ROOT / 'workflow/07_RUN_BIOGEOCHEMICAL_CYCLES.py'),
                       '--matrix', str(cls.matrix), '--metadata', str(cls.meta),
                       '--samples', str(cls.samples), '--unit', 'predicted_protein_records',
                       '--outdir', str(cls.output)]
        cls.before = [p.read_bytes() for p in (cls.matrix, cls.meta, cls.samples)]
        preview = subprocess.run(cls.command, cwd=cls.folder, capture_output=True, text=True, check=True)
        if cls.output.exists():
            raise AssertionError('Preview created outputs')
        if 'Preview only' not in preview.stdout:
            raise AssertionError('Preview status not reported')
        subprocess.run(cls.command + ['--execute'], cwd=cls.folder, capture_output=True, text=True, check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_missing_rows_are_unknown_not_zero_and_order_is_preserved(self):
        matrix = pd.read_csv(self.output / 'tables/matriz_97_reglas_ge1.tsv', sep='\t', index_col=0)
        self.assertEqual(list(matrix.columns), ['s3', 's1', 's2'])
        self.assertEqual(matrix.loc['N01_nifHDK', 's3'], 0)
        self.assertTrue(matrix.loc['N01_nifHDK', ['s1', 's2']].isna().all())
        self.assertTrue(matrix.loc['M01_MCR'].isna().all())
        self.assertEqual(self.before, [p.read_bytes() for p in (self.matrix, self.meta, self.samples)])

    def test_dimensions_integrity_and_hashes(self):
        diag = json.loads((self.output / 'tables/diagnostico_y_validacion.json').read_text())
        self.assertEqual((diag['n_muestras'], diag['n_reglas'], diag['n_procesos']), (3, 97, 65))
        self.assertEqual(diag['n_modulos_referencia'], 32)
        self.assertEqual(diag['comparaciones_independientes'], 8)
        self.assertIn('config/module_dependencies.json', diag['sha256'])
        status = json.loads((self.output / 'RUN_STATUS.json').read_text())
        self.assertEqual(status['status'], 'PASS')
        self.assertIn('tables/matriz_97_reglas_ge1.tsv', status['sha256_outputs'])

    def test_refuses_existing_output_directory(self):
        result = subprocess.run(self.command + ['--execute'], cwd=self.folder, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not empty', result.stderr)

    def test_verifier_detects_corrupted_summary(self):
        path = self.output / 'tables/resumen_reglas.tsv'
        original = path.read_bytes()
        try:
            frame = pd.read_csv(path, sep='\t')
            frame.loc[0, 'n_cumple'] += 1
            frame.to_csv(path, sep='\t', index=False)
            with self.assertRaises(AssertionError):
                verifier.verify(self.output / 'tables')
        finally:
            path.write_bytes(original)


if __name__ == '__main__':
    unittest.main()
