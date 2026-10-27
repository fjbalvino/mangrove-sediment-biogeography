"""Smoke-test real entry points without biological inputs or external jobs."""
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'tools'))
from repository_files import repository_files


class WorkflowPortabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='magenta portability ')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.repo = self.base / 'source with spaces'
        for path in repository_files(REPO):
            dest = self.repo / path.relative_to(REPO)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)
        self.env = dict(os.environ)
        for key in list(self.env):
            if key.startswith(('MAGENTA_', 'BIOGEOGRAPHY_', 'FUNCTIONAL_', 'ENVIRONMENTAL_')):
                del self.env[key]
        self.env['PYTHONDONTWRITEBYTECODE'] = '1'

    def run_command(self, *args, success=True):
        result = subprocess.run([str(v) for v in args], cwd=self.base, env=self.env,
                                capture_output=True, text=True, timeout=40)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def test_shell_previews_from_another_directory_leave_no_outputs(self):
        data = self.base / 'external data'
        outputs = self.base / 'external outputs'
        config = self.base / 'custom settings.env'
        config.write_text(f'export MAGENTA_DATA_DIR={shlex.quote(str(data))}\n'
                          f'export BIOGEOGRAPHY_RESULTS={shlex.quote(str(outputs))}\n')
        for number in ('01', '02', '03', '04', '06'):
            script, = self.repo.glob(f'workflow/{number}_*.sh')
            with self.subTest(workflow=number):
                result = self.run_command('bash', script, '--config', config)
                self.assertTrue(result.stdout.strip())
                self.assertFalse(data.exists())
                self.assertFalse(outputs.exists())
                self.assertFalse((self.repo / 'cache').exists())
        missing = self.run_command('bash', self.repo / 'workflow/01_RUN_TAXONOMIC_PREPROCESSING.sh',
                                   '--config', self.base / 'missing.env', success=False)
        self.assertIn('Configuration file not found', missing.stderr)

    def test_network_preview_does_not_require_inputs_or_create_outputs(self):
        output = self.base / 'network output'
        result = self.run_command(sys.executable, self.repo / 'workflow/05_RUN_FUNCTIONAL_NETWORK.py',
                                  '--counts', self.base / 'missing.tsv.gz',
                                  '--metadata', self.base / 'missing.csv',
                                  '--atlas-classes', self.base / 'missing.tsv',
                                  '--annotations', self.base / 'missing annotations',
                                  '--outdir', output, '--stop-after', 'G', '--dry-run')
        plan = json.loads(result.stdout)
        self.assertEqual(plan['stages'], ['prepare', 'B', 'A', 'C', 'D', 'E', 'F', 'G'])
        self.assertIn('PREVIEW ONLY', plan['status'])
        self.assertFalse(output.exists())

    def test_source_export_excludes_local_assets_and_detects_tampering(self):
        # A release extraction has no Git metadata; its checksum inventory
        # must remain authoritative even after users add local input/output files.
        for name in ('data/private.csv', 'results/summary.tsv', 'figures/panel.png', 'config/project.env'):
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('local-only asset\n')
        output = self.base / 'release.zip'
        self.run_command(sys.executable, self.repo / 'tools/build_code_release.py',
                         '--output', output, '--execute')
        with zipfile.ZipFile(output) as archive:
            names = archive.namelist()
            self.assertTrue(any(n.endswith('/README.md') for n in names))
            for forbidden in ('private.csv', 'summary.tsv', 'panel.png', 'project.env'):
                self.assertFalse(any(n.endswith('/' + forbidden) for n in names))
            self.assertIsNone(archive.testzip())
        readme = self.repo / 'README.md'
        readme.write_text(readme.read_text() + '\nUnrecorded edit.\n')
        result = self.run_command(sys.executable, self.repo / 'tools/validate_repository.py', success=False)
        self.assertIn('Manifest digest differs: README.md', result.stdout)


if __name__ == '__main__':
    unittest.main()
