"""Check that the static export keeps input matrices outside its payload."""
from pathlib import Path
import importlib.util
import json
import tempfile
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('exporter',ROOT/'workflow/09_EXPORT_ATLAS_SITE.py')
exporter=importlib.util.module_from_spec(spec);spec.loader.exec_module(exporter)


class AtlasExport(unittest.TestCase):
    def test_only_reviewed_site_files_and_links(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);atlas=root/'atlas';(atlas/'results').mkdir(parents=True)
            (atlas/'assets/vendor').mkdir(parents=True)
            (atlas/'assets/vendor/leaflet-LICENSE.txt').write_text('Synthetic license fixture')
            viewers={str(i):f'figures/cycle{i}/viewer.html' for i in range(8)}
            for name in viewers.values():
                path=atlas/name;path.parent.mkdir(parents=True);path.write_text('<html><body>synthetic view</body></html>')
            (atlas/'index.html').write_text(''.join(f'<a href="{p}">Cycle</a>' for p in viewers.values()))
            (atlas/'results/viewer_manifest.json').write_text(json.dumps({'profiles':16,'viewers':viewers}))
            (atlas/'results/secret_matrix.tsv').write_text('do not export')
            output=root/'site';exporter.export(atlas,output)
            with zipfile.ZipFile(root/'site.zip') as z:
                self.assertEqual(len([n for n in z.namelist() if n.endswith('.html')]),9)
                self.assertFalse(any('secret' in n or n.endswith('.tsv') for n in z.namelist()))
            (output/'index.html').write_text('<a href="../secret_matrix.tsv">bad</a>')
            with self.assertRaisesRegex(ValueError,'leaves website'):exporter.validate_site(output)


if __name__=='__main__':unittest.main()
