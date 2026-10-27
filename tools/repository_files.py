"""Enumerate release source files without including ignored local inputs/outputs."""
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def repository_files(root: Path = ROOT) -> "list[Path]":
    if (root / '.git').exists():
        names = subprocess.check_output(
            ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'],
            cwd=root,
        ).decode('utf-8').split('\0')
        return sorted({root / name for name in names if name and (root / name).is_file()})
    # An extracted release has no Git index. Its manifest supplies the inventory.
    manifest = root / 'MANIFEST.sha256'
    if not manifest.is_file():
        raise RuntimeError('A Git checkout or MANIFEST.sha256 source release is required.')
    files = [root / line.split('  ', 1)[1] for line in manifest.read_text().splitlines() if line]
    for path in files:
        if root.resolve() not in path.resolve().parents:
            raise ValueError(f'Manifest path leaves the source directory: {path}')
    return sorted(set(files + [manifest]))
