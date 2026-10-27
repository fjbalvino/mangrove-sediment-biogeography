#!/usr/bin/env python3
"""Validate a code-only checkout/release and its source-file checksums."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from repository_files import ROOT, repository_files

FORBIDDEN_SUFFIXES = {
    '.fastq', '.fq', '.fasta', '.fa', '.faa', '.fna', '.rds', '.rdata', '.rda',
    '.biom', '.k2d', '.npz', '.npy', '.parquet', '.sqlite', '.sqlite3', '.db',
    '.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg', '.pdf', '.docx', '.pptx',
    '.xlsx', '.xls', '.zip', '.gz', '.tar', '.bz2', '.xz', '.tif', '.tiff',
    '.geojson', '.gpkg', '.shp', '.shx', '.dbf', '.nc', '.newick', '.nwk',
    '.csv', '.log', '.env', '.pem', '.key',
}
OUTPUT_DIRECTORIES = {'data', 'results', 'figures', 'resultados', 'cache'}
REQUIRED = {'README.md', 'LICENSE', 'CITATION.cff', 'DATA_AVAILABILITY.md',
            'docs/provenance/SCRIPT_INDEX.tsv', 'docs/MANUSCRIPT_GUIDE.md',
            'workflow/00_README_HOW_TO_RUN.md', 'config/project.env.example',
            '.zenodo.json', 'MANIFEST.sha256'}
# Match host-specific endpoints and personal home paths; variable paths are valid.
PRIVATE_PATH = re.compile(r'(?i)(?:https?://[^\s/]*landa\w*|/home/[A-Za-z][\w.-]*/|/data/(?:Magenta|BD)(?:/|\b)|/data/' + r'bin\b)')


def validate(check_manifest: bool = True) -> list[str]:
    errors = []
    files = repository_files()
    names = {path.relative_to(ROOT).as_posix() for path in files}
    for name in sorted(REQUIRED - names):
        errors.append(f'Required source file is missing: {name}')
    for path in files:
        name = path.relative_to(ROOT).as_posix()
        if not path.is_file():
            errors.append(f'Manifest source is missing: {name}')
            continue
        if path.is_symlink():
            errors.append(f'Symlinks are not permitted in a source release: {name}')
            continue
        if path.stat().st_size > 20 * 1024 * 1024:
            errors.append(f'File exceeds 20 MiB: {name}')
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f'Data, figure, archive or local configuration in code release: {name}')
        if name.split('/')[0] in OUTPUT_DIRECTORIES and name not in {'data/README.md', 'results/README.md'}:
            errors.append(f'Generated or input file in code release: {name}')
        if path.suffix == '.tsv' and not (name.startswith('config/') or name == 'docs/provenance/SCRIPT_INDEX.tsv'):
            errors.append(f'Numerical tables belong in the data companion: {name}')
        try:
            content = path.read_text(encoding='utf-8')
            if '\0' in content:
                errors.append(f'Binary content in code release: {name}')
            if PRIVATE_PATH.search(content):
                errors.append(f'Institution-specific endpoint or path: {name}')
        except UnicodeError:
            errors.append(f'Non-text file in code release: {name}')
    index = ROOT / 'docs/provenance/SCRIPT_INDEX.tsv'
    if index.is_file():
        with index.open(newline='', encoding='utf-8') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                name = row.get('canonical_file', '')
                if name not in names:
                    errors.append(f'Indexed canonical file is missing: {name}')
    try:
        json.loads((ROOT / '.zenodo.json').read_text())
    except (OSError, ValueError) as exc:
        errors.append(f'Invalid software archive metadata: {exc}')
    if check_manifest and (ROOT / 'MANIFEST.sha256').is_file():
        expected = {}
        for line in (ROOT / 'MANIFEST.sha256').read_text().splitlines():
            try:
                sha, name = line.split('  ', 1)
                if not re.fullmatch('[a-f0-9]{64}', sha) or name in expected:
                    raise ValueError('invalid/duplicate digest entry')
                expected[name] = sha
            except ValueError:
                errors.append(f'Invalid manifest entry: {line}')
        actual_names = names - {'MANIFEST.sha256'}
        if set(expected) != actual_names:
            errors.append('Manifest inventory differs from release files; run make manifest.')
        for name in sorted(set(expected) & actual_names):
            if (ROOT / name).is_file() and hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected[name]:
                errors.append(f'Manifest digest differs: {name}')
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-manifest', action='store_true', help='Content checks while editing; not a release check.')
    args = parser.parse_args()
    errors = validate(check_manifest=not args.skip_manifest)
    print('REPOSITORY_VALIDATION_FAIL' if errors else f'REPOSITORY_VALIDATION_PASS: {len(repository_files())} source files')
    for error in errors:
        print(f'- {error}')
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
