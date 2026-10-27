#!/usr/bin/env python3
"""Build a verified code-only ZIP from the current checkout or source release."""
import argparse
import hashlib
from pathlib import Path
import zipfile
from repository_files import ROOT, repository_files
from validate_repository import validate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--execute', action='store_true', help='Write the archive; default only previews.')
    args = p.parse_args()
    errors = validate()
    if errors:
        p.error('\n'.join(errors))
    files = repository_files()
    print(f'{len(files)} verified source files -> {args.output}')
    if not args.execute:
        print('PREVIEW ONLY: add --execute to write the code archive.')
        return
    if args.output.exists():
        p.error('Output already exists; choose another destination.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, Path('mangrove-sediment-biogeography') / path.relative_to(ROOT))
    with zipfile.ZipFile(args.output) as archive:
        assert archive.testzip() is None
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix(args.output.suffix + '.sha256').write_text(f'{digest}  {args.output.name}\n')
    print(f'CODE_RELEASE_PASS: SHA-256 {digest}')


if __name__ == '__main__':
    main()
