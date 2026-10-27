#!/usr/bin/env python3
"""Create a deterministic SHA-256 manifest for code release files."""
import hashlib
from repository_files import ROOT, repository_files


def main():
    files = [path for path in repository_files() if path.name != 'MANIFEST.sha256']
    lines = [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT).as_posix()}" for path in files]
    output = ROOT / 'MANIFEST.sha256'
    output.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'Wrote {output.name} with {len(lines)} entries')


if __name__ == '__main__':
    main()
