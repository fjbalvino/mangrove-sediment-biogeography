"""Preview or execute the curated community-level biogeochemical evidence workflow."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BLOCK = ROOT / 'scripts/11_biogeochemical_cycles'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--samples', default=str(ROOT / 'config/samples/functional_accessions.txt'))
    parser.add_argument('--unit', required=True, choices=['predicted_protein_records'])
    parser.add_argument('--outdir', default=str(ROOT / 'results/biogeochemical_cycles'))
    parser.add_argument('--ko-labels', help='Optional external TSV containing KO and name')
    parser.add_argument('--plots', action='store_true', help='Generate figures locally after verification')
    parser.add_argument('--execute', action='store_true', help='Without this flag, print commands only')
    args = parser.parse_args()
    output = Path(args.outdir).resolve()
    if output.is_relative_to(ROOT) and not output.is_relative_to(ROOT / 'results'):
        parser.error('Outputs inside the repository must be under results/')
    commands = [[sys.executable, str(BLOCK / '11A_EVALUATE_BIOGEOCHEMICAL_RULES.py'),
                 '--matrix', args.matrix, '--metadata', args.metadata, '--samples', args.samples,
                 '--unit', args.unit, '--out', str(output / 'tables')]]
    if args.ko_labels:
        commands[0] += ['--ko-labels', args.ko_labels]
    commands.append([sys.executable, str(BLOCK / '11C_VERIFY_BIOGEOCHEMICAL_OUTPUTS.py'),
                     '--results', str(output / 'tables')])
    if args.plots:
        commands.append([sys.executable, str(BLOCK / '11B_PLOT_BIOGEOCHEMICAL_EVIDENCE.py'),
                         '--results', str(output / 'tables'), '--figures', str(output / 'figures')])
    for command in commands:
        print(shlex.join(command), flush=True)
    if not args.execute:
        print('Preview only; add --execute to run. No files were created.')
        return
    if output.exists() and any(output.iterdir()):
        parser.error('Output directory is not empty; use a new directory for each input release')
    for path in [args.matrix, args.metadata, args.samples] + ([args.ko_labels] if args.ko_labels else []):
        if not Path(path).is_file():
            parser.error('Missing input: ' + path)
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', TZ='UTC')
    status = {'status': 'FAIL', 'unit': args.unit, 'commands': commands}
    try:
        with (output / 'execution.log').open('w') as log:
            for command in commands:
                subprocess.run(command, env=environment, check=True, stdout=log, stderr=subprocess.STDOUT)
        status['status'] = 'PASS'
    finally:
        status['sha256_outputs'] = {
            str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'RUN_STATUS.json'
        }
        (output / 'RUN_STATUS.json').write_text(json.dumps(status, indent=2) + '\n')
        print('BIOGEOCHEMICAL_WORKFLOW_' + status['status'], flush=True)


if __name__ == '__main__':
    main()
