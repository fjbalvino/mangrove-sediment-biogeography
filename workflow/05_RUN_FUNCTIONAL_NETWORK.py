#!/usr/bin/env python3
"""Rerun the KO consensus network without imposing historical network counts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--skip-study-sensitivity', action='store_true')
    p.add_argument('--counts', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--atlas-classes', type=Path, required=True)
    p.add_argument('--annotations', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--start-at', choices=['prepare', 'B', 'A', 'C', 'D', 'E', 'F', 'G'], default='prepare')
    p.add_argument('--stop-after', choices=['B', 'A', 'C', 'D', 'E', 'F', 'G'], default='F')
    args = p.parse_args()
    repo = Path(__file__).resolve().parents[1]
    # Do not mix this frozen network with another cohort or annotation snapshot.
    expected = json.loads((repo / 'config/network_input_checksums.json').read_text())
    source_inputs = {
        'functional_KO_counts_original.tsv.gz': args.counts,
        'metadata_functional.csv': args.metadata,
        'KO_persistence_classes.tsv': args.atlas_classes,
        '137_KO_function_dictionary_used.tsv': args.annotations / '137_KO_function_dictionary_used.tsv',
        '144_KO_to_KEGG_MODULE_long.tsv': args.annotations / '144_KO_to_KEGG_MODULE_long.tsv',
    }
    input_verification = {}
    for name, path in source_inputs.items():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected[name]['sha256']:
            raise ValueError(f'Input checksum mismatch: {name}; verify the source data before rerunning.')
        input_verification[name] = {'path': str(path.resolve()), 'sha256': digest, 'status': 'PASS'}
    scripts = repo / 'scripts/10_functional_network'
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / 'source_input_verification.json').write_text(json.dumps(input_verification, indent=2) + '\n')
    logs = out / 'execution_logs'
    logs.mkdir(exist_ok=True)
    inp = out / 'inputs'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', TZ='UTC')
    history = out / 'commands.jsonl'

    def run(stage, command):
        started = datetime.now(timezone.utc).isoformat()
        print(f'{started} START {stage}', flush=True)
        log = logs / (stage + '_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '.log')
        with log.open('w') as f:
            proc = subprocess.run([str(v) for v in command], env=env, stdout=f, stderr=subprocess.STDOUT)
        with history.open('a') as f:
            f.write(json.dumps({'stage': stage, 'started_utc': started,
                                'finished_utc': datetime.now(timezone.utc).isoformat(),
                                'command': [str(v) for v in command],
                                'exit_code': proc.returncode, 'log': str(log)}) + '\n')
        print(f'END {stage}: exit={proc.returncode}; {log}', flush=True)
        if proc.returncode:
            print(log.read_text()[-6000:], flush=True)
            raise SystemExit(proc.returncode)

    run('runtime', [args.rscript, '--vanilla', '-e',
        'stopifnot(as.character(getRversion()) == "4.3.3"); '
        'stopifnot(as.character(packageVersion("igraph")) == "2.2.1"); '
        'stopifnot(as.character(packageVersion("data.table")) == "1.17.8"); '
        'f <- tempfile(fileext=".tsv.gz"); data.table::fwrite(data.frame(x=1), f); '
        'stopifnot(data.table::fread(f)$x == 1); unlink(f); '
        'sessionInfo()'])
    if args.start_at == 'prepare':
        run('prepare', [sys.executable, scripts / '10H_PREPARE_NETWORK_COHORT.py',
                       '--counts', args.counts, '--metadata', args.metadata,
                       '--atlas-classes', args.atlas_classes,
                       '--accessions', repo / 'config/samples/functional_accessions.txt',
                       '--outdir', inp])
    audit = json.loads((inp / 'cohort_audit.json').read_text())
    common = ['--base', out]

    def latest(number):
        found = sorted((out / 'resultados').glob('LATEST_' + number + '_*.txt'))
        paths = {f.read_text().strip() for f in found}
        if len(paths) != 1:
            raise ValueError(f'Ambiguous/missing stage {number}: {paths}')
        return Path(paths.pop())

    stages = ['B', 'A', 'C', 'D', 'E', 'F', 'G']
    start = stages.index('B' if args.start_at == 'prepare' else args.start_at)
    for stage in stages[start:stages.index(args.stop_after) + 1]:
        extra = []
        if stage == 'B':
            extra = ['--ko_matrix', args.counts, '--metadata', inp / 'metadata_191.tsv',
                     '--metadata_id_col', 'run_accession', '--missing_as_zero', 0,
                     '--strict_expected', 1, '--expected_n_samples', audit['n_samples'],
                     '--expected_n_kos', audit['n_KOs'], '--expected_n_assoc', audit['n_association_KOs'],
                     '--expected_n_network', 5000]
            for prefix, key in [('all', 'all_prevalence_classes'), ('assoc', 'association_prevalence_classes')]:
                for cls in ['Prevalent', 'Shell', 'Cloud']:
                    extra += ['--expected_' + prefix + '_' + cls.lower(), audit[key][cls]]
        elif stage == 'A':
            extra = ['--ko_table', latest('156') / 'tables/156_KO_network_filter_stats.tsv',
                     '--ko_dictionary', args.annotations / '137_KO_function_dictionary_used.tsv',
                     '--core_classes', inp / 'KO_core_classes.tsv',
                     '--kegg_snapshot', args.annotations / '144_KO_to_KEGG_MODULE_long.tsv',
                     '--allow_kegg_query', 0, '--expected_n_kos', 5000, '--strict_expected', 1]
        elif stage == 'C':
            extra = ['--run156', latest('156'), '--ko_core_classes', inp / 'KO_core_classes.tsv',
                     '--ko_dict', args.annotations / '137_KO_function_dictionary_used.tsv',
                     '--realm_col', 'REALM', '--country_col', 'Country', '--locality_col', 'Locality',
                     '--strict_expected', 1, '--enforce_reference', 0,
                     '--expected_n_samples', audit['n_samples'], '--expected_total_kos', audit['n_KOs'],
                     '--expected_bootstrap_median_n', audit['samples_per_balanced_subsample'],
                     '--seed', 123, '--n_boot', 500, '--max_kos', 1500,
                     '--max_n_per_stratum', 4, '--target_edges_per_bootstrap', 2000,
                     '--consensus_stability_min', .80, '--sign_consistency_min', 1,
                     '--consensus_mean_abs_min', .95, '--min_module_size', 10,
                     '--community_method', 'leiden', '--resolution', 1, '--n_starts', 50, '--partition_seed', 2001]
        elif stage == 'D':
            extra = ['--run171', latest('171'), '--strict_expected', 0, '--seed', 123, '--dpi', 300]
        elif stage == 'E':
            extra = ['--run172', latest('172'), '--run144_reference', latest('144'),
                     '--strict_expected', 0, '--make_figures', 0, '--dpi', 300]
        elif stage == 'F':
            extra = ['--run156', latest('156'), '--run171', latest('171'), '--run173', latest('173'),
                     '--strict_expected', 0, '--enforce_manuscript_reference', 0,
                     '--n_perm', 4999, '--n_rewire', 999, '--n_edge_boot', 500, '--workers', args.workers,
                     '--kegg_mapping_complete', 0,
                     '--seed', 123, '--dpi', 300]
        elif stage == 'G':
            extra = ['--run171', latest('171'), '--run172', latest('172'),
                     '--run173', latest('173'), '--run174', latest('174'),
                     '--outdir', out / 'main_figure', '--strict_expected', 1, '--dpi', 600]
        if stage == 'G' and args.start_at != 'G' and not args.skip_study_sensitivity:
            study = out / 'study_sensitivity'
            run('J_study_networks', [sys.executable, scripts / '10J_NETWORK_STUDY_SENSITIVITY.py',
                '--run156', latest('156'), '--run171', latest('171'),
                '--metadata', inp / 'metadata_191.tsv', '--outdir', study,
                '--rscript', args.rscript, '--workers', args.workers])
            run('J_study_partitions', [args.rscript, '--vanilla', scripts / 'lib/study_sensitivity.R',
                'evaluate', study, latest('173'), latest('171')])
            run('K_verify_study_sensitivity', [sys.executable, scripts / '10K_VERIFY_NETWORK_SENSITIVITY.py',
                '--study', study, '--run156', latest('156'), '--rscript', args.rscript])
        script, = scripts.glob('10' + stage + '_*.R')
        run(stage, [args.rscript, '--vanilla', script, *common, *extra])
    if args.stop_after in ['F', 'G']:
        run('verify', [sys.executable, scripts / '10I_VERIFY_NETWORK_OUTPUTS.py',
                       '--outdir', out])


if __name__ == '__main__':
    main()
