#!/usr/bin/env python3
"""Rebuild cycle evidence and gene–environment diagnostics from external inputs.

Preview is the default. Every input release needs a fresh run directory.
No upload, annotation, environmental download or metadata correction is performed.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / 'scripts/12_gene_environment/template'
CONFIG = ROOT / 'config/gene_biogeography'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(path):
    spec = importlib.util.spec_from_file_location('gene_models', path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def snapshot(root):
    files = [p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    return {p.relative_to(root).as_posix(): digest(p) for p in files}


def check_run(run):
    record = json.loads((run / 'run_manifest.json').read_text())
    for name, checksum in record['protected_files'].items():
        if not (run / name).is_file() or digest(run / name) != checksum:
            raise ValueError('Prepared input/source changed: ' + name + '; create a new run')
    return record


def prepare(args, run):
    import numpy as np
    import pandas as pd
    if run.exists() and any(run.iterdir()):
        raise ValueError('Use a new or empty run directory for every input release')
    matrix, metadata = Path(args.matrix).resolve(), Path(args.metadata).resolve()
    opener = gzip.open if matrix.suffix == '.gz' else open
    with opener(matrix, 'rt', encoding='utf-8-sig') as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    helper = module(TEMPLATE / 'analysis_nitrogen/scripts/12_gene_environment/12B_BUILD_NITROGEN_PILOT.py')
    helper.validate_ids(header[1:], 'matrix headers')
    ids = header[1:] if args.all_matrix_samples else [line for line in Path(args.samples).read_text().splitlines() if line.strip()]
    ids = helper.validate_ids(ids, 'selected accessions')
    if args.expected_profiles is not None and len(ids) != args.expected_profiles:
        raise ValueError('Selected accession count differs from --expected-profiles')
    meta = pd.read_csv(metadata, dtype={'run_accession': str, 'study_id': str, 'sample_accession': str})
    required = {'run_accession','study_id','sample_accession','lat','lon','Country','Locality','REALM','bio1','soc_mean','phh2o_mean'}
    if required - set(meta):
        raise ValueError('Missing metadata fields: ' + ', '.join(sorted(required - set(meta))))
    helper.validate_ids(meta.run_accession.tolist(), 'metadata run_accession')
    if not set(ids) <= set(meta.run_accession):
        raise ValueError('Selected IDs absent from metadata: ' + ', '.join(sorted(set(ids)-set(meta.run_accession))))
    if not set(ids) <= set(header[1:]):
        raise ValueError('Selected IDs absent from matrix')
    mat = pd.read_csv(matrix, sep='\t', index_col=0).loc[:, ids]
    helper.validate_ids(mat.index.tolist(), 'KO rows')
    import re
    if not all(re.fullmatch(r'K\d{5}', k) for k in mat.index):
        raise ValueError('Matrix must have KO identifiers K##### in its first column')
    values = mat.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values<0).any() or not np.equal(values,np.floor(values)).all() or (mat.sum()<=0).any():
        raise ValueError('Expected finite nonnegative predicted-protein record counts and positive profile totals; never round normalized data')
    meta = meta.set_index('run_accession').loc[ids].copy()
    if meta.study_id.isna().any() or meta.study_id.astype(str).str.strip().eq('').any():
        raise ValueError('Every selected run needs a nonempty study_id')
    for column, limit in [('lat',90),('lon',180)]:
        meta[column] = pd.to_numeric(meta[column], errors='raise')
        if not np.isfinite(meta[column]).all() or (meta[column].abs()>limit).any():
            raise ValueError('The geographic atlas requires finite coordinates for every selected run; audit missing coordinates or use workflow 07 for coordinate-free rule evidence')
    for column in ['bio1','soc_mean','phh2o_mean']:
        meta[column] = pd.to_numeric(meta[column], errors='raise')
        if np.isinf(meta[column]).any():
            raise ValueError('Infinite environmental value: ' + column)
    if args.soilgrids_storage == 'physical':
        # Analysis-local compatibility conversion only. Source metadata remain untouched.
        meta[['soc_mean','phh2o_mean']] *= 10
    conflicts = meta.groupby(['study_id','lat','lon'])[['bio1','soc_mean','phh2o_mean']].nunique()
    if (conflicts>1).any().any():
        raise ValueError('Conflicting environmental values at the same study/coordinate')
    shutil.copytree(TEMPLATE, run, dirs_exist_ok=True, ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    inputs = run / 'inputs'; inputs.mkdir(exist_ok=True)
    for name in ['operational_rules.tsv','curated_processes.tsv','unresolved_processes.json']:
        shutil.copy2(CONFIG/name, inputs/name)
    mat.to_csv(inputs/'functional_KO_counts_original.tsv.gz',sep='\t',index_label='KO',compression={'method':'gzip','mtime':0})
    meta.to_csv(inputs/'metadata_functional.csv',index_label='run_accession')
    (inputs/'functional_accessions.txt').write_text('\n'.join(ids)+'\n')
    for folder in ['outputs','analysis_nitrogen/inputs','analysis_nitrogen/config/gene_environment','analysis_remaining/curation','atlas/figures']:
        (run/folder).mkdir(parents=True,exist_ok=True)
    for name in ['metadata_functional.csv','functional_accessions.txt']:
        shutil.copy2(inputs/name,run/'analysis_nitrogen/inputs'/name)
    shutil.copy2(CONFIG/'nitrogen_reactions.json',run/'analysis_nitrogen/config/gene_environment/nitrogen_reactions.json')
    labels = pd.read_csv(CONFIG/'KO_labels.tsv',sep='\t').set_index('KO')
    labels['in_matrix'] = labels.index.isin(mat.index)
    panel = mat.reindex(labels.index)
    labels['n_ge1'] = panel.ge(1).sum(axis=1).where(labels.in_matrix)
    labels['n_ge2'] = panel.ge(2).sum(axis=1).where(labels.in_matrix)
    labels.to_csv(run/'outputs/KO_panel_inventory.tsv',sep='\t')
    inputs_record = {'matrix':{'path':str(matrix),'sha256':digest(matrix)},'metadata':{'path':str(metadata),'sha256':digest(metadata)}}
    if args.samples: inputs_record['accessions']={'path':str(Path(args.samples).resolve()),'sha256':digest(args.samples)}
    protected = snapshot(run)
    record = {'schema_version':1,'measurement':args.unit,'profiles':len(ids),'denominator_KO_rows':len(mat),'soilgrids_storage':args.soilgrids_storage,
              'environment_definition':'WorldClim BIO1 degrees C; SoilGrids 0–5 cm organic carbon and pH; encoded values divided by 10 for analysis',
              'original_inputs':inputs_record,'protected_files':protected,'python':platform.python_version(),
              'source_provenance':json.loads((CONFIG/'source_provenance.json').read_text()),
              'inference':'Exploratory only; no observed-gene p/q values are generated, including after calibration'}
    (run/'run_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def commands(stage, run, args):
    py = sys.executable
    ns = run/'analysis_nitrogen/scripts/12_gene_environment'
    rs = run/'analysis_remaining/scripts'
    if stage == 'prepare':
        return [[py,str(ns/'12B_BUILD_NITROGEN_PILOT.py'),'--metadata',str(run/'inputs/metadata_functional.csv'),'--matrix',str(run/'inputs/functional_KO_counts_original.tsv.gz'),'--samples',str(run/'inputs/functional_accessions.txt'),'--rules',str(run/'inputs/operational_rules.tsv'),'--reactions',str(run/'analysis_nitrogen/config/gene_environment/nitrogen_reactions.json'),'--outdir',str(run/'analysis_nitrogen/results')],
                [py,str(ns/'12D_DESIGN_SUPPORT.py')],[py,str(rs/'01_prepare.py')],[py,str(rs/'06_curate.py')],[py,str(run/'block_evidence.py')],[py,str(run/'atlas/scripts/00_integrate_analyses.py')],[py,str(run/'verify_outputs.py')]]
    if stage == 'atlas':
        if getattr(args,'atlas_cycle',None):
            return [[py,str(run/'atlas/scripts/08_build_atlas.py'),'--cycle',args.atlas_cycle],
                    [py,str(run/'verify_outputs.py'),'--viewers','--cycle',args.atlas_cycle]]
        return [[py,str(run/'atlas/scripts/08_build_atlas.py')],[py,str(run/'plot_supplementary_panels.py')],[py,str(run/'verify_outputs.py'),'--viewers']]
    if stage == 'calibration-plan':
        return [[py,str(ns/'12G_PREPARE_FINAL_CALIBRATION.py')],[py,str(rs/'03_prepare_calibration.py')]]
    if stage == 'calibration-check':
        return [['Rscript',str(ns/'12G_RUN_FINAL_CALIBRATION.R'),str(run/'analysis_nitrogen'),'validate'],[py,str(rs/'04_verify_models.py')]]
    if stage == 'calibrate':
        return [[py,str(ns/'run_final_calibration.py')],[py,str(rs/'05_calibrate.py')],[py,str(rs/'09_summarize.py')]]
    raise ValueError(stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['prepare','atlas','calibration-plan','calibration-check','calibrate'],default='prepare')
    parser.add_argument('--matrix');parser.add_argument('--metadata')
    cohort = parser.add_mutually_exclusive_group()
    cohort.add_argument('--samples',help='Explicit approved accession list; no historical manifest is applied automatically')
    cohort.add_argument('--all-matrix-samples',action='store_true',help='Explicitly include every matrix sample after metadata matching')
    parser.add_argument('--unit',choices=['predicted_protein_records'])
    parser.add_argument('--soilgrids-storage',choices=['encoded','physical'],default='encoded',help='encoded: SOC and pH stored x10; physical: SOC g/kg and pH units. BIO1 must be degrees C in either case')
    parser.add_argument('--expected-profiles',type=int)
    parser.add_argument('--run-dir',required=True,type=Path)
    parser.add_argument('--assets',type=Path,help='Local atlas/assets directory from the reference package; never uploaded')
    parser.add_argument('--atlas-cycle',choices=['N','C','CH4_C1','S','P','Fe','H2_electrons','O2_light'],help='Build a standalone cycle viewer; omit to build the complete atlas')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--execute',action='store_true')
    args = parser.parse_args();run=args.run_dir.resolve()
    if args.workers<1:parser.error('--workers must be positive')
    if ROOT in run.parents and ROOT/'results' not in run.parents:
        parser.error('Repository-local runs must be below results/')
    if args.stage=='prepare' and not (args.matrix and args.metadata and args.unit and (args.samples or args.all_matrix_samples)):
        parser.error('prepare requires --matrix, --metadata, --unit, and either --samples or --all-matrix-samples')
    todo=commands(args.stage,run,args)
    for command in todo: print(shlex.join(command),flush=True)
    if not args.execute:
        print('Preview only. No files created. Add --execute to run.');return
    if not __debug__:parser.error('Do not run with Python -O; scientific assertions must remain active')
    record=prepare(args,run) if args.stage=='prepare' else check_run(run)
    if args.stage=='atlas':
        for name in ['world_110m.geojson','natural_earth_relief.jpg','mangroves_USGS.jpg','nitrogen_reference_template.svg']:
            dest=run/'atlas/assets'/name
            if args.assets:
                original=args.assets/name
                if original.is_symlink():raise ValueError('Symlinked asset: '+name)
                shutil.copy2(original,dest)
            if not dest.is_file():raise ValueError('Missing local visual asset: '+name+'; supply --assets')
        (run/'atlas/assets/asset_checksums.json').write_text(json.dumps({p.name:digest(p) for p in (run/'atlas/assets').iterdir() if p.is_file()},indent=2)+'\n')
    if args.stage=='calibration-plan':
        for panel in ['analysis_nitrogen','analysis_remaining']:
            path=run/panel/'calibration_final'
            if path.exists() and any(path.iterdir()):raise ValueError('Calibration already prepared; keep it frozen or create a new run')
    if args.stage in {'calibration-check','calibrate'}:
        for panel,name in [('analysis_nitrogen','final_simulation_plan.json'),('analysis_remaining','simulation_plan.json')]:
            if not (run/panel/'calibration_final'/name).is_file():
                raise ValueError('Prepare both calibration plans before this stage')
    if args.stage=='calibrate':
        validation=run/'analysis_remaining/calibration_final/implementation_validation.json'
        if not validation.exists() or not json.loads(validation.read_text())['passed']:
            raise ValueError('Run calibration-check successfully before full simulation')
    environment=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',GENE_ATLAS_WORKERS=str(args.workers))
    status={'stage':args.stage,'status':'FAIL','commands':todo}
    try:
        with (run/(args.stage+'.log')).open('w') as log:
            for command in todo:
                subprocess.run(command,cwd=run,env=environment,check=True,stdout=log,stderr=subprocess.STDOUT)
        check_run(run)
        for item in record['original_inputs'].values():
            if Path(item['path']).exists() and digest(item['path']) != item['sha256']:
                raise ValueError('Original input changed during execution: '+item['path'])
        status['status']='PASS'
    finally:
        (run/(args.stage+'_status.json')).write_text(json.dumps(status,indent=2)+'\n')
        print('GENE_BIOGEOGRAPHY_'+args.stage.upper()+'_'+status['status'],flush=True)


if __name__=='__main__':
    main()
