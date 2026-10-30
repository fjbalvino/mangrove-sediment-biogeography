#!/usr/bin/env python3
"""Audit metadata support for gene-environment analyses without fitting effects.

Inputs remain unchanged. Coordinate groups are not assumed to be independent
biological sites. No environmental-unit conversions or missing-value imputations
are made. The matrix unit must be supplied, not inferred from integer values.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--metadata', required=True)
    p.add_argument('--matrix', required=True)
    p.add_argument('--samples', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--unit', required=True, choices=['predicted_protein_records'])
    a = p.parse_args()
    out = Path(a.outdir)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory must be new or empty')
    ids = Path(a.samples).read_text().split()
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError('Sample manifest must be nonempty and unique')
    meta = pd.read_csv(a.metadata, dtype={'run_accession': str, 'study_id': str,
                                         'sampling_site_id': str, 'sample_accession': str})
    required = {'run_accession', 'study_id', 'sampling_site_id', 'lat', 'lon', 'REALM'}
    if not required.issubset(meta.columns) or meta.run_accession.duplicated().any():
        raise ValueError('Missing metadata fields or duplicate accessions')
    if not set(ids).issubset(meta.run_accession):
        raise ValueError('Selected samples missing from metadata')
    meta = meta.set_index('run_accession').loc[ids].reset_index()
    if meta.study_id.isna().any():
        raise ValueError('All selected samples require a documented study label')
    matrix = pd.read_csv(a.matrix, sep='\t', index_col=0)
    if matrix.index.duplicated().any() or not set(ids).issubset(matrix.columns):
        raise ValueError('Duplicate KOs or missing matrix samples')
    matrix = matrix.loc[:, ids].apply(pd.to_numeric, errors='raise')
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError('Expected finite nonnegative predicted-protein record counts')
    totals = matrix.sum(axis=0)
    if (totals <= 0).any():
        raise ValueError('Zero-total profiles cannot support relative representation')
    meta['total_KO_assignments'] = meta.run_accession.map(totals)
    covariates = [f'bio{i}' for i in range(1, 20)] + [
        'soc_mean', 'phh2o_mean', 'clay_mean', 'sand_mean', 'silt_mean',
        'bdod_mean', 'cec_mean', 'cfvo_mean', 'human_footprint_2020',
        'gbif_occurrences_effective', 'gbif_richness_species_2km',
        'gbif_shannon_H_2km', 'gbif_pielou_J_2km', 'gbif_simpson_1minD_2km',
        'gbif_chao1_2km']
    records, conflicts, supports = [], [], []
    for var in covariates:
        if var not in meta:
            continue
        x = pd.to_numeric(meta[var], errors='coerce')
        work = meta.assign(value=x)
        nlevels = work.groupby('study_id').value.nunique()
        nsitelevels = work.groupby('sampling_site_id').value.nunique()
        records.append(dict(variable=var, n_profiles=int(x.notna().sum()),
            n_missing=int(x.isna().sum()), n_unique=int(x.nunique()),
            stored_min=None if x.dropna().empty else float(x.min()),
            stored_max=None if x.dropna().empty else float(x.max()),
            n_studies_with_data=int(work.loc[x.notna(), 'study_id'].nunique()),
            n_studies_with_within_variation=int((nlevels > 1).sum()),
            n_coordinate_groups_with_conflicting_values=int((nsitelevels > 1).sum())))
        for study, group in work.groupby('study_id'):
            supports.append(dict(variable=var, study_id=study,
                n_profiles=int(group.value.notna().sum()), n_unique=int(group.value.nunique()),
                n_coordinate_groups=int(group.loc[group.value.notna(), 'sampling_site_id'].nunique())))
        for site in nsitelevels[nsitelevels > 1].index:
            conflicts.append(dict(variable=var, coordinate_group=site,
                stored_values=';'.join(str(v) for v in sorted(work.loc[work.sampling_site_id == site, 'value'].dropna().unique())),
                reason='Review source, extraction and rounding before spatial aggregation'))
    designs = []
    for covs in [['bio1'], ['soc_mean'], ['bio1', 'soc_mean', 'phh2o_mean']]:
        if not set(covs).issubset(meta.columns):
            continue
        complete = meta.dropna(subset=covs + ['sampling_site_id'])
        # A descriptive support audit, not a declaration of biological replication.
        grouped = complete.groupby(['study_id', 'sampling_site_id'], as_index=False)[covs].first()
        for var in covs:
            if (complete.groupby(['study_id','sampling_site_id'])[var].nunique() > 1).any():
                raise ValueError('Target covariates conflict within study-coordinate group: ' + var)
        x = grouped[covs].to_numpy(dtype=float)
        sd = x.std(axis=0)
        if (sd == 0).any():
            continue
        x = (x - x.mean(axis=0)) / sd
        studies = pd.get_dummies(grouped.study_id).to_numpy(dtype=float)
        r0 = int(np.linalg.matrix_rank(studies))
        r1 = int(np.linalg.matrix_rank(np.column_stack([studies, x])))
        designs.append(dict(covariates=covs, n_profiles=len(complete),
            n_study_coordinate_groups=len(grouped), n_coordinate_groups=int(grouped.sampling_site_id.nunique()),
            n_studies=int(grouped.study_id.nunique()), rank_study=r0,
            rank_study_plus_covariates=r1, incremental_rank=r1-r0,
            nominal_residual_df=len(grouped)-r1,
            warning='Algebraic rank only; does not validate independence, power, units or per-gene fit'))
    dummy_study = pd.get_dummies(meta.study_id).to_numpy(dtype=float)
    dummy_realm = pd.get_dummies(meta.REALM).to_numpy(dtype=float)
    coordinates = pd.to_numeric(meta.lat, errors='coerce').notna() & pd.to_numeric(meta.lon, errors='coerce').notna()
    paired = meta[['sample_accession', 'run_accession']].dropna() if 'sample_accession' in meta else pd.DataFrame()
    repeated = paired[paired.sample_accession.duplicated(keep=False)] if len(paired) else paired
    report = dict(scope='Read-only audit; no gene-environment associations estimated',
        n_profiles=len(ids), matrix_unit=a.unit, n_KO_rows=len(matrix),
        n_KOs_detected=int((matrix.sum(axis=1)>0).sum()), n_zero_only_KOs=int((matrix.sum(axis=1)==0).sum()),
        n_studies=int(meta.study_id.nunique()), n_coordinate_groups=int(meta.sampling_site_id.nunique()),
        n_profiles_with_coordinates=int(coordinates.sum()),
        n_profiles_without_coordinates=int((~coordinates).sum()),
        total_KO_min=int(totals.min()), total_KO_max=int(totals.max()),
        repeated_BioSamples=int(repeated.sample_accession.nunique()) if len(repeated) else 0,
        rank_increment_realm_given_study=int(np.linalg.matrix_rank(np.column_stack([dummy_study,dummy_realm]))-np.linalg.matrix_rank(dummy_study)),
        design_checks=designs,
        absent_structured_columns=[c for c in ['salinity','DOC','TOC','redox','oxygen','sampling_date','depth_cm','sediment_temperature'] if c not in meta],
        unit_flags=[
            'Verify source encodings for each environmental variable; never infer units from column names.',
            'Do not interpret lack of structured columns as proof measurements were never made.'] +
            (['phh2o_mean exceeds 14; pH x 10 is a candidate encoding to verify.']
             if 'phh2o_mean' in meta and pd.to_numeric(meta.phh2o_mean,errors='coerce').max()>14 else []) +
            (['Human footprint conflicts within rounded coordinate groups; inspect extraction provenance.']
             if any(c['variable']=='human_footprint_2020' for c in conflicts) else []),
        sha256={name:digest(path) for name,path in [('metadata',a.metadata),('matrix',a.matrix),('samples',a.samples),('script',__file__)]},
        versions={'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__})
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out/'environment_variable_support.tsv', sep='\t', index=False)
    pd.DataFrame(supports).to_csv(out/'within_study_support.tsv', sep='\t', index=False)
    pd.DataFrame(conflicts,columns=['variable','coordinate_group','stored_values','reason']).to_csv(out/'coordinate_value_conflicts.tsv', sep='\t', index=False)
    repeated.to_csv(out/'repeated_biosamples.tsv', sep='\t', index=False)
    (out/'audit_summary.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({k:report[k] for k in ['scope','n_profiles','n_KO_rows','n_KOs_detected','n_studies','n_coordinate_groups','n_profiles_with_coordinates']}))


if __name__ == '__main__':
    main()
