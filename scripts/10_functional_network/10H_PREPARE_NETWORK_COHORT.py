#!/usr/bin/env python3
"""Prepare the network cohort and verify persistence classes against atlas output."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--accessions', type=Path, required=True)
    p.add_argument('--atlas-classes', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    args = p.parse_args()
    ids = args.accessions.read_text().splitlines()
    ids = [s.strip() for s in ids if s.strip() and not s.startswith('#')]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate accessions')
    meta = pd.read_csv(args.metadata, sep='\t' if args.metadata.suffix == '.tsv' else ',')
    if meta.run_accession.duplicated().any():
        raise ValueError('Duplicate metadata accessions')
    meta = meta.set_index('run_accession').loc[ids]
    if meta[['REALM', 'Country', 'Locality', 'core_locality']].isna().any().any():
        raise ValueError('Missing geographic stratum')
    raw = pd.read_csv(args.counts, sep='\t', index_col=0)
    if raw.index.duplicated().any() or raw.columns.duplicated().any():
        raise ValueError('Duplicate matrix identifiers')
    counts = raw.loc[:, ids]
    if not np.isfinite(counts.to_numpy()).all() or (counts < 0).any().any():
        raise ValueError('Nonfinite or negative counts')
    counts = counts.loc[counts.sum(axis=1) > 0]
    atlas_path = Path(__file__).resolve().parents[1] / '08_functional_atlas/08A_ANALYSE_FUNCTIONAL_ATLAS.py'
    spec = importlib.util.spec_from_file_location('functional_atlas', atlas_path)
    atlas = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(atlas)
    ann, _, support, _ = atlas.classification(counts, meta)
    ref = pd.read_csv(args.atlas_classes, sep='\t').set_index('KO').sort_index()
    computed = ann.set_index('KO').sort_index()
    if not ref.index.equals(computed.index):
        raise ValueError('Atlas KO set differs from current cohort')
    for col in ['class', 'samples_present', 'KO_assignment_count']:
        if not ref[col].equals(computed[col]):
            raise ValueError('Atlas mismatch: ' + col)
    presence = counts.T.gt(0).groupby(meta.REALM, sort=True).any().T
    n_realms = presence.sum(axis=1).reindex(ann.KO).to_numpy()
    mapping = {'Absolute core': 'absolute_global_core', 'Near-global': 'near_global',
               'Realm-core': 'realm_core_only', 'Periphery': 'multi_realm_non_core'}
    converted = ann['class'].map(mapping)
    converted.loc[(ann['class'] == 'Periphery') & (n_realms == 1)] = 'realm_specific'
    out = pd.DataFrame({
        'KO': ann.KO, 'functional_core_class': converted,
        'global_prevalence_n': ann.samples_present,
        'global_prevalence_prop': ann.sample_prevalence,
        'n_functional_samples': len(ids), 'total_abundance': ann.KO_assignment_count,
        'n_realms_detected': n_realms,
        'realms_detected': [';'.join(presence.columns[row]) for row in presence.reindex(ann.KO).to_numpy()],
        'realm_core_realms': ann.supporting_realms,
        'near_global_threshold': .95, 'min_realm_localities': 3,
        'min_mean_relative_abundance_realm_core': .0001,
        'presence_threshold': 0,
        'abundance_weighting': 'equal locality weights after pooling replicate counts',
    })
    for cls in mapping.values():
        out['is_' + cls] = out.functional_core_class.eq(cls)
    out['is_realm_specific'] = out.functional_core_class.eq('realm_specific')
    prev = counts.gt(0).sum(axis=1)
    abundance = counts.sum(axis=1)
    classes = pd.Series(np.where(prev / len(ids) >= .95, 'Prevalent',
                                np.where(prev / len(ids) >= .15, 'Shell', 'Cloud')),
                        index=counts.index)
    assoc = (prev >= 10) & (abundance >= 10)
    args.outdir.mkdir(parents=True, exist_ok=True)
    meta.reset_index().to_csv(args.outdir / 'metadata_191.tsv', sep='\t', index=False)
    out.to_csv(args.outdir / 'KO_core_classes.tsv', sep='\t', index=False)
    support.to_csv(args.outdir / 'realm_core_support.tsv', sep='\t', index=False)
    info = {
        'n_samples': len(ids), 'n_KOs': len(counts), 'n_association_KOs': int(assoc.sum()),
        'n_network_eligible_KOs': int(((prev >= 20) & (abundance >= 50)).sum()),
        'all_prevalence_classes': {k: int(v) for k, v in classes.value_counts().items()},
        'association_prevalence_classes': {k: int(v) for k, v in classes[assoc].value_counts().items()},
        'persistence_classes': {k: int(v) for k, v in out.functional_core_class.value_counts().items()},
        'samples_per_balanced_subsample': int(meta.groupby(['REALM', 'Country', 'Locality']).size().clip(upper=4).sum()),
        'excluded_raw_matrix_samples': sorted(set(raw.columns) - set(ids)),
        'atlas_classification_reproduced': True,
        'input_sha256': {str(f.resolve()): hashlib.sha256(f.read_bytes()).hexdigest()
                         for f in [args.counts, args.metadata, args.accessions, args.atlas_classes]},
    }
    (args.outdir / 'cohort_audit.json').write_text(json.dumps(info, indent=2) + '\n')
    print(json.dumps(info, indent=2))


if __name__ == '__main__':
    main()
