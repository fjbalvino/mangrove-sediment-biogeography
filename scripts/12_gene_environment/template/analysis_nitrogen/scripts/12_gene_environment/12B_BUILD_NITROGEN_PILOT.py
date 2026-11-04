#!/usr/bin/env python3
"""Nitrogen pilot: ID-safe observations, descriptive grouping, and fit diagnostics.

Protein-record assignments, not read abundances. Sources remain unchanged.
No p-values: exploratory slopes and leave-study-out ranges are diagnostics,
not confidence intervals or validated primary MaAsLin models.
"""
import argparse
import csv
import hashlib
import json
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def validate_ids(values, label):
    values = list(values)
    if not values or any(not isinstance(v, str) or not v or v.strip() != v for v in values):
        raise ValueError(label + ': empty or invalid identifier')
    if len(values) != len(set(values)):
        raise ValueError(label + ': duplicate identifier')
    return values


def load_inputs(metadata, matrix, manifest, expected):
    ids = validate_ids([line for line in Path(manifest).read_text().splitlines() if line.strip()], 'manifest')
    if expected is not None and len(ids) != expected:
        raise ValueError('Unexpected cohort size')
    meta = pd.read_csv(metadata, dtype={'run_accession': str, 'study_id': str, 'sample_accession': str})
    validate_ids(meta.run_accession, 'metadata')
    if set(meta.run_accession) != set(ids):
        raise ValueError('Metadata and manifest IDs must match exactly')
    # pandas mangles duplicate CSV headers. Inspect the original header first.
    import gzip
    opener = gzip.open if str(matrix).endswith('.gz') else open
    with opener(matrix, 'rt', encoding='utf-8-sig') as f:
        validate_ids(next(csv.reader(f, delimiter='\t'))[1:], 'matrix columns')
    mat = pd.read_csv(matrix, sep='\t', index_col=0)
    validate_ids(mat.index.tolist(), 'KO rows')
    if not set(ids).issubset(mat.columns):
        raise ValueError('Manifest IDs missing in matrix')
    mat = mat.loc[:, ids].apply(pd.to_numeric, errors='raise')
    a = mat.to_numpy(dtype=float)
    if not np.isfinite(a).all() or (a < 0).any() or not np.equal(a, np.floor(a)).all():
        raise ValueError('Invalid predicted-protein record counts')
    if (mat.sum() <= 0).any():
        raise ValueError('Zero-total profile')
    meta = meta.set_index('run_accession').loc[ids].copy()
    if meta.study_id.isna().any():
        raise ValueError('Missing study label')
    for c, limit in [('lat', 90), ('lon', 180)]:
        meta[c] = pd.to_numeric(meta[c], errors='raise')
        if not np.isfinite(meta[c]).all() or (meta[c].abs() > limit).any():
            raise ValueError('Pilot requires verified finite coordinates for all selected profiles')
    meta['coordinate_group'] = [f'coord_{x:.5f}_{y:.5f}' for x, y in zip(meta.lat, meta.lon)]
    meta['analysis_group'] = meta.study_id + '::' + meta.coordinate_group
    meta['biosample_group'] = [str(s) if pd.notna(s) and str(s).strip() else 'run:' + r
                              for r, s in zip(meta.index, meta.sample_accession)]
    return meta, mat


def hierarchical_summary(frame, meta):
    """Equal run weights within BioSample, then equal BioSample weights per study-point.

    This is a descriptive summary, never a merge of the input profiles or a claim
    that BioSamples/coordinate groups identify independent biological replicates.
    """
    a = frame.T.join(meta[['analysis_group', 'biosample_group']], validate='one_to_one')
    within = a.groupby(['analysis_group', 'biosample_group'], sort=True).mean()
    return within.groupby(level=0).mean()


def within_fit(x, y, study):
    """OLS after within-study centering; no inferential standard errors."""
    x, y, study = np.asarray(x, float), np.asarray(y, float), np.asarray(study)
    xc, yc = x.copy(), y.copy()
    for s in np.unique(study):
        sel = study == s
        xc[sel] -= x[sel].mean(axis=0)
        yc[sel] -= y[sel].mean()
    if not len(y) or not np.isfinite(xc).all() or not np.isfinite(yc).all() or np.linalg.matrix_rank(xc) != x.shape[1]:
        return None
    b = np.linalg.lstsq(xc, yc, rcond=None)[0]
    info = np.sum(xc * xc, axis=0)
    contrib = np.array([np.sum(xc[study == s] ** 2, axis=0) for s in np.unique(study)])
    # Frisch-Waugh-Lovell geometry: information for a partial coefficient is
    # the predictor residual after adjusting for every other predictor.
    # Marginal X^2 alone can conceal study dominance under collinearity.
    partial = np.empty_like(xc)
    for j in range(xc.shape[1]):
        other = np.delete(xc, j, axis=1)
        partial[:, j] = xc[:, j] - other @ np.linalg.lstsq(other, xc[:, j], rcond=None)[0]
    partial_info = np.sum(partial ** 2, axis=0)
    partial_contrib = np.array([np.sum(partial[study == s] ** 2, axis=0) for s in np.unique(study)])
    shares = partial_contrib / partial_info
    residual = yc - xc @ b
    df = len(y) - len(np.unique(study)) - x.shape[1]
    return dict(beta=b, df=df, residual=residual, xc=xc, yc=yc,
                n_informative=(contrib > 1e-12).sum(axis=0),
                max_study_information_share=contrib.max(axis=0) / info,
                partial_predictors=partial,
                conditional_variation_fraction=partial_info / info,
                n_studies_with_conditional_information=(shares > 1e-10).sum(axis=0),
                max_conditional_study_information_share=shares.max(axis=0),
                effective_information_studies=1 / np.sum(shares ** 2, axis=0),
                condition_number=float(np.linalg.cond(xc)))


def association_diagnostics(group_rel, group_detection, groups, labels, out):
    covs = ['bio1', 'soc_mean', 'phh2o_mean']
    complete = groups.dropna(subset=covs).copy()
    # Verified raw units; conversions only in analysis objects, never original metadata.
    complete['soc_mean'] /= 10.0
    complete['phh2o_mean'] /= 10.0
    scale = complete[covs].std(ddof=0)
    records, support = [], []
    for ko in group_rel.columns:
        for outcome, table in [('log2_positive_group_mean', group_rel),
                               ('mean_observed_detection', group_detection)]:
            y0 = table.reindex(complete.index)[ko]
            sel = y0.notna()
            if outcome.startswith('log2'):
                sel &= y0 > 0
            dat = complete.loc[sel]
            y = y0.loc[sel].to_numpy()
            if outcome.startswith('log2'):
                y = np.log2(y)
            x = (dat[covs] / scale).to_numpy()
            studies = dat.study_id.to_numpy()
            fit = within_fit(x, y, studies) if len(dat) > 0 else None
            status = 'descriptive_only' if fit is not None and fit['df'] > 0 else 'not_estimable'
            if fit is not None and np.ptp(fit['yc']) < 1e-12:
                status = 'no_within_study_response_variation'
            support.append(dict(KO=ko, outcome=outcome, n_groups=len(dat),
                                n_studies=dat.study_id.nunique(), status=status))
            if status != 'descriptive_only':
                continue
            loo = []
            for s in np.unique(studies):
                # Removing a study with no within-study predictor contrast leaves
                # the slopes unchanged and must not inflate stability fractions.
                if np.max(np.var(x[studies == s], axis=0)) < 1e-12:
                    continue
                take = studies != s
                f = within_fit(x[take], y[take], studies[take])
                if f is not None and f['df'] > 0:
                    loo.append(f['beta'])
            loo = np.array(loo)
            # Measurement sensitivity only, not a selected alternative model.
            lx = np.log2(dat.total_KO_assignments.to_numpy())
            adjusted = None
            if np.isfinite(lx).all() and lx.std() > 0:
                lx = (lx - lx.mean()) / lx.std()
                adjusted = within_fit(np.column_stack([x, lx]), y, studies)
            for i, c in enumerate(covs):
                b = float(fit['beta'][i])
                records.append(dict(KO=ko, gene=labels[ko], outcome=outcome,
                    predictor=c, beta_per_sd=b, predictor_sd=float(scale[c]),
                    n_groups=len(dat), n_studies=len(np.unique(studies)),
                    n_studies_with_predictor_variation=int(fit['n_informative'][i]),
                    max_study_information_share=float(fit['max_study_information_share'][i]),
                    conditional_variation_fraction=float(fit['conditional_variation_fraction'][i]),
                    n_studies_with_conditional_information=int(fit['n_studies_with_conditional_information'][i]),
                    max_conditional_study_information_share=float(fit['max_conditional_study_information_share'][i]),
                    effective_information_studies=float(fit['effective_information_studies'][i]),
                    nominal_residual_df=fit['df'], condition_number=fit['condition_number'],
                    loo_min=float(loo[:, i].min()) if len(loo) else np.nan,
                    loo_max=float(loo[:, i].max()) if len(loo) else np.nan,
                    loo_sign_agreement=float(np.mean(np.sign(loo[:, i]) == np.sign(b))) if len(loo) else np.nan,
                    loo_valid_fits=len(loo),
                    beta_with_annotation_yield_sensitivity=float(adjusted['beta'][i]) if adjusted else np.nan,
                    inference='Exploratory slopes only; LOO range is not a confidence interval; no p or q'))
    pd.DataFrame(records, columns=['KO','gene','outcome','predictor','beta_per_sd','predictor_sd','n_groups','n_studies','n_studies_with_predictor_variation','max_study_information_share','conditional_variation_fraction','n_studies_with_conditional_information','max_conditional_study_information_share','effective_information_studies','nominal_residual_df','condition_number','loo_min','loo_max','loo_sign_agreement','loo_valid_fits','beta_with_annotation_yield_sensitivity','inference']).to_csv(out/'environment_association_diagnostics.tsv', sep='\t', index=False)
    pd.DataFrame(support).to_csv(out/'per_KO_model_support.tsv', sep='\t', index=False)
    return dict(n_complete_groups=len(complete), n_complete_studies=complete.study_id.nunique(),
                n_diagnostic_coefficients=len(records), inference='No confirmatory tests performed',
                formula='group outcome ~ study fixed effects + standardized BIO1 + SOC + pH',
                unit='Descriptive study-coordinate group, equal BioSample weighting',
                p_values=False, confidence_intervals=False, selected_on_significance=False,
                sd_in_analysis_units=scale.to_dict())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['metadata', 'matrix', 'samples', 'rules', 'reactions', 'outdir']:
        p.add_argument('--' + key, required=True)
    p.add_argument('--expected-profiles', type=int, default=None)
    a = p.parse_args()
    out = Path(a.outdir)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use a new or empty output directory')
    meta, mat = load_inputs(a.metadata, a.matrix, a.samples, a.expected_profiles)
    reactions = json.loads(Path(a.reactions).read_text())
    rules = pd.read_csv(a.rules, sep='\t')
    nr = rules[rules.process_id.str.startswith('N')]
    required = set(re.findall(r'K\d{5}', ' '.join(nr.expression)))
    kos = list(dict.fromkeys(k for r in reactions for k in r['kos']))
    if set(kos) != required:
        raise ValueError('Reaction catalogue must cover exactly the frozen nitrogen-rule KOs')
    labels = {k:g for r in reactions for k, g in zip(r['kos'], r['genes'])}
    if any(len(r['kos']) != len(r['genes']) for r in reactions):
        raise ValueError('KO-label mismatch')
    totals = mat.sum()
    panel = mat.reindex(kos)
    rel = panel.div(totals, axis=1) * 1e6
    det = panel.ge(1).astype(float).mask(panel.isna())
    det2 = panel.ge(2).astype(float).mask(panel.isna())
    grouped = hierarchical_summary(rel, meta)
    gd = hierarchical_summary(det, meta)
    gd2 = hierarchical_summary(det2, meta)
    meta['total_KO_assignments'] = totals
    groups = meta.groupby('analysis_group')[['study_id', 'coordinate_group', 'lat', 'lon', 'Country', 'Locality', 'REALM']].first()
    for c in ['bio1', 'soc_mean', 'phh2o_mean']:
        conflicts = meta.groupby('analysis_group')[c].nunique()
        if (conflicts > 1).any():
            raise ValueError('Environmental conflict within study-coordinate group: ' + c)
        groups[c] = meta.groupby('analysis_group')[c].first()
    groups['n_runs'] = meta.groupby('analysis_group').size()
    groups['n_biosample_labels'] = meta.groupby('analysis_group').biosample_group.nunique()
    groups['total_KO_assignments'] = hierarchical_summary(pd.DataFrame([totals], index=['value']), meta)['value']
    out.mkdir(parents=True)
    meta.to_csv(out/'sample_analysis_manifest.tsv', sep='\t')
    groups.to_csv(out/'study_coordinate_groups.tsv', sep='\t')
    for name, table in [('nitrogen_raw_assignments', panel), ('nitrogen_per_million_KO_assignments', rel),
                        ('nitrogen_detection_ge1', det), ('nitrogen_detection_ge2', det2),
                        ('group_relative_representation', grouped), ('group_detection_ge1', gd), ('group_detection_ge2', gd2)]:
        table.to_csv(out/(name + '.tsv'), sep='\t', index_label='analysis_group' if name.startswith('group_') else 'KO')
    inv = pd.DataFrame([dict(KO=k, gene=labels[k], row_available=k in mat.index,
               detected_ge1=int(det.loc[k].sum()) if k in mat.index else np.nan,
               detected_ge2=int(det2.loc[k].sum()) if k in mat.index else np.nan,
               profiles=len(meta), median_per_million=float(rel.loc[k].median()) if k in mat.index else np.nan)
               for k in kos])
    inv.to_csv(out/'nitrogen_KO_inventory.tsv', sep='\t', index=False)
    reaction_rows = [dict(reaction_id=r['id'], label=r['label'], reaction=r['reaction'],
                          KO=k, gene=g, status=r['status'], note=r['note'], source=r['source'])
                     for r in reactions for k, g in zip(r['kos'], r['genes'])]
    pd.DataFrame(reaction_rows).to_csv(out/'nitrogen_reaction_gene_map.tsv', sep='\t', index=False)
    diag = association_diagnostics(grouped, gd, groups, labels, out)
    report = dict(n_profiles=len(meta), n_studies=meta.study_id.nunique(),
        n_coordinate_groups=meta.coordinate_group.nunique(), n_study_coordinate_groups=len(groups),
        n_nitrogen_KOs=len(kos), n_available_KOs=int(inv.row_available.sum()),
        unavailable_KOs=inv.loc[~inv.row_available, 'KO'].tolist(),
        metadata_sha256=digest(a.metadata), matrix_sha256=digest(a.matrix),
        samples_sha256=digest(a.samples), reactions_sha256=digest(a.reactions),
        scripts_sha256=digest(__file__),
        denominator=f'All {len(mat)} matrix KO rows before any panel filter; zero-only rows have no effect',
        group_summary='Mean runs within BioSample label, then mean labels within study-coordinate; descriptive only',
        units='KO assignments per million recovered KO assignments; not TPM, activity or flux',
        environmental_diagnostics=diag,
        versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__))
    (out/'pilot_summary.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
