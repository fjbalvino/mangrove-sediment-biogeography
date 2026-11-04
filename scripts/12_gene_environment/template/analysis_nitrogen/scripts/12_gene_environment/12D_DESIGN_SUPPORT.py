#!/usr/bin/env python3
"""Audit nitrogen estimands and design support; never fit gene associations.

Consumes the frozen 12B output, joins by IDs, and writes separate derived tables.
Study labels and coordinate bins are not certified independent biological units.
"""
import argparse
import base64
import hashlib
import html
import json
import platform
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

COVS = ['Temperature', 'Organic carbon', 'pH']
RAW = ['bio1', 'soc_mean', 'phh2o_mean']
FACTORS = [1., 0.1, 0.1]
UNITS = ['degrees C', 'g/kg', 'pH units']
TOL = 1e-10


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unique_index(frame, key, label):
    if frame[key].isna().any() or frame[key].duplicated().any():
        raise ValueError(f'{label}: missing or duplicated {key}')
    if not frame[key].map(lambda x:isinstance(x,str) and bool(x) and x.strip()==x).all():
        raise ValueError(f'{label}: invalid {key}')
    return frame.set_index(key)


def rank(x):
    x = np.asarray(x, dtype=float)
    if not x.size:
        return 0
    return int(np.linalg.matrix_rank(x, tol=TOL))


def within(frame, columns, block='study_id'):
    return frame[columns] - frame.groupby(block)[columns].transform('mean')


def standardized_within(frame, columns):
    w = within(frame, columns)
    # Scale before rank checks, removing units without amplifying zero columns.
    sd = frame[columns].std(ddof=0).replace(0, 1).fillna(1)
    return w / sd


def design_stats(frame, columns):
    w = standardized_within(frame, columns)
    r = rank(w)
    nblocks = frame.study_id.nunique()
    return dict(n_groups=len(frame), n_coordinate_groups=frame.coordinate_group.nunique(),
                n_archive_labels=nblocks, rank_study=nblocks,
                added_covariate_rank=r, full_design_rank=nblocks+r,
                nominal_residual_df=len(frame)-nblocks-r)


def predictor_support(frame, columns, context):
    """Design-only variance shares; these are NOT precision weights or effective n."""
    w = within(frame, columns)
    z = standardized_within(frame, columns)
    fullrank = rank(z)
    summary, contributions = [], []
    for c in columns:
        rawss = w[c].pow(2).groupby(frame.study_id).sum()
        shares = rawss / rawss.sum() if rawss.sum() > TOL else rawss * np.nan
        other = [x for x in columns if x != c]
        residual = z[c].to_numpy().copy()
        if other:
            residual -= z[other].to_numpy() @ np.linalg.lstsq(z[other], residual, rcond=None)[0]
        partialss = pd.Series(residual**2, index=frame.index).groupby(frame.study_id).sum()
        partialshares = partialss / partialss.sum() if partialss.sum() > TOL else partialss*np.nan
        ranges = frame.groupby('study_id')[c].agg(lambda a: a.max()-a.min())
        top = shares.idxmax() if shares.notna().any() else None
        d = dict(context=context, predictor=c, unit=UNITS[COVS.index(c)],
                 **design_stats(frame, columns),
                 studies_with_raw_predictor_variation=int((ranges>TOL).sum()),
                 unique_covariate_values=frame[c].nunique(), minimum=frame[c].min(), maximum=frame[c].max(),
                 conditional_rank_increment=fullrank-rank(z[other]) if other else fullrank,
                 largest_raw_within_variance_study=top,
                 largest_raw_within_variance_share=shares.max(),
                 largest_partial_design_share=partialshares.max(),
                 largest_partial_design_study=partialshares.idxmax() if partialshares.notna().any() else None,
                 variance_share_concentration=1/np.sum(shares**2) if shares.notna().any() else np.nan,
                 status='conditional', inference='Design support only; no association fitted')
        summary.append(d)
        for s in rawss.index:
            contributions.append(dict(context=context, predictor=c, study_id=s,
                n_groups=int(frame.study_id.eq(s).sum()), raw_within_range=ranges[s],
                raw_within_ss=rawss[s], raw_within_variance_share=shares[s],
                partial_design_share=partialshares[s]))
    return summary, contributions


def main():
    root=Path(__file__).resolve().parents[2];out=root/'design_review';out.mkdir(exist_ok=True)
    meta=pd.read_csv(root/'results/sample_analysis_manifest.tsv',sep='\t').set_index('run_accession')
    groups=pd.read_csv(root/'results/study_coordinate_groups.tsv',sep='\t').set_index('analysis_group')
    unit=meta[['analysis_group','study_id','biosample_group']].copy()
    nr=meta.groupby(['analysis_group','biosample_group']).biosample_group.transform('size')
    nb=meta.groupby('analysis_group').biosample_group.transform('nunique')
    unit['descriptive_run_weight_within_group']=1/(nr*nb)
    assert np.allclose(unit.groupby('analysis_group').descriptive_run_weight_within_group.sum(),1)
    unit.to_csv(out/'unit_design_manifest.tsv',sep='\t')
    complete=groups.dropna(subset=RAW).copy()
    for raw,name,factor in zip(RAW,COVS,FACTORS):complete[name]=complete[raw]*factor
    summaries,contributions=predictor_support(complete,COVS,'joint_complete_cases') if len(complete) else ([],[])
    pd.DataFrame(summaries).to_csv(out/'environmental_design_support.tsv',sep='\t',index=False)
    pd.DataFrame(contributions).to_csv(out/'environmental_information_by_study.tsv',sep='\t',index=False)
    groups[groups[RAW].isna().any(axis=1)].to_csv(out/'environmental_complete_case_exclusions.tsv',sep='\t')
    result={'profiles':len(meta),'groups':len(groups),'coordinates':meta.coordinate_group.nunique(),
            'complete_profiles':int(meta.analysis_group.isin(complete.index).sum()),'complete_groups':len(complete),'complete_studies':complete.study_id.nunique(),
            'design':design_stats(complete,COVS) if len(complete) else {},
            'interpretation':'Study contributions are shares of within-study predictor sums of squares, not variance explained in genes or effective sample sizes. Archive labels do not certify biological independence.'}
    (out/'design_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
