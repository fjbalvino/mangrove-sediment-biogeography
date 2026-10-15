#!/usr/bin/env python3
"""Describe Shell KO contributions using the exact primary balanced selections.

This is an exploratory decomposition of the adjusted realm sum of squares,
not pathway enrichment, a causal attribution, or a test of KO significance.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--dictionary', type=Path, required=True)
    p.add_argument('--tables', type=Path, required=True)
    a = p.parse_args()
    master = pd.read_csv(a.tables/'180_KO_master_classification.tsv', sep='\t')
    run = pd.read_csv(a.tables/'180_run_summary.tsv', sep='\t').set_index('item')['value']
    pseudocount = float(run['pseudocount'])
    draws = pd.read_csv(a.tables/'180_panel_c_prevalence_balanced_R2.tsv', sep='\t')
    draws = draws[(draws['design']=='3 x 6') & (draws['group']=='Shell')]
    meta = pd.read_csv(a.metadata).set_index('run_accession')
    counts = pd.read_csv(a.counts, sep='\t', index_col=0).loc[:, meta.index]
    kos = master.loc[master.prevalence_class=='Shell', 'KO'].tolist()
    assert len(draws) == int(run['balanced_iterations'])
    x = np.log(counts.loc[kos].T.to_numpy(dtype=float) + pseudocount)
    x -= x.mean(axis=1, keepdims=True)
    total = np.log10(counts.sum(axis=0).to_numpy())
    indexes = {k:i for i,k in enumerate(meta.index)}
    shares = []
    errors = []
    for row in draws.itertuples():
        ids = row.selected_run_accessions.split(';')
        ix = np.array([indexes[k] for k in ids])
        y = x[ix]
        y = y-y.mean(axis=0)
        group = pd.get_dummies(meta.loc[ids, 'REALM'], dtype=float).to_numpy()
        base = np.column_stack([np.ones(len(ix)), total[ix]])
        q0 = np.linalg.qr(base)[0]
        residual_group = group-q0@(q0.T@group)
        u, s, _ = np.linalg.svd(residual_group, full_matrices=False)
        assert (s > 1e-8).sum() == 2
        term = u[:, :2].T@y
        per_ko = (term*term).sum(axis=0)
        observed = per_ko.sum()/(y*y).sum()
        errors.append(abs(observed-row.adjusted_unique_R2))
        shares.append(per_ko/per_ko.sum())
    shares = np.asarray(shares)
    assert np.max(errors)<1e-9, np.max(errors)
    assert np.allclose(shares.sum(axis=1), 1)
    dictionary = pd.read_csv(a.dictionary, sep='\t')
    assert dictionary.KO.is_unique
    out = pd.DataFrame({'KO': kos,
                        'mean_adjusted_realm_SS_share':shares.mean(axis=0),
                        'median_adjusted_realm_SS_share':np.median(shares, axis=0),
                        'q05_share':np.quantile(shares,.05,axis=0),
                        'q95_share':np.quantile(shares,.95,axis=0)})
    out = out.merge(dictionary[['KO','KO_function']], on='KO', how='left')
    out['annotation_available'] = out.KO_function.notna()
    out['KO_function'] = out.KO_function.fillna('No description in frozen dictionary')
    out = out.sort_values('mean_adjusted_realm_SS_share', ascending=False)
    out['rank'] = np.arange(1,len(out)+1)
    out.to_csv(a.tables/'180_shell_KO_contributions.tsv', sep='\t', index=False)
    summary = {'scope':'Exploratory conditional SS decomposition within Shell CLR subcomposition',
               'n_draws':len(draws), 'n_KOs':len(kos),
               'n_KOs_with_description':int(out.annotation_available.sum()),
               'top20_mean_share':float(out.head(20).mean_adjusted_realm_SS_share.sum()),
               'KOs_needed_for_half_mean_share':int((out.mean_adjusted_realm_SS_share.cumsum()<.5).sum()+1),
               'maximum_R2_reconciliation_error':float(max(errors)),
               'dictionary_sha256':hashlib.sha256(a.dictionary.read_bytes()).hexdigest(),
               'pathway_enrichment_tested':False,
               'limitation':'Descriptions do not provide a complete validated pathway mapping; do not infer dominance of nitrogen, sulfur or methanogenesis from keyword matches.'}
    (a.tables/'180_functional_contribution_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    main()
