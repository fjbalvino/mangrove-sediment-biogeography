#!/usr/bin/env python3
"""Recompute original Figure 5 process scores and locality turnover.

Adapted from the recovered scripts 158 and 179. Uses the canonical cohort and
the exact balanced selections already used by 09A, rather than historical RDS
objects. Selected process-score R2 is not a contribution to total KO variance.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def regional_effects(y, realm, log_total):
    """Univariate eta2 and unique realm SS / total SS for each column."""
    y = y - y.mean(axis=0)
    total_ss = (y*y).sum(axis=0)
    assert (total_ss > 0).all()
    groups = pd.get_dummies(realm, dtype=float).to_numpy()
    group_means = groups.T @ y / groups.sum(axis=0)[:, None]
    between_ss = ((groups @ group_means)**2).sum(axis=0)
    baseline = np.column_stack([np.ones(len(y)), log_total])
    q = np.linalg.qr(baseline)[0]
    residual_groups = groups - q @ (q.T @ groups)
    u, s, _ = np.linalg.svd(residual_groups, full_matrices=False)
    rank = int((s > 1e-8).sum())
    assert rank == groups.shape[1]-1
    unique_ss = ((u[:, :rank].T @ y)**2).sum(axis=0)
    return between_ss / total_ss, unique_ss / total_ss


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['counts', 'metadata', 'accessions', 'catalog', 'tables', 'objects']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--n-bootstrap', type=int, default=500)
    p.add_argument('--seed', type=int, default=20261015)
    a = p.parse_args()
    a.objects.mkdir(parents=True, exist_ok=True)
    ids = a.accessions.read_text().split()
    assert len(ids) == len(set(ids)) == 191
    meta = pd.read_csv(a.metadata).set_index('run_accession').loc[ids]
    assert meta[['REALM', 'Country', 'core_locality']].notna().all().all()
    counts = pd.read_csv(a.counts, sep='\t', index_col=0).loc[:, ids]
    assert counts.index.is_unique
    assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy() >= 0).all()
    counts = counts.loc[counts.sum(axis=1)>0]
    master = pd.read_csv(a.tables/'180_KO_master_classification.tsv', sep='\t').set_index('KO')
    assert set(master.index) == set(counts.index) and len(master)==17482
    assert ((counts>0).sum(axis=1) == master.loc[counts.index, 'global_prevalence_n']).all()
    shell_ids = master.index[master.prevalence_class=='Shell'].tolist()
    assert len(shell_ids)==9886
    x = np.log(counts.loc[shell_ids].T.to_numpy(dtype=float)+1)
    x -= x.mean(axis=1, keepdims=True)
    assert np.abs(x.mean(axis=1)).max()<1e-12
    log_total = np.log10(counts.sum(axis=0).to_numpy())
    index = {k:i for i,k in enumerate(ids)}
    shell_index = {k:i for i,k in enumerate(shell_ids)}

    catalog = pd.read_csv(a.catalog, sep='\t')
    assert catalog.KO.is_unique and catalog.KO.isin(master.index).all()
    catalog = catalog.merge(master[['prevalence_class','persistence_class','global_prevalence_n']], left_on='KO', right_index=True, validate='one_to_one')
    catalog['included_in_score'] = catalog.prevalence_class.eq('Shell')
    catalog['exclusion_reason'] = np.where(catalog.included_in_score, '', 'Current prevalence is outside Shell')
    assert set(catalog.loc[~catalog.included_in_score,'KO'])=={'K00441'}
    catalog.to_csv(a.tables/'181_process_catalog_validation.tsv', sep='\t', index=False)
    used = catalog[catalog.included_in_score]
    processes = catalog.process_id.unique().tolist()
    labels = catalog.drop_duplicates('process_id').set_index('process_id').process_label.to_dict()
    scores = np.column_stack([x[:, [shell_index[k] for k in used.loc[used.process_id==g,'KO']]].mean(axis=1) for g in processes])
    z = (scores-scores.mean(axis=0))/scores.std(axis=0, ddof=1)
    sample_scores = pd.DataFrame(z, index=ids, columns=processes)
    sample_scores.index.name='run_accession'
    sample_scores.to_csv(a.tables/'181_sample_process_scores_z.tsv', sep='\t')
    pd.DataFrame(scores,index=ids,columns=processes).rename_axis('run_accession').to_csv(a.tables/'181_sample_process_scores_CLR.tsv', sep='\t')

    draws = pd.read_csv(a.tables/'180_panel_c_prevalence_balanced_R2.tsv',sep='\t')
    draws = draws[draws.group.eq('Shell')].copy()
    realms = ['Central Indo-Pacific','Tropical Atlantic','Western Indo-Pacific']
    iterations=[]; means=[]
    for row in draws.itertuples():
        selected=row.selected_run_accessions.split(';')
        m=meta.loc[selected]
        nloc=6 if row.design=='3 x 6' else 4
        assert len(selected)==len(set(selected))==3*nloc
        assert set(m.REALM)==set(realms)
        assert (m.groupby('REALM').core_locality.nunique()==nloc).all()
        ix=[index[s] for s in selected]
        raw, adjusted=regional_effects(z[ix],m.REALM,log_total[ix])
        for j,pid in enumerate(processes):
            iterations.append(dict(design=row.design,iteration=row.iteration,process_id=pid,R2=raw[j],adjusted_unique_R2=adjusted[j]))
        if row.design=='3 x 6':
            means.append(np.vstack([z[[index[k] for k in m.index[m.REALM.eq(r)]]].mean(axis=0) for r in realms]))
    effects=pd.DataFrame(iterations)
    assert len(means)==499 and len(effects)==2*499*6
    effects.to_csv(a.tables/'181_process_balanced_effects.tsv',sep='\t',index=False)
    summary=[]
    for (design,pid),g in effects.groupby(['design','process_id'],sort=False):
        r=dict(design=design,process_id=pid,process_label=labels[pid],n_KOs=int((used.process_id==pid).sum()),n_draws=len(g))
        for field in ['R2','adjusted_unique_R2']:
            r.update({field+'_median':float(g[field].median()),field+'_q05':float(g[field].quantile(.05)),field+'_q95':float(g[field].quantile(.95))})
        summary.append(r)
    summary=pd.DataFrame(summary)
    summary.to_csv(a.tables/'181_process_regional_summary.tsv',sep='\t',index=False)
    mean_z=np.mean(means,axis=0)
    heat=[dict(process_id=pid,process_label=labels[pid],realm=r,balanced_mean_process_z=mean_z[i,j]) for j,pid in enumerate(processes) for i,r in enumerate(realms)]
    pd.DataFrame(heat).to_csv(a.tables/'181_process_heatmap.tsv',sep='\t',index=False)
    draws[['design','iteration','selected_run_accessions']].to_csv(a.tables/'181_process_sampling_manifest.tsv',sep='\t',index=False)

    # One centroid per named locality. Country is retained to prevent collisions.
    # It does not split any canonical realm/core_locality unit in this cohort.
    assert (meta.groupby(['REALM','core_locality']).Country.nunique()==1).all()
    meta['locality_id']=meta.REALM+' | '+meta.Country+' | '+meta.core_locality
    units=sorted(meta.locality_id.unique())
    ix_units=[np.flatnonzero(meta.locality_id.to_numpy()==u) for u in units]
    centroids=np.vstack([x[ix].mean(axis=0) for ix in ix_units])
    global_centroid=centroids.mean(axis=0)
    turnover=np.sqrt(((centroids-global_centroid)**2).mean(axis=1))
    rng=np.random.default_rng(a.seed)
    boot=np.empty((a.n_bootstrap,len(units)))
    for b in range(a.n_bootstrap):
        cb=np.vstack([x[rng.choice(ix,len(ix),replace=True)].mean(axis=0) if len(ix)>1 else x[ix[0]] for ix in ix_units])
        boot[b]=np.sqrt(((cb-cb.mean(axis=0))**2).mean(axis=1))
    rows=[]
    for j,(u,ix) in enumerate(zip(units,ix_units)):
        m=meta.iloc[ix]
        coord=m.loc[m[['lat','lon']].notna().all(axis=1),['lat','lon']]
        lat=float(coord.lat.median()) if len(coord) else np.nan
        lon=float(coord.lon.median()) if len(coord) else np.nan
        rows.append(dict(locality_id=u,realm=m.REALM.iloc[0],country=m.Country.iloc[0],locality=m.core_locality.iloc[0],n_profiles=len(m),n_georeferenced_profiles=len(coord),n_coordinate_sites=len(coord.drop_duplicates()),lat=lat,lon=lon,mappable=bool(len(coord)),turnover_point_estimate=turnover[j],bootstrap_median=np.median(boot[:,j]),bootstrap_q05=np.quantile(boot[:,j],.05),bootstrap_q95=np.quantile(boot[:,j],.95)))
    locality=pd.DataFrame(rows)
    assert len(locality)==26 and locality.mappable.sum()==23 and locality.n_profiles.sum()==191
    assert locality.n_georeferenced_profiles.sum()==184
    locality.to_csv(a.tables/'182_shell_locality_turnover_summary.tsv',sep='\t',index=False)
    pd.DataFrame(boot,columns=units).rename_axis('bootstrap_iteration').to_csv(a.tables/'182_shell_locality_bootstrap.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
    meta.reset_index()[['run_accession','REALM','Country','core_locality','locality_id','lat','lon']].to_csv(a.tables/'182_sample_locality_manifest.tsv',sep='\t',index=False)
    np.savez_compressed(a.objects/'181_182_process_map_objects.npz',sample_ids=np.array(ids),shell_KOs=np.array(shell_ids),process_ids=np.array(processes),process_scores=scores,process_z=z,locality_ids=np.array(units),locality_CLR=centroids,global_centroid=global_centroid,turnover=turnover,turnover_bootstrap=boot)
    source_files=[a.counts,a.metadata,a.accessions,a.catalog,a.tables/'180_KO_master_classification.tsv',a.tables/'180_panel_c_prevalence_balanced_R2.tsv',Path(__file__)]
    report={'status':'PASS','n_profiles':191,'n_detected_KOs':17482,'n_shell_KOs':len(shell_ids),'n_catalog_KOs':len(catalog),'n_scored_KOs':len(used),'excluded_process_KO':'K00441','n_processes':len(processes),'n_primary_draws':len(means),'n_alternative_draws':499,'n_localities':len(units),'n_mapped_localities':int(locality.mappable.sum()),'n_mapped_profiles':int(locality.n_georeferenced_profiles.sum()),'unmapped_localities':locality.loc[~locality.mappable,'locality_id'].tolist(),'n_bootstrap':a.n_bootstrap,'seed':a.seed,'pseudocount':1,'process_score':'Equal-weight mean of selected KO coordinates in full Shell CLR, z-standardized across 191 profiles using sample SD','process_R2':'Univariate between-realm SS divided by total process-score SS; no process-level significance test','conditional_R2':'Unique realm SS after log10 total KO counts divided by total process-score SS; not df-adjusted','map_estimator':'RMS Aitchison distance from locality mean CLR to equally weighted centroid of all 26 localities','map_coordinates':'Median of complete profile coordinate pairs within each named locality; no invented coordinates or jitter','bootstrap_role':'Within-locality profile resampling and recomputation of global centroid; uncertainty only; point estimate uses all profiles','input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files}}
    (a.tables/'181_182_process_map_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
