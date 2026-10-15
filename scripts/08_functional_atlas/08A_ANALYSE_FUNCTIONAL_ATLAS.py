#!/usr/bin/env python3
"""Reproducible functional-atlas analyses; numpy, pandas and scipy only.

Inputs: original KO matrix, corrected sample metadata, paired rarefied taxonomic
counts. Exclusions are explicit; no iCAMP, read processing or implicit imputation.
Run with BLAS threads limited to one to avoid oversubscription.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, f as f_distribution
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261013
N_PERM = 9999
N_DRAWS = 1000
CLASSES = ['Absolute core', 'Near-global', 'Realm-core', 'Periphery']
EXCLUDED_AIR = ['SRR31095011', 'SRR31095534', 'SRR31095535', 'SRR31192804']


def log(message):
    print(time.strftime('[%H:%M:%S] ') + message, flush=True)


def save_json(path, obj):
    def encode(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)
    path.write_text(json.dumps(obj, indent=2, default=encode, ensure_ascii=False, allow_nan=False) + '\n')


def clr(counts):
    x = np.log1p(np.asarray(counts, dtype=float))
    return x - x.mean(axis=1, keepdims=True)


def pca(y):
    centered = y - y.mean(axis=0)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    return u * s, vt.T, s ** 2 / np.sum(s ** 2)


def reduced_coordinates(y):
    """An exact Euclidean embedding avoids permuting thousands of features."""
    z = y - y.mean(axis=0)
    gram = z @ z.T
    vals, vecs = np.linalg.eigh(gram)
    keep = vals > max(vals.max(), 1) * 1e-12
    return vecs[:, keep] * np.sqrt(vals[keep])


def orthogonal_basis(x):
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    tol = max(x.shape) * np.finfo(float).eps * max(s.max(), 1)
    return u[:, s > tol]


def permutation_model(y, groups, seed, nperm=N_PERM, covariates=None):
    """Euclidean PERMANOVA; reduced-model residual permutations if adjusted.

    F and SS equal multivariate linear-model partitioning in the Euclidean
    embedding. Conditional tests use Freedman-Lane residual permutations.
    R2 denotes the unique term SS / unadjusted total SS; partial_R2 also given.
    """
    y = reduced_coordinates(y)
    n = len(y)
    group_matrix = pd.get_dummies(pd.Series(groups), dtype=float).values
    base = np.ones((n, 1))
    if covariates is not None:
        covariates = np.asarray(covariates, dtype=float).reshape(n, -1)
        base = np.column_stack([base, covariates])
    q0 = orthogonal_basis(base)
    q1 = orthogonal_basis(np.column_stack([base, group_matrix]))
    qadd = orthogonal_basis(group_matrix - q0 @ (q0.T @ group_matrix))
    rank0, rank1 = q0.shape[1], q1.shape[1]
    df = rank1 - rank0
    if df == 0:
        return {'estimable': False, 'n': n, 'reason': 'No realm degrees of freedom after conditioning.'}
    # For near rank-deficient dummy columns, use the known model rank.
    residual_groups = group_matrix - q0 @ (q0.T @ group_matrix)
    u, _, _ = np.linalg.svd(residual_groups, full_matrices=False)
    qadd = u[:, :df]
    total = np.sum(y * y)
    sse = total - np.sum((q1.T @ y) ** 2)
    ss_term = np.sum((qadd.T @ y) ** 2)
    residual_df = n - rank1
    observed = (ss_term / df) / (sse / residual_df)
    e0 = y - q0 @ (q0.T @ y)
    e0_total = np.sum(e0 * e0)
    rng = np.random.default_rng(seed)
    permuted = []
    for start in range(0, nperm, 250):
        indexes = np.array([rng.permutation(n) for _ in range(min(250, nperm-start))])
        e = e0[indexes]
        ss1 = np.sum(np.einsum('ir,bic->brc', qadd, e, optimize=True) ** 2, axis=(1, 2))
        ssr = e0_total - np.sum(np.einsum('ir,bic->brc', q1, e, optimize=True) ** 2, axis=(1, 2))
        permuted.extend(((ss1 / df) / (ssr / residual_df)).tolist())
    null = np.asarray(permuted)
    return {'estimable': True, 'n': n, 'groups': len(np.unique(groups)), 'df': df,
            'residual_df': residual_df, 'pseudo_F': observed, 'R2': ss_term / total,
            'partial_R2': ss_term / (ss_term + sse), 'sum_of_squares': ss_term,
            'residual_sum_of_squares': sse, 'P': (np.count_nonzero(null >= observed-1e-10)+1)/(nperm+1),
            'permutations': nperm, 'seed': seed,
            'permutation_scheme': 'reduced-model residuals' if covariates is not None else 'site labels',
            'adjustment': 'mean within-site log10 KO annotation total' if covariates is not None else 'none'}


def geometric_median(y):
    init = y.mean(axis=0)
    def objective(x):
        delta = x - y
        norms = np.linalg.norm(delta, axis=1)
        gradient = (delta / np.maximum(norms[:, None], 1e-12)).sum(axis=0)
        return norms.sum(), gradient
    fit = minimize(objective, init, jac=True, method='BFGS', options={'gtol': 1e-8, 'maxiter': 2000})
    # BFGS can report loss of precision at a valid near-stationary solution.
    grad = np.linalg.norm(objective(fit.x)[1])
    if not fit.success and grad > 1e-4:
        raise RuntimeError('Spatial median failed to converge: ' + str(fit.message))
    return fit.x


def dispersion_test(y, groups, seed, nperm=N_PERM):
    """Spatial-median PERMDISP, sqrt(n/(n-1)) bias correction.

    Permutes within-group residuals of the one-way distance model, as in
    vegan::permutest.betadisper. Distances are in exact Euclidean CLR space.
    """
    y = reduced_coordinates(y)
    groups = np.asarray(groups)
    levels = np.unique(groups)
    n, k = len(y), len(levels)
    distances = np.zeros(n)
    means = {}
    for level in levels:
        mask = groups == level
        ng = mask.sum()
        if ng < 2:
            raise ValueError('Dispersion requires two or more sites per group.')
        center = geometric_median(y[mask])
        distances[mask] = np.linalg.norm(y[mask]-center, axis=1) * np.sqrt(ng/(ng-1))
        means[str(level)] = distances[mask].mean()
    gmat = pd.get_dummies(pd.Series(groups), dtype=float).values
    q = orthogonal_basis(gmat)
    fit = q @ (q.T @ distances)
    residuals = distances-fit
    observed = (np.sum((fit-distances.mean())**2)/(k-1))/(np.sum(residuals**2)/(n-k))
    rng = np.random.default_rng(seed)
    residual_ss = np.sum(residuals**2)
    exceed = 0
    for start in range(0, nperm, 500):
        ids = np.array([rng.permutation(n) for _ in range(min(500, nperm-start))])
        z = residuals[ids]
        mss = np.sum((z @ q)**2, axis=1) - n * z.mean(axis=1)**2
        fperm = (mss/(k-1))/((residual_ss-mss)/(n-k))
        exceed += np.count_nonzero(fperm >= observed - 1e-8)
    return {'n': n, 'df': k-1, 'residual_df': n-k, 'F': observed,
            'P': (exceed+1)/(nperm+1), 'parametric_P': f_distribution.sf(observed,k-1,n-k),
            'permutations': nperm, 'seed': seed, 'center': 'spatial median',
            'small_sample_bias_adjustment': True, 'mean_distance': means}, distances


def holm(values):
    values = np.asarray(values)
    ix = np.argsort(values)
    adjusted = np.empty(len(ix))
    adjusted[ix] = np.minimum(1, np.maximum.accumulate(values[ix]*(len(ix)-np.arange(len(ix)))))
    return adjusted


def classification(matrix, metadata):
    units = metadata.REALM.astype(str) + ' | ' + metadata.core_locality.astype(str)
    pooled = matrix.T.groupby(units, sort=True).sum()
    local_rel = pooled.div(pooled.sum(axis=1), axis=0)
    realm = metadata.REALM.groupby(units).first()
    support = pd.DataFrame(False, index=matrix.index, columns=sorted(realm.unique()))
    details = []
    for name in support.columns:
        indexes = realm.index[realm == name]
        nloc = len(indexes)
        present = pooled.loc[indexes].gt(0).sum()
        meanrel = local_rel.loc[indexes].mean()
        support[name] = (nloc >= 3) & present.eq(nloc) & meanrel.ge(0.0001)
        details.append(pd.DataFrame({'KO':matrix.index,'realm':name,'n_localities':nloc,
                                     'localities_present':present.values,'mean_relative_abundance':meanrel.values,
                                     'supports_realm_core':support[name].values}))
    prevalence_n = matrix.gt(0).sum(axis=1)
    labels = pd.Series('Periphery', index=matrix.index)
    labels.loc[support.any(axis=1)] = 'Realm-core'
    labels.loc[prevalence_n / matrix.shape[1] >= .95] = 'Near-global'
    labels.loc[prevalence_n == matrix.shape[1]] = 'Absolute core'
    assignments = matrix.sum(axis=1)
    ann = pd.DataFrame({'KO':matrix.index,'class':labels.values,'samples_present':prevalence_n.values,
                        'sample_prevalence':prevalence_n.values/matrix.shape[1], 'KO_assignment_count':assignments.values,
                        'supporting_realms':[';'.join(support.columns[row]) for row in support.values]})
    summary = ann.groupby('class').agg(KOs=('KO','size'), KO_assignments=('KO_assignment_count','sum')).reindex(CLASSES)
    summary['KO_richness_percent'] = 100*summary.KOs/matrix.shape[0]
    summary['KO_assignment_percent'] = 100*summary.KO_assignments/assignments.sum()
    return ann, summary, pd.concat(details, ignore_index=True), realm


def mantel(a, b, seed, nperm=N_PERM):
    """Pearson Mantel; upper-tail inference matches the original vegan test."""
    i,j = np.triu_indices(len(a),1)
    x = a[i,j]; x = x-x.mean()
    y = b[i,j]; y = y-y.mean()
    denom = np.linalg.norm(x)*np.linalg.norm(y)
    observed = np.dot(x,y)/denom
    rng = np.random.default_rng(seed)
    upper = 0; two_sided = 0
    for _ in range(nperm):
        p = rng.permutation(len(a))
        r = np.dot(x, b[p[i],p[j]]) / denom
        upper += r >= observed-1e-12
        two_sided += abs(r) >= abs(observed)-1e-12
    return {'n': len(a), 'pairs': len(i), 'r': observed, 'P':(upper+1)/(nperm+1),
            'P_two_sided':(two_sided+1)/(nperm+1), 'permutations':nperm, 'seed':seed}


def haversine(latlon):
    lat, lon = np.deg2rad(np.asarray(latlon).T)
    a = np.sin((lat[:,None]-lat)/2)**2 + np.cos(lat[:,None])*np.cos(lat)*np.sin((lon[:,None]-lon)/2)**2
    return 6371.0088 * 2*np.arcsin(np.sqrt(np.clip(a,0,1)))


def slope_data(tax_similarity, func_similarity, geo):
    i,j=np.triu_indices(len(geo),1)
    valid=geo[i,j]>.1
    i,j=i[valid],j[valid]
    x=np.log10(geo[i,j]); xc=x-x.mean()
    st=tax_similarity[i,j]; sf=func_similarity[i,j]
    bt=np.dot(xc,st)/np.dot(xc,xc); bf=np.dot(xc,sf)/np.dot(xc,xc)
    return {'slope_taxonomic':bt,'slope_functional':bf,'difference_tax_minus_function':bt-bf,
            'intercept_taxonomic':st.mean()-bt*x.mean(),'intercept_functional':sf.mean()-bf*x.mean(),
            'pairs':len(i)}, (i,j,x,xc,st,sf)


def slope_test(tax_similarity, func_similarity, geo, seed, nperm=N_PERM):
    result, arrays=slope_data(tax_similarity,func_similarity,geo)
    i,j,x,xc,st,sf=arrays
    scale=np.dot(xc,xc)
    delta=tax_similarity-func_similarity
    observed=result['difference_tax_minus_function']
    rng=np.random.default_rng(seed)
    nt=nf=nd=0
    for _ in range(nperm):
        p=rng.permutation(len(geo))
        bt=np.dot(xc,tax_similarity[p[i],p[j]])/scale
        bf=np.dot(xc,func_similarity[p[i],p[j]])/scale
        nd += abs(bt-bf) >= abs(observed)-1e-12
        nt += abs(bt) >= abs(result['slope_taxonomic'])-1e-12
        nf += abs(bf) >= abs(result['slope_functional'])-1e-12
    result.update({'P_difference_two_sided':(nd+1)/(nperm+1), 'P_taxonomic_slope_two_sided':(nt+1)/(nperm+1),
                   'P_functional_slope_two_sided':(nf+1)/(nperm+1),'n':len(geo),'permutations':nperm,'seed':seed})
    return result, pd.DataFrame({'log10_distance_km':x,'taxonomic_similarity':st,'functional_similarity':sf})


def group_mean(values, groups, levels):
    return np.array([values[np.asarray(groups)==level].mean(axis=0) for level in levels])


def simple_r2(y, labels):
    total=np.sum((y-y.mean(axis=0))**2)
    residual=sum(np.sum((y[np.asarray(labels)==g]-y[np.asarray(labels)==g].mean(axis=0))**2) for g in np.unique(labels))
    return (total-residual)/total


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,default=ROOT/'data/functional')
    parser.add_argument('--output',type=Path,default=ROOT/'results/functional')
    parser.add_argument('--sample-manifest', type=Path, default=ROOT/'config/samples/functional_accessions.txt')
    args=parser.parse_args()
    out=args.output; out.mkdir(parents=True,exist_ok=True)
    original=pd.read_csv(args.inputs/'functional_KO_counts_original.tsv.gz',sep='\t',index_col=0)
    meta=pd.read_csv(args.inputs/'metadata_functional.csv',low_memory=False).set_index('run_accession')
    omitted=sorted(set(original.columns)-set(meta.index))
    assert set(omitted).issubset(EXCLUDED_AIR), 'Unexpected sample exclusion; resolve explicitly.'
    ids=[x for x in original.columns if x in meta.index]
    expected={line.strip() for line in args.sample_manifest.read_text().splitlines() if line.strip()}
    if set(ids) != expected or len(ids) != len(expected):
        raise ValueError('Functional sample IDs do not match the declared accession manifest.')
    meta=meta.loc[ids]
    assert meta.index.is_unique and not meta[['REALM','core_locality']].isna().any().any()
    matrix=original[ids]
    empty=matrix.index[matrix.sum(axis=1)==0].tolist()
    matrix=matrix.loc[matrix.sum(axis=1)>0]
    assert matrix.index.is_unique and np.isfinite(matrix.values).all() and (matrix.values>=0).all()
    z=np.load(args.inputs/'taxonomic_matched_counts.npz',allow_pickle=False)
    tax=pd.DataFrame(z['counts'].T,index=z['samples'],columns=z['features']).loc[ids]
    assert (tax.sum(axis=1)==18145).all() and (tax.sum(axis=0)>0).all()
    meta['KO_assignment_total']=matrix.sum(axis=0)
    meta['KO_richness']=matrix.gt(0).sum(axis=0)
    meta['log10_KO_assignment_total']=np.log10(meta.KO_assignment_total)
    meta['functional_analysis_included']=True
    geo_ok=meta[['lat','lon']].notna().all(axis=1)&meta.lat.between(-90,90)&meta.lon.between(-180,180)
    assert meta.loc[geo_ok,'sampling_site_id'].notna().all()
    meta['geographic_analysis_included']=geo_ok
    meta.to_csv(out/'sample_inclusion_and_metadata.tsv',sep='\t')
    log('Classifying functional persistence in {} profiles'.format(len(ids)))
    ann,classes,support,locality_realms=classification(matrix,meta)
    ann.to_csv(out/'KO_persistence_classes.tsv',sep='\t',index=False)
    classes.to_csv(out/'functional_class_summary.tsv',sep='\t')
    support.to_csv(out/'realm_core_support.tsv',sep='\t',index=False)
    x=matrix.T.values.astype(float)
    func_clr=clr(x); tax_clr=clr(tax.values)
    scores,loadings,explained=pca(func_clr)
    # Orient axes deterministically, without changing distances or inference.
    if spearmanr(scores[:,0],meta.log10_KO_assignment_total).statistic < 0:
        scores[:,0]*=-1;loadings[:,0]*=-1
    top=np.argsort(-np.abs(loadings[:,0]))[:100]
    class_labels=ann.set_index('KO')['class'].loc[matrix.index].values
    loading_table=pd.DataFrame({'KO':matrix.index,'class':class_labels,'PC1_loading':loadings[:,0],
                                'PC2_loading':loadings[:,1],'top100_absolute_PC1':np.isin(np.arange(len(matrix)),top)})
    loading_table.to_csv(out/'PCA_KO_loadings.tsv',sep='\t',index=False)
    top_contrib=loading_table.iloc[top].assign(weight=lambda d:d.PC1_loading.abs()).groupby('class').weight.sum().reindex(CLASSES,fill_value=0)
    top_contrib=100*top_contrib/top_contrib.sum()
    rho=spearmanr(scores[:,0],meta.KO_assignment_total).statistic
    rho_rich=spearmanr(scores[:,0],meta.KO_richness).statistic
    yhell=np.sqrt(x/x.sum(axis=1,keepdims=True))
    hell_scores,_,hell_explained=pca(yhell)
    diag={'PC1_variance_percent':100*explained[0],'PC2_variance_percent':100*explained[1],
          'PC1_spearman_rho_KO_assignment_total':rho,'PC1_spearman_rho_KO_richness':rho_rich,
          'top100_absolute_PC1_loading_percent_by_class':top_contrib.to_dict(),
          'Hellinger_PC1_variance_percent':100*hell_explained[0],
          'Hellinger_PC2_variance_percent':100*hell_explained[1],
          'Hellinger_PC1_abs_rho_KO_assignment_total':abs(spearmanr(hell_scores[:,0],meta.KO_assignment_total).statistic)}
    save_json(out/'PCA_diagnostics.json',diag)
    pcdf=meta[['REALM','Country','core_locality','sampling_site_id','KO_assignment_total','KO_richness']].copy()
    pcdf['PC1']=scores[:,0];pcdf['PC2']=scores[:,1]
    pcdf['Hellinger_PC1']=hell_scores[:,0];pcdf['Hellinger_PC2']=hell_scores[:,1]
    pcdf.to_csv(out/'PCA_sample_scores.tsv',sep='\t')
    log('Testing regional differences on one centroid per coordinate-defined site')
    gm=meta.loc[geo_ok].copy()
    levels=sorted(gm.sampling_site_id.unique())
    gx=np.flatnonzero(geo_ok.values)
    groups=gm.sampling_site_id.values
    for col in ['REALM','Country']:
        assert gm.groupby('sampling_site_id')[col].nunique().max()==1, 'Conflicting '+col+' within a site'
    sm=gm.groupby('sampling_site_id').agg(REALM=('REALM','first'),Country=('Country','first'),study_id=('study_id',lambda x:';'.join(sorted(set(x)))),n_studies=('study_id','nunique'),
        lat=('lat','median'),lon=('lon','median'),samples=('REALM','size'),log10_KO_assignment_total=('log10_KO_assignment_total','mean')).loc[levels]
    sy=group_mean(func_clr[gx],groups,levels)
    sty=group_mean(tax_clr[gx],groups,levels)
    shy=group_mean(yhell[gx],groups,levels)
    counts=sm.REALM.value_counts()
    qualifying=sorted(counts.index[counts>=5])
    sm.to_csv(out/'site_metadata.tsv',sep='\t')
    results=[]; dispersions=[]; balance=[]
    comparisons=[('Three replicated realms',qualifying),('Tropical Atlantic vs Central Indo-Pacific',['Tropical Atlantic','Central Indo-Pacific']),
                 ('Tropical Atlantic vs Western Indo-Pacific',['Tropical Atlantic','Western Indo-Pacific']),
                 ('Central Indo-Pacific vs Western Indo-Pacific',['Central Indo-Pacific','Western Indo-Pacific'])]
    for c,(name,realms) in enumerate(comparisons):
        keep=sm.REALM.isin(realms).values
        labels=sm.loc[keep,'REALM'].values
        for method,y,covar in [('Aitchison',sy,None),('Aitchison adjusted for KO total',sy,sm.loc[keep,'log10_KO_assignment_total'].values),('Hellinger',shy,None)]:
            fit=permutation_model(y[keep],labels,SEED+100+c*10+len(results),covariates=covar)
            results.append({'comparison':name,'geometry':method,**fit})
        disp,distances=dispersion_test(sy[keep],labels,SEED+200+c)
        dispersions.append({'comparison':name,**{k:v for k,v in disp.items() if k!='mean_distance'}})
        pd.DataFrame({'site':sm.index[keep],'realm':labels,'distance_to_spatial_median':distances}).to_csv(out/('dispersion_distances_'+str(c)+'.tsv'),sep='\t',index=False)
        # Equal-site draws quantify effect-size stability without treating draws as new experiments.
        yy=reduced_coordinates(sy[keep]);rng=np.random.default_rng(SEED+300+c)
        pools=[np.flatnonzero(labels==label) for label in np.unique(labels)];nmin=min(map(len,pools))
        vals=[]
        for _ in range(N_DRAWS):
            take=np.concatenate([rng.choice(p,nmin,replace=False) for p in pools])
            vals.append(simple_r2(yy[take],labels[take]))
        balance.append({'comparison':name,'sites_per_realm':nmin,'draws':N_DRAWS,'median_R2':np.median(vals),
                        'R2_q025':np.quantile(vals,.025),'R2_q975':np.quantile(vals,.975),'min_R2':min(vals),'max_R2':max(vals)})
        log(name+' completed')
    perf=pd.DataFrame(results)
    perf['Holm_P_pairwise_within_geometry']=np.nan
    for geometry in perf.geometry.unique():
        mask=(perf.geometry==geometry)&(perf.comparison!='Three replicated realms')
        perf.loc[mask,'Holm_P_pairwise_within_geometry']=holm(perf.loc[mask,'P'])
    perf.to_csv(out/'functional_PERMANOVA.tsv',sep='\t',index=False)
    pd.DataFrame(dispersions).to_csv(out/'functional_PERMDISP.tsv',sep='\t',index=False)
    pd.DataFrame(balance).to_csv(out/'equal_site_effect_size_sensitivity.tsv',sep='\t',index=False)
    # Realm cannot be tested independently of study if each study occurs in only one realm.
    study_realms=gm.groupby('study_id').REALM.nunique()
    confounding={'studies_spanning_more_than_one_realm':int((study_realms>1).sum()),
                 'study_blocked_realm_test_estimable':bool((study_realms>1).any()),
                 'sites_containing_multiple_study_identifiers':int((sm.n_studies>1).sum()),
                 'reason':'Realm labels are invariant within study blocks; no cross-realm exchanges are available.'}
    save_json(out/'study_realm_confounding.json',confounding)
    coverage=gm[['sampling_site_id','study_id','REALM']].drop_duplicates()
    pd.crosstab(coverage.study_id,coverage.REALM).to_csv(out/'study_realm_site_coverage.tsv',sep='\t')
    pd.crosstab(sm.Country,sm.REALM).to_csv(out/'country_realm_site_coverage.tsv',sep='\t')
    log('Testing taxonomic-functional coupling and geographic distance-decay')
    da=squareform(pdist(tax_clr));db=squareform(pdist(func_clr))
    coupling={'metagenome_level':mantel(da,db,SEED+400),
              'site_centroid_level':mantel(squareform(pdist(sty)),squareform(pdist(sy)),SEED+401)}
    save_json(out/'taxonomic_functional_coupling.json',coupling)
    tr=tax.values/tax.values.sum(axis=1,keepdims=True)
    fr=x/x.sum(axis=1,keepdims=True)
    tsim=1-squareform(pdist(tr[gx],metric='braycurtis'))
    fsim=1-squareform(pdist(fr[gx],metric='braycurtis'))
    geo=haversine(gm[['lat','lon']].values)
    decay,pairs=slope_test(tsim,fsim,geo,SEED+500)
    pairs.to_csv(out/'metagenome_distance_decay_pairs.tsv.gz',sep='\t',index=False,compression='gzip')
    str_=group_mean(tr[gx],groups,levels);sfr=group_mean(fr[gx],groups,levels)
    sitegeo=haversine(sm[['lat','lon']].values)
    sdecay,spairs=slope_test(1-squareform(pdist(str_,metric='braycurtis')),1-squareform(pdist(sfr,metric='braycurtis')),sitegeo,SEED+501)
    spairs.to_csv(out/'site_distance_decay_pairs.tsv',sep='\t',index=False)
    save_json(out/'geographic_distance_decay.json',{'metagenome_level':decay,'site_equal_weight_profiles':sdecay})
    log('Drawing one paired metagenome per site, 1000 times')
    rng=np.random.default_rng(SEED+600)
    pools=[np.flatnonzero(groups==level) for level in levels]
    draws=[]
    for b in range(N_DRAWS):
        take=np.array([rng.choice(pool) for pool in pools])
        ds,_=slope_data(tsim[np.ix_(take,take)],fsim[np.ix_(take,take)],geo[np.ix_(take,take)])
        idx=gx[take];ii,jj=np.triu_indices(len(idx),1)
        ds['Mantel_r']=np.corrcoef(da[idx[ii],idx[jj]],db[idx[ii],idx[jj]])[0,1]
        draws.append({'draw':b+1,**ds})
    draws=pd.DataFrame(draws)
    draws.to_csv(out/'one_metagenome_per_site_draws.tsv',sep='\t',index=False)
    draw_summary={}
    for col in ['slope_taxonomic','slope_functional','difference_tax_minus_function','Mantel_r']:
        v=draws[col]
        draw_summary[col]={'median':v.median(),'q025':v.quantile(.025),'q975':v.quantile(.975),
                           'minimum':v.min(),'maximum':v.max(),'fraction_negative':(v<0).mean()}
    save_json(out/'one_metagenome_per_site_summary.json',draw_summary)
    ii,jj=np.triu_indices(len(meta),1)
    pd.DataFrame({'taxonomic_Aitchison':da[ii,jj],'functional_Aitchison':db[ii,jj]}).to_csv(out/'coupling_pairs.tsv.gz',sep='\t',index=False,compression='gzip')
    sitei,sitej=np.triu_indices(len(sm),1)
    pd.DataFrame({'taxonomic_Aitchison':squareform(pdist(sty))[sitei,sitej],
                  'functional_Aitchison':squareform(pdist(sy))[sitei,sitej]}).to_csv(out/'coupling_site_pairs.tsv',sep='\t',index=False)
    summary={'input_functional_profiles':len(original.columns),'excluded_air_profiles':omitted,'functional_profiles':len(meta),
             'input_KOs':len(original),'KOs_present_after_exclusion':len(matrix),'KOs_absent_after_air_exclusion':empty,
             'realm_profiles':meta.REALM.value_counts().to_dict(),'realm_localities':locality_realms.value_counts().to_dict(),
             'georeferenced_profiles':int(geo_ok.sum()),'coordinate_sites':len(sm),'realm_sites':counts.to_dict(),
             'geographic_exclusions':meta.index[~geo_ok].tolist(),'primary_realm_test_sites':int(sm.REALM.isin(qualifying).sum()),
             'eligible_test_realms':qualifying,'paired_taxonomic_features':tax.shape[1], 'taxonomy_rarefaction':18145,
             'classes':classes.reset_index().to_dict(orient='records'), 'PCA':diag,'coupling':coupling,'decay':{'metagenome':decay,'site':sdecay},
             'one_sample_per_site_sensitivity':draw_summary,'study_realm_confounding':confounding,
             'functional_unit':'Counts of assembled protein-cluster members assigned to KO terms; multiple KO annotations receive full counts.',
             'iCAMP_rerun':False,'permutations':N_PERM,'resampling_draws':N_DRAWS,'seed':SEED}
    save_json(out/'functional_results_summary.json',summary)
    import importlib.metadata
    save_json(out/'software_versions.json', {name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy']})
    manifest=[]
    for p in [Path(__file__),args.sample_manifest,*sorted(args.inputs.iterdir())]:
        if p.is_file():manifest.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    save_json(out/'input_provenance.json',manifest)
    log('FUNCTIONAL_ANALYSIS_PASS: '+str(out.resolve()))


if __name__=='__main__':
    main()
