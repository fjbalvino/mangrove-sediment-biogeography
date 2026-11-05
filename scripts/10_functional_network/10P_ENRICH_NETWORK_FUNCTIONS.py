#!/usr/bin/env python3
"""Conditional enrichment of the frozen repertoire-wide network annotation.

Uses SciPy exact hypergeometric probabilities and BY adjustment, scikit-learn
KMeans for covariate blocks, and common blockwise label permutations for maxT.
The graph is fixed. These are exploratory competitive gene-set tests, not
sample-level differential-abundance tests or evidence of biochemical activity.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.stats import hypergeom, false_discovery_control, binomtest
import sklearn
from sklearn.cluster import KMeans

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('broad_annotation', Path(__file__).with_name('10O_EXTEND_FUNCTIONAL_ANNOTATION.py'))
broad = importlib.util.module_from_spec(spec)
spec.loader.exec_module(broad)


def make_terms(annotation):
    terms = {}
    def add(key, name, kind, kos, source):
        if key in terms:
            raise ValueError('Duplicate term identifier: '+key)
        terms[key] = dict(term_id=key, term=name, kind=kind, KOs=set(kos), source=source)
    h = pd.read_csv(annotation/'KO_hierarchy_membership.tsv.gz', sep='\t').fillna('')
    for category, group in h[h.scope.eq('biological_pathway')].groupby('category'):
        add('CATEGORY:'+category[:5], category[6:], 'broad_category', group.KO,
            'https://www.kegg.jp/brite/ko00001')
    for term, group in h[h.scope.isin(['biological_pathway','protein_family'])].groupby('term'):
        match = re.search(r'\[(PATH|BR):(ko\d{5})\]', term)
        if not match:
            continue
        kind, tid = match.groups()
        if kind == 'PATH' and tid.startswith(('ko011','ko012','ko013')):
            continue
        name = re.sub(r'\s*\[(?:PATH|BR):.*\]', '', term[6:])
        add(('PATH:' if kind=='PATH' else 'BRITE:')+tid, name,
            'pathway' if kind=='PATH' else 'protein_family', group.KO,
            'https://www.kegg.jp/'+('pathway/' if kind=='PATH' else 'brite/')+tid)
    f = pd.read_csv(annotation/'KO_protein_family_detail.tsv.gz', sep='\t').fillna('')
    for (bid, path), group in f.groupby(['brite_id','family_path']):
        if not path:
            continue
        tid = 'DETAIL:'+bid+':'+hashlib.sha256(path.encode()).hexdigest()[:12]
        add(tid, path, 'detailed_family', group.KO, 'https://www.kegg.jp/brite/'+bid)
    m = pd.read_csv(annotation/'KO_module_logical_membership.tsv', sep='\t')
    for mid, group in m[~m.optional_only].groupby('module_id'):
        add('MODULE:'+mid, str(group.module_name.iloc[0]), 'biochemical_module', group.KO,
            'https://www.kegg.jp/entry/'+mid)
    return [terms[k] for k in sorted(terms)]


def match_blocks(frame, config):
    x = frame[config['covariates']].to_numpy(float)
    if not np.isfinite(x).all():
        raise ValueError('Missing/nonfinite matching covariates')
    med = np.median(x, axis=0)
    scale = np.median(np.abs(x-med), axis=0)
    scale = np.where(scale>0, scale, np.std(x, axis=0, ddof=1))
    keep = scale>0
    z = (x[:,keep]-med[keep])/scale[keep]
    k = min(max(1, len(frame)//config['block_target_size']), len(frame)//config['block_min_size'])
    if k <= 1 or not keep.any():
        return np.zeros(len(frame), dtype=int)
    labels = KMeans(n_clusters=k, n_init=config['kmeans_n_init'],
        max_iter=config['kmeans_max_iter'], random_state=config['seed'], algorithm='lloyd').fit_predict(z)
    while len(np.unique(labels)) > 1:
        levels, counts = np.unique(labels, return_counts=True)
        small = levels[counts<config['block_min_size']]
        if not len(small):
            break
        item = min(small, key=lambda c: (int((labels==c).sum()), int(c)))
        center = z[labels==item].mean(axis=0)
        others = levels[levels!=item]
        closest = min(others, key=lambda c: (float(np.linalg.norm(z[labels==c].mean(axis=0)-center)), int(c)))
        labels[labels==item] = closest
    # Stable numbering depends on member KO IDs, not arbitrary cluster indices.
    order = sorted(np.unique(labels), key=lambda c: tuple(frame.index[labels==c]))
    return np.array([order.index(c) for c in labels])


@lru_cache(maxsize=None)
def hypergeom_pmf(total, positives, draws):
    lo, hi = max(0, draws-(total-positives)), min(draws, positives)
    p = np.zeros(hi+1)
    p[lo:] = hypergeom.pmf(np.arange(lo,hi+1), total, positives, draws)
    return p/p.sum()


def exact_null(strata):
    """Independent conditional draws in disjoint strata; counts are summed."""
    distribution = np.ones(1)
    for total, positives, draws in strata:
        distribution = np.convolve(distribution, hypergeom_pmf(int(total),int(positives),int(draws)))
    distribution /= distribution.sum()
    x = np.arange(len(distribution))
    mean = float(x@distribution)
    variance = float(((x-mean)**2)@distribution)
    return distribution, mean, variance


def permutation_labels(labels, blocks, rng, count):
    result = np.tile(labels, (count,1))
    for block in np.unique(blocks):
        idx = np.flatnonzero(blocks==block)
        orders = np.argsort(rng.random((count,len(idx))), axis=1, kind='stable')
        result[:,idx] = labels[idx][orders]
    return result


def run_scenario(frame, terms, config, scenario, out):
    modules = sorted(set(frame.module_label)-{'outside'})
    labels = pd.Categorical(frame.module_label, categories=modules+['outside']).codes
    blocks = match_blocks(frame, config)
    frame = frame.copy(); frame['matching_block'] = blocks
    broad.write_tsv(frame.reset_index(), out/(scenario+'_KO_matching_blocks.tsv'))
    block_rows = []
    for b, group in frame.groupby('matching_block'):
        for cov in config['covariates']:
            block_rows.append(dict(block=int(b), n_KOs=len(group), n_labels=group.module_label.nunique(),
                covariate=cov, minimum=float(group[cov].min()), maximum=float(group[cov].max()),
                median=float(group[cov].median())))
    broad.write_tsv(pd.DataFrame(block_rows), out/(scenario+'_block_diagnostics.tsv'))
    A = np.array([[k in term['KOs'] for term in terms] for k in frame.index], dtype=np.int16)
    rows, means, sds = [], [], []
    for mi, module in enumerate(modules):
        observed = A[labels==mi].sum(axis=0)
        for ti, term in enumerate(terms):
            strata = [(int((blocks==b).sum()), int(A[blocks==b,ti].sum()), int(((blocks==b)&(labels==mi)).sum()))
                      for b in np.unique(blocks)]
            distribution, mean, var = exact_null(strata)
            n = int(observed[ti]); p = float(distribution[n:].sum())
            if not 0 < p <= 1+1e-12:
                raise ValueError('Invalid exact probability')
            sd = np.sqrt(var) if var>1e-14 else 0.
            z = (n-mean)/sd if sd else 0.
            rows.append(dict(scenario=scenario, statistical_module=module,
                term_id=term['term_id'], term=term['term'], kind=term['kind'],
                observed_KOs=n, module_KOs=int((labels==mi).sum()),
                percent_module=100*n/(labels==mi).sum(), background_KOs=len(frame),
                background_term_KOs=int(A[:,ti].sum()), expected_KOs=mean,
                conditional_fold=n/mean if mean>0 else np.nan,
                null_sd=sd, standardized_excess=z,
                null_count_p025=int(np.searchsorted(np.cumsum(distribution),.025)),
                null_count_p975=int(np.searchsorted(np.cumsum(distribution),.975)),
                p_exact=min(1.,p), exchangeability_variation=bool(sd), source=term['source']))
            means.append(mean);sds.append(sd)
    result = pd.DataFrame(rows)
    result['q_BY'] = false_discovery_control(result.p_exact.to_numpy(), method='by')
    mean = np.array(means).reshape(len(modules),len(terms))
    sd = np.array(sds).reshape(len(modules),len(terms))
    nonzero = sd>0
    safe_sd = np.where(nonzero,sd,1.)
    B = config['permutations']; rng = np.random.default_rng(config['seed']+1)
    maxima = []
    null_sums = np.zeros_like(mean); null_squares = np.zeros_like(mean)
    annotation = sparse.csr_matrix(A)
    for start in range(0,B,250):
        count = min(250,B-start)
        draws = permutation_labels(labels,blocks,rng,count)
        mx = np.zeros(count)
        for mi in range(len(modules)):
            counts = sparse.csr_matrix((draws==mi).astype(np.int16)).dot(annotation).toarray()
            zz = np.where(nonzero[mi],(counts-mean[mi])/safe_sd[mi],0.)
            mx = np.maximum(mx,zz.max(axis=1))
            null_sums[mi] += counts.sum(axis=0)
            null_squares[mi] += (counts.astype(float)**2).sum(axis=0)
        maxima.extend(mx)
    maxima = np.asarray(maxima)
    exceedances = [int(np.count_nonzero(maxima>=z-1e-12)) for z in result.standardized_excess]
    result['p_maxT_joint_null'] = [(1+value)/(B+1) for value in exceedances]
    intervals = [binomtest(value,B).proportion_ci(method='wilson') for value in exceedances]
    result['maxT_MC95_low'] = [(1+B*ci.low)/(B+1) for ci in intervals]
    result['maxT_MC95_high'] = [(1+B*ci.high)/(B+1) for ci in intervals]
    result['supported'] = result.q_BY.le(config['alpha']) & result.p_maxT_joint_null.le(config['alpha'])
    result['supported_with_MC_precision'] = result.q_BY.le(config['alpha']) & result.maxT_MC95_high.le(config['alpha'])
    result['empirical_null_mean'] = (null_sums/B).ravel()
    result['empirical_null_sd'] = np.sqrt(np.maximum(0,null_squares/B-(null_sums/B)**2)).ravel()
    result['n_hypotheses'] = len(result)
    broad.write_tsv(result, out/(scenario+'_enrichment.tsv'))
    broad.write_tsv(pd.DataFrame({'permutation':np.arange(1,B+1),'max_standardized_excess':maxima}),
                    out/(scenario+'_permutation_maxima.tsv.gz'))
    diagnostic = dict(n_KOs=len(frame), n_blocks=int(len(np.unique(blocks))),
        block_sizes=frame.groupby('matching_block').size().tolist(),
        n_sets=len(terms), hypotheses=len(result), raw_p05=int(result.p_exact.le(.05).sum()),
        BY_supported=int(result.q_BY.le(.05).sum()), jointly_supported=int(result.supported.sum()),
        max_abs_null_mean_difference=float(np.max(np.abs(null_sums/B-mean))),
        matching_method='MAD scaling; sklearn KMeans Lloyd with 50 starts; small clusters merged to nearest centroid',
        null='Exchangeable module labels within covariate blocks; preserves blockwise module sizes and complete annotation vectors')
    return result, diagnostic


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['annotation','nodes','candidates','edges','outdir']:
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--config',type=Path,default=ROOT/'config/network_functional_enrichment.json')
    args=p.parse_args();out=args.outdir
    if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be empty')
    out.mkdir(parents=True,exist_ok=True)
    config=json.loads(args.config.read_text());broad.write_json(out/'analysis_plan.json',config)
    nodes=pd.read_csv(args.nodes,sep='\t').set_index('KO')
    candidates=pd.read_csv(args.candidates,sep='\t').set_index('KO')
    edges=pd.read_csv(args.edges,sep='\t')
    if nodes.index.has_duplicates or candidates.index.has_duplicates:raise ValueError('Duplicate KO identifier')
    connected=set(edges.KO1)|set(edges.KO2)
    if not set(nodes.index)<=connected<=set(candidates.index):raise ValueError('Inconsistent KO universes')
    frame=candidates.copy();frame['module_label']=nodes.module_label.reindex(frame.index).fillna('outside')
    degree=pd.concat([edges.KO1,edges.KO2]).value_counts()
    strength=pd.concat([edges[['KO1','mean_abs_weight']].rename(columns={'KO1':'KO'}),
                        edges[['KO2','mean_abs_weight']].rename(columns={'KO2':'KO'})]).groupby('KO').mean_abs_weight.sum()
    frame['degree']=degree.reindex(frame.index).fillna(0)
    frame['strength_abs']=strength.reindex(frame.index).fillna(0)
    frame['log10_total_count']=np.log10(frame.total_count+1)
    frame=frame.sort_index()
    primary=frame.loc[frame.index.isin(nodes.index)]
    all_terms=make_terms(args.annotation)
    catalogue=[];eligible=[];major=set(primary.index)
    for term in all_terms:
        count=len(term['KOs']&major)
        ok=config['minimum_set_KOs_in_primary_universe']<=count<=len(primary)-config['minimum_nonset_KOs_in_primary_universe']
        catalogue.append({**term,'KOs':';'.join(sorted(term['KOs'])),'primary_universe_KOs':count,'eligible':ok})
        if ok:eligible.append(term)
    broad.write_tsv(pd.DataFrame(catalogue),out/'tested_term_catalogue.tsv.gz')
    # Eligibility is fixed from the union of modules, before inspecting any
    # module-specific count or P value; zero cells remain in both families.
    a,da=run_scenario(primary,eligible,config,'primary',out)
    b,db=run_scenario(frame.loc[frame.index.isin(connected)],eligible,config,'connected_sensitivity',out)
    comparison=a.merge(b[['statistical_module','term_id','q_BY','p_maxT_joint_null','supported','supported_with_MC_precision','maxT_MC95_low','maxT_MC95_high']],
        on=['statistical_module','term_id'],validate='one_to_one',suffixes=('','_sensitivity'))
    comparison['supported_both']=comparison.supported & comparison.supported_sensitivity
    broad.write_tsv(comparison,out/'enrichment_with_sensitivity.tsv')
    broad.write_tsv(comparison[comparison.supported],out/'supported_enrichments.tsv')
    summary=dict(primary=da,connected_sensitivity=db,
        supported_both=int(comparison.supported_both.sum()),
        primary_supported_unique_terms=int(comparison.loc[comparison.supported,'term_id'].nunique()),
        primary_supported_modules=sorted(comparison.loc[comparison.supported,'statistical_module'].unique()),
        scope='Conditional exploratory annotation enrichment; the same data selected the graph and its modules',
        original_candidates=len(candidates), original_edges=len(edges), original_connected_KOs=len(connected),
        versions=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__))
    broad.write_json(out/'summary.json',summary)
    inputs=[args.nodes,args.candidates,args.edges,args.config,*[args.annotation/f for f in
        ['KO_hierarchy_membership.tsv.gz','KO_protein_family_detail.tsv.gz','KO_module_logical_membership.tsv']]]
    broad.write_json(out/'provenance.json',dict(inputs={str(path):broad.digest(path) for path in inputs},
        source={str(path.relative_to(ROOT)):broad.digest(path) for path in
                [Path(__file__),Path(__file__).with_name('10O_EXTEND_FUNCTIONAL_ANNOTATION.py')]}))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
