#!/usr/bin/env python3
"""Calibrate reference-chain coherence under the observed sampling hierarchy.

No graph edges are added. Exact block sign flips require independent provenance
blocks and joint sign symmetry, not merely zero average covariance. Recovery
stress scenarios are distinct from statistical-null calibration.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import pandas as pd
import scipy
from scipy.stats import binomtest
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE / 'lib'))
import biogeochemical_context as context
spec = importlib.util.spec_from_file_location('network_pilot', HERE / '10L_AUDIT_PROCESS_CONTEXT.py')
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


def write_json(path, value):
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def write_table(frame, path):
    tmp = path.with_name(path.name + '.tmp')
    frame.to_csv(tmp, sep='\t', index=False, compression='gzip' if str(path).endswith('.gz') else None)
    tmp.replace(path)


def holm(p):
    p = np.asarray(p, float)
    if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any():
        raise ValueError('Invalid p-values')
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    result = np.empty_like(p)
    result[order] = adjusted
    return result


def exact_signs(n):
    if not 1 <= n <= 20:
        raise ValueError('Exact block enumeration is implemented for 1–20 informative labels; redesign for larger cohorts')
    return 2 * ((np.arange(2**n, dtype=np.uint32)[:, None] >> np.arange(n)) & 1).astype(float) - 1


def signflip_test(contributions, signs=None):
    c = np.asarray(contributions, float)
    if c.ndim == 1:
        c = c[:, None]
    if not np.isfinite(c).all():
        raise ValueError('Nonfinite block contributions')
    signs = exact_signs(len(c)) if signs is None else signs
    observed = c.mean(axis=0)
    null = signs @ c / len(c)
    # All transformations, including the identity, are enumerated exactly.
    tolerance = 100 * np.finfo(float).eps * np.maximum(1, np.abs(observed))
    p = np.mean(null >= observed - tolerance, axis=0)
    return observed, p


def aggregation(meta):
    """Equal BioSample weight within study-coordinate; equal runs within label."""
    groups = sorted(meta.analysis_group.unique())
    a = np.zeros((len(groups), len(meta)))
    for i, group in enumerate(groups):
        positions = np.flatnonzero(meta.analysis_group.to_numpy() == group)
        labels = meta.iloc[positions].biosample_group.to_numpy()
        for label in np.unique(labels):
            ids = positions[labels == label]
            a[i, ids] = 1 / (len(np.unique(labels)) * len(ids))
    assert np.allclose(a.sum(axis=1), 1)
    g = meta.groupby('analysis_group').agg(study_id=('study_id','first'),
        lat=('lat','first'), lon=('lon','first'), n_runs=('study_id','size'),
        n_biosamples=('biosample_group','nunique')).loc[groups]
    return a, g


def chain_pairs(chains, rules, feature_ids):
    index = {k: i for i, k in enumerate(feature_ids)}
    transitions, pair_set, membership = [], set(), []
    for chain in chains:
        steps, seen = [], set()
        for step_number, step in enumerate(chain['steps'], 1):
            expression = context.chain_expression(step, rules)
            expected = context.logic.kos(context.logic.parse(expression))
            if seen & expected:
                raise ValueError('Shared KO across chain steps needs explicit redesign: ' + chain['id'])
            seen |= expected
            available = sorted(expected & set(index))
            steps.append(available)
            membership.append(dict(chain=chain['id'], step=step_number, label=step['label'],
                expression=expression, available_KOs=';'.join(available),
                unavailable_KOs=';'.join(sorted(expected - set(index))),
                quantitative_meaning='Measured marker pairs; no enzyme or pathway abundance is inferred'))
        chain_transitions = []
        for j in range(len(steps)-1):
            pairs = sorted({tuple(sorted((index[k], index[l]))) for k in steps[j] for l in steps[j+1]})
            chain_transitions.append(pairs)
            pair_set.update(pairs)
        transitions.append(chain_transitions)
    pairs = sorted(pair_set)
    lookup = {p:i for i,p in enumerate(pairs)}
    weights = np.zeros((len(pairs), len(chains)))
    estimable = []
    for j, links in enumerate(transitions):
        ok = len(links) > 0 and all(len(link) > 0 for link in links)
        estimable.append(ok)
        if ok:
            for link in links:
                for pair in link:
                    weights[lookup[pair], j] += 1 / (len(links) * len(link))
    return np.asarray(pairs, int), weights, np.array(estimable), pd.DataFrame(membership)


def block_contributions(y, study_labels, pairs, weights, minimum_coordinates=2):
    y = np.asarray(y, float)
    labels = np.asarray(study_labels)
    studies = [s for s in sorted(set(labels)) if (labels == s).sum() >= minimum_coordinates]
    if not studies:
        raise ValueError('No within-study coordinate contrasts')
    covariances, variances = [], []
    for s in studies:
        x = y[labels == s]
        z = x - x.mean(axis=0)
        variances.append(np.mean(z*z, axis=0))
        covariances.append(np.mean(z[:,pairs[:,0]] * z[:,pairs[:,1]], axis=0))
    v = np.mean(variances, axis=0)
    denom = v[pairs[:,0]] + v[pairs[:,1]]
    contribution = np.zeros((len(studies), len(pairs)))
    np.divide(2*np.array(covariances), denom, out=contribution, where=denom > 1e-15)
    # A pair with two invariant markers has zero coherence, never evidence of
    # a strong association. Its prespecified weight is not redistributed.
    return contribution @ weights, studies, contribution.mean(axis=0), denom > 1e-15


def transformed(counts, reference_mask, feature_positions, pseudocount=1):
    logged = np.log(counts + pseudocount)
    return logged[:, feature_positions] - logged[:, reference_mask].mean(axis=1, keepdims=True)


def spatial_cholesky(groups, fraction):
    lat, lon = np.deg2rad(groups.lat.to_numpy()), np.deg2rad(groups.lon.to_numpy())
    xyz = np.column_stack([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)])
    result = {}
    for study in sorted(groups.study_id.unique()):
        ix = np.flatnonzero(groups.study_id.to_numpy() == study)
        distance = cdist(xyz[ix], xyz[ix])
        positive = distance[distance > 1e-12]
        scale = np.median(positive) if len(positive) else 1
        kernel = fraction*np.exp(-distance/scale) + (1-fraction)*np.eye(len(ix))
        result[study] = (ix, np.linalg.cholesky(kernel))
    return result


def simulate_gaussian(rng, scenario, means, sd, groups, kernels, target, cfg):
    z = rng.normal(size=means.shape)
    if scenario == 'gaussian_spatial':
        for ix, chol in kernels.values():
            z[ix] = chol @ z[ix]
    elif scenario == 't3_unequal_precision':
        z = rng.standard_t(3, size=means.shape) / np.sqrt(3)
        z /= np.sqrt(groups.n_biosamples.to_numpy())[:, None]
    elif scenario == 'gaussian_positive_control':
        rho = cfg['positive_control_correlation']
        factor = rng.normal(size=(len(groups), 1))
        z[:, target] = np.sqrt(rho)*factor + np.sqrt(1-rho)*z[:, target]
    return means + sd*z


def simulate_counts(rng, scenario, base_rates, run_study, run_group, run_bio,
                    depths, n_groups, n_bio, target, reference_mask, feature_positions, a, cfg):
    p = base_rates.shape[1]
    z = rng.normal(size=(n_groups, p))
    if scenario == 'counts_positive_control':
        rho = cfg['positive_control_correlation']
        factor = rng.normal(size=(n_groups, 1))
        z[:, target] = np.sqrt(rho)*factor + np.sqrt(1-rho)*z[:, target]
    cs, bs = cfg['coordinate_log_sd'], cfg['biosample_log_sd']
    log_multiplier = cs*z[run_group] + bs*rng.normal(size=(n_bio,p))[run_bio] - .5*(cs*cs+bs*bs)
    expected = base_rates[run_study] * np.exp(log_multiplier)
    if scenario == 'counts_extra_dropout':
        expected *= (rng.random((n_groups,p)) >= cfg['dropout_probability'])[run_group]
    yield_per_run = np.full(len(depths), np.median(depths)) if scenario == 'counts_equal_recovery' else depths
    counts = rng.poisson(expected * yield_per_run[:,None])
    y = a @ transformed(counts, reference_mask, feature_positions)
    return y, float(np.mean(counts == 0)), float(np.median(counts.sum(axis=1)))


def analyse(y, groups, pairs, weights, estimable, chain_ids, minimum_coordinates=2):
    c, studies, rho, available = block_contributions(y, groups.study_id, pairs, weights, minimum_coordinates)
    score, p = signflip_test(c)
    p[~estimable] = 1
    adjusted = holm(p)
    table = pd.DataFrame(dict(chain=chain_ids, score=score, p_assumption_conditional=p,
        p_holm_assumption_conditional=adjusted, study_labels=len(studies),
        positive_block_contributions=(c>1e-12).sum(axis=0),
        negative_block_contributions=(c<-1e-12).sum(axis=0),
        nonzero_block_contributions=(np.abs(c)>1e-12).sum(axis=0),
        sign_configurations=2**len(studies), measurable_chain=estimable))
    return table, c, studies, rho


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ('matrix','metadata','samples','reference-edges','reference-candidates','outdir'):
        ap.add_argument('--'+key, type=Path, required=True)
    ap.add_argument('--config', type=Path, default=ROOT/'config/network_calibration.json')
    ap.add_argument('--curation', type=Path, default=ROOT/'config/gene_biogeography')
    ap.add_argument('--replicates', type=int, help='Explicit smoke-test override; incomplete campaigns cannot pass the release gate')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    plan = json.loads(args.config.read_text())
    reps = args.replicates or plan['simulation_replicates']
    if reps < 1 or not 1 <= args.workers <= 8:
        ap.error('Positive replicate count and 1–8 workers required')
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()) and not args.resume:
        raise ValueError('Use an empty output directory or an identity-checked --resume')
    out.mkdir(parents=True,exist_ok=True)
    fingerprints = {k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in vars(args).items()
                    if isinstance(v,Path) and v.is_file()}
    for name in ('operational_rules.tsv','network_reference_chains.json'):
        fingerprints[name] = hashlib.sha256((args.curation/name).read_bytes()).hexdigest()
    fingerprints['source'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    identity = dict(inputs=fingerprints,replicates=reps,seed=plan['seed'])
    if (out/'run_identity.json').exists():
        if json.loads((out/'run_identity.json').read_text()) != identity:
            raise ValueError('Resume identity mismatch')
    else:
        write_json(out/'run_identity.json',identity)
    write_json(out/'calibration_plan.json',plan)
    counts = context.read_matrix(args.matrix,args.samples)
    meta = pd.read_csv(args.metadata)
    if meta.run_accession.isna().any() or meta.run_accession.duplicated().any():
        raise ValueError('Metadata accessions must be nonmissing and unique')
    meta = meta.set_index('run_accession').loc[counts.columns].copy()
    if meta.study_id.isna().any() or not np.isfinite(meta[['lat','lon']]).all().all():
        raise ValueError('Study labels and coordinates required')
    if (meta.lat.abs()>90).any() or (meta.lon.abs()>180).any():
        raise ValueError('Coordinates outside geographic bounds')
    meta['coordinate_group'] = [f'coord_{la:.5f}_{lo:.5f}' for la,lo in zip(meta.lat,meta.lon)]
    meta['analysis_group'] = meta.study_id+'::'+meta.coordinate_group
    meta['biosample_group'] = [str(b) if pd.notna(b) and str(b).strip() else 'run:'+r
                             for r,b in zip(meta.index,meta.sample_accession)]
    a, groups = aggregation(meta)
    rules = pd.read_csv(args.curation/'operational_rules.tsv',sep='\t').set_index('rule_id')
    chains = json.loads((args.curation/'network_reference_chains.json').read_text())['chains']
    chain_ids = [c['id'] for c in chains]
    edges = pd.read_csv(args.reference_edges,sep='\t')
    candidates = pd.read_csv(args.reference_candidates,sep='\t')
    if candidates.KO.duplicated().any() or not set(candidates.KO).issubset(counts.index):
        raise ValueError('Invalid reference candidates')
    if not (set(edges.KO1)|set(edges.KO2)).issubset(candidates.KO):
        raise ValueError('Historical edges outside supplied candidate universe')
    candidate_pairs = len(candidates)*(len(candidates)-1)//2
    queried = set()
    for chain in chains:
        for step in chain['steps']:
            queried |= context.logic.kos(context.logic.parse(context.chain_expression(step,rules)))
    features = sorted((queried | set(edges.KO1) | set(edges.KO2)) & set(counts.index))
    positions = counts.index.get_indexer(features)
    reference = ((counts.gt(0).sum(axis=1)>=10) & (counts.sum(axis=1)>=10)).to_numpy()
    raw = counts.T.to_numpy(float)
    pairs, weights, estimable, membership = chain_pairs(chains,rules,features)
    write_table(membership,out/'chain_marker_contract.tsv')
    write_table(groups.reset_index(),out/'analysis_groups.tsv')
    y = a @ transformed(raw, reference, positions)
    observed, contributions, studies, pair_values = analyse(y,groups,pairs,weights,estimable,chain_ids)
    write_table(pd.DataFrame(contributions,index=studies,columns=chain_ids).reset_index(names='study_id'),out/'observed_block_contributions.tsv')
    write_table(pd.DataFrame(dict(KO1=np.array(features)[pairs[:,0]],KO2=np.array(features)[pairs[:,1]],within_study_rho=pair_values)),out/'adjacent_marker_pair_estimates.tsv')
    sensitivity = []
    for pc in (.5,1,2):
        yy = a @ transformed(raw, reference, positions,pc)
        t,_,_,_ = analyse(yy,groups,pairs,weights,estimable,chain_ids)
        t['scenario']=f'pseudocount_{pc:g}';sensitivity.append(t)
    yy = a @ transformed(raw,counts.sum(axis=1).to_numpy()>0,positions)
    t,_,_,_ = analyse(yy,groups,pairs,weights,estimable,chain_ids)
    t['scenario']='all_detected_reference';sensitivity.append(t)
    for study in studies:
        keep = groups.study_id.ne(study).to_numpy()
        t,_,_,_ = analyse(y[keep],groups.loc[keep],pairs,weights,estimable,chain_ids)
        t['scenario']='omit_'+study;sensitivity.append(t)
    t,_,_,_ = analyse(y,groups,pairs,weights,estimable,chain_ids,3)
    t['scenario']='at_least_three_coordinates';sensitivity.append(t)
    write_table(pd.concat(sensitivity,ignore_index=True),out/'observed_sensitivity.tsv')
    # Fixed parameters, derived before inspecting observed conditional P values.
    mean_by_study = pd.DataFrame(y,index=groups.index).groupby(groups.study_id).mean()
    means = mean_by_study.loc[groups.study_id].to_numpy()
    residual = y-means
    sd = np.maximum(np.sqrt(np.mean(residual*residual,axis=0)),.05)
    kernels = spatial_cholesky(groups,plan['generative_parameters']['spatial_kernel_fraction'])
    target_kos = set(membership.loc[membership.chain=='N_denitrification_Nap','available_KOs'].str.split(';').explode())
    target = [features.index(k) for k in sorted(target_kos)]
    target_raw = counts.index.get_indexer(sorted(target_kos))
    depths = raw.sum(axis=1)
    rates = a @ (raw/depths[:,None])
    rates_by_study = pd.DataFrame(rates,index=groups.index).groupby(groups.study_id).mean()
    global_rate = rates_by_study.mean(axis=0).to_numpy()
    cfg=plan['generative_parameters']
    base_rates = cfg['study_profile_weight']*rates_by_study.to_numpy()+cfg['global_profile_weight']*global_rate
    run_study = rates_by_study.index.get_indexer(meta.study_id)
    run_group = groups.index.get_indexer(meta.analysis_group)
    bio_labels = meta.analysis_group+'::'+meta.biosample_group
    _,run_bio = np.unique(bio_labels,return_inverse=True)
    edge_left = np.array([features.index(k) for k in edges.KO1]);edge_right = np.array([features.index(k) for k in edges.KO2])
    write_json(out/'design_and_resolution.json',dict(profiles=len(meta),study_labels=meta.study_id.nunique(),
        coordinate_groups=meta.coordinate_group.nunique(),study_coordinate_groups=len(groups),
        informative_labels=len(studies),informative_coordinate_groups=int(groups.study_id.isin(studies).sum()),
        independent_blocks_assumed_not_verified=studies,exact_configurations=2**len(studies),
        minimum_one_sided_p=1/(2**len(studies)),historical_candidates=len(candidates),historical_candidate_pairs=candidate_pairs,
        first_Holm_threshold_for_candidate_pairs=plan['alpha']/candidate_pairs,
        chain_family=len(chains),first_Holm_threshold_for_chains=plan['alpha']/len(chains),
        frozen_chain_membership=True,new_edges_selected=False))
    checkpoint=out/'simulation_replicates.tsv.gz'
    rows=pd.read_csv(checkpoint,sep='\t').to_dict('records') if checkpoint.exists() else []
    complete={(r['scenario'],int(r['replicate'])) for r in rows}
    start=time.monotonic()
    def one_simulation(job):
        scenario_number, rep = job
        scenario = plan['scenarios'][scenario_number]
        name = scenario['id']
        rng=np.random.default_rng(np.random.SeedSequence([plan['seed'],scenario_number,rep]))
        if name.startswith('counts_'):
            sim,zero_fraction,yield_median=simulate_counts(rng,name,base_rates,run_study,run_group,run_bio,
                depths,len(groups),len(set(run_bio)),target_raw,reference,positions,a,cfg)
        else:
            sim=simulate_gaussian(rng,name,means,sd,groups,kernels,target,cfg)
            zero_fraction,yield_median=np.nan,np.nan
        tt,_,_,_=analyse(sim,groups,pairs,weights,estimable,chain_ids)
        reject=tt.p_holm_assumption_conditional.to_numpy()<=plan['alpha']
        global_rho,_,_=pilot.pair_rho(sim,edge_left,edge_right,groups.study_id)
        within_rho,_,_=pilot.pair_rho(sim,edge_left,edge_right,groups.study_id,True)
        row=dict(scenario=name,kind=scenario['kind'],replicate=rep,any_chain_rejected=bool(reject.any()),
            rejected_chains=int(reject.sum()),target_rejected=bool(reject[chain_ids.index('N_denitrification_Nap')]),
            global_fixed_edges_above_095=int((global_rho>=.95).sum()),
            within_fixed_edges_above_095=int((within_rho>=.95).sum()),
            zero_fraction=zero_fraction,median_annotation_yield=yield_median)
        for j,cid in enumerate(chain_ids):
            row['reject__'+cid]=bool(reject[j])
            row['score__'+cid]=float(tt.score.iloc[j])
        return row
    jobs=[(i,rep) for i,scenario in enumerate(plan['scenarios']) for rep in range(reps)
          if (scenario['id'],rep) not in complete]
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for row in executor.map(one_simulation,jobs):
            rows.append(row)
            rep=int(row['replicate'])
            if (rep+1)%50==0 or rep==reps-1:
                write_table(pd.DataFrame(rows),checkpoint)
                print(f"{row['scenario']}: {rep+1}/{reps}; elapsed {time.monotonic()-start:.1f}s",flush=True)
    simulations=pd.DataFrame(rows)
    summary=[]
    for scenario in plan['scenarios']:
        d=simulations[simulations.scenario==scenario['id']]
        metrics={'family_rejection':'any_chain_rejected','target_detection':'target_rejected'}
        metrics.update({'chain:'+cid:'reject__'+cid for cid in chain_ids})
        for metric,col in metrics.items():
            k=int(d[col].sum());n=len(d);ci=binomtest(k,n).proportion_ci(method='wilson')
            summary.append(dict(scenario=scenario['id'],kind=scenario['kind'],metric=metric,rejections=k,
                repetitions=n,rate=k/n,lower_95=ci.low,upper_95=ci.high))
    summary=pd.DataFrame(summary);write_table(summary,out/'calibration_rates.tsv')
    null_family=summary[(summary.kind=='statistical_null')&(summary.metric=='family_rejection')]
    stress_family=summary[(summary.kind=='measurement_stress')&(summary.metric=='family_rejection')]
    completed=reps>=plan['simulation_replicates']
    statistical_pass=completed and bool((null_family.upper_95<=plan['simulation_acceptance_upper_95_limit']).all())
    stress_pass=completed and bool((stress_family.upper_95<=plan['simulation_acceptance_upper_95_limit']).all())
    status='CONDITIONAL_METHOD_ONLY' if statistical_pass and stress_pass else 'NOT_RELEASED_FOR_MANUSCRIPT_INFERENCE'
    observed['inference_status']=status
    observed['interpretation']='Conditional diagnostics only; provenance independence and real recovery model remain unverified'
    write_table(observed,out/'observed_chain_diagnostics.tsv')
    report=dict(status=status,simulation_campaign_complete=completed,statistical_null_gate=statistical_pass,
        measurement_stress_gate=stress_pass,simulated_datasets=len(simulations),replicates_per_scenario=reps,
        family_size=len(chains),alpha=plan['alpha'],diagnostic_upper_limit=plan['simulation_acceptance_upper_95_limit'],
        scope='Eight frozen chains, exact provenance-block sign flips. No full discovery-network topology/FDR calibration.',
        claims_released=False,source_files_modified_inputs=False,
        runtime=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
    write_json(out/'calibration_summary.json',report)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
