#!/usr/bin/env python3
"""Stress-test the existing complete consensus-edge selector with known generators."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return m
cal=load_module('selector_calibration',HERE/'10M_CALIBRATE_PROCESS_COHERENCE.py')
engine=load_module('selector_engine',HERE/'10J_NETWORK_STUDY_SENSITIVITY.py')


class RecordingGenerator:
    """Retain the Poisson draw while reusing the unchanged calibration generator."""
    def __init__(self,seed):self.rng=np.random.default_rng(seed);self.counts=None
    def normal(self,*a,**k):return self.rng.normal(*a,**k)
    def random(self,*a,**k):return self.rng.random(*a,**k)
    def poisson(self,*a,**k):
        self.counts=self.rng.poisson(*a,**k)
        return self.counts


def run_case(job):
    return engine.infer_case(job)


def prepare_case(root,name,clr,candidates,sample_ids,design):
    case_root=root/'cases'/name
    (case_root/'inputs').mkdir(parents=True,exist_ok=True)
    (case_root/'designs').mkdir(exist_ok=True)
    np.save(case_root/'inputs/CLR_candidates.npy',clr)
    pd.DataFrame({'KO':candidates}).to_csv(case_root/'inputs/candidates.tsv',sep='\t',index=False)
    (case_root/'inputs/sample_ids.json').write_text(json.dumps(sample_ids))
    shutil.copy2(design,case_root/'designs'/f'{name}.tsv.gz')
    return str(case_root),name


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('matrix','metadata','samples','reference-edges','reference-candidates','reference-design','reference-membership','outdir'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args();root=args.outdir.resolve()
    if root.exists() and any(root.iterdir()):raise ValueError('Use an empty stress-test directory')
    root.mkdir(parents=True,exist_ok=True)
    plan=json.loads((ROOT/'config/network_selector_stress.json').read_text())
    cfg=json.loads((ROOT/'config/network_calibration.json').read_text())['generative_parameters']
    cal.write_json(root/'plan.json',plan)
    fingerprints={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in vars(args).items() if isinstance(v,Path) and v.is_file()}
    for source in (Path(__file__),HERE/'10J_NETWORK_STUDY_SENSITIVITY.py',HERE/'10M_CALIBRATE_PROCESS_COHERENCE.py',ROOT/'config/network_calibration.json',ROOT/'config/network_selector_stress.json'):
        fingerprints[str(source.relative_to(ROOT))]=hashlib.sha256(source.read_bytes()).hexdigest()
    cal.write_json(root/'input_and_source_hashes.json',fingerprints)
    counts=cal.context.read_matrix(args.matrix,args.samples)
    reference,candidates,y=cal.pilot.candidate_selection(counts,plan['gene_candidate_cap'])
    original_candidates=pd.read_csv(args.reference_candidates,sep='\t')
    if candidates!=original_candidates.KO.tolist():raise ValueError('Historical candidates not reproduced')
    np.testing.assert_allclose(y.loc[candidates].var(axis=1),original_candidates.recalculated_CLR_variance,rtol=1e-10,atol=1e-12)
    baseline=prepare_case(root,'baseline',y.loc[candidates].T.to_numpy(),candidates,counts.columns.tolist(),args.reference_design)
    print('Replaying the full historical selector before generating controls',flush=True)
    baseline_result=run_case(baseline)
    baseline_root=Path(baseline[0]);idx={k:i for i,k in enumerate(candidates)}
    membership=pd.read_csv(args.reference_membership,sep='\t')
    draws=np.load(baseline_root/'networks/baseline_edge_indices.npz')
    for b,d in membership.groupby('bootstrap_id'):
        first=d.KO1.map(idx).to_numpy();second=d.KO2.map(idx).to_numpy()
        lo=np.minimum(first,second);hi=np.maximum(first,second)
        np.testing.assert_array_equal(np.sort(hi*(hi-1)//2+lo),draws[str(b)])
    old=pd.read_csv(args.reference_edges,sep='\t').set_index(['KO1','KO2']).sort_index()
    new=pd.read_csv(baseline_root/'networks/baseline_edges.tsv',sep='\t').set_index(['KO1','KO2']).sort_index()
    assert old.index.equals(new.index)
    delta=float(np.max(np.abs(old.mean_abs_weight-new.mean_abs_weight)))
    assert delta<1e-12
    cal.write_json(root/'baseline_verification.json',dict(status='PASS',edge_selections=len(membership),
        subsamples=len(draws.files),consensus_edges=len(old),maximum_weight_difference=delta,engine_reused='10J.infer_case'))
    print('Baseline verified:',baseline_result,flush=True)
    meta=pd.read_csv(args.metadata)
    if meta.run_accession.isna().any() or meta.run_accession.duplicated().any():raise ValueError('Invalid metadata identifiers')
    meta=meta.set_index('run_accession').loc[counts.columns].copy()
    meta['coordinate_group']=[f'coord_{la:.5f}_{lo:.5f}' for la,lo in zip(meta.lat,meta.lon)]
    meta['analysis_group']=meta.study_id+'::'+meta.coordinate_group
    meta['biosample_group']=[str(b) if pd.notna(b) and str(b).strip() else 'run:'+r for r,b in zip(meta.index,meta.sample_accession)]
    a,groups=cal.aggregation(meta);raw=counts.T.to_numpy(float);depths=raw.sum(axis=1)
    group_rates=a@(raw/depths[:,None])
    study_rates=pd.DataFrame(group_rates,index=groups.index).groupby(groups.study_id).mean()
    base_rates=np.tile(study_rates.mean(axis=0).to_numpy(),(len(study_rates),1))
    run_study=study_rates.index.get_indexer(meta.study_id);run_group=groups.index.get_indexer(meta.analysis_group)
    _,run_bio=np.unique(meta.analysis_group+'::'+meta.biosample_group,return_inverse=True)
    mask=counts.index.isin(reference);jobs=[];manifest=[]
    for condition_number,condition in enumerate(plan['conditions']):
        for rep in range(plan['replicates_per_condition']):
            name=f'{condition}_{rep:02d}'
            seed=[plan['seed'],100+condition_number,rep]
            rng=RecordingGenerator(np.random.SeedSequence(seed))
            cal.simulate_counts(rng,condition,base_rates,run_study,run_group,run_bio,depths,len(groups),
                len(set(run_bio)),[],mask,np.array([0]),a,cfg)
            simulated=pd.DataFrame(rng.counts.T,index=counts.index,columns=counts.columns)
            ref,selected,clr=cal.pilot.candidate_selection(simulated,plan['gene_candidate_cap'])
            jobs.append(prepare_case(root,name,clr.loc[selected].T.to_numpy(),selected,counts.columns.tolist(),args.reference_design))
            manifest.append(dict(case=name,condition=condition,replicate=rep,seed=json.dumps(seed),
                reference_KOs=len(ref),candidate_KOs=len(selected),zero_fraction=float(np.mean(rng.counts==0))))
    cal.write_table(pd.DataFrame(manifest),root/'case_manifest.tsv')
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(run_case,jobs):
            results.append(result);print(json.dumps(result),flush=True)
            cal.write_table(pd.DataFrame(results),root/'selector_results.tsv')
    combined=pd.DataFrame(manifest).merge(pd.DataFrame(results),on='case',validate='one_to_one')
    cal.write_table(combined,root/'selector_stress_summary.tsv')
    cal.write_json(root/'completion.json',dict(status='COMPLETE',controls=len(results),
        selector_subsamples_per_control=baseline_result['n_subsamples'],baseline_replayed=True,
        candidate_selection_repeated=True,community_clustering_executed=False,
        interpretation=plan['interpretation']))


if __name__=='__main__':main()
