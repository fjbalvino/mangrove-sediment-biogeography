#!/usr/bin/env python3
"""Reinfer the fixed-candidate network after removing each metadata study unit.

500 locality-capped samples per omission. The CLR reference and candidate IDs
remain fixed, so this is conditional sensitivity, not independent validation.
The numerical engine must reproduce every edge from the primary R run first.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import argparse, concurrent.futures, hashlib, json, subprocess, time
from pathlib import Path
import numpy as np
import pandas as pd


def infer_case(job):
    root,case=job; root=Path(root); started=time.monotonic()
    X=np.load(root/'inputs/CLR_candidates.npy',mmap_mode='r')
    names=pd.read_csv(root/'inputs/candidates.tsv',sep='\t').KO.to_numpy()
    sample_ids=json.loads((root/'inputs/sample_ids.json').read_text())
    lookup={v:i for i,v in enumerate(sample_ids)}
    design=pd.read_csv(root/'designs'/f'{case}.tsv.gz',sep='\t')
    jj,ii=np.tril_indices(len(names),k=-1) # R upper.tri column-major order
    n_pairs=len(ii); counts=np.zeros(n_pairs,dtype=np.int32)
    sums=np.zeros(n_pairs); positive=np.zeros(n_pairs,dtype=np.int32)
    selected_runs=[]
    for b,d in design.groupby('bootstrap_id',sort=True):
        rows=np.array([lookup[v] for v in d.run_accession])
        z=np.array(X[rows],copy=True); z-=z.mean(axis=0)
        cov=z.T@z/(len(rows)-1); variance=np.diag(cov)
        den=variance[:,None]+variance[None,:]
        with np.errstate(divide='ignore',invalid='ignore'):
            values=(2*cov[ii,jj])/den[ii,jj]
        absolute=np.abs(values); valid=np.flatnonzero(np.isfinite(absolute))
        k=min(2000,len(valid)); cutoff=np.partition(absolute[valid],len(valid)-k)[len(valid)-k]
        eligible=valid[absolute[valid]>=cutoff]
        ix=eligible[np.lexsort((eligible,-absolute[eligible]))[:k]]
        counts[ix]+=1; sums[ix]+=absolute[ix]; positive[ix]+=(values[ix]>0)
        if case=='baseline': selected_runs.append((int(b),np.sort(ix)))
    seen=np.flatnonzero(counts)
    out=pd.DataFrame({'KO1':names[ii[seen]],'KO2':names[jj[seen]],
      'n_boot_present':counts[seen], 'edge_stability':counts[seen]/design.bootstrap_id.nunique(),
      'mean_abs_weight':sums[seen]/counts[seen],
      'sign_consistency':np.maximum(positive[seen],counts[seen]-positive[seen])/counts[seen]})
    out=out.sort_values(['edge_stability','sign_consistency','mean_abs_weight','KO1','KO2'],ascending=[False,False,False,True,True])
    edges=out.query('edge_stability>=0.8 and sign_consistency>=1 and mean_abs_weight>=0.95')
    dest=root/'networks';dest.mkdir(exist_ok=True)
    out.to_csv(dest/f'{case}_edge_pool.tsv.gz',sep='\t',index=False)
    edges.to_csv(dest/f'{case}_edges.tsv',sep='\t',index=False)
    if case=='baseline':
        np.savez_compressed(dest/'baseline_edge_indices.npz',**{str(b):v for b,v in selected_runs})
    result={'case':case,'n_subsamples':int(design.bootstrap_id.nunique()),
      'n_subsample':int(design.groupby('bootstrap_id').size().iloc[0]),
      'n_edges':len(edges),'n_connected':len(set(edges.KO1)|set(edges.KO2)),
      'seconds':time.monotonic()-started}
    (dest/f'{case}_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run156',type=Path,required=True);p.add_argument('--run171',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--rscript',default='Rscript');p.add_argument('--workers',type=int,default=4)
    p.add_argument('--design-only',action='store_true');p.add_argument('--resume',action='store_true')
    a=p.parse_args();root=a.outdir.resolve();(root/'inputs').mkdir(parents=True,exist_ok=True)
    t171=a.run171/'tables';t156=a.run156/'tables'
    candidates=pd.read_csv(t171/'171_candidate_KO_selection.tsv',sep='\t')
    clr=pd.read_csv(t156/'156_KO_CLR_network_filtered.tsv.gz',sep='\t').set_index('KO')
    meta=pd.read_csv(a.metadata,sep='\t'); meta=meta.set_index('run_accession').loc[clr.columns].rename_axis('run_accession').reset_index()
    assert len(meta)==191 and candidates.KO.nunique()==1500
    X=clr.loc[candidates.KO].to_numpy(dtype=np.float64).T
    np.save(root/'inputs/CLR_candidates.npy',X)
    candidates[['KO']].to_csv(root/'inputs/candidates.tsv',sep='\t',index=False)
    meta.to_csv(root/'inputs/metadata.tsv',sep='\t',index=False)
    (root/'inputs/sample_ids.json').write_text(json.dumps(clr.columns.to_list()))
    helper=Path(__file__).parent/'lib/study_sensitivity.R'
    subprocess.run([a.rscript,str(helper),'design',str(root)],check=True)
    manifest=pd.read_csv(root/'study_units.tsv',sep='\t')
    if a.design_only:return
    baseline=infer_case((str(root),'baseline'))
    ref=pd.read_csv(t171/'171_bootstrap_edge_membership.tsv.gz',sep='\t')
    design=pd.read_csv(root/'designs/baseline.tsv.gz',sep='\t')
    original_design=pd.read_csv(t171/'171_bootstrap_design.tsv.gz',sep='\t')
    assert design[['bootstrap_id','run_accession']].equals(original_design[['bootstrap_id','run_accession']])
    idx={v:i for i,v in enumerate(candidates.KO)}
    comparisons=np.load(root/'networks/baseline_edge_indices.npz')
    for b,d in ref.groupby('bootstrap_id'):
        first=d.KO1.map(idx).to_numpy();second=d.KO2.map(idx).to_numpy()
        lo=np.minimum(first,second);hi=np.maximum(first,second)
        assert np.array_equal(np.sort(hi*(hi-1)//2+lo),comparisons[str(b)]),f'Edge replay failed in draw {b}'
    original=pd.read_csv(t171/'171_consensus_edges.tsv.gz',sep='\t').set_index(['KO1','KO2']).sort_index()
    replay=pd.read_csv(root/'networks/baseline_edges.tsv',sep='\t').set_index(['KO1','KO2']).sort_index()
    assert original.index.equals(replay.index)
    delta=float(np.max(np.abs(original.mean_abs_weight-replay.mean_abs_weight)))
    assert delta<1e-12
    validation={'status':'PASS','n_exact_edge_sets':500,'edges_checked':len(ref),
      'same_sampling_design':True,'max_weight_difference':delta,'reference_run171':str(a.run171.resolve())}
    (root/'numerical_engine_verification.json').write_text(json.dumps(validation,indent=2)+'\n')
    print('NUMERICAL ENGINE VERIFIED',json.dumps(validation),flush=True)
    jobs=[(str(root),v) for v in manifest.case if not (a.resume and (root/'networks'/f'{v}_summary.json').exists())]
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as ex:
        for result in ex.map(infer_case,jobs):
            print('STUDY',json.dumps(result),flush=True)
    summary=pd.DataFrame([json.loads((root/'networks'/f'{v}_summary.json').read_text()) for v in manifest.case])
    summary=manifest.merge(summary,on='case',validate='one_to_one')
    summary.to_csv(root/'study_network_summary.tsv',sep='\t',index=False)
    print('ALL STUDY NETWORKS COMPLETE',len(summary),flush=True)

if __name__=='__main__':main()
