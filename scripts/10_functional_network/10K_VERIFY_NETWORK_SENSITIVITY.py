#!/usr/bin/env python3
"""Verify every study omission and independently reproduce influential-case edges."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True);p.add_argument('--run156',type=Path,required=True)
    p.add_argument('--rscript',default='Rscript');p.add_argument('--reuse-r-check',action='store_true')
    a=p.parse_args();root=a.study.resolve();scripts=Path(__file__).resolve().parent
    meta=pd.read_csv(root/'inputs/metadata.tsv',sep='\t');manifest=pd.read_csv(root/'study_units.tsv',sep='\t');checks=[]
    for u in manifest.itertuples():
        d=pd.read_csv(root/'designs'/f'{u.case}.tsv.gz',sep='\t')
        excluded=set(meta.loc[meta.study_id.eq(u.study_id),'run_accession'])
        assert not excluded & set(d.run_accession)
        assert d.bootstrap_id.nunique()==500
        assert d.groupby(['bootstrap_id','BALANCE_STRATUM']).size().le(4).all()
        assert set(d.run_accession)<=set(meta.run_accession)-excluded
        checks.append({'study_id':u.study_id,'case':u.case,'n_removed':len(excluded),
          'n_subsamples':500,'no_omitted_accessions_present':True,'stratum_cap_valid':True})
    (root/'study_design_verification.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n')
    results=pd.read_csv(root/'study_exclusion_results.tsv',sep='\t')
    u=results.query("unit_type=='BioProject'").nsmallest(1,'ARI_reference_major').iloc[0]
    if not a.reuse_r_check:
        with (root/'influential_omission_R_check.log').open('w') as f:
            subprocess.run([a.rscript,'--vanilla',str(scripts/'lib/verify_study_edges.R'),
              str(scripts/'10C_FUNCTIONAL_KO_CONSENSUS_NETWORK.R'),str(a.run156.resolve()),str(root),u.case],
              stdout=f,stderr=subprocess.STDOUT,check=True)
    r=pd.read_csv(root/'influential_omission_R_check.tsv.gz',sep='\t')
    D=pd.read_csv(root/'designs'/f'{u.case}.tsv.gz',sep='\t');Y=np.load(root/'inputs/CLR_candidates.npy')
    names=pd.read_csv(root/'inputs/candidates.tsv',sep='\t').KO.to_list();ids=json.loads((root/'inputs/sample_ids.json').read_text())
    sid={v:i for i,v in enumerate(ids)};kid={v:i for i,v in enumerate(names)};jj,ii=np.tril_indices(len(names),k=-1)
    for b,e in r.groupby('bootstrap_id'):
        rows=[sid[v] for v in D.loc[D.bootstrap_id.eq(b),'run_accession']]
        z=Y[rows].copy();z-=z.mean(axis=0);cv=z.T@z/(len(rows)-1);vv=np.diag(cv)
        with np.errstate(divide='ignore',invalid='ignore'):v=2*cv[ii,jj]/(vv[ii]+vv[jj])
        ab=np.abs(v);valid=np.flatnonzero(np.isfinite(ab));cutoff=np.partition(ab[valid],len(valid)-2000)[len(valid)-2000]
        eligible=valid[ab[valid]>=cutoff];ix=eligible[np.lexsort((eligible,-ab[eligible]))[:2000]]
        left=e.KO1.map(kid).to_numpy();right=e.KO2.map(kid).to_numpy();lo=np.minimum(left,right);hi=np.maximum(left,right)
        assert np.array_equal(np.sort(hi*(hi-1)//2+lo),np.sort(ix)),int(b)
    result={'status':'PASS','case':u.case,'BioProject':u.study_id,'n_independent_R_subsamples':int(r.bootstrap_id.nunique()),
      'n_edge_selections_compared':len(r),'all_edge_sets_identical':True,'omission_cohort_n':int(len(meta)-u.n_removed)}
    assert result['n_independent_R_subsamples']==50 and len(r)==100000
    (root/'influential_omission_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'study_designs_verified':len(checks),'influential_case':result},indent=2))

if __name__=='__main__':main()
