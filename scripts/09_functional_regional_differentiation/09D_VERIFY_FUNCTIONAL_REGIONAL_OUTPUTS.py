#!/usr/bin/env python3
"""Verify cohorts, paired sampling, complete diagnostics and numerical summaries."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tables', type=Path, required=True)
    p.add_argument('--accessions', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--counts', type=Path, required=True)
    a = p.parse_args()
    def read(name):
        return pd.read_csv(a.tables/name, sep='\t')
    allowed = a.accessions.read_text().split()
    assert len(allowed) == len(set(allowed)) == 191
    meta = pd.read_csv(a.metadata).set_index('run_accession').loc[allowed]
    assert meta.index.is_unique and meta[['REALM','core_locality','study_id']].notna().all().all()
    assert sum(meta.Country == 'Guadeloupe') == 23 and 'TX1_3' in meta.index
    counts = pd.read_csv(a.counts, sep='\t', index_col=0)
    excluded = {'SRR31095011','SRR31095534','SRR31095535'}
    assert set(counts.columns)-set(allowed) == excluded
    counts = counts.loc[:, allowed]
    counts = counts.loc[counts.sum(axis=1)>0]
    assert counts.shape == (17482,191)
    master = read('180_KO_master_classification.tsv').set_index('KO').loc[counts.index]
    prev_n = (counts>0).sum(axis=1)
    assert (master.global_prevalence_n==prev_n).all()
    labels = np.where(prev_n/191 >= .95, 'Prevalent', np.where(prev_n/191 < .15,'Cloud','Shell'))
    assert (master.prevalence_class==labels).all()
    assert master.prevalence_class.value_counts().to_dict()=={'Shell':9886,'Cloud':6468,'Prevalent':1128}
    assert read('180_atlas_classification_reconciliation.tsv')['pass'].all()
    assert len(read('180_balanced_diagnostic_errors.tsv'))==0
    run = read('180_run_summary.tsv').set_index('item')['value']
    assert int(run['balanced_iterations'])==499 and int(run['balanced_permutations'])==499
    prev=read('180_panel_c_prevalence_balanced_R2.tsv')
    inter=read('180_panel_b_intersection_balanced_R2.tsv')
    realms={'Central Indo-Pacific','Tropical Atlantic','Western Indo-Pacific'}
    for frame in [prev,inter]:
        for (design,it),z in frame.groupby(['design','iteration']):
            nlocal=6 if '6' in design else 4
            assert z.selected_run_accessions.nunique()==1
            ids=z.selected_run_accessions.iloc[0].split(';')
            m=meta.loc[ids]
            assert len(ids)==len(set(ids))==3*nlocal
            assert set(m.REALM)==realms
            assert (m.groupby('REALM').size()==nlocal).all()
            assert m.groupby(['REALM','core_locality']).size().max()==1
        assert frame.diagnostic_iteration.all()
        assert frame[['R2','F','adjusted_unique_R2','adjusted_partial_R2','adjusted_P','permanova_P','permdisp_P']].notna().all().all()
        assert frame.diagnostic_error.isna().all()
        for c in ['adjusted_P','permanova_P','permdisp_P']:
            assert frame[c].between(.002,1).all()
            assert np.allclose(frame[c]*500,np.round(frame[c]*500))
        assert (frame.groupby(['design','group']).size()==499).all()
    psel=prev[(prev.design=='3 x 6')&(prev.group=='Prevalent')].set_index('iteration')
    isel=inter[inter.group=='Prevalent'].set_index('iteration')
    assert (psel.selected_run_accessions==isel.selected_run_accessions).all()
    assert np.allclose(psel.R2,isel.R2)
    assert np.allclose(psel.adjusted_unique_R2,isel.adjusted_unique_R2)
    rows=[]
    for label,frame in [('prevalence',prev),('intersection',inter)]:
        for (design,group),z in frame.groupby(['design','group']):
            item={'comparison':label,'design':design,'group':group,'n_draws':len(z),'n_KOs':int(z.n_KOs.iloc[0])}
            for field in ['R2','adjusted_unique_R2','adjusted_partial_R2']:
                item.update({field+'_median':float(z[field].median()),field+'_q05':float(z[field].quantile(.05)),field+'_q95':float(z[field].quantile(.95))})
            for field in ['permanova_P','adjusted_P','permdisp_P']:
                item[field+'_nominal_count']=int((z[field]<=.05).sum())
                item[field+'_nominal_percent']=float(100*(z[field]<=.05).mean())
            rows.append(item)
    pd.DataFrame(rows).to_csv(a.tables/'180_verified_main_statistics.tsv',sep='\t',index=False)
    contrasts=[]
    for design,z in prev.groupby('design'):
        for field in ['R2','adjusted_unique_R2']:
            w=z.pivot(index='iteration',columns='group',values=field)
            d=w.Shell-w[['Prevalent','Cloud']].max(axis=1)
            contrasts.append({'design':design,'metric':field,'median':float(d.median()),'q05':float(d.quantile(.05)),'q95':float(d.quantile(.95)),'positive_count':int((d>0).sum()),'n':len(d)})
    pd.DataFrame(contrasts).to_csv(a.tables/'180_verified_paired_contrasts.tsv',sep='\t',index=False)
    raw=read('180_panel_d_robustness_raw.tsv')
    intervals=[]
    for analysis,z in raw.groupby('analysis'):
        v=z.delta_R2.dropna();deletion=analysis.startswith('Leave-one-')
        intervals.append({'analysis':analysis,'median':float(v.median()),'lower':float(v.min() if deletion else v.quantile(.05)),
                          'upper':float(v.max() if deletion else v.quantile(.95)),'n':len(v),'positive_count':int((v>0).sum()),
                          'interval':'range across deletion medians' if deletion else '5th–95th percentiles across paired draws'})
    pd.DataFrame(intervals).to_csv(a.tables/'180_verified_sensitivity_intervals.tsv',sep='\t',index=False)
    persistence=master.groupby('persistence_class').agg(n_KOs=('prevalence_class','size'),assignments=('total_abundance','sum')).reset_index()
    persistence['richness_percent']=100*persistence.n_KOs/len(master)
    persistence['assignment_percent']=100*persistence.assignments/persistence.assignments.sum()
    persistence.to_csv(a.tables/'180_verified_persistence_totals.tsv',sep='\t',index=False)
    result={'status':'PASS','n_profiles':191,'n_nonzero_KOs':17482,'n_eligible_profiles':int(meta.REALM.isin(realms).sum()),
            'all_main_draws_one_profile_per_named_locality':True,'main_draws_paired_across_KO_groups':True,
            'intersection_and_prevalence_primary_selections_match':True,'n_complete_main_diagnostics':len(prev)+len(inter),
            'permutations_per_diagnostic':499,'all_44_study_ids_nested_in_realm':bool((meta.groupby('study_id').REALM.nunique()==1).all()),
            'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.accessions,a.metadata,a.counts]},
            'main_statistics':rows,'paired_contrasts':contrasts,'sensitivity_intervals':intervals,'persistence':persistence.to_dict(orient='records')}
    (a.tables/'180_independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['main_statistics','paired_contrasts','sensitivity_intervals','persistence']}))


if __name__=='__main__':main()
