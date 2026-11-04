"""Summarize prespecified calibration without selecting on real gene effects."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import norm

root=Path(__file__).resolve().parents[1];out=root/'calibration_final'
plan=json.loads((out/'simulation_plan.json').read_text())
raw=pd.read_csv(out/'chunk_sufficient_statistics.tsv.gz',sep='\t',keep_default_na=False,na_values=['NA',''])
keys=['scenario','signal','predictor','method','component','gene']
m=raw.groupby(keys,sort=True).sum(numeric_only=True).reset_index()
assert len(m)==len(plan['family_genes'])*120
assert set(m.loc[m.signal=='null','replicates'])=={1000}
assert set(m.loc[m.signal=='mixed','replicates'])=={500}
m['truth']=m.truth/m.replicates
assert np.allclose(m.loc[m.signal=='null','truth'],0)
assert np.isfinite(m.select_dtypes('number')).all().all()
def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-half,mid+half

m['valid_fraction']=m.valid/m.replicates
m['mean_estimate']=m.b/m.valid
m['bias']=m.mean_estimate-m.truth
m['empirical_sd']=np.sqrt(np.maximum(0,(m.b2-m.b*m.b/m.valid)/(m.valid-1)))
m['bias_mcse']=m.empirical_sd/np.sqrt(m.valid)
m['mean_se']=m.se/m.valid;m['mean_df']=m.df/m.valid
m['mean_interval_width']=m.width/m.valid
m['rejection_rate']=m.reject/m.valid
m['rejection_mcse']=np.sqrt(m.rejection_rate*(1-m.rejection_rate)/m.valid)
m['rejection_low'],m['rejection_high']=wilson(m.reject,m.valid)
m['coverage_rate']=m.coverage/m.valid
m['coverage_low'],m['coverage_high']=wilson(m.coverage,m.valid)
m['cell_pass']=(m.valid_fraction>=.95)&(m.rejection_high<=.075)&(m.coverage_low>=.925)
m.loc[m.signal!='null','cell_pass']=False # Only the prespecified all-null cells define this gate.
m.to_csv(out/'final_performance_metrics.tsv',sep='\t',index=False)

f=pd.read_csv(out/'replicate_FDR.tsv',sep='\t',keep_default_na=False,na_values=['NA',''])
assert len(f)==15000 and not f.duplicated(['scenario','signal','replicate','method']).any()
assert ((f.false<=f.discoveries)&(f.true+f.false==f.discoveries)&(f.valid<=plan['family_tests'])).all()
assert np.allclose(f.FDP,f.false/np.maximum(1,f.discoveries))
assert set(f.loc[f.signal=='mixed','n_nonnull'])=={len(plan['signal_genes'])*6}
fr=[]
for (scenario,signal,method),a in f.groupby(['scenario','signal','method']):
    n=len(a);estimate=a.FDP.mean();se=a.FDP.std(ddof=1)/np.sqrt(n)
    if signal=='null' or (a.FDP==0).all():
        lo,hi=wilson((a.false>0).sum(),n)
        ci='Wilson for any false discovery (equals FDR under global null; conservative upper bound when mixed FDP is uniformly zero)'
    else:
        lo=max(0,estimate-1.96*se);hi=min(1,estimate+1.96*se)
        ci='Normal Monte Carlo interval for mean false-discovery proportion'
    fr.append(dict(scenario=scenario,signal=signal,method=method,replicates=n,FDR=estimate,FDR_mcse=se,FDR_low=lo,FDR_high=hi,
                   interval_method=ci,power=a.power.mean(),power_mcse=a.power.std(ddof=1)/np.sqrt(n),mean_valid_fraction=a.valid.mean()/plan['family_tests'],
                   mean_discoveries=a.discoveries.mean(),pass_FDR=hi<=.075))
fm=pd.DataFrame(fr);fm.to_csv(out/'final_FDR_metrics.tsv',sep='\t',index=False)

rows=[]
for (method,component),a in m[m.signal=='null'].groupby(['method','component']):
    rows.append(dict(method=method,component=component,cells=len(a),valid_min=a.valid_fraction.min(),valid_max=a.valid_fraction.max(),
         typeI_min=a.rejection_rate.min(),typeI_max=a.rejection_rate.max(),coverage_min=a.coverage_rate.min(),coverage_max=a.coverage_rate.max(),
         cells_pass=int(a.cell_pass.sum()),all_cells_pass=bool(a.cell_pass.all())))
table=pd.DataFrame(rows);table.to_csv(out/'Table_S41_source.tsv',sep='\t',index=False)
decisions={}
for method in ['conventional','CR2']:
    a=m[(m.method==method)&(m.signal=='null')];b=fm[fm.method==method]
    decisions[method]={'accepted':bool(a.cell_pass.all() and b.pass_FDR.all()),'null_cells':len(a),'null_cells_passing':int(a.cell_pass.sum()),
                       'minimum_valid_fraction':float(a.valid_fraction.min()),'maximum_null_typeI':float(a.rejection_rate.max()),
                       'minimum_null_coverage':float(a.coverage_rate.min()),'maximum_FDR':float(b.FDR.max()),
                       'mixed_BH_power_range':[float(b[b.signal=='mixed'].power.min()),float(b[b.signal=='mixed'].power.max())]}
assert not m[m.signal=='null'].empty
inputs_ok={p:hashlib.sha256((out/p).read_bytes()).hexdigest()==h for p,h in plan['input_hashes'].items()}
assert all(inputs_ok.values())
summary={'status':'Final calibration completed','synthetic_joint_datasets':7500,'genes':len(plan['family_genes']),'family_tests':plan['family_tests'],
 'null_replicates_per_scenario':1000,'mixed_replicates_per_scenario':500,'scenarios':5,
 'performance_rows':len(m),'FDR_replicate_rows':len(f),'seed':plan['seed'],'decisions':decisions,
 'input_hashes_unchanged':inputs_ok,'plan_sha256':hashlib.sha256((out/'simulation_plan.json').read_bytes()).hexdigest(),
 'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/'scripts').glob('0[2345]*'))},
 'interpretation':'No automatic certification from low rejection: validity, coverage, type-I error and multiplicity gates assessed jointly. Calibration is conditional on working DGPs, not a claim of biological independence or full compositional robustness.'}
(out/'final_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# These are descriptive views of the fixed global family, not new BH families.
inventory=pd.read_csv(root/'results/KO_inventory.tsv',sep='\t')
cycles=json.loads((root/'curation/cycles.json').read_text())
cycle_rows=[]
for cycle in cycles:
    genes=inventory.loc[inventory.cycles.fillna('').str.split(';').map(lambda a:cycle in a)&inventory.row_available,'KO']
    for method in ['conventional','CR2']:
        for component in ['detection','positive']:
            a=m[(m.signal=='null')&(m.method==method)&(m.component==component)&m.gene.isin(genes)]
            assert len(a)==len(genes)*15 and len(a)>0
            cycle_rows.append(dict(cycle=cycle,method=method,component=component,unique_available_KOs=len(genes),
                null_cells=len(a),null_cells_passing=int(a.cell_pass.sum()),valid_fraction_min=a.valid_fraction.min(),
                valid_fraction_max=a.valid_fraction.max(),typeI_min=a.rejection_rate.min(),typeI_max=a.rejection_rate.max(),
                coverage_min=a.coverage_rate.min(),coverage_max=a.coverage_rate.max(),all_null_cells_pass=bool(a.cell_pass.all())))
pd.DataFrame(cycle_rows).to_csv(out/'per_cycle_null_calibration.tsv',sep='\t',index=False)
print(json.dumps(summary,indent=2))
