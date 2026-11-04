"""ID-safe extension of the frozen nitrogen design to every remaining rule."""
from pathlib import Path
import importlib.util,json,re,hashlib,shutil
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
b=module('base',ROOT/'scripts/12B_BUILD_NITROGEN_PILOT.py')
ev=module('rules',SOURCE/'evaluate.py')
meta,mat=b.load_inputs(SOURCE/'inputs/metadata_functional.csv',SOURCE/'inputs/functional_KO_counts_original.tsv.gz',SOURCE/'inputs/functional_accessions.txt',None)
rules=pd.read_csv(SOURCE/'inputs/operational_rules.tsv',sep='\t');rules=rules[~rules.process_id.str.startswith('N')].copy()
processes=pd.read_csv(SOURCE/'inputs/curated_processes.tsv',sep='\t');processes=processes[processes.cycle!='N'].copy()
ks=sorted(set(re.findall(r'K\d{5}',' '.join(rules.expression))))
inv=pd.read_csv(SOURCE/'outputs/KO_panel_inventory.tsv',sep='\t').set_index('KO').loc[ks].copy()
inv['gene']=inv.name.str.split(';').str[0].str.slice(0,70)
cycles={p.id:p.cycle for p in processes.itertuples()}
rules['cycle']=rules.process_id.map(cycles)
inv['cycles']=[';'.join(sorted(set(rules.loc[rules.expression.str.contains(k),'cycle']))) for k in ks]
inv['row_available']=[k in mat.index for k in ks]
panel=mat.reindex(ks);rel=panel.div(mat.sum(),axis=1)*1e6
det=panel.ge(1).astype(float).mask(panel.isna());det2=panel.ge(2).astype(float).mask(panel.isna())
out=ROOT/'results';out.mkdir(exist_ok=True)
for name,df in [('raw_assignments',panel),('relative_representation',rel),('detection_ge1',det),('detection_ge2',det2),('group_relative_representation',b.hierarchical_summary(rel,meta)),('group_detection_ge1',b.hierarchical_summary(det,meta)),('group_detection_ge2',b.hierarchical_summary(det2,meta))]:
    df.to_csv(out/(name+'.tsv.gz'),sep='\t',index_label='analysis_group' if name.startswith('group_') else 'KO')
inv['n_ge1']=det.sum(axis=1,min_count=1);inv['n_ge2']=det2.sum(axis=1,min_count=1)
inv.to_csv(out/'KO_inventory.tsv',sep='\t')
rules.to_csv(ROOT/'curation/operational_rules.tsv',sep='\t',index=False)
processes.to_csv(ROOT/'curation/process_scope_original.tsv',sep='\t',index=False)
ev.self_tests()
ruleframes=[]
for threshold in [1,2]:
    values={k:(mat.loc[k].to_numpy()>=threshold).astype('uint8') if k in mat.index else np.full(len(meta),2,dtype='uint8') for k in ks}
    df=pd.DataFrame({r.rule_id:ev.evaluate(ev.parse(r.expression),values,len(meta)) for r in rules.itertuples()},index=meta.index)
    df.to_csv(out/f'rule_status_ge{threshold}.tsv',sep='\t')
    ruleframes.append(df)
# Rebuild derived grouping and weights from the current ID-matched inputs.
meta['total_KO_assignments']=mat.sum()
groups=meta.groupby('analysis_group')[['study_id','coordinate_group','lat','lon','Country','Locality','REALM']].first()
for c in ['bio1','soc_mean','phh2o_mean']:
    assert (meta.groupby('analysis_group')[c].nunique()<=1).all()
    groups[c]=meta.groupby('analysis_group')[c].first()
groups['n_runs']=meta.groupby('analysis_group').size()
groups['n_biosample_labels']=meta.groupby('analysis_group').biosample_group.nunique()
groups['total_KO_assignments']=b.hierarchical_summary(pd.DataFrame([mat.sum()],index=['value']),meta)['value']
meta.to_csv(out/'sample_analysis_manifest.tsv',sep='\t')
groups.to_csv(out/'study_coordinate_groups.tsv',sep='\t')
unit=meta[['analysis_group','study_id','biosample_group']].copy()
nb=meta.groupby('analysis_group').biosample_group.transform('nunique')
nr=meta.groupby(['analysis_group','biosample_group']).biosample_group.transform('size')
unit['descriptive_run_weight_within_group']=1/(nb*nr)
assert np.allclose(unit.groupby('analysis_group').descriptive_run_weight_within_group.sum(),1)
unit.to_csv(out/'unit_design_manifest.tsv',sep='\t')
gr=pd.read_csv(out/'group_relative_representation.tsv.gz',sep='\t',index_col=0)
gd=pd.read_csv(out/'group_detection_ge1.tsv.gz',sep='\t',index_col=0)
diagnostics=b.association_diagnostics(gr,gd,groups,inv.gene.to_dict(),out)
summary=[]
for cycle,rs in rules.groupby('cycle',sort=False):
    kk=sorted(set(re.findall(r'K\d{5}',' '.join(rs.expression))))
    summary.append(dict(cycle=cycle,processes=int((processes.cycle==cycle).sum()),rules=len(rs),KOs=len(kk),available=int(inv.loc[kk,'row_available'].sum()),unavailable=int((~inv.loc[kk,'row_available']).sum()),unresolved_processes=int((~processes.loc[processes.cycle==cycle,'id'].isin(rs.process_id)).sum())))
pd.DataFrame(summary).to_csv(out/'cycle_scope.tsv',sep='\t',index=False)
report={'profiles':len(meta),'remaining_processes':len(processes),'rules':len(rules),'unique_KOs':len(ks),'available_unique_KOs':int(inv.row_available.sum()),'unavailable_unique_KOs':int((~inv.row_available).sum()),'unresolved_processes':int((~processes.id.isin(rules.process_id)).sum()),'diagnostics':diagnostics,'original_input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE/'inputs/metadata_functional.csv',SOURCE/'inputs/functional_KO_counts_original.tsv.gz']},'denominator_KO_rows':len(mat),'states':{'0':'observed rule not met','1':'required evidence co-detected within sample','2':'indeterminate because an unavailable row is needed'},'interpretation':'Community-level genetic potential; co-detection does not establish one complete genome or realized flux.'}
(out/'extension_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
