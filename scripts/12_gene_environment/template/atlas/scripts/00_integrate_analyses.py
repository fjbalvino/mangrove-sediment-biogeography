"""Integrate two input-matched analyses by KO and accession. No model fitting or new test family."""
from pathlib import Path
import json,re,importlib.util
import pandas as pd,numpy as np
r=Path(__file__).resolve().parents[2];a=r/'atlas';n=r/'analysis_nitrogen';rem=r/'analysis_remaining'
# Write matrices atomically, preserving gzip format for repeatable local rebuilds.
import uuid,shutil
for folder in ['results','curation']:(a/folder).mkdir(parents=True,exist_ok=True)
_csv_original=pd.DataFrame.to_csv
def _atomic_csv(self,path,*args,**kwargs):
    path=Path(path);temporary=path.with_name('.'+uuid.uuid4().hex+'-'+path.name)
    if path.name.endswith('.gz'):kwargs.setdefault('compression','gzip')
    _csv_original(self,temporary,*args,**kwargs);temporary.replace(path)
pd.DataFrame.to_csv=_atomic_csv
for name in ['sample_analysis_manifest.tsv','study_coordinate_groups.tsv','unit_design_manifest.tsv']:
    shutil.copy2(rem/'results'/name,a/'results'/name)
rr=pd.read_csv(r/'inputs/operational_rules.tsv',sep='\t');nreact=json.loads((n/'config/gene_environment/nitrogen_reactions.json').read_text())
scope=pd.read_csv(r/'inputs/curated_processes.tsv',sep='\t')
ps=json.loads((rem/'curation/processes_english.json').read_text());nps=[]
for p in nreact:
 p=dict(p);p['cycle']='N';p['rules']=[{'id':row.rule_id,'expression':row.expression,'kos':sorted(set(re.findall(r'K\d{5}',row.expression))),'source':row.source,'kind':row.evidence_kind,'primary':bool(row.primary)} for row in rr[rr.rule_id.isin(p['rules'])].itertuples()];nps.append(p)
ps=nps+ps;(a/'curation/processes_english.json').write_text(json.dumps(ps,ensure_ascii=False,indent=2)+'\n')
cycles=json.loads((rem/'curation/cycles.json').read_text());cycles={'N':{'name':'Nitrogen','slug':'nitrogen','symbol':'N','intro':'Follow nitrogen fixation, nitrogen-oxide transformations and ammonium-producing pathways. Shared enzyme families and unresolved anammox markers retain their limits.'},**cycles};(a/'curation/cycles.json').write_text(json.dumps(cycles,ensure_ascii=False,indent=2)+'\n')
inv=pd.read_csv(rem/'results/KO_inventory.tsv',sep='\t').set_index('KO');ni=pd.read_csv(n/'results/nitrogen_KO_inventory.tsv',sep='\t').set_index('KO');full=pd.read_csv(r/'outputs/KO_panel_inventory.tsv',sep='\t').set_index('KO')
shared_kos=inv.index.intersection(ni.index)
for ko,row in ni.iterrows():
 if ko in inv.index:inv.loc[ko,'cycles']='N;'+inv.loc[ko,'cycles']
 else:
  inv.loc[ko]=None
  for col in inv.columns:
   if col in full:inv.loc[ko,col]=full.loc[ko,col]
  inv.loc[ko,'cycles']='N';inv.loc[ko,'row_available']=row.row_available;inv.loc[ko,'n_ge1']=row.detected_ge1;inv.loc[ko,'n_ge2']=row.detected_ge2
 inv.loc[ko,'gene']=row.gene
inv['row_available']=inv['row_available'].map(lambda v:str(v).lower() in ['true','1','1.0'])
inv.to_csv(a/'results/KO_inventory.tsv',sep='\t')
for dest,ns in [('raw_assignments','nitrogen_raw_assignments.tsv'),('relative_representation','nitrogen_per_million_KO_assignments.tsv'),('detection_ge1','nitrogen_detection_ge1.tsv'),('detection_ge2','nitrogen_detection_ge2.tsv'),('group_relative_representation','group_relative_representation.tsv'),('group_detection_ge1','group_detection_ge1.tsv'),('group_detection_ge2','group_detection_ge2.tsv')]:
 x=pd.read_csv(rem/'results'/(dest+'.tsv.gz'),sep='\t',index_col=0);y=pd.read_csv(n/'results'/ns,sep='\t',index_col=0)
 if dest.startswith('group_'):x=x.T;y=y.T
 assert set(x.columns)==set(y.columns);y=y[x.columns]
 shared=x.index.intersection(y.index);assert np.allclose(x.loc[shared],y.loc[shared],equal_nan=True),dest
 z=pd.concat([y.loc[~y.index.isin(x.index)],x]);z=z.T if dest.startswith('group_') else z
 z.to_csv(a/'results'/(dest+'.tsv.gz'),sep='\t')
for filename,keys in [('environment_association_diagnostics.tsv',['KO','outcome','predictor']),('per_KO_model_support.tsv',['KO','outcome'])]:
 x=pd.read_csv(rem/'results'/filename,sep='\t');y=pd.read_csv(n/'results'/filename,sep='\t')
 # Shared KOs must describe identical fits before keeping one record. Gene
 # display labels can differ by panel; KO IDs, support and numerical results cannot.
 for label,frame in [('remaining',x),('nitrogen',y)]:
  if frame[keys].isna().any().any() or frame.duplicated(keys).any():
   raise ValueError(f'{filename}: missing or duplicate diagnostic keys in {label}')
 if set(x.columns)!=set(y.columns):
  raise ValueError(f'{filename}: diagnostic schemas differ between panels')
 compare_columns=[column for column in x.columns if column not in keys+['gene']]
 left=x[x.KO.isin(shared_kos)].set_index(keys)[compare_columns].sort_index()
 right=y[y.KO.isin(shared_kos)].set_index(keys)[compare_columns].sort_index()
 try:
  pd.testing.assert_frame_equal(left,right,check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-10)
 except AssertionError as exc:
  raise ValueError(f'{filename}: conflicting diagnostics for shared KOs') from exc
 z=pd.concat([x,y],ignore_index=True).drop_duplicates(keys)
 if 'gene' in z:
  shared=z.KO.isin(shared_kos)
  z.loc[shared,'gene']=z.loc[shared,'KO'].map(inv.gene)
 z.to_csv(a/'results'/filename,sep='\t',index=False)
# Rule evidence uses the unchanged complete count matrix and original sample IDs.
spec=importlib.util.spec_from_file_location('rules',r/'evaluate.py');ev=importlib.util.module_from_spec(spec);spec.loader.exec_module(ev)
mat=pd.read_csv(r/'inputs/functional_KO_counts_original.tsv.gz',sep='\t',index_col=0);ids=pd.read_csv(a/'results/sample_analysis_manifest.tsv',sep='\t').run_accession.tolist();mat=mat[ids];ks=sorted(set(re.findall(r'K\d{5}',' '.join(rr.expression))))
for th in [1,2]:
 vals={k:(mat.loc[k].to_numpy()>=th).astype('uint8') if k in mat.index else np.full(len(ids),2,dtype='uint8') for k in ks};df=pd.DataFrame({q.rule_id:ev.evaluate(ev.parse(q.expression),vals,len(ids)) for q in rr.itertuples()},index=ids);df.to_csv(a/'results'/f'rule_status_ge{th}.tsv',sep='\t')
samples=pd.read_csv(a/'results/sample_analysis_manifest.tsv',sep='\t').set_index('run_accession').loc[ids]
summary={'profiles':len(ids),'panels':len(cycles),'curated_processes':len(scope),'display_reaction_groups':len(ps),'rules':len(rr),'unique_KOs':len(inv),'available_unique_KOs':int(inv.row_available.sum()),'nitrogen_KOs':len(ni),'nitrogen_available_KOs':int(ni.row_available.sum()),'coordinate_groups':samples.coordinate_group.nunique(),'study_coordinate_groups':samples.analysis_group.nunique(),'scope':'Presentation integration of two input-matched analyses. No pooled statistical family, recalibration or new environmental significance tests.'};(a/'results/integration_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

# Complete supporting tables, including nitrogen, share the same frozen rules.
rr.to_csv(a/'curation/operational_rules.tsv',sep='\t',index=False)
scope.to_csv(a/'curation/process_scope_original.tsv',sep='\t',index=False)
rmap=[]
for p in ps:
    for rule in p['rules']:
        for k in rule['kos']:
            rmap.append({'cycle':p['cycle'],'process':p['id'],'label':p['label'],'reaction':p['reaction'],'interpretation':p['note'],'status':p['status'],'rule':rule['id'],'expression':rule['expression'],'KO':k,'name':inv.loc[k,'name'],'source':p['source']})
pd.DataFrame(rmap).to_csv(a/'results/reaction_gene_map.tsv',sep='\t',index=False)
scopes=[]
for cycle in cycles:
    pp=scope[scope.cycle==cycle];rs=rr[rr.process_id.isin(pp.id)]
    kk=sorted(set(re.findall(r'K\d{5}',' '.join(rs.expression))))
    scopes.append({'cycle':cycle,'processes':len(pp),'display_reaction_groups':sum(p['cycle']==cycle for p in ps),'rules':len(rs),'KOs':len(kk),'available':int(inv.loc[kk,'row_available'].sum()),'unavailable':int((~inv.loc[kk,'row_available']).sum()),'unresolved_processes':int((~pp.id.isin(rs.process_id)).sum())})
pd.DataFrame(scopes).to_csv(a/'results/cycle_scope.tsv',sep='\t',index=False)
rulesummary=[]
rcycle=rr.set_index('rule_id').process_id.map(scope.set_index('id').cycle)
for th in [1,2]:
    states=pd.read_csv(a/'results'/f'rule_status_ge{th}.tsv',sep='\t',index_col=0)
    for rid in states:
        z=states[rid]
        # Count each coordinate once when any of its profiles meets the rule.
        # These are descriptive reporting bins, not independent field sites.
        coordinates=samples.loc[z.index,'coordinate_group']
        rulesummary.append({'rule':rid,'cycle':rcycle[rid],'threshold':th,'met':int((z==1).sum()),'not_met':int((z==0).sum()),'indeterminate':int((z==2).sum()),'profiles':len(z),'met_coordinate_groups':int(coordinates[z==1].nunique()),'total_coordinate_groups':int(coordinates.nunique())})
pd.DataFrame(rulesummary).to_csv(a/'results/rule_evidence_summary.tsv',sep='\t',index=False)
diag=pd.read_csv(a/'results/environment_association_diagnostics.tsv',sep='\t')
rows=[]
for cycle in cycles:
    kk=inv.index[inv.cycles.str.split(';').map(lambda q:cycle in q)]
    d=diag[(diag.KO.isin(kk))&(diag.outcome=='log2_positive_group_mean')]
    row={'cycle':cycle,'positive_models':d.KO.nunique()}
    for predictor,label in [('bio1','temperature'),('soc_mean','carbon'),('phh2o_mean','pH')]:
        q=d[d.predictor==predictor];row[label+'_sign_change_or_zero']=int(((q.beta_per_sd*q.loo_min<=0)|(q.beta_per_sd*q.loo_max<=0)).sum())
    rows.append(row)
pd.DataFrame(rows).to_csv(a/'results/cycle_sensitivity_summary.tsv',sep='\t',index=False)
# Replace the older scope-only report with the integrated presentation scope.
(a/'results/extension_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
