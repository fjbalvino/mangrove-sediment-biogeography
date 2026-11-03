import os
"""Fill only pre-existing empty environmental cells; preserve every other string cell."""
from pathlib import Path
import csv,json,hashlib,math
import pandas as pd,numpy as np
from id_contracts import index_by_id, transfer_by_id, equivalent
D=Path(os.environ['RECOVERY_WORKDIR']).resolve()
EXCLUDE={'SRR31131347','SRR31131364'}
GROUPS={'climate':([f'bio{i}' for i in range(1,20)],'climate30_diagnostics.csv'),
'human':(['human_footprint_2020'],'human_no_osm_diagnostics.csv'),
'soil':([p+'_mean' for p in ['clay','soc','sand','silt','bdod','cec','cfvo','phh2o']],'soil_verified_diagnostics.csv')}
ENV=json.loads((D/'outputs/environmental_columns.json').read_text())
def empty(x):return str(x).strip().lower() in {'','nan','na','n/a','null','none'}
def read(p):
 with p.open(newline='',encoding='utf-8-sig') as f:
  r=csv.DictReader(f);return r.fieldnames,list(r)
def key(r):return (round(float(r['lat']),8),round(float(r['lon']),8))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,cols,rows):
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
NAMES=['Metadata_reinicio_401_reconciliada.csv','metadata_functional.csv']
EXPECTED={NAMES[0]:int(os.environ.get('EXPECTED_MASTER_ROWS','401')),NAMES[1]:int(os.environ.get('EXPECTED_FUNCTIONAL_ROWS','191'))}
original={name:read(D/'inputs'/name) for name in NAMES}
for name,(columns,rows) in original.items():
 index_by_id(rows,name,EXPECTED[name])
 if len(columns)!=len(set(columns)) or any(c not in columns for c in ['run_accession','lat','lon']+ENV):
  raise ValueError(f'{name}: invalid existing column contract')
# Check identity, membership and existing cross-table values before any output write.
transfer_by_id(original[NAMES[0]][1],original[NAMES[1]][1],ENV,EXCLUDE)
lookup={};checks={}
for group,(cols,name) in GROUPS.items():
 diag=pd.read_csv(D/'outputs'/name);assert not diag[['lat','lon','role']].duplicated().any()
 control=diag[diag.role=='control'];matches={}
 for c in cols:
  pair=control[[c,c+'_historical']].dropna();matches[c]={'n':len(pair),'exact':int(np.isclose(pair[c],pair[c+'_historical'],rtol=0,atol=1e-5).sum())}
 checks[group]=matches
 for c,comparison in matches.items():
  if comparison['n']<3 or comparison['exact']/comparison['n']<0.5:
   raise ValueError(f'Historical compatibility check failed for {group}/{c}: {comparison}')
 lookup[group]={}
 for _,r in diag[diag.role=='target'].iterrows():
  if not all(c in r and pd.notna(r[c]) and math.isfinite(float(r[c])) for c in cols):continue
  assert pd.isna(r.get('error',np.nan)),r.to_dict()
  assert r.distance_m<=2000.001
  if group=='soil':
   assert r.phh2o_mean>0 and r.bdod_mean>0
   assert abs(r.clay_mean+r.sand_mean+r.silt_mean-1000)<=2
   assert 0<r.phh2o_mean<=140
  lookup[group][key(r)]={c:format(float(r[c]),'.12g') for c in cols}
  lookup[group][key(r)]['distance_m']=float(r.distance_m)
summary={};audit=[];prepared={};master_cell_audit={}
for name in NAMES:
 inp=D/'inputs'/name;cols,before=original[name]
 counts={g:0 for g in GROUPS};touched=set()
 if name==NAMES[0]:
  after=[r.copy() for r in before]
  # Geographic lookup belongs only to extraction at the master sample's point.
  # It never joins the two sample databases.
  for r in after:
   run=r['run_accession']
   if run in EXCLUDE or empty(r['lat']) or empty(r['lon']):continue
   for g,(fields,source) in GROUPS.items():
    measurement=lookup[g].get(key(r))
    if measurement is None:continue
    for c in fields:
     if empty(r[c]):
      entry={'file':name,'run_accession':run,'field':c,'previous':r[c],'value':measurement[c],'source_diagnostic':source,'anchor_distance_m':measurement['distance_m']}
      audit.append(entry);master_cell_audit[(run,c)]=entry
      r[c]=measurement[c];counts[g]+=1;touched.add(run)
 else:
  # The only transfer between metadata databases is one-to-one by run_accession.
  after=transfer_by_id(prepared[NAMES[0]][1],before,ENV,EXCLUDE)
  before_by_id=index_by_id(before,name)
  for r in after:
   run=r['run_accession'];old=before_by_id[run]
   for c in ENV:
    if r[c]==old[c]:continue
    source=master_cell_audit.get((run,c),{})
    audit.append({'file':name,'run_accession':run,'field':c,'previous':old[c],'value':r[c],'source_diagnostic':source.get('source_diagnostic','master_by_run_accession'),'anchor_distance_m':source.get('anchor_distance_m','')})
    for g,(fields,_) in GROUPS.items():
     if c in fields:counts[g]+=1
    touched.add(run)
 if [r['run_accession'] for r in after]!=[r['run_accession'] for r in before]:
  raise ValueError(f'{name}: row membership or order changed')
 index_by_id(after,name,EXPECTED[name])
 for old,r in zip(before,after):
  if old.keys()!=r.keys():raise ValueError(f'{name}: schema changed')
  for c in cols:
   if old[c]!=r[c] and not (c in ENV and empty(old[c]) and not empty(r[c]) and r['run_accession'] not in EXCLUDE):
    raise ValueError(f'{name}: unauthorized cell change for {r["run_accession"]}/{c}')
 prepared[name]=(cols,after)
 summary[name]={'rows':len(after),'unique_run_accession':len(index_by_id(after,name)),'duplicate_ids':0,'columns':len(cols),'filled_cells':sum(counts.values()),'filled_by_group':counts,'samples_changed':len(touched),'before_sha256':sha(inp),'missing_before':{c:sum(empty(r[c]) for r in before) for c in ENV},'missing_after':{c:sum(empty(r[c]) for r in after) for c in ENV}}
master=index_by_id(prepared[NAMES[0]][1],'updated master')
for r in prepared[NAMES[1]][1]:
 for c in ['lat','lon']+ENV:
  if not equivalent(r[c],master[r['run_accession']][c]):
   raise ValueError(f'Cross-table conflict after ID transfer: {r["run_accession"]}/{c}')
# All identity and content checks have passed before writing either table.
for name,(cols,after) in prepared.items():
 out=D/'outputs'/name;save(out,cols,after)
 c2,r2=read(out)
 if c2!=cols or r2!=after:raise RuntimeError(f'{name}: serialization changed data')
 summary[name]['after_sha256']=sha(out)
save(D/'outputs/changed_cells.csv',['file','run_accession','field','previous','value','source_diagnostic','anchor_distance_m'],audit)
pending=[]
for r in master.values():
 for c in ENV:
  if not empty(r[c]):continue
  if r['run_accession'] in EXCLUDE:reason='excluded_from_recovery_no_verified_coordinates'
  elif c.startswith('gbif_'):
   reason='undefined_diversity_historical_zero_occurrences' if not empty(r['gbif_occurrences_effective']) and float(r['gbif_occurrences_effective'])==0 else 'GBIF_HTTP_503_service_unavailable'
  else:reason='no_valid_raster_cell_within_2000m'
  pending.append({'run_accession':r['run_accession'],'field':c,'reason':reason})
save(D/'outputs/pending_cells.csv',['run_accession','field','reason'],pending)
summary['primary_key']='run_accession';summary['metadata_join']='one_to_one_by_run_accession';summary['controls']=checks;summary['extraction_exclusions']=sorted(EXCLUDE);summary['validation']='PASS: only blank existing environmental cells changed; schemas, rows, all other cells and both exclusion rows preserved; shared cohort matches'
(D/'outputs/update_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps({k:{a:v for a,v in x.items() if a not in ['missing_before','missing_after']} for k,x in summary.items() if k.endswith('.csv')},indent=2))
