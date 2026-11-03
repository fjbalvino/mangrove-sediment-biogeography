import os
"""Fill only pre-existing empty environmental cells; preserve every other string cell."""
from pathlib import Path
import csv,json,hashlib,math
import pandas as pd,numpy as np
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
summary={};audit=[]
for name in ['Metadata_reinicio_401_reconciliada.csv','metadata_functional.csv']:
 inp=D/'inputs'/name;cols,before=read(inp);after=[r.copy() for r in before];counts={g:0 for g in GROUPS};touched=set()
 assert len({r['run_accession'] for r in before})==len(before)
 for r in after:
  if r['run_accession'] in EXCLUDE or empty(r['lat']) or empty(r['lon']):continue
  for g,(fields,source) in GROUPS.items():
   donor=lookup[g].get(key(r))
   if donor is None:continue
   for c in fields:
    assert c in cols
    if empty(r[c]):
     audit.append({'file':name,'run_accession':r['run_accession'],'field':c,'previous':r[c],'value':donor[c],'source_diagnostic':source,'anchor_distance_m':donor['distance_m']})
     r[c]=donor[c];counts[g]+=1;touched.add(r['run_accession'])
 for a,b in zip(before,after):
  assert a.keys()==b.keys()
  for c in cols:
   if a[c]!=b[c]:assert c in ENV and empty(a[c]) and not empty(b[c]) and a['run_accession'] not in EXCLUDE
 out=D/'outputs'/name;save(out,cols,after)
 c2,r2=read(out);assert c2==cols and r2==after
 summary[name]={'rows':len(after),'columns':len(cols),'filled_cells':sum(counts.values()),'filled_by_group':counts,'samples_changed':len(touched),'before_sha256':sha(inp),'after_sha256':sha(out),'missing_before':{c:sum(empty(r[c]) for r in before) for c in ENV},'missing_after':{c:sum(empty(r[c]) for r in after) for c in ENV}}
master={r['run_accession']:r for r in read(D/'outputs/Metadata_reinicio_401_reconciliada.csv')[1]}
for r in read(D/'outputs/metadata_functional.csv')[1]:
 for c in ['lat','lon']+ENV:
  a,b=r[c],master[r['run_accession']][c]
  assert (empty(a) and empty(b)) or (not empty(a) and not empty(b) and math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-8)),(r['run_accession'],c)
save(D/'outputs/changed_cells.csv',list(audit[0]),audit)
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
summary['controls']=checks;summary['extraction_exclusions']=sorted(EXCLUDE);summary['validation']='PASS: only blank existing environmental cells changed; schemas, rows, all other cells and both exclusion rows preserved; shared cohort matches'
(D/'outputs/update_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps({k:{a:v for a,v in x.items() if a not in ['missing_before','missing_after']} for k,x in summary.items() if k.endswith('.csv')},indent=2))
