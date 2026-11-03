import os
"""Re-query all coordinates with missing GBIF fields; never turn API failure into zero."""
from pathlib import Path
import sys,importlib.util,json,time,concurrent.futures as cf,urllib.parse,urllib.request
from collections import Counter
import pandas as pd,numpy as np
D=Path(os.environ['RECOVERY_WORKDIR']).resolve()
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
spec=importlib.util.spec_from_file_location('gbif',Path(__file__).resolve().parent.parent/'02D_EXTRACT_GBIF_PLANT_DIVERSITY_2KM.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
from metadata_pipeline_common import haversine_m
C=D/'cache/gbif';C.mkdir(exist_ok=True)
limiter=mod.RateLimiter(1)
def query(p,params,page):
 f=C/f'{p.lat:.6f}_{p.lon:.6f}_page{page}.json'
 if f.exists():return json.loads(f.read_text())
 u=mod.GBIF_OCCURRENCE_URL+'?'+urllib.parse.urlencode(params);limiter.acquire()
 try:
  req=urllib.request.Request(u,headers={'User-Agent':'MAGENTA-metadata-recovery/2026','Accept':'application/json'})
  with urllib.request.urlopen(req,timeout=20) as r:payload=json.load(r)
  f.write_text(json.dumps(payload));return payload
 except Exception as e:
  f.with_suffix('.error.json').write_text(json.dumps({'url':u,'error':repr(e),'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}));raise

def job(p):
 r=dict(lat=p.lat,lon=p.lon,example_run=p.run_accession,role='target',previous_zero=bool(p.gbif_occurrences_effective==0))
 params={'hasCoordinate':'true','hasGeospatialIssue':'false','occurrenceStatus':'PRESENT','kingdomKey':6,'geoDistance':f'{p.lat:.6f},{p.lon:.6f},2000m','limit':300,'offset':0};counts=Counter();seen=set()
 try:
  page=0
  while True:
   j=query(p,params,page);items=j.get('results',[]);total=int(j.get('count',0))
   if total>100000:raise ValueError('GBIF query exceeds API paging limit; download required')
   for x in items:
    key=x.get('key');sk=x.get('speciesKey')
    if key in seen or not sk:continue
    seen.add(key)
    if haversine_m(p.lat,p.lon,float(x['decimalLatitude']),float(x['decimalLongitude']))>2000.001:raise ValueError('GBIF returned record outside requested 2km radius')
    counts[int(sk)]+=1
   params['offset']+=len(items);page+=1
   if j.get('endOfRecords') or params['offset']>=total:break
   if not items:raise ValueError('incomplete GBIF pagination')
  r['reported_count']=total;r['species_counts_n']=sum(counts.values())
  if r['previous_zero']:
   r['error']='historical_zero_block_preserved; query_success';return r
  metric=mod.diversity_metrics(np.array(list(counts.values()),dtype=int))
  r['gbif_occurrences_effective']=total
  for k,c in {'richness':'gbif_richness_species_2km','shannon':'gbif_shannon_H_2km','pielou':'gbif_pielou_J_2km','simpson':'gbif_simpson_1minD_2km','chao1':'gbif_chao1_2km'}.items():r[c]=metric[k]
 except Exception as e:r['error']=repr(e)
 return r
m=pd.read_csv(D/'inputs/eligible399.csv');cols=[c for c in m if c.startswith('gbif_')];todo=m[m[cols].isna().any(axis=1)].drop_duplicates(['lat','lon']);out=[]
with cf.ThreadPoolExecutor(max_workers=2) as ex:
 for r in ex.map(job,[p for _,p in todo.iterrows()]):
  out.append(r);pd.DataFrame(out).to_csv(D/'outputs/gbif_diagnostics.csv',index=False);print(r['example_run'],r.get('reported_count'),r.get('error'),flush=True)
