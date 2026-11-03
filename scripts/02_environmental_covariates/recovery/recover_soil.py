import os
"""Retrieve native SoilGrids subsets; verify depth against historical controls."""
from pathlib import Path
import sys, math, json, time, concurrent.futures as cf, hashlib, threading, uuid, os
from collections import defaultdict
import requests,pandas as pd,numpy as np,rasterio
from pyproj import Transformer
sys.path.insert(0,str(Path(__file__).resolve().parent))
from recover_rasters import at
D=Path(os.environ['RECOVERY_WORKDIR']).resolve();C=D/'cache/soil_wcs';C.mkdir(exist_ok=True)
PROPS=['clay','soc','sand','silt','bdod','cec','cfvo','phh2o']
T=Transformer.from_crs('EPSG:4326','+proj=igh +datum=WGS84 +units=m +no_defs',always_xy=True)
LOCKS=defaultdict(threading.Lock)
def fetch(prop,lat,lon):
 x,y=T.transform(lon,lat)
 key=(prop,math.floor(x/250),math.floor(y/250))
 with LOCKS[key]:return _fetch(prop,lat,lon)
def _fetch(prop,lat,lon):
 x,y=T.transform(lon,lat);x0=math.floor(x/250)*250-2500;y0=math.floor(y/250)*250-2500
 key=f'{prop}_0-5cm_mean_{x0}_{y0}'
 f=C/(key+'.tif')
 p={'map':f'/map/{prop}.map','SERVICE':'WCS','VERSION':'2.0.1','REQUEST':'GetCoverage','COVERAGEID':f'{prop}_0-5cm_mean','FORMAT':'image/tiff','SUBSET':[f'X({x0},{x0+5250})',f'Y({y0},{y0+5250})'],'SUBSETTINGCRS':'http://www.opengis.net/def/crs/EPSG/0/152160'}
 if f.exists():
  try:
   with rasterio.open(f) as ss:assert ss.count==1
  except Exception:f.unlink()
 if not f.exists():
  for attempt in range(3):
   try:
    r=requests.get('https://maps.isric.org/mapserv',params=p,timeout=40);r.raise_for_status()
    with rasterio.MemoryFile(r.content) as mf:
     with mf.open() as s:assert s.count==1
    tmp=f.with_name(f.name+'.'+uuid.uuid4().hex+'.tmp');tmp.write_bytes(r.content);os.replace(tmp,f);f.with_suffix('.json').write_text(json.dumps({'url':r.url,'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'sha256':hashlib.sha256(r.content).hexdigest(),'headers':dict(r.headers)},indent=2));break
   except Exception:
    if attempt==2:raise
    time.sleep(1+attempt)
 return f

def job(p):
 r=dict(lat=p.lat,lon=p.lon,role=p.role,example_run=p.run_accession,depth='0-5cm',encoding='native_integer')
 try:
  with rasterio.open(fetch('clay',p.lat,p.lon)) as clay, rasterio.open(fetch('phh2o',p.lat,p.lon)) as ph:
   assert clay.transform==ph.transform and clay.shape==ph.shape
   data=clay.read(1); phdata=ph.read(1); profile=clay.profile.copy();profile.update(nodata=-32768)
   # WCS fills ocean with untagged zeros. pH<=0 cannot be a valid prediction.
   data[(phdata<=0)|(data<0)]=-32768
   with rasterio.MemoryFile() as mf:
    with mf.open(**profile) as masksrc:
     masksrc.write(data,1);a=at(masksrc,p.lat,p.lon)
  if not a:r['error']='no_valid_cell_within_2000m';return r
  r.update({k:v for k,v in a.items() if k!='value'})
  for prop in PROPS:
   with rasterio.open(fetch(prop,p.lat,p.lon)) as s:
    xx,yy=T.transform(a['anchor_lon'],a['anchor_lat']);v=next(s.sample([(xx,yy)],masked=True))[0]
    if not np.ma.is_masked(v) and np.isfinite(v):r[prop+'_mean']=float(v)
   if p.role=='control':r[prop+'_mean_historical']=float(p[prop+'_mean'])
 except Exception as e:r['error']=str(e)[:300]
 return r
if __name__=='__main__':
 m=pd.read_csv(D/'inputs/eligible399.csv');cols=[p+'_mean' for p in PROPS]
 todo=m[m[cols].isna().any(axis=1)].drop_duplicates(['lat','lon']).assign(role='target')
 ctrl=m[m[cols].notna().all(axis=1)].drop_duplicates(['lat','lon']).iloc[:20].assign(role='control')
 items=[p for _,p in pd.concat([ctrl,todo]).iterrows()];out=[]
 with cf.ThreadPoolExecutor(max_workers=4) as ex:
  for r in ex.map(job,items):
   out.append(r);pd.DataFrame(out).to_csv(D/'outputs/soil_verified_diagnostics.csv',index=False);print(r['example_run'],r.get('clay_mean'),r.get('clay_mean_historical'),r.get('error'),flush=True)
