import os
"""Audited raster extraction; output diagnostics only, never overwrite metadata."""
from pathlib import Path
import sys,json, math, os
import numpy as np,pandas as pd,rasterio
from rasterio.warp import transform
from rasterio.windows import Window
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from metadata_pipeline_common import nearest_valid_raster_cell,haversine_m
D=Path(os.environ['RECOVERY_WORKDIR']).resolve()
m=pd.read_csv(D/'inputs/eligible399.csv')
def at(src,lat,lon):
 x,y=transform('EPSG:4326',src.crs,[lon],[lat]); row,col=src.index(x[0],y[0])
 if 0<=row<src.height and 0<=col<src.width:
  v=src.read(1,window=Window(col,row,1,1),masked=True)
  if not np.ma.is_masked(v[0,0]) and np.isfinite(v[0,0]):
   ax,ay=src.xy(row,col);al,at=transform(src.crs,'EPSG:4326',[ax],[ay])
   return dict(value=float(v[0,0]),row=row,col=col,anchor_lat=at[0],anchor_lon=al[0],distance_m=haversine_m(lat,lon,at[0],al[0]),method='sample_cell')
 return nearest_valid_raster_cell(src,lat,lon,2000)
def run(kind):
 if kind!='human':raise ValueError('Use recover_worldclim30.py for climate')
 cols=['human_footprint_2020'] if kind=='human' else [f'bio{i}' for i in range(1,20)]
 todo=m[m[cols].isna().any(axis=1)].drop_duplicates(['lat','lon'])
 controls=m[m[cols].notna().all(axis=1)].drop_duplicates(['lat','lon']).iloc[:20]
 target=pd.concat([todo.assign(role='target'),controls.assign(role='control')])
 sources={cols[0]:str(D/'cache/hii_no_osm_2020-01-01.tif')} if kind=='human' else {c:str(D/'cache/wc2.1_2.5m_bio'/f'wc2.1_2.5m_bio_{c[3:]}.tif') for c in cols}
 opened={c:rasterio.open(p) for c,p in sources.items()}; ref=opened[cols[0]];records=[]
 try:
  for _,p in target.iterrows():
   r=dict(lat=p.lat,lon=p.lon,role=p.role,example_run=p.run_accession)
   a=at(ref,p.lat,p.lon)
   if a:
    r.update({k:v for k,v in a.items() if k!='value'})
    for c,s in opened.items():
     v=s.read(1,window=Window(a['col'],a['row'],1,1),masked=True)[0,0]
     if not np.ma.is_masked(v) and np.isfinite(v):r[c]=float(v)
     if p.role=='control':r[c+'_historical']=float(p[c])
   else:r['error']='no_valid_cell_within_2000m'
   records.append(r);print(kind,r['example_run'],r.get(cols[0]),r.get('distance_m'),r.get('error'),flush=True)
   pd.DataFrame(records).to_csv(D/f'outputs/{kind}_no_osm_diagnostics.csv',index=False)
 finally:
  for s in opened.values():s.close()
if __name__=='__main__':
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif',GDAL_HTTP_TIMEOUT='30',GDAL_HTTP_CAINFO=os.environ.get('SSL_CERT_FILE','/etc/ssl/certs/ca-certificates.crt')):
  run(sys.argv[1])
