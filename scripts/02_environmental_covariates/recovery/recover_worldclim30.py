import os
from pathlib import Path
import pandas as pd,numpy as np,rasterio
from rasterio.windows import Window
from recover_rasters import at
D=Path(os.environ['RECOVERY_WORKDIR']).resolve();m=pd.read_csv(D/'inputs/eligible399.csv');cols=[f'bio{i}' for i in range(1,20)]
todo=m[m[cols].isna().any(axis=1)].drop_duplicates(['lat','lon']).assign(role='target')
ctrl=m[m[cols].notna().all(axis=1)].drop_duplicates(['lat','lon']).iloc[:20].assign(role='control')
out=[]
for _,p in pd.concat([todo,ctrl]).iterrows():
 tile=int((90-p.lat)//30)*12+int((p.lon+180)//30)+1;f=D/f'cache/tile_{tile}_wc2.1_30s_bio.tif'
 if not f.exists():
  if p.role=='target':raise FileNotFoundError(f)
  continue
 r=dict(lat=p.lat,lon=p.lon,role=p.role,example_run=p.run_accession,source_file=f.name)
 with rasterio.open(f) as s:
  assert s.count==19 and abs(s.res[0]-1/120)<1e-10
  a=at(s,p.lat,p.lon)
  if a:
   r.update({k:v for k,v in a.items() if k!='value'})
   values=s.read(window=Window(a['col'],a['row'],1,1),masked=True)[:,0,0]
   for c,v in zip(cols,values):
    if not np.ma.is_masked(v) and np.isfinite(v):r[c]=float(v)
    if p.role=='control':r[c+'_historical']=float(p[c])
  else:r['error']='no_valid_cell_within_2000m'
 out.append(r);print(r['example_run'],r.get('bio1'),r.get('bio1_historical'),r.get('distance_m'),r.get('error'),flush=True)
pd.DataFrame(out).to_csv(D/'outputs/climate30_diagnostics.csv',index=False)
