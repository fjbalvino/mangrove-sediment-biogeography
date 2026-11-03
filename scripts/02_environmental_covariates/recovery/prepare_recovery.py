from pathlib import Path
import pandas as pd, json, shutil, os, argparse
a=argparse.ArgumentParser();a.add_argument('--master',required=True);a.add_argument('--functional',required=True);a=a.parse_args()
D=Path(os.environ['RECOVERY_WORKDIR']).resolve();D.mkdir(parents=True,exist_ok=True)
M=Path(a.master);F=Path(a.functional)
if M.name!='Metadata_reinicio_401_reconciliada.csv' or F.name!='metadata_functional.csv':raise ValueError('Use the canonical filenames')
(D/'inputs').mkdir(exist_ok=True); (D/'outputs').mkdir(exist_ok=True); (D/'cache').mkdir(exist_ok=True)
for p in [M,F]:
 if p.resolve()!=(D/'inputs'/p.name).resolve():shutil.copyfile(p,D/'inputs'/p.name)
m=pd.read_csv(M);f=pd.read_csv(F)
cols=[c for c in m if c.startswith('bio') and c[3:].isdigit() or c.startswith('gbif_') or c.endswith('_mean') or c=='human_footprint_2020']
e=m.loc[~m.run_accession.isin(['SRR31131347','SRR31131364']) & m.lat.notna() & m.lon.notna()].copy()
e.to_csv(D/'inputs/eligible399.csv',index=False)
t=e.loc[e[cols].isna().any(axis=1)]
t[['run_accession','lat','lon','Country','Locality']+cols].to_csv(D/'outputs/recovery_targets.csv',index=False)
print('targets',len(t),'coords',len(t[['lat','lon']].drop_duplicates()))
for label,cc in {'climate':[c for c in cols if c.startswith('bio')],'human':['human_footprint_2020'],'gbif':[c for c in cols if c.startswith('gbif')],'soil':[c for c in cols if c.endswith('_mean')]}.items():
 sub=e[e[cc].isna().any(axis=1)]
 print(label,len(sub),'coords',len(sub[['lat','lon']].drop_duplicates()))
 print(sub[['run_accession','lat','lon','Country']].to_string(index=False))
(D/'outputs/environmental_columns.json').write_text(json.dumps(cols))
