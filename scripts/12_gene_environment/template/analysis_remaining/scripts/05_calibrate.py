"""Joint current-panel KO conditional calibration; resumable, provenance guarded chunks."""
from pathlib import Path
import os,sys,json,hashlib,importlib.util,itertools,time,warnings
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np,pandas as pd
from scipy.special import ndtr,expit
from scipy.stats import t
r=Path(__file__).resolve().parents[1];o=r/'calibration_final';o.mkdir(exist_ok=True)
s=importlib.util.spec_from_file_location('models',r/'scripts/02_models.py');model=importlib.util.module_from_spec(s);s.loader.exec_module(model)
plan=json.loads((o/'simulation_plan.json').read_text());validation=json.loads((o/'implementation_validation.json').read_text());assert validation['passed']
g=pd.read_csv(o/'design.tsv',sep='\t');u=pd.read_csv(o/'run_weights.tsv',sep='\t')
alpha=pd.read_csv(o/'detection_intercepts.tsv',sep='\t',index_col=0);genes=alpha.columns.tolist();alpha=alpha.to_numpy();mu=pd.read_csv(o/'positive_intercepts.tsv',sep='\t',index_col=0).to_numpy()
K=len(genes);N=len(g);U=len(u);X=g[['x1','x2','x3']].to_numpy();study=g.study_id.to_numpy();recovery=g.recovery.to_numpy()
levels,si=np.unique(study,return_inverse=True);rg=pd.Index(g.analysis_group).get_indexer(u.analysis_group)
W=np.zeros((N,U));W[rg,np.arange(U)]=u.descriptive_run_weight_within_group;assert np.max(abs(W.sum(1)-1))<1e-12
def spatial(lat,lon):
    a=np.deg2rad(lat);b=np.deg2rad(lon);h=np.sin((a[:,None]-a)/2)**2+np.cos(a)[:,None]*np.cos(a)*np.sin((b[:,None]-b)/2)**2
    distance=12742*np.arcsin(np.sqrt(np.clip(h,0,1)))
    return np.linalg.cholesky(.5*np.eye(len(a))+.5*np.exp(-distance/50))
Lg=spatial(g.lat.to_numpy(),g.lon.to_numpy());Lr=spatial(g.lat.to_numpy()[rg],g.lon.to_numpy()[rg])
signalk=np.isin(genes,plan['signal_genes']);truth_alt=np.zeros((K,2,2,3));truth_alt[signalk,0,:,:]=np.log(2);truth_alt[signalk,1,:,:]=.5
idx=pd.DataFrame(itertools.product(genes,['detection','positive'],['conventional','CR2'],['x1','x2','x3']),columns=['gene','component','method','predictor'])
names=['valid','b','b2','se','se2','df','coverage','reject','width','truth']
def one(scenario,signal,replicate):
    seed=plan['seed']+plan['scenarios'].index(scenario)*10000+(5000 if signal=='mixed' else 0)+replicate
    rng=np.random.default_rng(seed)
    def ng(n):return np.sqrt(.4)*rng.standard_normal((n,1))+np.sqrt(.6)*rng.standard_normal((n,K))
    latent=ng(U)
    if scenario=='study_dependence':latent=np.sqrt(.4)*ng(len(levels))[si[rg]]+np.sqrt(.3)*ng(N)[rg]+np.sqrt(.3)*latent
    if scenario=='spatial_50km':latent=Lr@latent
    bd=np.zeros((3,K));bp=bd.copy()
    if signal=='mixed':bd[:,signalk]=np.log(2);bp[:,signalk]=.5
    prob=expit(alpha+.75*recovery[:,None]+X@bd);fraction=W@(ndtr(latent)<prob[rg])
    err=ng(N)
    if scenario=='study_dependence':
        for j in range(3):err+=X[:,j,None]*ng(len(levels))[si]
    if scenario=='dominant_heterogeneity':err*=np.where(study==plan['dominant_study'],4.,1.)[:,None]
    if scenario=='spatial_50km':err=Lg@err
    if scenario=='heavy_tails':err=t.ppf(np.clip(ndtr(err),1e-12,1-1e-12),3)/np.sqrt(3)
    positive=mu+X@bp+err;positive[fraction==0]=np.nan
    a=np.full((K,2,2,3,4),np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        for k in range(K):
            a[k,0]=model.fit_component(fraction[:,k],'detection',X,study,recovery)
            a[k,1]=model.fit_component(positive[:,k],'positive',X,study,recovery)
        b,se,df,p=np.moveaxis(a,-1,0);width=2*t.ppf(.975,df)*se
    truth=truth_alt if signal=='mixed' else np.zeros_like(truth_alt)
    valid=np.isfinite(a).all(-1)&(se>1e-10)&(df>0)&np.isfinite(width)
    p[~valid]=1;coverage=abs(b-truth)<=width/2
    values=np.stack([valid,b,b*b,se,se*se,df,coverage,p<.05,width,truth],axis=-1).astype(float)
    values[~valid,1:9]=0
    fdr=[]
    for j,method in enumerate(['conventional','CR2']):
        pp=p[:,:,j,:].ravel();tr=(truth[:,:,j,:]!=0).ravel();order=np.argsort(pp,kind='stable');ok=pp[order]<=.05*np.arange(1,len(pp)+1)/len(pp)
        reject=np.zeros(len(pp),bool)
        if ok.any():reject[order[:np.flatnonzero(ok)[-1]+1]]=True
        false=int((reject&~tr).sum());true=int((reject&tr).sum());count=int(reject.sum())
        fdr.append(dict(scenario=scenario,signal=signal,replicate=replicate,method=method,discoveries=count,false=false,true=true,n_nonnull=int(tr.sum()),FDP=false/max(1,count),power=true/tr.sum() if tr.any() else np.nan,valid=int(valid[:,:,j,:].sum())))
    return values.reshape(-1,10),fdr
def chunk(task):
    scenario,signal,number=task;dest=o/'chunks'/f'{scenario}_{signal}_{number:03d}.npz'
    if dest.exists() and dest.with_suffix('.tsv').exists():return dest.name
    count=25 if signal=='null' else 25;start=number*count
    sums=None;rows=[]
    for i in range(start+1,start+count+1):
        v,f=one(scenario,signal,i);sums=v if sums is None else sums+v;rows+=f
    tmp=dest.with_suffix('.tmp.npz');np.savez_compressed(tmp,sums=sums,replicates=count);tmp.replace(dest)
    pd.DataFrame(rows).to_csv(dest.with_suffix('.tsv'),sep='\t',index=False)
    return dest.name
if __name__=='__main__':
    for name,digest in plan['input_hashes'].items():
        assert hashlib.sha256((o/name).read_bytes()).hexdigest()==digest, 'Prepared simulation input changed: '+name
    if '--benchmark' in sys.argv:
        tm=time.time();v,f=one('independent','mixed',1);print('Seconds per joint dataset',time.time()-tm,'metrics',v.shape,'FDR',f);sys.exit()
    paths=[r/'scripts/02_models.py',Path(__file__),o/'simulation_plan.json',o/'implementation_validation.json']
    provenance={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    guard=o/'execution_provenance.json'
    if guard.exists():assert json.loads(guard.read_text())==provenance,'Changed simulation code/plan: use a fresh output directory'
    else:guard.write_text(json.dumps(provenance,indent=2)+'\n')
    (o/'chunks').mkdir(exist_ok=True)
    tasks=[(s,signal,c) for s in plan['scenarios'] for signal in ['null','mixed'] for c in range(40 if signal=='null' else 20)]
    start=time.time()
    with ProcessPoolExecutor(max_workers=int(os.environ.get('GENE_ATLAS_WORKERS','2'))) as pool:
        futures=[pool.submit(chunk,q) for q in tasks]
        for i,f in enumerate(as_completed(futures)):print(f'{i+1}/{len(tasks)} {f.result()}',flush=True)
    statistics=[];fdr=[]
    for scenario,signal,num in tasks:
        path=o/'chunks'/f'{scenario}_{signal}_{num:03d}.npz';part=np.load(path)
        df=pd.concat([idx,pd.DataFrame(part['sums'],columns=names)],axis=1);df.insert(0,'signal',signal);df.insert(0,'scenario',scenario);df['replicates']=int(part['replicates']);statistics.append(df)
        fdr.append(pd.read_csv(path.with_suffix('.tsv'),sep='\t'))
    pd.concat(statistics).to_csv(o/'chunk_sufficient_statistics.tsv.gz',sep='\t',index=False)
    pd.concat(fdr).to_csv(o/'replicate_FDR.tsv',sep='\t',index=False)
    (o/'execution_summary.json').write_text(json.dumps({'joint_datasets':7500,'genes':K,'family_tests':K*6,'planted_nonnulls':len(plan['signal_genes'])*6,'elapsed_minutes':(time.time()-start)/60,'chunks':len(tasks)},indent=2)+'\n')
    print('COMPLETE',flush=True)
