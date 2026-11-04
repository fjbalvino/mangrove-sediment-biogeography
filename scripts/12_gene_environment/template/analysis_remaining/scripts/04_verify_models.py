from pathlib import Path
import importlib.util,json,subprocess,time
import numpy as np,pandas as pd
from scipy.special import expit
r=Path(__file__).resolve().parents[1];o=r/'calibration_final'
s=importlib.util.spec_from_file_location('models',r/'scripts/02_models.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
g=pd.read_csv(o/'design.tsv',sep='\t');x=g[['x1','x2','x3']].to_numpy();ss=g.study_id.to_numpy();rec=g.recovery.to_numpy()
u=pd.read_csv(o/'run_weights.tsv',sep='\t');rg=pd.Index(g.analysis_group).get_indexer(u.analysis_group);W=np.zeros((len(g),len(u)));W[rg,np.arange(len(u))]=u.descriptive_run_weight_within_group
al=pd.read_csv(o/'detection_intercepts.tsv',sep='\t',index_col=0).to_numpy();rng=np.random.default_rng(2026110402)
cases=[];out=[];tm=time.time()
for j in range(240):
    component='positive' if j%2 else 'detection'
    p=expit(al[:,j%al.shape[1]]+.75*rec+(x@np.array([.3,-.2,.4]) if j%3 else 0))
    f=W@(rng.random(len(u))<p[rg])
    y=(x@np.array([.3,-.2,.4])+rng.standard_t(3,len(g))) if component=='positive' else f
    if component=='positive':y[f==0]=np.nan
    if j in [238,239]:y[:]=0 # deliberate degenerate boundary
    z=m.fit_component(y,component,x,ss,rec)
    for a,method in enumerate(['conventional','CR2']):
        for k in range(3):out.append(dict(case=j,component=component,method=method,predictor=k+1,beta=z[a,k,0],se=z[a,k,1],df=z[a,k,2],p=z[a,k,3]))
    df=pd.DataFrame({'case':j,'component':component,'y':y});cases.append(df)
pd.concat(cases).to_csv(o/'verification_cases.tsv',sep='\t',index=False);pd.DataFrame(out).to_csv(o/'verification_python.tsv',sep='\t',index=False)
print('Python 240 fits seconds',time.time()-tm,flush=True)
subprocess.run(['Rscript',str(r/'scripts/04_verify_models.R'),str(r)],check=True)
a=pd.read_csv(o/'verification_python.tsv',sep='\t');b=pd.read_csv(o/'verification_R.tsv',sep='\t')
c=a.merge(b,on=['case','component','method','predictor'],suffixes=('_py','_R'),validate='one_to_one')
stats={};passed=True
for key in ['beta','se','df','p']:
    aa=c[key+'_py'].to_numpy();bb=c[key+'_R'].to_numpy();finite=np.isfinite(aa)&np.isfinite(bb)
    equalmask=np.array_equal(np.isfinite(aa),np.isfinite(bb));err=abs(aa[finite]-bb[finite]);rel=err/(1+abs(bb[finite]));good=equalmask and (rel.max(initial=0)<1e-5)
    stats[key]={'matched_finite_pattern':equalmask,'max_absolute_error':float(err.max(initial=0)),'max_scaled_error':float(rel.max(initial=0)),'passed':bool(good)};passed &=good
result={'passed':bool(passed),'fits':240,'coefficient_comparisons':len(c),'tolerance':'abs(Python-R)/(1+abs(R)) < 1e-5 and identical finite patterns','metrics':stats,'R':subprocess.check_output(['Rscript','--version'],text=True).strip(),'clubSandwich':subprocess.check_output(['Rscript','-e','cat(as.character(packageVersion("clubSandwich")))'],text=True).strip(),'epsilon':1e-11}
(o/'implementation_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert passed
