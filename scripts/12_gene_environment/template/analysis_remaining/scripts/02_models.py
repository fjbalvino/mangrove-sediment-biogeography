"""Numerical fixed-effect fits and the official CR2 working-model construction.

This implementation must pass the accompanying R/clubSandwich cross-check
before use. No biological data p-values are exported.
"""
import numpy as np
from scipy.special import expit,xlogy
from scipy.stats import t

def matrix_power(a,p):
    v,q=np.linalg.eigh((a+a.T)/2)
    # Exact default tolerance in clubSandwich 0.5.10.
    keep=v>1e-12
    return (q[:,keep]*(v[keep]**p))@q[:,keep].T

def fit_component(y,component,x,study,recovery,details=False):
    ans=np.full((2,3,4),np.nan)
    mask=np.isfinite(y)
    for s in np.unique(study):
        sel=(study==s)&mask
        if sel.sum()<2 or (component=='detection' and (np.all(y[sel]==0) or np.all(y[sel]==1))):mask[sel]=False
    yy=y[mask];xx=x[mask];ss=study[mask];levels=np.unique(ss)
    z=np.c_[(ss[:,None]==levels).astype(float),xx]
    if component=='detection':z=np.c_[z,recovery[mask]]
    n,p=z.shape;df=n-p;ix=np.arange(len(levels),len(levels)+3)
    if len(levels)<2 or df<=1 or np.linalg.matrix_rank(z)<p:return ans
    try:
        if component=='positive':
            inv=np.linalg.inv(z.T@z);coef=inv@z.T@yy;resid=yy-z@coef
            v=np.ones(n);u=z;bread=inv;disp=resid@resid/df
        else:
            mu=(yy+.5)/2;eta=np.log(mu/(1-mu))
            olddev=2*np.sum(xlogy(yy,yy/mu)+xlogy(1-yy,(1-yy)/(1-mu)))
            converged=False
            for iteration in range(60):
                w=mu*(1-mu)
                if w.min()<1e-14:return ans
                zw=z*np.sqrt(w)[:,None]
                target=(eta+(yy-mu)/w)*np.sqrt(w)
                coef=np.linalg.lstsq(zw,target,rcond=1e-12)[0]
                eta=z@coef;mu=expit(eta)
                if np.any(mu<=0)|np.any(mu>=1):return ans
                dev=2*np.sum(xlogy(yy,yy/mu)+xlogy(1-yy,(1-yy)/(1-mu)))
                if abs(dev-olddev)/(.1+abs(dev))<1e-11:converged=True;break
                olddev=dev
            if not converged or mu.min()<1e-8 or mu.max()>1-1e-8:return ans
            v=mu*(1-mu);u=v[:,None]*z
            # R's fitted-object bread and Pearson dispersion use last IRLS weights.
            bread=np.linalg.inv(z.T@(w[:,None]*z));resid=yy-mu
            disp=np.sum(w*(resid/v)**2)/df
        beta=coef[ix];se=np.sqrt(np.diag(bread)[ix]*disp)
        ans[0]=np.c_[beta,se,np.full(3,df),2*t.sf(abs(beta/se),df)]
        # This is the original-scale clubSandwich CR2 adjustment, not a
        # substituted symmetric adjustment on a whitened regression.
        m=matrix_power(z.T@(v[:,None]*z),-1)
        scores=[];hs=[];gdiag=[]
        mb=np.linalg.cholesky(bread)
        for s in levels:
            ids=np.flatnonzero(ss==s);vv=v[ids]
            ih=np.eye(len(ids))-u[ids]@m@z[ids].T
            gv=np.sqrt(vv)[:,None]*(ih*vv[None,:])*np.sqrt(vv)[None,:]
            adj=np.sqrt(vv)[:,None]*matrix_power(gv,-.5)*np.sqrt(vv)[None,:]
            e=bread[ix]@z[ids].T@adj
            scores.append(e@resid[ids])
            gdiag.append(np.sum(e*e*v[ids][None,:],axis=1))
            hs.append(e@u[ids]@mb)
        scores=np.asarray(scores);hs=np.asarray(hs);gdiag=np.asarray(gdiag)
        se=np.sqrt(np.sum(scores*scores,axis=0));dfs=[]
        for j in range(3):
            gram=-(hs[:,j]@hs[:,j].T)
            gram[np.diag_indices(len(levels))]+=gdiag[:,j]
            dfs.append(np.trace(gram)**2/np.sum(gram*gram))
        dfs=np.asarray(dfs)
        ans[1]=np.c_[beta,se,dfs,2*t.sf(abs(beta/se),dfs)]
        return ans
    except (ValueError,np.linalg.LinAlgError,FloatingPointError):return ans
