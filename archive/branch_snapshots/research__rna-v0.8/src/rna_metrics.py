"""Quantitative effects only: no correlation of treated expression with baseline."""
from __future__ import annotations
import numpy as np

def row_correlation(a,b,center=True):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    if center:a=a-a.mean(axis=1,keepdims=True);b=b-b.mean(axis=1,keepdims=True)
    den=np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1))
    return np.divide(np.sum(a*b,axis=1),den,out=np.zeros(len(a)),where=den>1e-14)

def effect_metrics(prediction,truth):
    p=np.asarray(prediction,dtype=float);t=np.asarray(truth,dtype=float)
    if p.ndim!=2 or p.shape!=t.shape or len(p)<2:raise ValueError('Need matching >=2 condition effect matrices')
    if not np.isfinite(p).all() or not np.isfinite(t).all():raise ValueError('Missing values are not zero effects')
    e=p-t; bias=e.mean(axis=0);cmse=float(np.mean((e-bias)**2));bmse=float(np.mean(bias**2))
    pc=p-p.mean(axis=0);tc=t-t.mean(axis=0)
    pc=pc-pc.mean(axis=1,keepdims=True);tc=tc-tc.mean(axis=1,keepdims=True)
    pn=np.linalg.norm(pc,axis=1);tn=np.linalg.norm(tc,axis=1)
    sim=(pc@tc.T)/np.maximum(pn[:,None]*tn[None,:],1e-14)
    diag=np.diag(sim)
    better=(sim>diag[:,None]+1e-12).sum(axis=1);ties=np.isclose(sim,diag[:,None],atol=1e-12,rtol=0).sum(axis=1)
    retrieval=float(np.mean(np.where(better==0,1/np.maximum(ties,1),0)))
    return {'mse':float(np.mean(e*e)),'centered_mse':cmse,'bias_mse':bmse,
            'mae':float(np.mean(np.abs(e))),'effect_pearson':float(row_correlation(p,t).mean()),
            'effect_cosine':float(row_correlation(p,t,False).mean()),
            'centered_effect_pearson':float(row_correlation(p-p.mean(axis=0),t-t.mean(axis=0)).mean()),
            'centered_retrieval_top1':retrieval,'query_count':len(p),'gene_count':p.shape[1]}
