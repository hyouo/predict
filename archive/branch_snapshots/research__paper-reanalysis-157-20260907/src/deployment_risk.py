"""Blend target-anchor LOO risk with source-only deployment-shaped validation.

The target query outcomes are NOT arguments. Each observed source is used as a
pseudo-target, and models are rebuilt from the remaining source contexts. The
query intervention catalog is known and is measured in source contexts.

This transfers a risk estimate across contexts. It requires a working assumption
that source context behavior informs target context behavior, not a guarantee of
causal transportability or calibrated biological uncertainty.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from .covariance import build_kernels,Fit
from .aggregation import bank


def stack_gram(C,penalty_ratio=.01):
    C=np.asarray(C,float);m=len(C);prior=np.ones(m)/m
    penalty=penalty_ratio*max(float(np.trace(C)/m),1e-12)
    fun=lambda w:float(w@C@w+penalty*np.sum((w-prior)**2))
    jac=lambda w:2*(C@w+penalty*(w-prior))
    result=minimize(fun,prior,jac=jac,bounds=[(0,1)]*m,method='SLSQP',
        constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(m)},
        options={'ftol':1e-10,'maxiter':500})
    if not result.success:raise RuntimeError('Risk stack failed '+result.message)
    w=np.maximum(result.x,0);return w/w.sum()

def fit_deployment(source,cues,anchors,y_anchor,queries,kernels=None,
                   families=('cue','joint','output_disagreement','output_marginal'),
                   weights=(.25,.5,1.)):
    source=np.asarray(source,float);a=np.asarray(anchors,int);q=np.asarray(queries,int)
    if source.shape[0]<3 or np.intersect1d(a,q).size:raise ValueError('Need >=3 sources and disjoint query set')
    k=build_kernels(source,cues) if kernels is None else kernels
    pp,ee,pars=bank(source,cues,k,a,np.asarray(y_anchor,float),families)
    E=ee.reshape(len(pp),-1);Ct=E@E.T/E.shape[1]
    sourcegrams=[]
    for j in range(len(source)):
        other=np.delete(source,j,axis=0)
        kj=build_kernels(other,cues)
        pj,_,pnames=bank(other,cues,kj,a,source[j,a].copy(),families)
        if pars!=pnames:raise AssertionError('Model slots not aligned across source folds')
        err=pj[:,q,:]-source[j,q][None,:,:]
        e=err.reshape(len(pp),-1)
        sourcegrams.append(e@e.T/e.shape[1])
    Cs=np.mean(sourcegrams,axis=0)
    result={}
    for alpha in weights:
        # Risk-mixing weights, not claimed to equal effective sample counts.
        C=(1-alpha)*Ct+alpha*Cs
        w=stack_gram(C)
        p=np.einsum('m,mng->ng',w,pp)
        result[f'meta_{alpha}']=Fit(p,{'meta_risk_weight':alpha,'weights':w.tolist(),'components':pars,
            'source_pseudoquery_risk':float(w@Cs@w),'anchor_oof_training_risk':float(w@Ct@w),
            'source_context_count':len(source),'query_count':len(q),
            'target_query_outcomes_read':False},np.empty(0))
    return result
