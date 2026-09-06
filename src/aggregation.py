"""Output-specific covariance transfer and aggregation of fixed-parameter predictors.

This is classical kernel shrinkage/model averaging, not a new mechanistic theorem.
Aggregation resamples complete conditions, not individual proteins as independent
experiments. Bootstrap weights are an algorithm, NOT confidence intervals.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from .covariance import Fit,LAMBDAS,MIXES,normalize_kernel

def output_kernels(source,kind='disagreement'):
    s=np.asarray(source,float);S,N,G=s.shape
    if S<2:raise ValueError('>=2 source contexts required')
    d=s-s.mean(0,keepdims=True) if kind=='disagreement' else s.copy()
    d=d-d.mean(1,keepdims=True)
    return np.stack([normalize_kernel(d[:,:,g].T@d[:,:,g]) for g in range(G)])

def path(K,anchors,y,lambdas=LAMBDAS):
    """Batched KRR by output. K may be NxN (shared) or GxNxN (output-specific)."""
    a=np.asarray(anchors,int);y=np.asarray(y,float);K=np.asarray(K,float)
    n,g=y.shape
    if n<3 or len(set(a))!=n:raise ValueError('invalid anchors')
    if K.ndim==2:K=np.broadcast_to(K,(g,)+K.shape)
    if K.shape[0]!=g:raise ValueError('outcome count mismatch')
    k=K[:,a,:][:,:,a];km=k.mean(2);tm=k.mean((1,2))
    kc=k-km[:,:,None]-km[:,None,:]+tm[:,None,None]
    ev,U=np.linalg.eigh((kc+kc.transpose(0,2,1))/2);ev=np.maximum(ev,0)
    ym=y.mean(0);yc=(y-ym).T
    uy=np.einsum('gij,gi->gj',U,yc)
    coefs=np.einsum('gij,gj,lgj->lgi',U,uy,1/(ev[None,:,:]+np.array(lambdas)[:,None,None]))
    cross=K[:,:,a];cross=cross-cross.mean(2)[:,:,None]-km[:,None,:]+tm[:,None,None]
    pred=ym+np.einsum('gni,lgi->lng',cross,coefs)
    h=1/n+np.einsum('gij,lgj->lgi',U*U,ev[None,:,:]/(ev[None,:,:]+np.array(lambdas)[:,None,None]))
    # preserve SIGNED leave-condition-out errors for stacking.
    err=(y[None,:,:]-pred[:,a,:])/np.maximum(1-h.transpose(0,2,1),1e-12)
    return pred,err

def bank(source,cues,k,anchors,y_anchor,families=('cue','joint','output_disagreement','output_marginal')):
    b=source.mean(0);r=y_anchor-b[anchors]
    specs=[]
    if 'cue' in families:specs.append(('cue',k['cue_original']))
    if 'joint' in families:
        specs.extend((f'joint_{w}',(1-w)*k['source_rbf']+w*k['cue']) for w in MIXES)
    for family,kind in [('output_disagreement','disagreement'),('output_marginal','marginal')]:
        if family in families:
            kg=output_kernels(source,kind)
            specs.extend((f'{family}_{w}',(1-w)*k['cue'][None,:,:]+w*kg) for w in MIXES[1:])
    pp=[];ee=[];params=[]
    for name,K in specs:
        p,e=path(K,anchors,r)
        pp.extend(p);ee.extend(e)
        params.extend({'kernel':name,'lambda':float(l)} for l in LAMBDAS)
    return b+np.asarray(pp),np.asarray(ee),params

def select(pp,ee,params):
    scores=np.mean(ee*ee,axis=(1,2));m=scores.min()
    ids=np.flatnonzero(scores<=m+1e-10*(1+abs(m)))
    idx=sorted(ids,key=lambda j:(-params[j]['lambda'],j))[0]
    return Fit(pp[idx],{**params[idx],'anchor_loo_mse':float(scores[idx])},ee[idx]**2)

def aggregate(pp,ee,params,mode='bootstrap',seed=801,replicates=256):
    m,n,g=ee.shape;loss=np.mean(ee*ee,axis=2)
    if mode=='bootstrap':
        rng=np.random.default_rng(seed)
        # Bayesian bootstrap conditions, keep all proteins from one condition together.
        boot=rng.dirichlet(np.ones(n),replicates)
        scores=boot@loss.T
        # When CV cannot distinguish penalties, use the declared deterministic
        # strong-penalty prior, not round-off. This matches hard selection.
        priority=np.array(sorted(range(m),key=lambda j:(-params[j]['lambda'],j)))
        minimum=scores.min(1)
        eligible=scores[:,priority] <= minimum[:,None]+1e-10*(1+np.abs(minimum[:,None]))
        winners=priority[np.argmax(eligible,axis=1)]
        weights=np.bincount(winners,minlength=m)/replicates
    elif mode=='stacking':
        E=ee.reshape(m,-1);C=E@E.T/(n*g)
        penalty=.01*max(float(np.trace(C)/m),1e-12)
        prior=np.ones(m)/m
        f=lambda w:float(w@C@w+penalty*np.sum((w-prior)**2))
        j=lambda w:2*(C@w+penalty*(w-prior))
        res=minimize(f,prior,jac=j,bounds=[(0.,1.)]*m,
                     constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(m)},
                     method='SLSQP',options={'ftol':1e-10,'maxiter':500})
        if not res.success:raise RuntimeError('Stacking optimization failed: '+res.message)
        weights=np.maximum(res.x,0);weights/=weights.sum()
    elif mode=='uniform':weights=np.ones(m)/m
    else:raise ValueError(mode)
    pred=np.einsum('m,mng->ng',weights,pp)
    err=np.einsum('m,mng->ng',weights,ee)
    return Fit(pred,{'method':mode,'weights':weights.tolist(),'components':params,
                    'anchor_oof_training_mse':float(np.mean(err*err)),
                    'note':'meta-training fit, NOT unbiased estimate of stack performance'},err*err)

def fit_stage2(source,cues,a,y,k):
    out={}
    for fams,name in [(('cue','joint'),'old_bank'),
                      (('cue','output_disagreement'),'output_disagreement'),
                      (('cue','output_marginal'),'output_marginal'),
                      (('cue','joint','output_disagreement','output_marginal'),'expanded')]:
        pp,ee,pars=bank(source,cues,k,a,y,fams)
        out[name+'_select']=select(pp,ee,pars)
        if name in ['old_bank','expanded']:
            for mode in ['bootstrap','stacking']:
                out[name+'_'+mode]=aggregate(pp,ee,pars,mode)
    return out
