"""Two-pass readout-blocked equivalent of the expanded model bank.

First pass accumulates complete-condition losses and the candidate-error Gram
matrix. Second pass writes only the combined prediction. This avoids retaining
all candidate x condition x gene predictions at once. Inputs can be NumPy arrays
or compatible read-only memory maps; it does not fetch/preprocess RNA data.

The global source/cue kernels must be built from ALL source readouts, not rebuilt
inside each block. Changing that would change the statistical model.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from .aggregation import bank
from .covariance import Fit,build_kernels


def _weights(loss,C,params,mode,seed=801,replicates=256):
    m,n=loss.shape
    priority=np.array(sorted(range(m),key=lambda j:(-params[j]['lambda'],j)))
    if mode=='bootstrap':
        boot=np.random.default_rng(seed).dirichlet(np.ones(n),replicates)
        score=boot@loss.T;v=score.min(1)
        eligible=score[:,priority]<=v[:,None]+1e-10*(1+np.abs(v[:,None]))
        winner=priority[np.argmax(eligible,axis=1)]
        return np.bincount(winner,minlength=m)/replicates
    if mode=='stacking':
        penalty=.01*max(float(np.trace(C)/m),1e-12);prior=np.ones(m)/m
        res=minimize(lambda w:float(w@C@w+penalty*np.sum((w-prior)**2)),prior,
            jac=lambda w:2*(C@w+penalty*(w-prior)),bounds=[(0.,1.)]*m,
            constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(m)},
            method='SLSQP',options={'ftol':1e-10,'maxiter':500})
        if not res.success:raise RuntimeError('Stacking optimization failed: '+res.message)
        w=np.maximum(res.x,0);return w/w.sum()
    if mode=='uniform':return np.ones(m)/m
    raise ValueError('mode must be stacking, bootstrap or uniform')


def fit_streaming(source,cues,anchors,y_anchor,*,block_size=128,mode='stacking',kernels=None,out=None):
    """Fit without target-query labels, in two passes over output blocks.

    `out` may be a writable NumPy memory map of shape (n_conditions,n_readouts).
    Working kernels still scale quadratically with the condition count; this is
    readout-blocking, not a claimed solution to million-condition scaling.
    """
    s=np.asarray(source);x=np.asarray(cues);a=np.asarray(anchors,int);y=np.asarray(y_anchor)
    if s.ndim!=3 or y.shape!=(len(a),s.shape[2]) or not isinstance(block_size,int) or block_size<1:
        raise ValueError('Aligned source/anchor responses and positive integer block size required')
    if not np.isfinite(s).all() or not np.isfinite(y).all():raise ValueError('finite data required')
    G=s.shape[2]
    if not G:raise ValueError('at least one readout required')
    k=build_kernels(s,x) if kernels is None else kernels
    C=None;loss=None;params=None
    for first in range(0,G,block_size):
        sl=slice(first,min(first+block_size,G))
        pp,e,pars=bank(s[:,:,sl],x,k,a,y[:,sl])
        if params is None:
            params=pars;M=len(pars);C=np.zeros((M,M));loss=np.zeros((M,len(a)))
        elif params!=pars:raise RuntimeError('Model bank changed across readout blocks')
        E=e.reshape(len(pars),-1);C+=E@E.T;loss+=np.sum(e*e,axis=2)
        del pp,e,E
    C/=len(a)*G;loss/=G;w=_weights(loss,C,params,mode)
    if out is None:out=np.empty((s.shape[1],G),float)
    if out.shape!=(s.shape[1],G) or not np.issubdtype(out.dtype,np.floating):
        raise ValueError('Output must have matching shape and floating dtype')
    loo=np.empty(y.shape,float)
    for first in range(0,G,block_size):
        sl=slice(first,min(first+block_size,G))
        pp,e,pars=bank(s[:,:,sl],x,k,a,y[:,sl])
        out[:,sl]=np.einsum('m,mng->ng',w,pp)
        loo[:,sl]=np.einsum('m,mng->ng',w,e)**2
    return Fit(out,{'method':mode,'weights':w.tolist(),'components':params,'block_size':block_size,
       'anchor_oof_training_mse':float(w@C@w),'note':'OOF meta-training loss, NOT independent performance; mathematically equivalent blocked computation.'},loo)
