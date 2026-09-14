"""GLS/GP residual fitting with an explicit, estimated control-error covariance.

This follows exploratory protocol 005. It does not modify the frozen 71-candidate
RNA implementation. Conditional Gaussian assumptions are not biological facts.
"""
from __future__ import annotations
import numpy as np
from src.rna_models import (candidate_bank, shared_kernels, kernel_blocks, maps_for,
                           prior_block, fit_weights, _finite_shape)
from src.reference_risk import reference_correction


def correlated_smoother(kernel, anchors, regularization, intercept, noise):
    """Exact latent-response LOO and final maps, including a flat GLS intercept.

    Measurement-error cross-covariance is NOT used as latent query signal.
    The ordinary I-Q/diag(Q) shortcut is invalid for non-scalar noise.
    """
    k=np.asarray(kernel,dtype=float);ix=np.asarray(anchors,dtype=int)
    v=np.asarray(noise,dtype=float);n=len(ix)
    if (k.ndim!=3 or k.shape[1]!=k.shape[2] or v.shape!=(len(k),n,n)
        or n<2 or len(set(ix.tolist()))!=n or np.any(ix<0) or np.any(ix>=k.shape[1])):
        raise ValueError('Incompatible kernel, anchor or covariance shape')
    if regularization<=0 or not np.isfinite(regularization):raise ValueError('Invalid nugget')
    if not np.isfinite(k).all() or not np.isfinite(v).all():raise ValueError('Nonfinite covariance')
    if not np.allclose(v,v.transpose(0,2,1),atol=1e-10,rtol=1e-10):raise ValueError('Noise must be symmetric')
    if np.min(np.linalg.eigvalsh(v)) < -1e-8:raise ValueError('Noise must be positive semidefinite')
    if intercept not in ('none','shrink','free'):raise ValueError('Invalid intercept mode')
    if intercept=='shrink':k=k+1.
    ka=k[:,:,ix];kaa=ka[:,ix,:]
    matrix=kaa+v+regularization*np.eye(n)
    b=np.linalg.solve(matrix,np.broadcast_to(np.eye(n),matrix.shape))
    if intercept=='free':
        z=b.sum(axis=2);den=z.sum(axis=1)
        if np.any(den<=0):raise ValueError('Invalid GLS intercept precision')
        q=b-z[:,:,None]*z[:,None,:]/den[:,None,None]
        final=ka@q+(z/den[:,None])[:,None,:]
    else:q=b;final=ka@q
    s=final[:,ix,:];diagq=np.diagonal(q,axis1=1,axis2=2)
    if np.any(diagq<=1e-14):raise ValueError('Degenerate leave-one-out precision')
    h=s-(np.diagonal(s,axis1=1,axis2=2)/diagq)[:,:,None]*q
    h[:,np.arange(n),np.arange(n)]=0.
    return h,final


def fit_correlated_bank(source, anchors, reference_c, reference_d, factors,
                        *, noise_mode='full',block_size=128):
    source,yc,yd,ix=_finite_shape(source,reference_c,reference_d,anchors)
    f=np.asarray(factors,dtype=float)
    if noise_mode not in ('full','diagonal'):raise ValueError('Invalid covariance mode')
    if (f.ndim!=3 or f.shape[1:]!=yc.shape or len(f)<1 or not np.isfinite(f).all()
        or not np.allclose(f.sum(axis=0),yc-yd,rtol=1e-8,atol=1e-9)):
        raise ValueError('Invalid or inconsistent reference factors')
    if isinstance(block_size,bool) or not isinstance(block_size,int) or block_size<1:raise ValueError('Invalid block size')
    bank=candidate_bank();n,ng=yc.shape;nc=source.shape[1];m=len(bank)
    shared=shared_kernels(source);mean=source.mean(axis=0);y=(yc+yd)/2
    gram=np.zeros((m,m));corr=np.zeros(m)
    def prepare(start,stop):
        kernels=kernel_blocks(source,start,stop,shared)
        # Under exchangeable independent C/D references, mean-reference
        # measurement covariance is one quarter of Var(C-D).
        block=f[:,:,start:stop]
        noise=.25*np.einsum('uig,ujg->gij',block,block,optimize=True)
        if noise_mode=='diagonal':noise=np.einsum('gi,ij->gij',np.diagonal(noise,axis1=1,axis2=2),np.eye(n))
        return kernels,noise,block
    def mapping(can,kernels,noise,b):
        if can.kernel in ('fixed','offset'):return maps_for(can,kernels,ix,b,nc)
        return correlated_smoother(kernels[can.kernel],ix,can.regularization,can.intercept,noise)
    for start in range(0,ng,block_size):
        stop=min(start+block_size,ng);b=stop-start;kernels,noise,fb=prepare(start,stop)
        err=np.empty((m,n,b))
        for j,can in enumerate(bank):
            h,_=mapping(can,kernels,noise,b);pr=prior_block(can,mean,start,stop)
            resid=y[:,start:stop]-pr[ix]
            err[j]=np.einsum('gij,jg->ig',h,resid)-resid
            corr[j]+=reference_correction(h,fb)
        flat=err.reshape(m,-1);gram+=flat@flat.T
    gram/=n*ng;corr/=n*ng
    choices={};optim={}
    for group,mask in [('simple',np.array([not c.kernel.startswith('mixed') for c in bank])),
                       ('mixed',np.ones(m,dtype=bool))]:
        ids=np.flatnonzero(mask)
        for adjusted in (False,True):
            label=group+'_stack_'+('adjusted' if adjusted else 'cv')
            w,rec=fit_weights(gram[np.ix_(ids,ids)],corr[ids],adjusted=adjusted)
            full=np.zeros(m);full[ids]=w;choices[label]=full;optim[label]=rec
    preds={k:np.zeros((nc,ng)) for k in choices}
    for start in range(0,ng,block_size):
        stop=min(start+block_size,ng);b=stop-start;kernels,noise,_=prepare(start,stop)
        for j,can in enumerate(bank):
            needed=[k for k,w in choices.items() if w[j]>1e-12]
            if not needed:continue
            _,final=mapping(can,kernels,noise,b);pr=prior_block(can,mean,start,stop)
            pred=pr+np.einsum('gij,jg->ig',final,y[:,start:stop]-pr[ix])
            for label in needed:preds[label][:,start:stop]+=choices[label][j]*pred
    return preds,{'noise_mode':noise_mode,'candidate_names':[c.name for c in bank],
                  'ordinary_cv':np.diag(gram).tolist(),'adjusted_cv':(np.diag(gram)+.5*corr).tolist(),
                  'reference_correction_before_half':corr.tolist(),'gram':gram.tolist(),
                  'weights':{k:v.tolist() for k,v in choices.items()},'optimization':optim,
                  'target_test_outcomes_received':False,'reference_units':len(f),
                  'scientific_status':'exploratory; covariance is an estimated working model, not validated biological uncertainty'}
