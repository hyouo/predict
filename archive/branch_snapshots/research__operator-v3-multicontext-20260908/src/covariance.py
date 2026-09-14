"""Source-only kernels and residual KRR. No biological/SOTA guarantee.

Target outcomes can enter only as y_anchor. Catalog-level source and intervention
features are allowed, including source measurements of queried interventions.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

LAMBDAS = np.array([0.001,0.01,0.1,1.,10.,100.,1000.])
MIXES = (0.,.25,.5,.75,1.)

@dataclass
class Fit:
    prediction: np.ndarray
    parameters: dict
    loo_losses: np.ndarray

def normalize_kernel(k):
    k=np.asarray(k,dtype=float)
    if k.ndim!=2 or k.shape[0]!=k.shape[1] or not np.isfinite(k).all():
        raise ValueError('finite square kernel required')
    # Uniform catalog centering is source/descriptor-only, NOT target-outcome based.
    c=k-k.mean(0)[None,:]-k.mean(1)[:,None]+k.mean()
    tr=float(np.mean(np.diag(c)))
    if tr<=1e-12: return np.zeros_like(c)
    return (c+c.T)/(2*tr)

def rbf(x):
    x=np.asarray(x,dtype=float)
    x=x-x.mean(0)
    dd=np.maximum(0.,np.sum(x*x,1)[:,None]+np.sum(x*x,1)[None,:]-2*x@x.T)
    pos=dd[np.triu_indices(len(x),1)];pos=pos[pos>1e-12]
    bw=np.median(pos) if len(pos) else 1.
    return np.exp(-dd/(2*bw))

def cue_kernel(x):
    x=np.asarray(x,float)
    pairs=np.stack([x[:,i]*x[:,j] for i in range(x.shape[1]) for j in range(i+1,x.shape[1])],axis=1)
    return rbf(np.concatenate([x,pairs],axis=1))

def build_kernels(source,cues):
    s=np.asarray(source,float);x=np.asarray(cues,float)
    if s.ndim!=3 or s.shape[0]<2 or len(x)!=s.shape[1] or not np.isfinite(s).all():
        raise ValueError('Need >=2 source contexts with finite aligned outcomes')
    n=s.shape[1]
    d=s-s.mean(0,keepdims=True)
    d=d-d.mean(1,keepdims=True)  # remove each source-output offset
    z=d.transpose(1,0,2).reshape(n,-1)
    f=s.transpose(1,0,2).reshape(n,-1)
    # Groupwise ANOVA kernel; all administered cues, no inferred exact drug targets.
    kl=normalize_kernel(rbf(x[:,:5]));ki=normalize_kernel(rbf(x[:,5:]))
    return {'cue_original':cue_kernel(x),
            'cue':normalize_kernel(cue_kernel(x)),
            'source_rbf':normalize_kernel(rbf(f)),
            'disagreement':normalize_kernel(z@z.T),
            'marginal_linear':normalize_kernel(f@f.T),
            'anova_additive':normalize_kernel(kl+ki),
            'anova_interaction':normalize_kernel(kl*ki)}

def kernel_path(K, anchors, residual, lambdas=LAMBDAS):
    """Exact full-condition LOO with free output intercepts. Lambda not scaled by n."""
    K=np.asarray(K,float);a=np.asarray(anchors,int);r=np.asarray(residual,float)
    ll=np.asarray(lambdas,float);n=len(a)
    if n<3 or len(set(a))!=n or r.ndim!=2 or r.shape[0]!=n:
        raise ValueError('Need >=3 distinct anchors and a 2D outcome matrix')
    if np.any(a<0) or np.any(a>=len(K)) or not np.isfinite(r).all() or np.any(ll<=0):
        raise ValueError('invalid indices/data/penalty')
    k=K[np.ix_(a,a)];km=k.mean(1);tm=k.mean()
    kc=k-km[:,None]-km[None,:]+tm
    ev,U=np.linalg.eigh((kc+kc.T)/2);ev=np.maximum(ev,0.)
    ym=r.mean(0);yc=r-ym
    inv=1/(ev[None,:]+ll[:,None])
    coef=np.einsum('ij,lj,jg->lig',U,inv,U.T@yc)
    cross=K[:,a]-K[:,a].mean(1)[:,None]-km[None,:]+tm
    pred=ym+np.einsum('ni,lig->lng',cross,coef)
    leverage=1/n+np.einsum('ij,lj->li',U*U,ev[None,:]*inv)
    error=(r[None,:,:]-pred[:,a,:])/np.maximum(1-leverage[:,:,None],1e-12)
    return pred,error*error

def tune(kernels,anchors,residual,lambdas=LAMBDAS,per_output=False):
    """Target-anchor-only choice; stable numeric tie rule, no target queries accepted."""
    preds=[];losses=[];pars=[]
    for name,K in kernels:
        pp,ee=kernel_path(K,anchors,residual,lambdas)
        for j,lam in enumerate(lambdas):
            preds.append(pp[j]);losses.append(ee[j].mean(0))
            pars.append({'kernel':name,'lambda':float(lam)})
    pp=np.asarray(preds);ee=np.asarray(losses)
    if per_output:
        minima=ee.min(0);elig=ee<=minima[None,:]+1e-10*(1+np.abs(minima[None,:]))
        selected=[]
        for g in range(ee.shape[1]):
            ids=np.flatnonzero(elig[:,g])
            selected.append(sorted(ids,key=lambda j:(-pars[j]['lambda'],j))[0])
        selected=np.array(selected)
        p=np.stack([pp[j,:,g] for g,j in enumerate(selected)],axis=1)
        params={'per_output':True,'selected':[pars[j] for j in selected],
                'anchor_loo_mse':float(ee[selected,np.arange(ee.shape[1])].mean())}
    else:
        risk=ee.mean(1);best=risk.min()
        ids=np.flatnonzero(risk<=best+1e-10*(1+abs(best)))
        j=sorted(ids,key=lambda j:(-pars[j]['lambda'],j))[0]
        p=pp[j];params={**pars[j],'per_output':False,'anchor_loo_mse':float(risk[j])}
    return Fit(p,params,ee)

def candidates(k, family):
    if family=='cue':return [('cue_original',k['cue_original'])]
    if family=='joint':return [(f'joint_{w}',(1-w)*k['source_rbf']+w*k['cue']) for w in MIXES]
    if family=='disagreement':return [(f'disagreement_{w}',(1-w)*k['cue']+w*k['disagreement']) for w in MIXES]
    if family=='marginal':return [(f'marginal_{w}',(1-w)*k['cue']+w*k['marginal_linear']) for w in MIXES]
    if family=='anova':return [(f'anova_{w}',(1-w)*k['anova_additive']+w*k['anova_interaction']) for w in MIXES]
    raise ValueError(family)

def fit_all(source,cues,anchors,y_anchor,kernels=None):
    source=np.asarray(source,float);a=np.asarray(anchors,int);y=np.asarray(y_anchor,float)
    b=source.mean(0);r=y-b[a]
    k=build_kernels(source,cues) if kernels is None else kernels
    out={'source_mean':Fit(b.copy(),{},np.empty(0)),
         'source_offset':Fit(b+r.mean(0),{},np.empty(0))}
    for family in ['cue','joint','disagreement','marginal','anova']:
        z=tune(candidates(k,family),a,r)
        out[family]=Fit(b+z.prediction,z.parameters,z.loo_losses)
    for family in ['cue','disagreement']:
        z=tune(candidates(k,family),a,r,per_output=True)
        out[family+'_per_output']=Fit(b+z.prediction,z.parameters,z.loo_losses)
    return out
