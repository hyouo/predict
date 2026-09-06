"""Registered low-dimensional response transfer baselines for OP3.

All learning is train-only. Per-gene slopes are ordinary shrinkage regression;
KRR is classical and no mechanism/novelty claim is attached to these estimators.
"""
from __future__ import annotations
import hashlib
import numpy as np
from scipy.optimize import minimize
from .op3_data import SOURCES, TARGETS

LAMBDAS = (.01, .1, 1., 10., 100.)
GATE_LAMBDAS = (0., .01, .1, 1., 10., 100.)


def _check(x,y):
    if x.ndim!=3 or y.shape!=x.shape[:2] or x.shape[0]<3:
        raise ValueError('Need at least 3 anchors, aligned gene/source axes')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Missing values cannot be treated as zeros')


def source_kernels(source):
    if source.ndim!=3 or not np.isfinite(source).all():raise ValueError('Invalid sources')
    x=source.transpose(1,0,2).reshape(source.shape[1],-1)
    x=x-x.mean(0)
    raw=x@x.T/x.shape[1]
    scale=max(float(np.trace(raw)/len(raw)),1e-12)
    linear=raw/scale
    d=np.maximum(np.diag(raw)[:,None]+np.diag(raw)[None,:]-2*raw,0.)
    bandwidth=max(float(np.median(d[np.triu_indices(len(d),1)])),1e-12)
    rbf=np.exp(-d/bandwidth)
    return {'linear':linear, 'rbf':rbf, 'joint':.5*linear+.5*rbf}


def krr_path(kernel, anchors, y, prior, lambdas=LAMBDAS):
    """Exact unpenalized-intercept KRR and analytic leave-drug-out predictions."""
    a=np.asarray(anchors,int);y=np.asarray(y,float);prior=np.asarray(prior,float)
    n=len(a)
    if n<3 or len(set(a))!=n:raise ValueError('Invalid anchor identities')
    if y.shape!=(n,prior.shape[1]) or kernel.shape!=(len(prior),len(prior)):
        raise ValueError('Kernel/outcome mismatch')
    pp=[];ll=[]
    residual=y-prior[a]
    for lam in lambdas:
        inv=np.linalg.inv(kernel[np.ix_(a,a)]+lam*np.eye(n))
        v=inv.sum(1);total=v.sum()
        c=inv-np.outer(v,v)/total
        alpha=c@residual
        intercept=v@residual/total
        pp.append(prior+intercept+kernel[:,a]@alpha)
        ll.append(y-alpha/np.diag(c)[:,None])
    return np.array(pp),np.array(ll)


def slope_fit(x,y,query,lam=1.,kind='gene'):
    """Per-gene 2-source slope with a pooled-slope shrinkage center and gene intercept."""
    _check(x,y)
    xm=x.mean(0);ym=y.mean(0);xc=x-xm;yc=y-ym
    cov=np.einsum('ngs,ngt->gst',xc,xc)
    rhs=np.einsum('ngs,ng->gs',xc,yc)
    glob=cov.sum(0);grhs=rhs.sum(0)
    # Stable common slope; fixed tiny ridge is not tuned on public outcomes.
    penalty=1e-3*max(float(np.trace(glob)/glob.shape[0]),1e-12)
    shared=np.linalg.solve(glob+penalty*np.eye(glob.shape[0]),grhs)
    if kind=='shared':beta=np.broadcast_to(shared,rhs.shape)
    elif kind=='gene':
        scale=np.maximum(np.trace(cov,axis1=1,axis2=2)/cov.shape[-1],1e-8)
        reg=lam*scale
        beta=np.linalg.solve(cov+reg[:,None,None]*np.eye(cov.shape[-1]),
                             (rhs+reg[:,None]*shared)[...,None])[...,0]
    else:raise ValueError(kind)
    return ym+np.einsum('ngs,gs->ng',query-xm,beta)


def slope_path(source,a,y,lambdas=LAMBDAS):
    x=source.transpose(1,2,0);n=len(a)
    pp=[];loo=[];spec=[]
    for kind,lam in [('shared',1.)]+[('gene',l) for l in lambdas]:
        pp.append(slope_fit(x[a],y,x,lam,kind))
        lo=[]
        for j in range(n):
            ix=np.arange(n)!=j
            lo.append(slope_fit(x[a[ix]],y[ix],x[a[j:j+1]],lam,kind)[0])
        loo.append(lo);spec.append({'family':kind+'_slope','lambda':lam})
    return np.array(pp),np.array(loo),spec


def stable_best(loss,spec):
    scores=np.mean(loss*loss,axis=(1,2));minimum=scores.min()
    eligible=np.flatnonzero(scores<=minimum+1e-10*(1+abs(minimum)))
    return sorted(eligible,key=lambda j:(-spec[j].get('lambda',0),j))[0]


def stack_weights(errors):
    m,n,g=errors.shape;e=errors.reshape(m,-1);c=e@e.T/(n*g)
    pen=.01*max(float(np.trace(c)/m),1e-12);w0=np.ones(m)/m
    fit=minimize(lambda w:float(w@c@w+pen*np.sum((w-w0)**2)),w0,
                 jac=lambda w:2*(c@w+pen*(w-w0)),bounds=[(0.,1.)]*m,
                 constraints={'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(m)},
                 method='SLSQP',options={'ftol':1e-10,'maxiter':500})
    if not fit.success:raise RuntimeError(fit.message)
    w=np.maximum(fit.x,0);return w/w.sum()


def fewshot(source,a,y_fit,y_validation):
    """Only target anchors are accepted; no final reference or query outcome argument."""
    if y_fit.shape!=y_validation.shape or len(a)!=len(y_fit):raise ValueError('Anchor mismatch')
    base=source.mean(0);zero=np.zeros_like(base);k=source_kernels(source)
    pp=[base,base+(y_fit-base[a]).mean(0)]
    lo=[base[a],base[a]+((y_fit-base[a]).sum(0)-(y_fit-base[a]))/(len(a)-1)]
    specs=[{'family':'source_mean','lambda':0.},{'family':'source_offset','lambda':0.}]
    sp,sl,ss=slope_path(source,a,y_fit);pp.extend(sp);lo.extend(sl);specs.extend(ss)
    for style,prior in [('direct',zero),('residual',base)]:
        for name,kernel in k.items():
            p,l=krr_path(kernel,a,y_fit,prior)
            pp.extend(p);lo.extend(l)
            specs.extend({'family':style+'_krr','kernel':name,'lambda':v} for v in LAMBDAS)
    pp=np.array(pp);lo=np.array(lo)
    out={'source_offset':pp[1], 'shared_slopes':pp[2]};pars={}
    for mode,truth in [('same',y_fit),('separated',y_validation)]:
        error=truth[None]-lo
        for fam in ('gene_slope','direct_krr','residual_krr'):
            ids=[j for j,s in enumerate(specs) if s['family']==fam]
            best=ids[stable_best(error[ids],[specs[j] for j in ids])]
            name=fam+'_'+mode;out[name]=pp[best]
            pars[name]={**specs[best],'validation_mse':float(np.mean(error[best]**2))}
        w=stack_weights(error)
        name='stack_'+mode;out[name]=np.einsum('m,mng->ng',w,pp)
        pars[name]={'weights':w.tolist(),'components':specs,
                    'validation_mse':float(np.mean(np.einsum('m,mng->ng',w,error)**2)),
                    'warning':'meta-training score, not unbiased performance estimate'}
    return out,pars


def _gate_features(effect, source_base, target_base, degree):
    diff=np.clip(target_base-source_base,-4.,4.)/4.
    terms=[np.ones_like(diff)] + ([diff,diff*diff] if degree==2 else [])
    return effect[:,:,None]*np.stack(terms,axis=1)[None,:,:]


def _gate_fit(features, truth, lam):
    f=features.reshape(-1,features.shape[-1]);y=truth.ravel()
    cov=f.T@f;rhs=f.T@y;prior=np.zeros(cov.shape[0]);prior[0]=1.
    scale=max(float(np.trace(cov)/len(prior)),1e-12)
    return np.linalg.solve(cov+(lam*scale+1e-10)*np.eye(len(prior)),rhs+lam*scale*prior)


def zero_shot(bundle):
    """Reciprocal source tasks only; no target treatment values are accessed."""
    s=bundle.source;sv=bundle.source_validation
    valid=np.ones(len(bundle.compounds),bool)
    for row in bundle.missing_source:valid[bundle.compounds.index(row['compound'])]=False
    folds=np.array([int(hashlib.sha256(p.encode()).hexdigest()[:8],16)%5 for p in bundle.compounds])
    perm=np.random.default_rng(321).permutation(s.shape[-1])
    outcomes={};params={}
    for degree,shuffle,name in [(0,False,'zs_calibrated'),(2,False,'zs_baseline_gate'),
                                 (2,True,'zs_shuffled_gate')]:
        fx=[];ya=[];yb=[];fo=[]
        for a,b in [(0,1),(1,0)]:
            tb=bundle.baseline[SOURCES[b]]
            if shuffle:tb=tb[perm]
            fx.append(_gate_features(s[a,valid],bundle.baseline[SOURCES[a]],tb,degree))
            ya.append(s[b,valid]);yb.append(sv[b,valid]);fo.append(folds[valid])
        fx=np.concatenate(fx);ya=np.concatenate(ya);yb=np.concatenate(yb);fo=np.concatenate(fo)
        losses=[]
        for lam in GATE_LAMBDAS:
            err=[]
            for fold in range(5):
                train=fo!=fold;val=~train
                coef=_gate_fit(fx[train],ya[train],lam)
                err.append(np.mean((np.einsum('ngk,k->ng',fx[val],coef)-yb[val])**2))
            losses.append(float(np.mean(err)))
        chosen=min(range(len(losses)),key=lambda j:(round(losses[j],10),-GATE_LAMBDAS[j]))
        coef=_gate_fit(fx,ya,GATE_LAMBDAS[chosen])
        params[name]={'coefficients':coef.tolist(),'lambda':GATE_LAMBDAS[chosen],
                      'source_validation_losses':losses,'degree':degree,'shuffle':shuffle}
        for ct in TARGETS:
            tb=bundle.baseline[ct]
            if shuffle:tb=tb[perm]
            pred=np.stack([np.einsum('ngk,k->ng',_gate_features(s[a],bundle.baseline[SOURCES[a]],tb,degree),coef)
                           for a in range(2)]).mean(0)
            outcomes[(ct,name)]=pred
    for ct in TARGETS:
        outcomes[(ct,'zero')]=np.zeros_like(s[0]);outcomes[(ct,'source_mean')]=s.mean(0)
    return outcomes,params
