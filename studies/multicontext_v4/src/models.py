"""Shared response and support-centred context interactions.

Classical reduced-rank/ridge ingredients; no claim of a new universal estimator.
Only supplied training effects are consumed. NaN denotes an unobserved condition,
never a biological zero. Gene axis is common and finite for observed rows.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import solve


def validate(y, baseline):
    y=np.asarray(y,dtype=float);b=np.asarray(baseline,dtype=float)
    if y.ndim!=3 or b.shape!=(y.shape[0],y.shape[-1]):
        raise ValueError('Expected contexts x compounds x genes and contexts x genes')
    if np.isinf(y).any():raise ValueError('Infinite values are not missing profiles')
    m=np.isfinite(y).all(-1)
    if not np.array_equal(m,np.isfinite(y).any(-1)):
        raise ValueError('A profile must be fully observed or fully missing')
    if not np.isfinite(b).all() or not m.any():raise ValueError('Invalid baselines or no observed profiles')
    return y,b,m


def pca_fit(x,rank,scale=False):
    x=np.asarray(x,float);mean=x.mean(0);xc=x-mean
    sd=np.sqrt(np.mean(xc**2,0)) if scale else np.ones(x.shape[1])
    # bounded scaling; absent features do not alter positive-variance scale
    positive=sd[sd>1e-10];floor=max(float(np.median(positive))*.1,1e-8) if len(positive) else 1.
    sd=np.maximum(sd,floor)
    xc=xc/sd
    u,s,v=np.linalg.svd(xc,full_matrices=False)
    r=min(rank,int(np.sum(s>max(float(s[0])*1e-9,1e-10))))
    if r<1:return mean,sd,np.zeros((x.shape[1],1)),1.
    basis=v[:r].T
    norm=max(float(np.sqrt(np.mean((xc@basis)**2))),1e-8)
    return mean,sd,basis,norm


def project(x,pca):
    mu,sd,v,norm=pca;return ((np.asarray(x,float)-mu)/sd)@v/norm


def consensus(y):
    finite=np.isfinite(y).all(-1);count=finite.sum(0)
    out=np.nansum(y,axis=0)/np.maximum(count[:,None],1)
    out[count==0]=np.nan
    return out,count,finite


def ridge(x,y,penalty,weights=None):
    if penalty<=0:raise ValueError('Positive ridge required')
    weights=np.ones(len(x)) if weights is None else np.asarray(weights,float)
    a=x.T@(weights[:,None]*x)
    positive=np.diag(a);positive=positive[positive>1e-12]
    s=float(np.mean(positive)) if len(positive) else 1.
    return solve(a+penalty*s*np.eye(a.shape[0]),x.T@(weights[:,None]*y),assume_a='pos')


@dataclass
class Predictor:
    kind: str = 'interaction_centered'
    penalty: float = 1.
    strength: float = 1.
    rank_drug: int = 16
    rank_context: int = 6

    def fit(self,y,baseline):
        y,b,m=validate(y,baseline)
        self.y=y.copy();self.b=b.copy();self.mask=m
        self.f,self.count,_=consensus(y)
        self.zpca=pca_fit(b,self.rank_context,True)
        self.z=project(b,self.zpca)
        valid=self.count>=2
        self.upca=pca_fit(self.f[valid],self.rank_drug,False)
        self.u=project(np.nan_to_num(self.f),self.upca)
        # source means for a drug respect actual observed source support
        self.zbar=np.einsum('cp,ck->pk',m.astype(float),self.z)/np.maximum(self.count[:,None],1)
        if self.kind in ('zero','source_mean','context_rbf','context_linear'):
            return self
        # Low-dimensional common response correction trained with leave-context-out input.
        xx=[];tt=[];ww=[]
        for c in range(len(y)):
            ok=m[c]&valid
            fm=(self.f[ok]*self.count[ok,None]-y[c,ok])/(self.count[ok,None]-1)
            u=project(fm,self.upca)
            xx.append(np.c_[np.ones(len(u)),u]);tt.append(y[c,ok]);ww.append(np.ones(len(u))/max(len(u),1))
        x=np.concatenate(xx);t=np.concatenate(tt);w=np.concatenate(ww)
        if self.kind=='scalar':
            fs=[]
            for c in range(len(y)):
                ok=m[c]&valid;fs.append((self.f[ok]*self.count[ok,None]-y[c,ok])/(self.count[ok,None]-1))
            fs=np.concatenate(fs)
            self.alpha=float(np.clip(np.sum(w[:,None]*fs*t)/max(float(np.sum(w[:,None]*fs*fs)),1e-12),0,2))
            return self
        self.trunk=ridge(x,t,self.penalty,w)
        if self.kind=='shared_ridge':return self
        # Train context-specific deviations, with the exact available-context mean removed.
        # Outcome moments here are from training contexts only, not unseen targets.
        fi=[];ti=[];wi=[]
        for c in range(len(y)):
            ok=m[c]&valid
            z=(self.z[c]-self.zbar[ok]) if self.kind!='interaction_uncentered' else np.broadcast_to(self.z[c],(ok.sum(),len(self.z[c])))
            if self.kind=='context_additive':
                a=z
            else:
                du=np.c_[np.ones(ok.sum()),self.u[ok]]
                a=np.einsum('nk,nr->nkr',z,du).reshape(len(z),-1)
            fi.append(a);ti.append(y[c,ok]-self.f[ok]);wi.append(np.ones(ok.sum())/max(int(ok.sum()),1))
        self.interaction=ridge(np.concatenate(fi),np.concatenate(ti),self.penalty,np.concatenate(wi))
        return self

    def predict(self,baseline):
        bt=np.asarray(baseline,float)
        if bt.ndim!=2 or bt.shape[1]!=self.b.shape[1] or not np.isfinite(bt).all():raise ValueError('Invalid target baseline')
        zt=project(bt,self.zpca);pred=[]
        if self.kind.startswith('context_') and self.kind!='context_additive':
            dist=np.mean((zt[:,None]-self.z[None,:])**2,axis=-1)
            if self.kind=='context_rbf':w=np.exp(-dist/(2*self.penalty))
            else:
                k=self.z@self.z.T/self.z.shape[1];kt=zt@self.z.T/self.z.shape[1]
                w=solve(k+self.penalty*np.eye(len(k)),kt.T,assume_a='pos').T
            for i in range(len(bt)):
                p=[]
                for drug in range(self.y.shape[1]):
                    mask=self.mask[:,drug];weights=w[i,mask]
                    if not mask.any():p.append(np.full(self.y.shape[-1],np.nan));continue
                    if self.kind=='context_rbf':
                        weights=weights/max(float(weights.sum()),1e-12)
                        p.append(self.strength*(weights@self.y[mask,drug]))
                    else:
                        km=self.z[mask]@self.z[mask].T/self.z.shape[1]
                        ktm=zt[i]@self.z[mask].T/self.z.shape[1]
                        weights=solve(km+self.penalty*np.eye(mask.sum()),ktm,assume_a='pos')
                        p.append(self.strength*(self.f[drug]+weights@(self.y[mask,drug]-self.f[drug])))
                pred.append(np.asarray(p))
            return np.asarray(pred)
        if self.kind=='zero':common=np.zeros_like(self.f)
        elif self.kind=='source_mean':common=self.f.copy()
        elif self.kind=='scalar':common=self.alpha*self.f
        else:common=np.c_[np.ones(len(self.u)),self.u]@self.trunk
        for zz in zt:
            out=common.copy()
            if self.kind in ('interaction_centered','interaction_uncentered','context_additive'):
                z=zz-self.zbar if self.kind!='interaction_uncentered' else np.broadcast_to(zz,self.zbar.shape)
                a=z if self.kind=='context_additive' else np.einsum('pk,pr->pkr',z,np.c_[np.ones(len(self.u)),self.u]).reshape(len(z),-1)
                out=out+self.strength*(a@self.interaction)
            out[self.count==0]=np.nan
            pred.append(out)
        return np.asarray(pred)
