"""Support-aware mean calibration and matched biological kernel priors.

All operations below consume permitted source effects/baselines only. The
constant-mean solver is classical ordinary kriging, not a new theorem.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import solve


def membership(gmt, symbols, min_size=5, max_size=150):
    """Deduplicate human pathway membership on the supplied measured gene axis."""
    symbols=list(map(str,symbols))
    if len(symbols)!=len(set(symbols)):raise ValueError('Duplicate gene symbols')
    lookup={s:i for i,s in enumerate(symbols)};cols=[];records=[];seen=set()
    with open(gmt) as stream:
        for line in stream:
            f=line.rstrip('\n').split('\t')
            if len(f)<3 or not f[1].startswith('R-HSA-'):continue
            ids=tuple(sorted({lookup[s] for s in f[2:] if s in lookup}))
            if not min_size<=len(ids)<=max_size or ids in seen:continue
            seen.add(ids);a=np.zeros(len(symbols));a[list(ids)]=1
            cols.append(a);records.append({'pathway':f[0],'stable_id':f[1],'landmark_size':len(ids)})
    if not cols:raise ValueError('No eligible real pathway annotations')
    return np.stack(cols,axis=1),records


def shuffled_membership(a,seed):
    """Exact degree-stratified relabelling preserves all set sizes and overlaps."""
    a=np.asarray(a);degree=a.sum(1);r=np.random.default_rng(seed);perm=np.arange(len(a))
    for d in np.unique(degree):
        ix=np.flatnonzero(degree==d);perm[ix]=r.permutation(ix)
    return a[perm],perm


def pathway_features(x,a):
    a=np.asarray(a,float);deg=np.maximum(a.sum(1),1);size=np.maximum(a.sum(0),1)
    return x@(a/np.sqrt(deg[:,None]*size[None,:]))


class Geometry:
    """Fit feature scales on source baselines; optional curated/prior-only kernel."""
    def __init__(self,mix=0.,prior=None):self.mix=mix;self.prior=prior
    def fit(self,b):
        b=np.asarray(b,float)
        if b.ndim!=2 or len(b)<2 or not np.isfinite(b).all():raise ValueError('Invalid source baselines')
        self.mean=b.mean(0);sd=np.sqrt(np.mean((b-self.mean)**2,0));pos=sd[sd>1e-10]
        self.scale=np.maximum(sd,max(float(np.median(pos))*.1,1e-8) if len(pos) else 1.)
        self.x=(b-self.mean)/self.scale
        self.k0=self.x@self.x.T/self.x.shape[1]
        if self.mix:
            if self.prior is None or len(self.prior)!=self.x.shape[1]:raise ValueError('Prior gene-axis mismatch')
            self.u=pathway_features(self.x,self.prior);k=self.u@self.u.T/max(self.u.shape[1],1)
            self.rescale=float(np.trace(self.k0)/max(np.trace(k),1e-12))
            self.k=(1-self.mix)*self.k0+self.mix*self.rescale*k
        else:self.k=self.k0.copy()
        return self
    def cross(self,q):
        q=np.asarray(q,float)
        if q.ndim!=2 or q.shape[1]!=len(self.mean) or not np.isfinite(q).all():raise ValueError('Invalid target baseline')
        x=(q-self.mean)/self.scale;k=x@self.x.T/self.x.shape[1]
        if self.mix:
            u=pathway_features(x,self.prior)
            k=(1-self.mix)*k+self.mix*self.rescale*(u@self.u.T/max(u.shape[1],1))
        return k


def affine_weights(k,kt,penalty,noise=None):
    """Posterior mean weights with an unknown unpenalized constant intercept."""
    k=np.asarray(k,float);kt=np.asarray(kt,float)
    if penalty<=0 or not np.isfinite(penalty):raise ValueError('Positive penalty required')
    n=len(k);v=np.ones(n) if noise is None else np.asarray(noise,float)
    if k.shape!=(n,n) or kt.ndim!=2 or kt.shape[1]!=n or n<1 or v.shape!=(n,) or not np.isfinite(v).all() or np.any(v<=0):raise ValueError('Invalid kernel/noise shapes')
    c=k+penalty*np.diag(v)
    rhs=np.column_stack([kt.T,np.ones(n)])
    ans=solve(c,rhs,assume_a='pos');base=ans[:,:-1].T;ci1=ans[:,-1]
    return base+(1-base.sum(1))[:,None]*ci1[None,:]/ci1.sum()


class SupportTransport:
    def __init__(self,solver='affine',mix=0.,prior=None,noise_gamma=0.):
        self.solver=solver;self.mix=mix;self.prior=prior;self.noise_gamma=noise_gamma
    def fit(self,y,baseline,repeat_variability=None):
        y=np.asarray(y,float);b=np.asarray(baseline,float)
        if y.ndim!=3 or b.shape!=(y.shape[0],y.shape[2]) or np.isinf(y).any():raise ValueError('Invalid effect tensor')
        self.mask=np.isfinite(y).all(-1)
        if not np.array_equal(self.mask,np.isfinite(y).any(-1)):raise ValueError('Partially missing profiles')
        self.y=y;self.count=self.mask.sum(0)
        self.mean=np.nansum(y,0)/np.maximum(self.count[:,None],1);self.mean[self.count==0]=np.nan
        self.geometry=Geometry(self.mix,self.prior).fit(b)
        self.groups={}
        for j in range(y.shape[1]):self.groups.setdefault(self.mask[:,j].tobytes(),[]).append(j)
        self.noise=None
        if self.noise_gamma:
            if repeat_variability is None:raise ValueError('Source variability proxy required')
            v=np.asarray(repeat_variability,float)
            if v.shape!=self.mask.shape:raise ValueError('Source variability proxy shape mismatch')
            good=np.isfinite(v)&(v>1e-12)&self.mask
            med=float(np.median(v[good])) if good.any() else 1.
            normalized=np.ones_like(v);normalized[good]=np.clip(v[good]/med,.1,10)
            self.noise=(1-self.noise_gamma)+self.noise_gamma*normalized
        return self
    def components(self,query_baseline,penalty):
        kt=self.geometry.cross(query_baseline);k=self.geometry.k
        delta=np.full((len(query_baseline),self.y.shape[1],self.y.shape[2]),np.nan)
        for code,js in self.groups.items():
            mask=self.mask[:,js[0]]
            if not mask.any():continue
            km=k[np.ix_(mask,mask)];kmq=kt[:,mask]
            if self.noise is None:
                if self.solver=='affine':w=affine_weights(km,kmq,penalty)
                elif self.solver=='legacy':w=solve(km+penalty*np.eye(mask.sum()),kmq.T,assume_a='pos').T
                else:raise ValueError('Unknown solver')
                yy=self.y[mask][:,js];resid=yy-self.mean[js][None]
                delta[:,js]=np.einsum('tc,cpg->tpg',w,resid,optimize=True)
            else:
                for j in js:
                    w=affine_weights(km,kmq,penalty,self.noise[mask,j]);delta[:,j]=w@(self.y[mask,j]-self.mean[j])
        return self.mean.copy(),delta
    def predict(self,query_baseline,penalty,common=1.,correction=1.):
        mu,d=self.components(query_baseline,penalty)
        return common*mu[None]+correction*d
