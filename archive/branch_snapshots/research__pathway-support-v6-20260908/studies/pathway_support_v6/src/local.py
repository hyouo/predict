"""Gene-specific biological neighbourhood kernels with the matched affine solver."""
import numpy as np
from .kernels import Geometry,SupportTransport


def neighbourhood(a):
    a=np.asarray(a,float);deg=a.sum(1);sizes=np.maximum(a.sum(0),1)
    w=(a/sizes[None])@a.T/np.maximum(deg[:,None],1)
    w[deg==0]=1/len(a)
    np.testing.assert_allclose(w.sum(1),1,atol=1e-12)
    return w


class LocalTransport(SupportTransport):
    def __init__(self,weights,mix=.5):
        super().__init__();self.weights=np.asarray(weights,float);self.mix_local=mix
    def fit(self,y,baseline,repeat_variability=None):
        super().fit(y,baseline,repeat_variability=None)
        if self.weights.shape!=(self.y.shape[-1],self.y.shape[-1]) or not np.isfinite(self.weights).all() or np.any(self.weights<0):raise ValueError('Invalid gene-local weights')
        x=self.geometry.x
        self.kgene=(1-self.mix_local)*self.geometry.k0[None]+self.mix_local*np.einsum('ch,gh,dh->gcd',x,self.weights,x,optimize=True)
        return self
    def components(self,query_baseline,penalty):
        if penalty<=0 or not np.isfinite(penalty):raise ValueError('Positive penalty required')
        kt=self.geometry.cross(query_baseline)
        q=(np.asarray(query_baseline,float)-self.geometry.mean)/self.geometry.scale
        kgq=(1-self.mix_local)*kt[None]+self.mix_local*np.einsum('th,gh,ch->gtc',q,self.weights,self.geometry.x,optimize=True)
        delta=np.full((len(q),self.y.shape[1],self.y.shape[2]),np.nan)
        for code,js in self.groups.items():
            mask=np.flatnonzero(self.mask[:,js[0]])
            if not len(mask):continue
            k=self.kgene[:,mask][:,:,mask]+penalty*np.eye(len(mask))[None]
            qq=kgq[:,:,mask];rhs=np.concatenate([qq.transpose(0,2,1),np.ones((len(k),len(mask),1))],axis=2)
            ans=np.linalg.solve(k,rhs);base=ans[:,:,:-1].transpose(0,2,1);ci1=ans[:,:,-1]
            w=base+(1-base.sum(2))[:,:,None]*ci1[:,None,:]/ci1.sum(1)[:,None,None]
            resid=self.y[mask][:,js]-self.mean[js][None]
            delta[:,js]=np.einsum('gtc,cpg->tpg',w,resid,optimize=True)
        return self.mean.copy(),delta
