"""Post-score attribution: drug-shared own-baseline slope and residual kriging.

The additive-only model cannot predict centred drug-by-context interactions.
This is a regularized fixed-effect baseline, not a mechanistic model.
"""
import numpy as np
from .kernels import SupportTransport

class DriftTransport(SupportTransport):
    def __init__(self,residual=False):
        super().__init__();self.residual=residual
    def fit(self,y,baseline,repeat_variability=None):
        super().fit(y,baseline,repeat_variability)
        x=self.geometry.x
        self.bbar=np.einsum('cp,cg->pg',self.mask.astype(float),x)/np.maximum(self.count[:,None],1)
        xx=np.zeros(x.shape[1]);xy=xx.copy()
        for c in range(len(y)):
            ok=self.mask[c]&(self.count>=2)
            d=x[c]-self.bbar[ok];e=self.y[c,ok]-self.mean[ok]
            xx+=np.sum(d*d,0);xy+=np.sum(d*e,0)
        positive=xx[xx>1e-12];lam=float(np.median(positive)) if len(positive) else 1.
        self.slope=xy/(xx+lam)
        if self.residual:
            self.inner=SupportTransport().fit(self.y-x[:,None,:]*self.slope,baseline)
        return self
    def components(self,query_baseline,penalty):
        self.geometry.cross(query_baseline) # validates the common gene axis
        q=(np.asarray(query_baseline,float)-self.geometry.mean)/self.geometry.scale
        delta=(q[:,None,:]-self.bbar[None])*self.slope
        if self.residual:
            _,correction=self.inner.components(query_baseline,penalty)
            delta=delta+correction
        delta[:,self.count==0]=np.nan
        return self.mean.copy(),delta
