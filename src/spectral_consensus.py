"""Source-only empirical spectral shrinkage; protocol 006.

The two views are biological contexts, not assumed technical replicates.
The estimated nuisance includes context-specific effects. No target input exists.
"""
from __future__ import annotations
import numpy as np


def spectral_consensus(source, *, ridge: float=1.0, gene_weighted: bool=True):
    x=np.asarray(source,dtype=float)
    if x.ndim!=3 or x.shape[0]!=2 or x.shape[1]<3 or x.shape[2]<1 or not np.isfinite(x).all():
        raise ValueError('Two finite source-context matrices are required')
    if not np.isfinite(ridge) or ridge<=0:raise ValueError('Ridge must be positive')
    p,g=x.shape[1:];mean=x.mean(axis=0);template=mean.mean(axis=0)
    center=x-x.mean(axis=1,keepdims=True)
    m=center.mean(axis=0);d=(center[0]-center[1])/2
    n=d@d.T/g;scale=max(float(np.trace(n)/p),1e-12)
    v=n+ridge*scale*np.eye(p)
    eig,u=np.linalg.eigh((v+v.T)/2)
    if np.min(eig)<=0:raise ValueError('Nonpositive regularized nuisance covariance')
    sqrt=(u*np.sqrt(eig))@u.T;whitener=(u/np.sqrt(eig))@u.T
    z=whitener@m;c=z@z.T/g
    lam,basis=np.linalg.eigh((c+c.T)/2)
    shrink=np.maximum(1-1/np.maximum(lam,1e-12),0)
    filtered=sqrt@((basis*shrink)@(basis.T@z))
    alpha=np.clip(2*np.sum(x[0]*x[1],axis=0)/np.maximum(np.sum(x[0]**2+x[1]**2,axis=0),1e-8),0,1)
    if gene_weighted:filtered=filtered*alpha[None,:]
    pred=filtered+alpha[None,:]*template[None,:]
    if not np.isfinite(pred).all():raise FloatingPointError('Nonfinite spectral prediction')
    return pred,{'ridge':float(ridge),'gene_weighted':bool(gene_weighted),'positive_modes':int((shrink>0).sum()),
                 'effective_modes':float(shrink.sum()),'eigenvalues':lam.tolist(),'shrinkage':shrink.tolist(),
                 'nuisance_diagonal_mean':float(np.trace(n)/p),'source_contexts':2,
                 'target_treatment_rows_received':0,'target_baseline_rows_received':0,
                 'nuisance_interpretation':'context disagreement plus measurement; not pure technical noise'}
