"""Evaluation contrasts that remove exactly the registered plate-reference nuisance."""
from __future__ import annotations
import numpy as np


def reference_design(obs, row_ids, compounds):
    rows=obs.iloc[np.asarray(row_ids,int)]
    plates=sorted(rows.plate.unique())
    h=np.zeros((len(compounds),len(plates)))
    for i,p in enumerate(compounds):
        sub=rows[rows.perturbagen==p]
        if not len(sub):raise ValueError('Query lacks target replicates')
        for j,plate in enumerate(plates):h[i,j]=np.mean(sub.plate==plate)
    return h


def orthogonal_projector(h):
    h=np.asarray(h,float)
    if h.ndim!=2 or not np.isfinite(h).all():raise ValueError('Invalid nuisance design')
    u,s,_=np.linalg.svd(h,full_matrices=False)
    tol=max(h.shape)*np.finfo(float).eps*(s[0] if len(s) else 0.)
    rank=int(np.sum(s>tol));q=u[:,:rank]
    if rank>=len(h):raise ValueError('No estimable residual contrasts')
    return np.eye(len(h))-q@q.T,rank


def contrast_metrics(y,p,h):
    """No expression-dependent gene/contrast selection. H uses public metadata only."""
    proj,rank=orthogonal_projector(h)
    yy=proj@np.asarray(y,float);pp=proj@np.asarray(p,float)
    denominator=(len(h)-rank)*y.shape[1]
    norm=np.linalg.norm(yy)*np.linalg.norm(pp)
    return {'plate_orthogonal_mse':float(np.sum((pp-yy)**2)/denominator),
            'plate_orthogonal_cosine':float(np.sum(pp*yy)/norm) if norm>1e-12 else 0.,
            'reference_design_rank':rank,'contrast_degrees_of_freedom':len(h)-rank}
