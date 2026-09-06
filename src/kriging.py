"""Universal-kriging acquisition, explicitly matching an unpenalized intercept.

K and noise are working statistical priors, NOT validated biological uncertainty.
Prediction and design may use full source catalog, not target query outcomes.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import cho_factor,cho_solve

def uk_covariance(K,anchors,noise=0.1):
    k=np.asarray(K,float);a=np.asarray(anchors,int);n=len(k)
    if len(a)<1 or len(set(a))!=len(a):raise ValueError('Need >=1 distinct anchor')
    if noise<=0:raise ValueError('Noise must be positive')
    v=k[np.ix_(a,a)]+float(noise)*np.eye(len(a))
    cf=cho_factor((v+v.T)/2,lower=True)
    inv1=cho_solve(cf,np.ones(len(a)))
    denom=float(inv1.sum())
    if denom<=0:raise ValueError('Nonpositive intercept information')
    A=cho_solve(cf,k[a,:])
    h=1-k[:,a]@inv1
    c=k-k[:,a]@A+np.outer(h,h)/denom
    return (c+c.T)/2

def uk_predict(K,anchors,y_anchor,noise):
    k=np.asarray(K,float);a=np.asarray(anchors,int);y=np.asarray(y_anchor,float)
    v=k[np.ix_(a,a)]+noise*np.eye(len(a))
    cf=cho_factor(v,lower=True);iv1=cho_solve(cf,np.ones(len(a)))
    mean=iv1@y/iv1.sum()
    return mean+k[:,a]@cho_solve(cf,y-mean)

def _utility(cross,denom,contrast_weight):
    # w=0: absolute response risk. w=1: centered query contrasts only.
    # Uniform query weights; query endpoints are never measured as anchors.
    z=cross-cross.mean(0,keepdims=True)
    return ((1-contrast_weight)*np.mean(cross*cross,0)+contrast_weight*np.mean(z*z,0))/denom

def acquire(K,pool,queries,budget,noise=.1,criterion='universal',contrast_weight=0.):
    k=np.asarray(K,float);pool=np.asarray(pool,int);q=np.asarray(queries,int)
    if len(set(pool))!=len(pool) or len(set(q))!=len(q) or np.intersect1d(pool,q).size:
        raise ValueError('Need distinct, disjoint candidate and query sets')
    if not 0<budget<=len(pool) or not len(q) or not 0<=contrast_weight<=1:raise ValueError('invalid design')
    selected=[];cov=None;history=[]
    if criterion=='proper': cov=k+1.  # v0.6 comparator: N(0,1) intercept prior
    for j in range(budget):
        remain=np.array([i for i in pool if i not in selected],int)
        if not selected and criterion=='universal':
            # One observation estimates the free intercept; predictive error = f(q)-f(i)-eps(i).
            vd=np.diag(k)[q,None]+np.diag(k)[None,remain]-2*k[np.ix_(q,remain)]+noise
            risk=(1-contrast_weight)*vd.mean(0)
            # After centering queries, first-anchor-independent part does not select an index.
            # For pure contrast design, break by absolute response risk rather than arbitrary index.
            if contrast_weight==1.: risk=vd.mean(0)
            chosen=int(remain[np.argmin(risk)]);gain=float('nan')
        else:
            if criterion=='universal': cov=uk_covariance(k,selected,noise)
            cross=cov[np.ix_(q,remain)]
            denom=np.maximum(np.diag(cov)[remain]+noise,1e-12)
            score=_utility(cross,denom,contrast_weight)
            chosen=int(remain[np.argmax(score)]);gain=float(score.max())
        selected.append(chosen);history.append({'chosen':chosen,'gain_working_risk':gain})
        if criterion=='proper':
            vec=cov[:,chosen].copy();cov=cov-np.outer(vec,vec)/(cov[chosen,chosen]+noise)
    return np.array(selected,int),history
