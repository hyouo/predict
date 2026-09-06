"""Anchor-pattern weighted source risk (developmental negative-transfer check)."""
import numpy as np
from .covariance import Fit,build_kernels
from .aggregation import bank
from .deployment_risk import stack_gram

def source_weights(source,anchors,y,beta=1.,uniform_fraction=.25):
    s=np.asarray(source,float)[:,np.asarray(anchors,int),:]
    s=s-s.mean(1,keepdims=True);yc=y-y.mean(0,keepdims=True)
    d=np.mean((s-yc[None,:,:])**2,axis=(1,2))
    sd=np.array([np.mean((s[i]-s[j])**2) for i in range(len(s)) for j in range(i)])
    bw=max(float(np.median(sd[sd>1e-12])) if np.any(sd>1e-12) else 1.,1e-12)
    logw=-beta*d/bw;w=np.exp(logw-logw.max());w/=w.sum()
    w=(1-uniform_fraction)*w+uniform_fraction/len(w)
    return w,{'source_anchor_distances':d.tolist(),'bandwidth':bw,'beta':beta,'effective_source_count':float(1/np.sum(w*w))}

def fit_weighted(source,cues,a,y,q,betas=(1.,4.),alpha=.5):
    k=build_kernels(source,cues);pp,ee,pars=bank(source,cues,k,a,y)
    E=ee.reshape(len(pp),-1);Ct=E@E.T/E.shape[1];grams=[]
    for j in range(len(source)):
        other=np.delete(source,j,axis=0);ko=build_kernels(other,cues)
        pj,_,ps=bank(other,cues,ko,a,source[j,a])
        if pars!=ps:raise AssertionError('alignment')
        er=(pj[:,q]-source[j,q]).reshape(len(pp),-1);grams.append(er@er.T/er.shape[1])
    out={}
    for beta in betas:
        sw,diag=source_weights(source,a,y,beta)
        Cs=np.einsum('s,sij->ij',sw,np.array(grams));w=stack_gram((1-alpha)*Ct+alpha*Cs)
        out[f'weighted_meta_{beta}']=Fit(np.einsum('m,mng->ng',w,pp),
            {'source_weights':sw.tolist(),**diag,'meta_weight':alpha,'weights':w.tolist(),'components':pars},np.empty(0))
    return out
