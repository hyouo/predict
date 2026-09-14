"""Small quantitative-RNA baselines and paired-reference CV correction.

Every candidate is affine in the target anchor effects conditional on source
inputs. No target evaluation outcome is an argument to any fit function.
This is an experimental estimator, not a claim of a new statistical principle.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import minimize

@dataclass(frozen=True)
class Candidate:
    name: str
    kernel: str='fixed'
    regularization: float=1.0
    intercept: str='none'
    prior: str='source'
    strength: float=0.0

def candidate_bank():
    result=[Candidate('zero',prior='zero'),Candidate('source'),Candidate('template',prior='template'),
            Candidate('anchor_mean','offset',prior='zero',strength=1.0)]
    for t in (.25,.5,1.0): result.append(Candidate(f'source_offset_{t}','offset',strength=t))
    for kind in ('shared_linear','shared_rbf','gene_linear','mixed_linear','mixed_rbf'):
        for intercept in ('none','shrink','free'):
            for lam in (.03,.3,3.,30.):
                result.append(Candidate(f'{kind}/{intercept}/{lam}',kind,lam,intercept))
    for lam in (.03,.3,3.,30.):
        result.append(Candidate(f'direct_ridge/{lam}','shared_linear',lam,'free','zero'))
    return result

def _finite_shape(source, yc, yd, indices):
    source=np.asarray(source,dtype=float);yc=np.asarray(yc,dtype=float);yd=np.asarray(yd,dtype=float)
    ix=np.asarray(indices)
    if source.ndim!=3 or source.shape[0]<1 or source.shape[1]<2: raise ValueError('Source must be context x condition x gene')
    if yc.ndim!=2 or yc.shape!=yd.shape or yc.shape[1]!=source.shape[2] or len(yc)<3: raise ValueError('Paired anchor shape invalid')
    if ix.ndim!=1 or not np.issubdtype(ix.dtype,np.integer) or len(ix)!=len(yc) or len(set(ix.tolist()))!=len(ix): raise ValueError('Anchor IDs invalid')
    if np.any(ix<0) or np.any(ix>=source.shape[1]): raise ValueError('Anchor IDs out of bounds')
    if not all(np.isfinite(z).all() for z in (source,yc,yd)): raise ValueError('Missing/nonfinite observations must be explicitly resolved before fitting')
    return source,yc,yd,ix.astype(int)

def shared_kernels(source):
    s=source.mean(axis=0); z=s-s.mean(axis=0)
    k=z@z.T/z.shape[1]; scale=float(np.trace(k)/len(k))
    k=k/max(scale,1e-12)
    dist=np.maximum(np.diag(k)[:,None]+np.diag(k)[None,:]-2*k,0.)
    positive=dist[np.triu_indices(len(k),1)];med=float(np.median(positive))
    return {'shared_linear':k, 'shared_rbf':np.exp(-dist/max(med,1e-12))}

def kernel_blocks(source, start, stop, shared):
    z=source[:,:,start:stop];z=z-z.mean(axis=1,keepdims=True)
    gene=np.einsum('cpg,cqg->gpq',z,z,optimize=True)/len(source)
    scale=np.trace(gene,axis1=1,axis2=2)/source.shape[1]
    gene=gene/np.maximum(scale[:,None,None],1e-8)
    b=stop-start
    lin=np.broadcast_to(shared['shared_linear'],(b,)+shared['shared_linear'].shape)
    rbf=np.broadcast_to(shared['shared_rbf'],(b,)+shared['shared_rbf'].shape)
    return {'gene_linear':gene,'shared_linear':lin,'shared_rbf':rbf,
            'mixed_linear':.5*(gene+lin),'mixed_rbf':.5*(gene+rbf)}

def smoother(kernel, anchors, regularization, intercept):
    """Return exact leave-one-anchor-out and final query linear maps.

    Batch dimension is gene. 'free' is a flat intercept prior; 'shrink' adds
    a unit constant kernel; 'none' shrinks residuals directly to the source.
    """
    k=np.asarray(kernel,dtype=float);ix=np.asarray(anchors,dtype=int)
    if regularization<=0 or not np.isfinite(regularization): raise ValueError('Positive regularization required')
    if intercept not in ('none','shrink','free'): raise ValueError('Invalid intercept mode')
    if intercept=='shrink':k=k+1.
    aa=k[:,ix][:,:,ix]; n=len(ix);ident=np.eye(n)
    inv=np.linalg.solve(aa+regularization*ident,np.broadcast_to(ident,aa.shape))
    if intercept=='free':
        v=inv.sum(axis=2);norm=v.sum(axis=1)
        if np.any(norm<=0):raise ValueError('Invalid intercept precision')
        q=inv-v[:,:,None]*v[:,None,:]/norm[:,None,None]
        final=k[:,:,ix]@q+np.broadcast_to((v/norm[:,None])[:,None,:],(len(k),k.shape[1],n))
    else:
        q=inv;final=k[:,:,ix]@q
    diag=np.diagonal(q,axis1=1,axis2=2)
    if np.any(diag<=1e-14):raise ValueError('Degenerate LOO precision')
    loo=ident[None,:,:]-q/diag[:,:,None]
    loo[:,np.arange(n),np.arange(n)]=0.
    return loo,final

def prior_block(candidate, source_mean, start, stop):
    z=source_mean[:,start:stop]
    if candidate.prior=='zero':return np.zeros_like(z)
    if candidate.prior=='template':return np.broadcast_to(z.mean(axis=0),z.shape)
    return z

def maps_for(candidate,kernels,anchors,n_genes,n_conditions):
    n=len(anchors)
    if candidate.kernel=='fixed':
        return np.zeros((n_genes,n,n)),np.zeros((n_genes,n_conditions,n))
    if candidate.kernel=='offset':
        loo=candidate.strength*(np.ones((n,n))-np.eye(n))/(n-1)
        final=np.full((n_conditions,n),candidate.strength/n)
        return np.broadcast_to(loo,(n_genes,n,n)),np.broadcast_to(final,(n_genes,n_conditions,n))
    return smoother(kernels[candidate.kernel],anchors,candidate.regularization,candidate.intercept)

def fit_weights(gram, correction, *, adjusted: bool):
    """Convex simplex stacking; the reference correction is LINEAR in weights.

    Do not PSD-project a 'corrected Gram': that would change the estimand.
    Selection on estimated risk still has ordinary selection uncertainty.
    """
    g=(np.asarray(gram)+np.asarray(gram).T)/2;c=np.asarray(correction) if adjusted else np.zeros(len(g))
    if not np.isfinite(g).all() or not np.isfinite(c).all():raise ValueError('Nonfinite risk')
    ridge=1e-7*max(float(np.diag(g).mean()),1e-8)
    g=g+ridge*np.eye(len(g)); scale=max(float(np.diag(g).mean()),1e-8)
    g=g/scale;c=c/scale
    init=np.zeros(len(g));init[int(np.argmin(np.diag(g)+.5*c))]=1
    result=minimize(lambda w:float(w@g@w+.5*w@c),init,
                    jac=lambda w:2*g@w+.5*c,method='SLSQP',bounds=[(0,1)]*len(g),
                    constraints=[{'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(len(w))}],
                    options={'maxiter':1500,'ftol':1e-12})
    if not result.success or abs(result.x.sum()-1)>1e-6:raise RuntimeError('Stacking failed: '+result.message)
    w=np.clip(result.x,0,1);w=w/w.sum()
    return w,{'success':bool(result.success),'iterations':int(result.nit),'objective':float(result.fun),'weight_ridge':ridge}

def fit_rna_bank(source, anchor_indices, reference_c, reference_d, *, block_size=256, gene_permutation=None):
    source,yc,yd,ix=_finite_shape(source,reference_c,reference_d,anchor_indices)
    if isinstance(block_size,bool) or not isinstance(block_size,int) or block_size<1:raise ValueError('Invalid block size')
    ncond=source.shape[1];ng=source.shape[2];n=len(ix);bank=candidate_bank();m=len(bank)
    sh=shared_kernels(source);mean=source.mean(axis=0);ybar=(yc+yd)/2;difference=yc-yd
    gsum=np.zeros((m,m));csum=np.zeros(m);crosssum=np.zeros(m)
    # Per-compound risk is retained for resampling/diagnostics, not extra biological replicates.
    ordinary_rows=np.zeros((m,n));correct_rows=np.zeros((m,n))
    for start in range(0,ng,block_size):
        stop=min(start+block_size,ng);b=stop-start
        kernels=kernel_blocks(source,start,stop,sh)
        if gene_permutation is not None:
            perm=np.asarray(gene_permutation)
            if sorted(perm.tolist())!=list(range(ng)):raise ValueError('Invalid gene permutation')
            z=source[:,:,perm[start:stop]];z=z-z.mean(axis=1,keepdims=True)
            kg=np.einsum('cpg,cqg->gpq',z,z,optimize=True)/len(source)
            kg/=np.maximum((np.trace(kg,axis1=1,axis2=2)/ncond)[:,None,None],1e-8)
            kernels['gene_linear']=kg;kernels['mixed_linear']=.5*(kg+kernels['shared_linear']);kernels['mixed_rbf']=.5*(kg+kernels['shared_rbf'])
        errors=np.empty((m,n,b));diffpred=np.empty_like(errors)
        for j,can in enumerate(bank):
            loo,_=maps_for(can,kernels,ix,b,ncond);pr=prior_block(can,mean,start,stop)
            residual=ybar[:,start:stop]-pr[ix]
            errors[j]=np.einsum('gij,jg->ig',loo,residual)-residual
            diffpred[j]=np.einsum('gij,jg->ig',loo,difference[:,start:stop])
        flat=errors.reshape(m,-1);gsum+=flat@flat.T
        co=np.sum(diffpred*difference[None,:,start:stop],axis=2)
        csum+=co.sum(axis=1); ordinary_rows+=np.sum(errors**2,axis=2);correct_rows+=co
        crosssum+=np.sum((diffpred+difference[None,:,start:stop])**2,axis=(1,2))
    gram=gsum/(n*ng);corr=csum/(n*ng);ordinary_rows/=ng;correct_rows/=ng
    risk=np.diag(gram);adj=risk+.5*corr;cross=risk+.25*crosssum/(n*ng)
    choices={}
    def selected(label,mask,adjusted=False):
        ids=np.flatnonzero(mask);j=ids[np.argmin((adj if adjusted else risk)[ids])]
        w=np.zeros(m);w[j]=1;choices[label]=w
    for label in ('zero','source','template','anchor_mean','source_offset_1.0'):
        selected(label,[c.name==label for c in bank])
    for label in ('shared_linear','shared_rbf','gene_linear','direct_ridge'):
        mask=[(c.kernel==label and c.prior=='source') if label!='direct_ridge' else c.name.startswith('direct_ridge/') for c in bank]
        selected(label+'_cv',mask);selected(label+'_adjusted',mask,True)
    simple=np.array([not c.kernel.startswith('mixed') for c in bank])
    selected('simple_select_cv',simple);selected('simple_select_adjusted',simple,True)
    diagnostics={}
    for label,mask,adjusted in [('simple_stack_cv',simple,False),('simple_stack_adjusted',simple,True),
                                ('mixed_stack_cv',np.ones(m,dtype=bool),False),('mixed_stack_adjusted',np.ones(m,dtype=bool),True)]:
        ids=np.flatnonzero(mask);w0,rec=fit_weights(gram[np.ix_(ids,ids)],corr[ids],adjusted=adjusted)
        w=np.zeros(m);w[ids]=w0;choices[label]=w;diagnostics[label]=rec
    predictions={k:np.zeros((ncond,ng)) for k in choices}
    for start in range(0,ng,block_size):
        stop=min(start+block_size,ng);b=stop-start;kernels=kernel_blocks(source,start,stop,sh)
        if gene_permutation is not None:
            z=source[:,:,np.asarray(gene_permutation)[start:stop]];z-=z.mean(axis=1,keepdims=True)
            kg=np.einsum('cpg,cqg->gpq',z,z,optimize=True)/len(source);kg/=np.maximum((np.trace(kg,axis1=1,axis2=2)/ncond)[:,None,None],1e-8)
            kernels['gene_linear']=kg;kernels['mixed_linear']=.5*(kg+kernels['shared_linear']);kernels['mixed_rbf']=.5*(kg+kernels['shared_rbf'])
        for j,can in enumerate(bank):
            needed=[k for k,w in choices.items() if w[j]>1e-12]
            if not needed:continue
            _,final=maps_for(can,kernels,ix,b,ncond);pr=prior_block(can,mean,start,stop)
            pred=pr+np.einsum('gij,jg->ig',final,ybar[:,start:stop]-pr[ix])
            for k in needed:predictions[k][:,start:stop]+=choices[k][j]*pred
    record={'candidates':[asdict(c) for c in bank],'ordinary_cv':risk.tolist(),'adjusted_cv':adj.tolist(),
            'cross_reference_cv':cross.tolist(),'control_correction':corr.tolist(),'gram':gram.tolist(),
            'weights':{k:v.tolist() for k,v in choices.items()},'optimization':diagnostics,
            'ordinary_risk_by_anchor':ordinary_rows.tolist(),'correction_by_anchor':correct_rows.tolist(),
            'target_anchor_count':n,'target_test_outcomes_received':False,
            'note':'Correction assumes exchangeable independent reference errors conditional on fixed features; selected-model unbiasedness is NOT guaranteed.'}
    return predictions,record

def zero_shot(source, source_baselines, target_baseline):
    """Small source-only baselines. Gate cannot create an effect from no source signal."""
    s=np.asarray(source,dtype=float);b=np.asarray(source_baselines,dtype=float);bt=np.asarray(target_baseline,dtype=float)
    if s.shape[0]!=2 or b.shape!=(2,s.shape[2]) or bt.shape!=(s.shape[2],):raise ValueError('Two source contexts required')
    if not all(np.isfinite(z).all() for z in (s,b,bt)):raise ValueError('Nonfinite inputs')
    def features(x,src_b,tgt_b):
        return np.stack([x,x*np.tanh((tgt_b-src_b)[None,:]/2),x*np.tanh((tgt_b+src_b-14)[None,:]/4)],axis=-1)
    gram=np.zeros((3,3));rhs=np.zeros(3);xy=0.;xx=0.
    for a,t in ((0,1),(1,0)):
        x=features(s[a],b[a],b[t]).reshape(-1,3); y=s[t].reshape(-1)
        gram+=x.T@x;rhs+=x.T@y;xy+=float(s[a].ravel()@y);xx+=float(s[a].ravel()@s[a].ravel())
    coef=np.linalg.solve(gram+.01*np.trace(gram)/3*np.eye(3),rhs)
    alpha=float(np.clip(xy/max(xx,1e-12),0,1.5))
    z=s.mean(axis=0);x=features(z,b.mean(axis=0),bt)
    gene_alpha=2*np.sum(s[0]*s[1],axis=0)/np.maximum(np.sum(s[0]**2+s[1]**2,axis=0),1e-8)
    gene_alpha=np.clip(gene_alpha,0,1)
    return {'zero':np.zeros_like(z),'source':z,'template':np.broadcast_to(z.mean(axis=0),z.shape).copy(),
            'source_scalar':alpha*z,'source_consensus_shrink':gene_alpha*z,'baseline_gate':x@coef},\
           {'scalar':alpha,'gate_coefficients':coef.tolist(),'fit_contexts':['T','NK'],'target_treatment_rows_used':0}
