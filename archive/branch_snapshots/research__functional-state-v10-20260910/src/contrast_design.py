"""Costed linear-Gaussian design for controls, new treatments and repeats.

Unlike selecting rows of an effect matrix as if independent, represent baseline
and effects as latent variables. A control observes baseline; a treated well
observes baseline + its effect. Correlated batch/donor effects must be represented
as extra latent variables (or a joint observation covariance), not ignored.

Exact risk reduction is conditional on specified Gaussian/linear assumptions.
No default 'biological variance', surrogate guarantee or synthetic performance.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Action:
    name: str
    design: np.ndarray
    noise_variance: float
    cost: float = 1.
    biological_unit: str = ''


def validate_covariance(C):
    c=np.asarray(C,float)
    if c.ndim!=2 or c.shape[0]!=c.shape[1] or not np.isfinite(c).all():raise ValueError('finite square covariance required')
    if not np.allclose(c,c.T,atol=1e-9):raise ValueError('covariance must be symmetric')
    if np.linalg.eigvalsh(c).min()<-1e-8:raise ValueError('covariance must be positive semidefinite')
    return (c+c.T)/2

def _checked_action(action, dimension):
    h=np.asarray(action.design,float)
    if h.shape!=(dimension,) or not np.isfinite(h).all():
        raise ValueError('finite action vector of latent dimension required')
    if not np.isfinite([action.noise_variance,action.cost]).all() or action.noise_variance<=0 or action.cost<=0:
        raise ValueError('explicit finite positive noise variance and cost required')
    return h

def _checked_contrasts(L,weights,dimension):
    L=np.asarray(L,float)
    if L.ndim!=2 or L.shape[1]!=dimension or not len(L) or not np.isfinite(L).all():
        raise ValueError('nonempty finite contrast matrix of latent dimension required')
    w=np.ones(len(L))/len(L) if weights is None else np.asarray(weights,float)
    if w.shape!=(len(L),) or not np.isfinite(w).all() or np.any(w<0) or w.sum()<=0:
        raise ValueError('finite nonnegative nonzero contrast weights required')
    return L,w

def posterior_update(mean,covariance,action,observed_value):
    c=validate_covariance(covariance);m=np.asarray(mean,float)
    if m.shape!=(len(c),) or not np.isfinite(m).all():raise ValueError('finite latent mean with matching dimension required')
    h=_checked_action(action,len(c))
    if not np.isfinite(observed_value):raise ValueError('nonfinite observation')
    u=c@h;v=float(h@u+action.noise_variance)
    new_m=m+u*(observed_value-h@m)/v
    new_c=c-np.outer(u,u)/v
    return new_m,(new_c+new_c.T)/2

def risk(covariance,L,weights=None):
    c=validate_covariance(covariance);L,w=_checked_contrasts(L,weights,len(c))
    return float(np.sum(w*np.einsum('ij,jk,ik->i',L,c,L)))

def action_gain(covariance,L,action,weights=None):
    c=validate_covariance(covariance);h=_checked_action(action,len(c));L,w=_checked_contrasts(L,weights,len(c))
    v=float(h@c@h+action.noise_variance);z=L@c@h
    gain=float(np.sum(w*z*z)/v)
    return {'absolute_risk_reduction':gain,'risk_reduction_per_cost':gain/action.cost}

def control_and_treatment_actions(n_conditions,control_variance,treatment_variances,costs=None):
    v=np.asarray(treatment_variances,float)
    if n_conditions<1 or v.shape!=(n_conditions,) or not np.isfinite(v).all() or not np.isfinite(control_variance) or np.any(v<=0) or control_variance<=0:raise ValueError('explicit positive variances required')
    costs=np.ones(n_conditions+1) if costs is None else np.asarray(costs,float)
    if costs.shape!=(n_conditions+1,) or not np.isfinite(costs).all() or np.any(costs<=0):raise ValueError('invalid costs')
    h=np.zeros(n_conditions+1);h[0]=1
    result=[Action('control',h.copy(),control_variance,float(costs[0]))]
    for j in range(n_conditions):
        z=h.copy();z[j+1]=1
        result.append(Action(f'treatment_{j}',z,float(v[j]),float(costs[j+1])))
    return result

def continuous_replication_allocation(weights,treated_sd,control_sd,budget,costs=None):
    """Minimize sum_j w_j*(sigma_j^2/n_j + sigma_c^2/n_c), no prior shrinkage.

Continuous allocation, not rounded lab instructions. Variances must come from
appropriate independent biological units. A common control is not duplicated.
"""
    w=np.asarray(weights,float);s=np.asarray(treated_sd,float)
    if w.ndim!=1 or not len(w) or w.shape!=s.shape or not np.isfinite(w).all() or not np.isfinite(s).all() or not np.isfinite([control_sd,budget]).all() or np.any(w<=0) or np.any(s<=0) or control_sd<=0 or budget<=0:raise ValueError('positive aligned inputs required')
    a=np.r_[w.sum()*control_sd**2,w*s*s]
    cost=np.ones(len(a)) if costs is None else np.asarray(costs,float)
    if cost.shape!=a.shape or not np.isfinite(cost).all() or np.any(cost<=0):raise ValueError('positive cost required')
    n=budget*np.sqrt(a/cost)/np.sum(np.sqrt(a*cost))
    return {'control':float(n[0]),'treatments':n[1:],'predicted_measurement_risk':float(np.sum(a/n)), 'continuous':True}

@dataclass(frozen=True)
class VectorAction:
    """One experimental unit with several jointly observed readouts.

    Cost is charged ONCE for the whole observation vector, not per gene.
    The noise covariance must be explicitly specified and positive definite.
    Correlation with past observations needs a joint latent/covariance model;
    a string biological_unit label alone cannot establish independence.
    """
    name: str
    design: np.ndarray
    noise_covariance: np.ndarray
    cost: float = 1.
    biological_unit: str = ''


def _checked_vector_action(action,dimension):
    H=np.asarray(action.design,float);R=validate_covariance(action.noise_covariance)
    if H.ndim!=2 or H.shape[1]!=dimension or not len(H) or not np.isfinite(H).all():
        raise ValueError('finite nonempty observation matrix of latent dimension required')
    if R.shape!=(len(H),len(H)) or not np.isfinite(action.cost) or action.cost<=0:
        raise ValueError('matching observation covariance and finite positive total cost required')
    try:np.linalg.cholesky(R)
    except np.linalg.LinAlgError as exc:raise ValueError('observation noise covariance must be positive definite') from exc
    return H,R


def vector_posterior_update(mean,covariance,action,observed_values):
    C=validate_covariance(covariance);m=np.asarray(mean,float)
    if m.shape!=(len(C),) or not np.isfinite(m).all():raise ValueError('invalid latent mean')
    H,R=_checked_vector_action(action,len(C));y=np.asarray(observed_values,float)
    if y.shape!=(len(H),) or not np.isfinite(y).all():raise ValueError('invalid observation vector')
    U=C@H.T;S=H@U+R
    m1=m+U@np.linalg.solve(S,y-H@m)
    C1=C-U@np.linalg.solve(S,U.T)
    return m1,(C1+C1.T)/2


def vector_action_gain(covariance,L,action,weights=None):
    C=validate_covariance(covariance);L,w=_checked_contrasts(L,weights,len(C))
    H,R=_checked_vector_action(action,len(C));U=C@H.T;S=H@U+R
    V=L@U;Q=np.linalg.solve(S,V.T).T
    gain=float(np.sum(w*np.sum(V*Q,axis=1)))
    return {'absolute_risk_reduction':gain,'risk_reduction_per_cost':gain/action.cost}


def multiplex_actions(n_conditions,n_readouts,control_noise,treated_noise,costs=None):
    """Independent new wells, latent order: [baseline,effects] for each readout.

    `treated_noise`: one (G,G) covariance explicitly assumed common across
    conditions, or (P,G,G) covariances. Extra donor/batch latents must instead be
    supplied directly through VectorAction. This is a dense reference API.
    """
    if n_conditions<1 or n_readouts<1:raise ValueError('positive condition and readout counts required')
    G=n_readouts;P=n_conditions;R0=validate_covariance(control_noise)
    if R0.shape!=(G,G):raise ValueError('control covariance/readout mismatch')
    Rt=np.asarray(treated_noise,float)
    if Rt.shape==(G,G):Rt=np.broadcast_to(Rt,(P,G,G))
    if Rt.shape!=(P,G,G):raise ValueError('treated covariance/readout mismatch')
    cost=np.ones(P+1) if costs is None else np.asarray(costs,float)
    if cost.shape!=(P+1,) or not np.isfinite(cost).all() or np.any(cost<=0):raise ValueError('explicit valid well costs required')
    H=np.zeros((G,G*(P+1)))
    for g in range(G):H[g,g*(P+1)]=1.
    acts=[VectorAction('control_well',H.copy(),R0,float(cost[0]))]
    for j in range(P):
        h=H.copy()
        for g in range(G):h[g,g*(P+1)+j+1]=1.
        acts.append(VectorAction(f'treatment_{j}_well',h,Rt[j].copy(),float(cost[j+1])))
    for a in acts:_checked_vector_action(a,G*(P+1))
    return acts
