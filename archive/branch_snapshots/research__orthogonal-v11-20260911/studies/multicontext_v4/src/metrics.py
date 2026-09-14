from __future__ import annotations
import itertools
import numpy as np


def score_cells(pred,truth,contexts,drugs,min_sources,source_counts):
 rows=[];bydrug=[]
 for c,ctx in enumerate(contexts):
  ok=np.isfinite(truth[c]).all(-1)&np.isfinite(pred[c]).all(-1)&(source_counts>=min_sources)
  p=pred[c,ok];t=truth[c,ok]
  if len(t)==0:raise ValueError('No comparable held conditions')
  e=p-t; bias=np.mean(e,0)
  pc=p-p.mean(0);tc=t-t.mean(0)
  pn=np.linalg.norm(pc,axis=1);tn=np.linalg.norm(tc,axis=1)
  sims=pc@tc.T/np.maximum(pn[:,None]*tn[None,:],1e-14)
  diag=np.diag(sims);ties=np.isclose(sims,sims.max(1,keepdims=True),atol=1e-12)
  hit=float(np.mean(ties[np.arange(len(p)),np.arange(len(p))]/ties.sum(1)))
  mse=float(np.mean(e**2));null=float(np.mean(t**2))
  rows.append({'context':ctx,'n_compounds':int(ok.sum()),'mse':mse,'zero_mse':null,'skill':1-mse/null,'centered_mse':float(np.mean((e-bias)**2)),'bias_mse':float(np.mean(bias**2)),'rms_prediction':float(np.sqrt(np.mean(p*p))),'rms_observation':float(np.sqrt(null)),'cosine_drug_centered':float(np.mean(diag)),'retrieval_top1':hit})
  for j,d in enumerate(np.asarray(drugs)[ok]):bydrug.append({'context':ctx,'compound':d,'mse':float(np.mean(e[j]**2)),'zero_mse':float(np.mean(t[j]**2))})
 return rows,bydrug


def score_pairs(pred,truth,contexts,source_counts,min_sources=3):
 rows=[]
 for i,j in itertools.combinations(range(len(contexts)),2):
  ok=np.isfinite(truth[i]).all(-1)&np.isfinite(truth[j]).all(-1)&np.isfinite(pred[i]).all(-1)&np.isfinite(pred[j]).all(-1)&(source_counts>=min_sources)
  if ok.sum()<10:continue
  t=truth[i,ok]-truth[j,ok];p=pred[i,ok]-pred[j,ok];e=p-t;n=float(np.mean(t*t));m=float(np.mean(e*e))
  ep=e-e.mean(0);tp=t-t.mean(0)
  rows.append({'context1':contexts[i],'context2':contexts[j],'n_common_drugs':int(ok.sum()),'contrast_mse':m,'zero_contrast_mse':n,'contrast_skill':1-m/n,'centered_contrast_mse':float(np.mean(ep*ep)),'zero_centered_contrast_mse':float(np.mean(tp*tp)),'contrast_rms':float(np.sqrt(np.mean(p*p))),'observed_contrast_rms':float(np.sqrt(n))})
 return rows
