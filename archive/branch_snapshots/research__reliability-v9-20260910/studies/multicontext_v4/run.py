"""Fit/validate/freeze, then score test contexts in a separate command."""
from __future__ import annotations
import argparse,json,hashlib,time,traceback,platform,sys
from pathlib import Path
import numpy as np,pandas as pd
from src.data import Study,metadata_plan,sha
from src.models import Predictor
from src.metrics import score_cells,score_pairs

KINDS=['zero','source_mean','scalar','shared_ridge','context_rbf','context_linear','context_additive','interaction_uncentered','interaction_centered']
PENALTIES=[.01,.1,1.,10.,100.]
STRENGTHS=[0.,.25,.5,1.]

def hashes():return {str(p):sha(p) for p in sorted(Path('.').glob('src/*.py'))}|{'run.py':sha('run.py')}
def write(path,obj):Path(path).write_text(json.dumps(obj,indent=2,default=str))

def fit(args):
 out=Path(args.out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True)
 try:
  started=time.time();d=Study(args.data,args.plan);plan=d.plan
  yt,_=d.effects(plan['train']);bt=d.baselines(plan['train'])
  yv,_=d.effects(plan['validation']);bv=d.baselines(plan['validation'])
  btest=d.baselines(plan['test']);validation=[];selected={};cache={}
  for kind in KINDS:
   penalties=[1.] if kind in ['zero','source_mean','scalar'] else PENALTIES
   strengths=STRENGTHS if kind.startswith('interaction_') or kind in ('context_additive','context_rbf','context_linear') else [1.]
   for penalty in penalties:
    m=Predictor(kind=kind,penalty=penalty).fit(yt,bt)
    for st in strengths:
     m.strength=st;p=m.predict(bv)
     cr=[]
     for c in range(len(yv)):
      ok=np.isfinite(yv[c]).all(-1)&np.isfinite(p[c]).all(-1)&(m.count>=3)
      if not ok.any():raise ValueError('Validation context has no supported profiles')
      cr.append({'mse':float(np.mean((p[c,ok]-yv[c,ok])**2))})
     pr=score_pairs(p,yv,plan['validation'],m.count)
     val={'kind':kind,'penalty':penalty,'strength':st,'validation_mse':float(np.mean([v['mse'] for v in cr])),'validation_contrast_mse':float(np.mean([v['contrast_mse'] for v in pr])) if pr else None}
     validation.append(val)
   cand=[v for v in validation if v['kind']==kind]
   selected[kind]=min(cand,key=lambda v:(v['validation_mse'],v['strength'],-v['penalty']))
   if kind.startswith('interaction_'):
    selected[kind+'_contrast_selected']=min(cand,key=lambda v:(v['validation_contrast_mse'] if v['validation_contrast_mse'] is not None else float('inf'),v['strength'],-v['penalty']))
   print('SELECTED',kind,selected[kind],flush=True)
  pd.DataFrame(validation).to_csv(out/'validation_grid.csv',index=False);write(out/'selected.json',selected)
  del yt,yv,bt,bv,m
  contexts=plan['train']+plan['validation'];y,records=d.effects(contexts);b=d.baselines(contexts)
  preds={};models={}
  for name,v in selected.items():
   key=(v['kind'],v['penalty'])
   if key not in cache:
    cache.clear()
    m=None
    cache[key]=Predictor(kind=v['kind'],penalty=v['penalty']).fit(y,b)
   m=cache[key];m.strength=v['strength'];p=m.predict(btest);preds[name]=p
   coeff={}
   for field in ['f','count','z','u','zbar','trunk','interaction','alpha']:
    if hasattr(m,field):coeff[field]=np.asarray(getattr(m,field))
   for field in ['zpca','upca']:
    for k,a in enumerate(getattr(m,field)):coeff[field+str(k)]=np.asarray(a)
   np.savez_compressed(out/(name+'_model.npz'),**coeff)
   models[name]=sha(out/(name+'_model.npz'))
  # Every method has identical source support; no hidden target treatment read.
  np.savez_compressed(out/'predictions.npz',**preds)
  np.savez_compressed(out/'source_inputs.npz',y=y,baseline=b,target_baseline=btest,source_counts=m.count,contexts=np.asarray(contexts),target_contexts=np.asarray(plan['test']),drugs=np.asarray(plan['drugs']),genes=np.asarray(plan['genes']))
  write(out/'fit_access.json',d.log);write(out/'source_records.json',records)
  freeze={'complete':True,'test_scored':False,'time_unix':time.time(),'elapsed_seconds':time.time()-started,'python':sys.version,'platform':platform.platform(),'code_hashes':hashes(),'plan_sha256':sha(args.plan),'prediction_sha256':sha(out/'predictions.npz'),'model_hashes':models,'test_contexts':plan['test'],'source_contexts':contexts,'selected':selected,'fit_access_sha256':sha(out/'fit_access.json'),'target_treatment_calibrations':0}
  write(out/'FREEZE.json',freeze)
 except Exception as e:
  write(out/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise

def score(args):
 root=Path(args.run);out=Path(args.out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True)
 try:
  freeze=json.loads((root/'FREEZE.json').read_text())
  if sha(root/'predictions.npz')!=freeze['prediction_sha256']:raise ValueError('Predictions changed after freeze')
  if sha(args.plan)!=freeze['plan_sha256']:raise ValueError('Plan changed')
  d=Study(args.data,args.plan);d.allow_score(root/'FREEZE.json');truth,records=d.effects(d.plan['test'],role='score')
  preds=np.load(root/'predictions.npz',allow_pickle=False);inputs=np.load(root/'source_inputs.npz',allow_pickle=False);counts=inputs['source_counts']
  rows=[];drugrows=[];pairs=[]
  for name in preds.files:
   r,b=score_cells(preds[name],truth,d.plan['test'],d.drugs,3,counts);ps=score_pairs(preds[name],truth,d.plan['test'],counts)
   rows.extend(dict(model=name,**v) for v in r);drugrows.extend(dict(model=name,**v) for v in b);pairs.extend(dict(model=name,**v) for v in ps)
  frame=pd.DataFrame(rows);frame.to_csv(out/'by_context.csv',index=False);pd.DataFrame(drugrows).to_csv(out/'by_drug.csv',index=False);pd.DataFrame(pairs).to_csv(out/'paired_contrasts.csv',index=False)
  frame.groupby('model').mean(numeric_only=True).reset_index().to_csv(out/'macro.csv',index=False)
  if pairs:pd.DataFrame(pairs).groupby('model').mean(numeric_only=True).reset_index().to_csv(out/'contrast_macro.csv',index=False)
  np.savez_compressed(out/'truth.npz',truth=truth,contexts=np.asarray(d.plan['test']),drugs=np.asarray(d.drugs),genes=np.asarray(d.genes))
  write(out/'score_access.json',d.log);write(out/'target_records.json',records);write(out/'COMPLETE.json',{'complete':True,'freeze_sha256':sha(root/'FREEZE.json'),'truth_sha256':sha(out/'truth.npz'),'time_unix':time.time()})
  print(frame.groupby('model').mse.mean().sort_values().to_string())
 except Exception as e:write(out/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise

def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 s=sp.add_parser('plan');s.add_argument('--data',required=True);s.add_argument('--out',required=True)
 for cmd in ['fit','score']:
  s=sp.add_parser(cmd);s.add_argument('--data',required=True);s.add_argument('--plan',required=True);s.add_argument('--out',required=True)
  if cmd=='score':s.add_argument('--run',required=True)
 a=p.parse_args()
 if a.cmd=='plan':print(json.dumps(metadata_plan(a.data,a.out),default=str)[:1000])
 elif a.cmd=='fit':fit(a)
 else:score(a)
if __name__=='__main__':main()
