"""Full-feature direct/residual ridge; extra comparator, locked before test scoring."""
from pathlib import Path
import sys,json,time,traceback
import numpy as np,pandas as pd
from scipy.linalg import solve
from src.data import Study,sha
from src.models import consensus

class FullRidge:
 def fit(self,y,penalty,residual):
  f,n,m=consensus(y);valid=n>=2;x=[];t=[];w=[]
  for c in range(len(y)):
   ok=m[c]&valid
   a=(n[ok,None]*f[ok]-y[c,ok])/(n[ok,None]-1)
   x.append(a);t.append(y[c,ok]-a if residual else y[c,ok]);w.append(np.repeat(1/max(ok.sum(),1),ok.sum()))
  x=np.concatenate(x);t=np.concatenate(t);w=np.concatenate(w);self.mx=np.average(x,axis=0,weights=w);self.my=np.average(t,axis=0,weights=w)
  x=x-self.mx;t=t-self.my;gram=x.T@(w[:,None]*x);rhs=x.T@(w[:,None]*t);di=np.diag(gram);scale=np.mean(di[di>1e-12]) if np.any(di>1e-12) else 1.
  self.coef=solve(gram+penalty*scale*np.eye(len(gram)),rhs,assume_a='pos');self.residual=residual;self.f=f;self.n=n
  return self
 def predict(self,n_targets):
  out=(np.nan_to_num(self.f)-self.mx)@self.coef+self.my
  if self.residual:out+=self.f
  out[self.n==0]=np.nan
  return np.repeat(out[None],n_targets,axis=0)

def main():
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--data',required=True);a.add_argument('--plan',required=True);a.add_argument('--out',required=True);args=a.parse_args();out=Path(args.out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True)
 try:
  d=Study(args.data,args.plan);yt,_=d.effects(d.plan['train']);yv,_=d.effects(d.plan['validation']);rows=[];chosen={};penalties=[.01,.1,1.,10.,100.]
  for residual in [False,True]:
   name='full_residual_ridge' if residual else 'full_direct_ridge'
   for lam in penalties:
    m=FullRidge().fit(yt,lam,residual);p=m.predict(len(yv));rr=[]
    for c in range(len(yv)):
     ok=np.isfinite(yv[c]).all(-1)&(m.n>=3);rr.append(float(np.mean((p[c,ok]-yv[c,ok])**2)))
    rows.append(dict(model=name,penalty=lam,residual=residual,validation_mse=float(np.mean(rr))))
   chosen[name]=min([v for v in rows if v['model']==name],key=lambda v:(v['validation_mse'],-v['penalty']))
   print(chosen[name],flush=True)
  pd.DataFrame(rows).to_csv(out/'validation.csv',index=False);(out/'selected.json').write_text(json.dumps(chosen,indent=2))
  y,_=d.effects(d.plan['train']+d.plan['validation']);preds={}
  for name,v in chosen.items():
   m=FullRidge().fit(y,v['penalty'],v['residual']);preds[name]=m.predict(len(d.plan['test']));np.savez_compressed(out/(name+'_model.npz'),coef=m.coef,mx=m.mx,my=m.my)
  np.savez_compressed(out/'predictions.npz',**preds);(out/'access.json').write_text(json.dumps(d.log))
  f=dict(complete=True,test_scored=False,time_unix=time.time(),source_sha256=d.plan['source_sha256'],plan_sha256=sha(args.plan),prediction_sha256=sha(out/'predictions.npz'),code_sha256=sha(__file__),chosen=chosen)
  (out/'FREEZE.json').write_text(json.dumps(f,indent=2))
 except Exception as e:(out/'FAILED.json').write_text(json.dumps(dict(error=str(e),traceback=traceback.format_exc())));raise
if __name__=='__main__':main()
