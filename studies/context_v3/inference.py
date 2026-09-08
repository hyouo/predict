"""Portable state-conditioned inference. No target treatment measurements required.

Expects the same gene IDs and log2(1+CPM) effect/baseline space as the study.
The drug response must already be observed in the allowed source contexts.
"""
from pathlib import Path
import argparse,json,math
import numpy as np,torch
import hashlib

def sha(path):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

class Decoder(torch.nn.Module):
 def __init__(self,dim,rank):
  super().__init__();self.net=torch.nn.Sequential(torch.nn.Linear(dim,64),torch.nn.GELU(),torch.nn.Linear(64,64),torch.nn.GELU(),torch.nn.Linear(64,rank+1))
 def forward(self,h,features,base):
  d=self.net(features)
  return base+torch.einsum('cpr,cgr->cpg',h,d[:,:,:-1])/math.sqrt(h.shape[-1])+d[:,:,-1][:,None,:]

def predict(bundle,source_mean,target_baseline,genes):
 with np.load(bundle,allow_pickle=False) as z:a={k:z[k] for k in z.files}
 original=[str(x) for x in genes];saved=a['genes'].astype(str).tolist()
 if len(set(original))!=len(original) or set(original)!=set(saved):raise ValueError('Gene IDs must be unique and match the frozen panel exactly')
 x=np.asarray(source_mean,float);b=np.asarray(target_baseline,float)
 if x.ndim!=2 or b.ndim!=2 or x.shape[1]!=len(saved) or b.shape[1]!=len(saved):raise ValueError('Invalid dimensions')
 ix=np.array([original.index(g) for g in saved]);x=x[:,ix];b=b[:,ix]
 if not np.isfinite(x).all() or not np.isfinite(b).all():raise ValueError('Missing/invalid effect or baseline')
 if x.ndim!=2 or b.ndim!=2 or x.shape[1]!=len(saved) or b.shape[1]!=len(saved):raise ValueError('Invalid dimensions')
 h=(x-a['hm'])@a['hv']/a['hs'];base=a['beta']*x+np.column_stack([np.ones(len(x)),h])@a['ridge_weights']
 Z=np.clip((b-a['baseline_mean'])/a['baseline_sd'],-5,5);Zb=np.clip((b-a['bm'])@a['bv']/a['bs'],-5,5)
 C,G=b.shape;gv=a['hv']*np.sqrt(G)
 dyn=np.concatenate([Z[:,:,None],np.broadcast_to(Zb[:,None,:],(C,G,Zb.shape[1])),Z[:,:,None]*gv[None,:,:8],Zb[:,None,:8]*gv[None,:,:8]],axis=-1)
 fixed=a['static_features'];features=np.concatenate([np.broadcast_to(fixed,(C,*fixed.shape)),dyn],axis=-1).astype(np.float32)
 def t(v):return torch.from_numpy(np.array(v,dtype=np.float32))
 torch.set_num_threads(1);outputs=[]
 for i in range(3):
  m=Decoder(features.shape[-1],h.shape[-1]);prefix=f'seed{i}__';state={k[len(prefix):]:torch.from_numpy(v) for k,v in a.items() if k.startswith(prefix)};m.load_state_dict(state);m.eval()
  with torch.no_grad():outputs.append(m(t(np.broadcast_to(h,(C,*h.shape))),t(features),t(np.broadcast_to(base,(C,*base.shape)))).numpy())
 # Float64 averaging matches the pre-scoring frozen ensemble exactly.
 return np.mean(np.stack(outputs).astype(float),axis=0),np.array(saved)

def main():
 p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.output.exists():raise FileExistsError('Prediction outputs are never overwritten')
 with np.load(a.input,allow_pickle=False) as z:
  pred,g=predict(a.model,z['source_mean'],z['target_baseline'],z['genes']);ids=z['drug_ids'].copy();ct=z['contexts'].copy()
 a.output.parent.mkdir(parents=True,exist_ok=True)
 np.savez_compressed(a.output,prediction=pred,genes=g,drug_ids=ids,contexts=ct)
 print(json.dumps({'shape':list(pred.shape),'model_sha256':sha(a.model),'prediction_sha256':sha(a.output),'target_treatment_used':False}))
if __name__=='__main__':main()
