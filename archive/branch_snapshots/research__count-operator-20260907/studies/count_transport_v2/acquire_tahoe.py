"""One fixed public Tahoe pseudobulk asset -> metadata-locked portable subset.

Does not train models or compute target-test performance. Biological effects
are never consulted for the cohort, role assignment, or gene-panel selection.
Source-control-only gene panel; reserved test records exported separately.
"""
from pathlib import Path
import urllib.request,json,hashlib,time,datetime,os
import h5py,numpy as np,pandas as pd

REV='6c54c8eb0321cceff4f888c54e199077d055e20b'
EXPECTED='09ac7c8f63a77dc730156c78e361273c3369959640120e5c54241f1a6f82d7b2'
SIZE=4473155146
URL=f'https://huggingface.co/datasets/theislab/chem-perturbridge/resolve/{REV}/tahoe/tahoe_processed.h5ad?download=true'
OUT=Path('tahoe_subset');OUT.mkdir(exist_ok=False)

def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def digest(path):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rankid(x):return hashlib.sha256(('count-operator-20260907|'+str(x)).encode()).hexdigest()
def col(n):
 if isinstance(n,h5py.Dataset):return n.asstr()[:] if h5py.check_string_dtype(n.dtype) else n[:]
 if 'codes' in n:
  c=col(n['categories']);return np.array([c[i] if i>=0 else None for i in n['codes'][:]])
 raise ValueError('Unsupported metadata encoding')

def download():
 tmp=Path('tahoe.fixed.partial');h=hashlib.sha256();size=0;t=time.monotonic()
 req=urllib.request.Request(URL,headers={'User-Agent':'hyouo-predict-fixed-public-data'})
 with urllib.request.urlopen(req,timeout=60) as r,open(tmp,'wb') as f:
  while b:=r.read(4*1024*1024):
   size+=len(b)
   if size>SIZE or time.monotonic()-t>600:raise RuntimeError('Download size or time cap')
   f.write(b);h.update(b)
 if size!=SIZE or h.hexdigest()!=EXPECTED:raise ValueError('Upstream checksum/size mismatch')
 p=Path('tahoe.fixed.h5ad');tmp.rename(p)
 dump(OUT/'retrieval.json',{'source':URL,'revision':REV,'bytes':size,'sha256':h.hexdigest(),'verified_utc':stamp(),'seconds':time.monotonic()-t})
 print('Verified public file',size,flush=True);return p

def choose_metadata(o):
 o=o.copy();o['control']=o.is_control.astype(str).str.lower().isin(['true','1','yes']);o['idx']=np.arange(len(o))
 o['time']=pd.to_numeric(o.pert_time_h,errors='coerce');o['dose']=pd.to_numeric(o.pert_dose_uM,errors='coerce');o['exposure']=pd.to_numeric(o.psbulk_counts,errors='coerce')
 non=o[(~o.control)&np.isfinite(o.time)&np.isfinite(o.dose)&(o.exposure>0)]
 tm=float(non.time.value_counts().sort_index().idxmax());n=non[np.isclose(non.time,tm)];ctrl=o[o.control&np.isclose(o.time,tm)&(o.exposure>0)]
 eligible=[]
 for c,grp in n.groupby('cell_type',observed=True):
  if grp.perturbagen.nunique()>=40 and len(ctrl[ctrl.cell_type==c])>=4:eligible.append(str(c))
 lines=sorted(eligible,key=rankid)[:48]
 if len(lines)<12:raise ValueError('Insufficient contexts for fixed multicontext study')
 ntest=min(8,max(2,len(lines)//6));nval=min(4,max(2,len(lines)//10))
 split={c:('test_sealed' if i<ntest else 'validation' if i<ntest+nval else 'train') for i,c in enumerate(lines)}
 training={c for c in lines if split[c]=='train'}
 choices=[]
 for drug,gr in n[n.cell_type.isin(training)].groupby('perturbagen',observed=True):
  tab=gr.groupby('dose',observed=True).cell_type.nunique();best=sorted(tab.items(),key=lambda v:(-v[1],v[0]))[0]
  if best[1]>=max(2,int(np.ceil(.8*len(training)))):choices.append((str(drug),float(best[0]),int(best[1])))
 choices=sorted(choices,key=lambda x:(-x[2],rankid(x[0])))[:192]
 if len(choices)<40:raise ValueError('Insufficient source-covered chemicals')
 conditions=[{'id':f'{d}|{dose:.8g}uM|{tm:.8g}h','drug':d,'dose':dose,'time':tm,'train_coverage':cov} for d,dose,cov in choices]
 controls={};treatments={};missing={}
 for c in lines:
  ci=ctrl[ctrl.cell_type==c].copy();ci['rk']=ci.sample_id.astype(str).map(rankid)
  batches=[list(g.sort_values('rk').idx) for _,g in ci.groupby('plate',sort=True,observed=True)]
  inter=[v[k] for k in range(max(map(len,batches))) for v in batches if k<len(v)]
  nb=min(8,max(2,len(inter)//2));controls[c]={'baseline':list(map(int,inter[:nb])),'reserved_audit':list(map(int,inter[nb:]))}
  treatments[c]={};missing[c]=[]
  for cond in conditions:
   ix=n[(n.cell_type==c)&(n.perturbagen==cond['drug'])&np.isclose(n.dose,cond['dose'])].sort_values('sample_id').idx.to_numpy()[:2]
   if len(ix):treatments[c][cond['id']]=list(map(int,ix))
   else:missing[c].append(cond['id'])
 plan={'created_utc_before_count_selection':stamp(),'data_sha256':EXPECTED,'metadata_only_selection':True,'time_h':tm,'cell_types':lines,'roles':split,'conditions':conditions,'controls':controls,'treatments':treatments,'missing_pairs':missing,'max_treatment_records_per_pair':2,'max_baseline_controls_per_context':8,'gene_panel_policy':'top 5000 average TRAINING-context baseline-count proportions, stable index tie-break; all remaining recorded molecules in OTHER','query_target_values_used_for_selection':False,'neural_training_performed':False,'biological_independence_of_records_not_asserted':True}
 return o,plan

def main():
 data=download();start=time.monotonic()
 with h5py.File(data,'r',rdcc_nbytes=128*1024*1024) as f:
  o=pd.DataFrame({k:col(v) for k,v in f['obs'].items()});x=f['X'];shape=tuple(map(int,x.attrs['shape']))
  if x.attrs.get('encoding-type') not in ['csr_matrix',b'csr_matrix']:raise ValueError('Expected CSR counts')
  key=f['var'].attrs.get('_index','_index');key=key.decode() if isinstance(key,bytes) else key;genes=col(f['var'][key]).astype(str)
  symbols=col(f['var/symbol']).astype(str) if 'symbol' in f['var'] else genes.copy()
  if len(set(genes))!=len(genes):raise ValueError('Gene identifiers are not unique')
  if 'sample_id' not in o:
   ki=f['obs'].attrs.get('_index','_index');ki=ki.decode() if isinstance(ki,bytes) else ki;o['sample_id']=col(f['obs'][ki])
  o,plan=choose_metadata(o);dump(OUT/'selection_pre_counts.json',plan)
  o.to_csv(OUT/'metadata_all.csv.gz',index=False);ptr=x['indptr'][:];vals=x['data'];inds=x['indices'];n=o.exposure.to_numpy(float)
  print('Metadata locked',shape,'contexts',len(plan['cell_types']),'conditions',len(plan['conditions']),flush=True)
  def readrow(i):
   lo,hi=map(int,ptr[i:i+2]);v=vals[lo:hi];j=inds[lo:hi]
   if not np.isfinite(v).all() or np.any(v<0) or np.any(v!=np.floor(v)):raise ValueError('Non-count observation')
   if not np.isfinite(n[i]) or n[i]!=np.floor(n[i]) or np.sum(v,dtype=float)>n[i]+.5:raise ValueError('Invalid exposure')
   return np.asarray(j,int),np.asarray(v,np.int64)
  relevance=np.zeros(shape[1]);ns=0
  for c in plan['cell_types']:
   if plan['roles'][c]!='train':continue
   agg=np.zeros(shape[1]);exposure=0
   for i in plan['controls'][c]['baseline']:
    j,v=readrow(i);np.add.at(agg,j,v);exposure+=n[i]
   relevance+=agg/exposure;ns+=1
  relevance/=ns;panel=np.sort(np.lexsort((np.arange(shape[1]),-relevance))[:min(5000,shape[1])]);mapping=np.full(shape[1],-1,int);mapping[panel]=np.arange(len(panel))
  np.savez_compressed(OUT/'gene_panel.npz',indices=panel,genes=genes[panel],symbols=symbols[panel],source_control_relevance=relevance[panel])
  def extract(ix):
   a=np.zeros((len(ix),len(panel)+1),np.int64)
   for k,i in enumerate(ix):
    j,v=readrow(i);p=mapping[j];ok=p>=0;np.add.at(a[k],p[ok],v[ok]);a[k,-1]=int(n[i])-int(a[k,:-1].sum())
   return a
  inventory=[]
  for ordinal,c in enumerate(plan['cell_types']):
   role=plan['roles'][c];dest=OUT/role;dest.mkdir(exist_ok=True);bi=plan['controls'][c]['baseline'];ai=plan['controls'][c]['reserved_audit'];tr=[];ci=[]
   for k,cond in enumerate(plan['conditions']):
    ids=plan['treatments'][c].get(cond['id'],[]);tr.extend(ids);ci.extend([k]*len(ids))
   path=dest/(f'context_{ordinal:03d}.npz')
   np.savez_compressed(path,baseline_counts=extract(bi),baseline_rows=np.array(bi),audit_control_counts=extract(ai),audit_control_rows=np.array(ai),treatment_counts=extract(tr),treatment_rows=np.array(tr),condition_indices=np.array(ci),treatment_cells=o.iloc[tr].psbulk_cells.to_numpy(),treatment_plates=o.iloc[tr].plate.to_numpy(dtype=str),genes=np.concatenate([genes[panel],['__OTHER__']]))
   inventory.append({'cell_type':c,'role':role,'relative_path':str(path.relative_to(OUT)),'sha256':digest(path),'treatment_records':len(tr),'conditions_observed':len(set(ci)),'baseline_records':len(bi),'reserved_controls':len(ai),'bytes':path.stat().st_size})
   print('Saved',c,role,len(tr),flush=True)
  dump(OUT/'inventory.json',{'created_utc':stamp(),'shape_original':shape,'n_genes_selected':len(panel),'other_bin':True,'profiles':inventory,'selection_sha256':digest(OUT/'selection_pre_counts.json'),'gene_panel_sha256':digest(OUT/'gene_panel.npz'),'training_control_only_panel':True,'source_file_sha256':EXPECTED,'extraction_seconds':time.monotonic()-start,'no_prediction_or_scoring_performed':True})
 (OUT/'ATTRIBUTION.md').write_text('Derived from Tahoe-100M (Tahoe Bio), upstream CC0 1.0; source https://huggingface.co/datasets/tahoebio/Tahoe-100M . Standardized pseudobulk collection: theislab/Chem-PerturBridge, curation/packaging CC BY 4.0. Preserve upstream provenance, notices and citations. This is a metadata-selected subset, not the full Tahoe-100M cell matrix. Raw count profiles are aggregates, not individual cells.\n')
 data.unlink();print('COMPLETE: verified raw cache removed, only fixed public subset retained',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  dump(OUT/'FAILED.json',{'timestamp':stamp(),'error':repr(exc),'no_model_performance_claim':True});raise
