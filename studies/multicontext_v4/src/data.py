"""Metadata-locked L1000 reader. It keeps supplied Level-3 units unchanged."""
from __future__ import annotations
import hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA='34d198df9eddac5b535a0408794caaf41073fa014c74429d1050e483576fb581'

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def order_key(x):return hashlib.sha256(('multicontext-v4-20260908:'+str(x)).encode()).hexdigest()

def metadata_plan(root,out):
 root=Path(root);out=Path(out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True)
 obs=pd.read_csv(root/'slice_obs.csv.gz',index_col=0,dtype={'perturbagen':str,'cell_type':str,'plate':str,'well':str},low_memory=False);var=pd.read_csv(root/'slice_var.csv.gz',index_col=0)
 required=['cell_type','plate','well','is_control','perturbagen','original_row']
 if any(x not in obs for x in required):raise ValueError('Missing required metadata')
 obs['is_control']=obs.is_control.astype(str).str.lower().isin(['true','1'])
 if obs.index.duplicated().any() or obs.original_row.duplicated().any() or var.index.duplicated().any():raise ValueError('Duplicate identifiers')
 refs={};eligible=[];rejected=[]
 for (c,p),o in obs.groupby(['cell_type','plate'],sort=True):
  ctl=o[o.is_control].sort_values(['well','original_row'])
  if len(ctl)<4:
   rejected.append({'cell_type':c,'plate':p,'reason':'fewer_than_4_controls','n_controls':len(ctl)});continue
  cut=len(ctl)//2
  refs[c+'|'+p]={'A':ctl.original_row.iloc[:cut].astype(int).tolist(),'B':ctl.original_row.iloc[cut:].astype(int).tolist()}
  eligible.extend(o[~o.is_control].original_row.astype(int).tolist())
 tr=obs[obs.original_row.isin(eligible)&~obs.is_control]
 counts=tr.groupby(['cell_type','perturbagen']).size()
 valid=counts[counts>=2].reset_index(name='n_records')
 cells=valid.groupby('cell_type').size();cells=cells[cells>=50]
 contexts=sorted(cells.index,key=order_key)
 if len(contexts)<8:raise ValueError(f'Only {len(contexts)} eligible contexts; protocol requires >=8')
 n=min(6,max(2,math.ceil(.2*len(contexts))))
 test=contexts[:n];val=contexts[n:2*n];train=contexts[2*n:]
 if len(train)<4:raise ValueError('Insufficient training contexts')
 valid=valid[valid.cell_type.isin(contexts)]
 drugs=sorted(valid.perturbagen.unique())
 plan={'source_sha256':SOURCE_SHA,'contexts':contexts,'train':train,'validation':val,'test':test,'drugs':drugs,'genes':var.index.astype(str).tolist(),'refs':refs,'profile_counts':valid.to_dict('records'),'metadata_only':True,'held_target_treatments_used':False,'context_names_field':'cell_type','test_outcomes_scored':False}
 (out/'plan.json').write_text(json.dumps(plan,indent=2));valid.to_csv(out/'eligible_profiles.csv',index=False)
 pd.DataFrame(rejected).to_csv(out/'excluded_plates.csv',index=False)
 obs.to_csv(out/'metadata.csv.gz');var.to_csv(out/'gene_metadata.csv.gz')
 return plan

class Study:
 def __init__(self,root,plan):
  self.root=Path(root);self.plan=json.loads(Path(plan).read_text()) if isinstance(plan,(str,Path)) else plan
  self.obs=pd.read_csv(self.root/'slice_obs.csv.gz',index_col=0,dtype={'perturbagen':str,'cell_type':str,'plate':str,'well':str},low_memory=False)
  self.obs['is_control']=self.obs.is_control.astype(str).str.lower().isin(['true','1'])
  self.obs=self.obs.set_index('original_row',drop=False)
  self.drugs=self.plan['drugs'];self.genes=self.plan['genes'];self.parts=[];self.log=[];self.scoring=False
  manifest=json.loads((self.root/'acquisition.json').read_text())
  if manifest['sha256']!=SOURCE_SHA:raise ValueError('Wrong data source')
  for p in manifest['parts']:
   path=self.root/Path(p['file']).name
   if sha(path)!=p['sha256']:raise ValueError('Part checksum mismatch')
   z=np.load(path,allow_pickle=False);rr=z['original_rows'];x=z['X']
   if x.shape!=(len(rr),len(self.genes)) or not np.isfinite(x).all():raise ValueError('Invalid numeric part')
   self.parts.append((rr,x))
  self.lookup={int(r):(i,j) for i,(rr,_) in enumerate(self.parts) for j,r in enumerate(rr)}
  if set(self.lookup)!=set(self.obs.index):raise ValueError('Metadata/matrix row mismatch')
  self.valid={(r['cell_type'],r['perturbagen']) for r in self.plan['profile_counts']}
 def values(self,rows,role):
  rows=list(map(int,rows));meta=self.obs.loc[rows]
  if not self.scoring:
   excluded=meta.cell_type.isin(self.plan['test'])
   if (excluded&~meta.is_control).any():raise PermissionError('Held test treatment before joint freeze')
   b=set(r for k,v in self.plan['refs'].items() if k.split('|')[0] in self.plan['test'] for r in v['B'])
   if set(rows)&b:raise PermissionError('Held test scoring reference before freeze')
  self.log.append({'role':role,'rows':rows})
  return np.stack([self.parts[self.lookup[r][0]][1][self.lookup[r][1]] for r in rows]).astype(float)
 def baselines(self,contexts):
  b=[]
  for c in contexts:
   arrays=[self.values(v['A'],'baseline_A').mean(0) for k,v in self.plan['refs'].items() if k.split('|')[0]==c]
   if not arrays:raise ValueError('No usable controls')
   b.append(np.mean(arrays,0))
  return np.asarray(b)
 def effects(self,contexts,role='fit',control='B'):
  if role=='fit' and set(contexts)&set(self.plan['test']):raise PermissionError('Test context passed to fit')
  if role=='score' and not self.scoring:raise PermissionError('Scoring not enabled')
  if role not in ('fit','score'):raise ValueError('Unknown purpose')
  y=np.full((len(contexts),len(self.drugs),len(self.genes)),np.nan);records=[]
  for i,c in enumerate(contexts):
   o=self.obs[(self.obs.cell_type==c)&~self.obs.is_control]
   refs={k.split('|')[1]:self.values(v[control],'reference_'+control).mean(0) for k,v in self.plan['refs'].items() if k.split('|')[0]==c}
   for j,p in enumerate(self.drugs):
    if (c,p) not in self.valid:continue
    oo=o[(o.perturbagen==p)&o.plate.isin(refs)]
    x=self.values(oo.original_row,role+'_treatment')
    eff=x-np.stack([refs[plate] for plate in oo.plate])
    y[i,j]=eff.mean(0)
    records.append({'cell_type':c,'compound':p,'n_observations':len(oo),'n_plates':oo.plate.nunique(),'rows':oo.original_row.astype(int).tolist()})
  return y,records
 def allow_score(self,freeze_path):
  if not Path(freeze_path).is_file():raise ValueError('Freeze record missing')
  self.scoring=True
