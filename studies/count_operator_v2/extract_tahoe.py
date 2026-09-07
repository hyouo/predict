"""Pinned public Tahoe extraction; panel uses training-cell controls only."""
from __future__ import annotations
import argparse, datetime, hashlib, json, pathlib, urllib.request
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
REV='6c54c8eb0321cceff4f888c54e199077d055e20b'
HASH='09ac7c8f63a77dc730156c78e361273c3369959640120e5c54241f1a6f82d7b2'
SIZE=4473155146
URL=f'https://huggingface.co/datasets/theislab/chem-perturbridge/resolve/{REV}/tahoe/tahoe_processed.h5ad?download=true'
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def col(node):
    if isinstance(node,h5py.Dataset):
        return node.asstr()[:] if h5py.check_string_dtype(node.dtype) else node[:]
    if 'categories' in node:
        c=col(node['categories']);return np.array([c[i] if i>=0 else None for i in node['codes'][:]],object)
    if 'values' in node:
        a=np.asarray(col(node['values']),object)
        if 'mask' in node:a[node['mask'][:]]=None
        return a
    raise ValueError('Unsupported metadata encoding '+node.name)
def acquire(p):
    if p.exists():
        if p.stat().st_size!=SIZE or sha(p)!=HASH:raise ValueError('Existing asset mismatch')
        return
    tmp=p.with_suffix('.part');h=hashlib.sha256();n=0
    try:
        with urllib.request.urlopen(URL,timeout=120) as r,tmp.open('wb') as f:
            while b:=r.read(8*1024*1024):
                n+=len(b)
                if n>4500000000:raise ValueError('Download size limit')
                h.update(b);f.write(b)
        if n!=SIZE or h.hexdigest()!=HASH:raise ValueError('Asset size/SHA mismatch')
        tmp.rename(p)
    except Exception:
        tmp.unlink(missing_ok=True);raise

def extract(p,out):
    out.mkdir(parents=True,exist_ok=False)
    start=datetime.datetime.now(datetime.timezone.utc).isoformat()
    if p.stat().st_size!=SIZE or sha(p)!=HASH:raise ValueError('Data verification failed')
    with h5py.File(p) as f:
        obs=pd.DataFrame({k:col(v) for k,v in f['obs'].items()})
        var=pd.DataFrame({k:col(v) for k,v in f['var'].items()})
        genes=var['_index'].astype(str).to_numpy()
        if len(set(genes))!=len(genes):raise ValueError('Duplicate genes')
        cells=sorted(obs.cell_type.unique(),key=lambda c:hashlib.sha256(('predict-tahoe-v2-20260907|'+c).encode()).hexdigest())
        if len(cells)!=50:raise ValueError('Expected 50 contexts')
        ci={c:i for i,c in enumerate(cells)}
        x=f['X'];ptr=x['indptr'][:];shape=tuple(x.attrs['shape'])
        if shape!=(67018,61483):raise ValueError('Unexpected matrix shape')
        def rows(lo,hi):
            a,b=int(ptr[lo]),int(ptr[hi])
            return csr_matrix((x['data'][a:b],x['indices'][a:b],ptr[lo:hi+1]-a),shape=(hi-lo,shape[1]))
        fullctrl=np.zeros((50,2,shape[1]),np.int64);mappings=[]
        for (cell,plate),g in obs[obs.is_control].groupby(['cell_type','plate'],sort=True):
            g=g.sort_values('sample_id')
            if len(g)!=2:raise ValueError('Exactly two controls per cell/plate required')
            for role,(i,row) in enumerate(g.iterrows()):
                v=rows(i,i+1).toarray()[0]
                if np.any(v<0) or np.any(v!=np.floor(v)):raise ValueError('Expected nonnegative integer counts')
                fullctrl[ci[cell],role]+=v.astype(np.int64)
                mappings.append({'row':int(i),'cell_type':cell,'plate':plate,'role':'AB'[role],'sample_id':row.sample_id,'count_sum':int(v.sum())})
        train=fullctrl[:32,0].astype(float)
        cp=train/train.sum(1,keepdims=True)*1e6
        ok=np.where(cp.mean(0)>=1)[0];vv=np.var(np.log2(cp+1),axis=0)
        panel=ok[np.lexsort((genes[ok],-vv[ok]))[:3000]]
        if len(panel)!=3000:raise ValueError('Insufficient expressed genes')
        pd.DataFrame({'gene':genes[panel],'source_variance':vv[panel],'source_mean_cpm':cp.mean(0)[panel]}).to_csv(out/'gene_panel.csv',index=False)
        cp=None;train=None
        ctrl=np.concatenate((fullctrl[:,:,panel],(fullctrl.sum(2)-fullctrl[:,:,panel].sum(2))[:,:,None]),2)
        np.savez_compressed(out/'controls.npz',counts=ctrl,cells=np.array(cells),genes=np.r_[genes[panel],'__OTHER__'])
        fullctrl=None
        mask=(~obs.is_control)&np.isclose(obs.pert_dose_uM,5)&np.isclose(obs.pert_time_h,24)
        presence=[set(obs.loc[mask&(obs.cell_type==c),'perturbagen']) for c in cells[:32]]
        drugs=sorted(set.intersection(*presence))
        if len(drugs)!=379:raise ValueError('Unexpected training drug coverage')
        di={d:i for i,d in enumerate(drugs)}
        values=np.zeros((50,len(drugs),3001),np.int64)
        libraries=np.zeros((50,len(drugs)),np.int64);nrec=libraries.copy();ncells=libraries.copy()
        rowsout=[];max_total_difference=0
        for lo in range(0,len(obs),512):
            hi=min(lo+512,len(obs));local=obs.iloc[lo:hi]
            selected=np.where(mask.iloc[lo:hi].to_numpy()&local.perturbagen.isin(drugs).to_numpy())[0]
            if not len(selected):continue
            block=rows(lo,hi)[selected]
            if np.any(block.data<0) or np.any(block.data!=np.floor(block.data)):raise ValueError('Invalid raw count entries')
            totals=np.asarray(block.sum(1)).ravel().astype(np.int64)
            kept=block[:,panel].toarray().astype(np.int64);packed=np.c_[kept,totals-kept.sum(1)]
            if np.any(packed<0):raise ValueError('OTHER bin invalid')
            for j,relative in enumerate(selected):
                i=lo+int(relative);r=obs.iloc[i];c,d=ci[r.cell_type],di[r.perturbagen]
                values[c,d]+=packed[j];libraries[c,d]+=totals[j];nrec[c,d]+=1;ncells[c,d]+=int(r.psbulk_cells)
                max_total_difference=max(max_total_difference,abs(int(totals[j])-int(r.psbulk_counts)))
                rowsout.append({'row':i,'cell_type':r.cell_type,'drug':r.perturbagen,'plate':r.plate,'well':r.well,'sample_id':r.sample_id,'count_sum':int(totals[j]),'psbulk_counts':int(r.psbulk_counts),'psbulk_cells':int(r.psbulk_cells)})
            if lo%8192==0:print('Rows extracted',hi,flush=True)
        if not np.array_equal(values.sum(2),libraries):raise ValueError('Count conservation failed')
        for role,ix in [('train',range(32)),('validation',range(32,41)),('test',range(41,50))]:
            ix=list(ix)
            if np.any(nrec[ix]==0):raise ValueError('Missing target condition; do not impute')
            dest=out/('SEALED/test_counts.npz' if role=='test' else role+'_counts.npz')
            dest.parent.mkdir(exist_ok=True)
            np.savez_compressed(dest,counts=values[ix],libraries=libraries[ix],n_records=nrec[ix],n_cells=ncells[ix],cells=np.array(cells)[ix],drugs=np.array(drugs),genes=np.r_[genes[panel],'__OTHER__'])
        pd.DataFrame(rowsout).to_csv(out/'treatment_rows.csv.gz',index=False)
        pd.DataFrame(mappings).to_csv(out/'control_rows.csv',index=False)
        info={'input_sha256':HASH,'input_bytes':SIZE,'started_utc':start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_revision':REV,'raw_shape':[int(v) for v in shape],'train_cells':cells[:32],'validation_cells':cells[32:41],'test_cells':cells[41:],'drug_ids':drugs,'gene_selection':'TRAIN A controls only; log2(CPM+1) cross-cell variance, mean CPM >=1','n_genes':3000,'includes_OTHER':True,'counts_conserved':True,'max_recorded_total_difference':max_total_difference,'test_treatments_used_for_panel_or_model':False,'expression_extraction_not_model_training':True,'files':{str(q.relative_to(out)):{'sha256':sha(q),'bytes':q.stat().st_size} for q in out.rglob('*') if q.is_file()}}
        (out/'extraction_manifest.json').write_text(json.dumps(info,indent=2))
        print(json.dumps({k:v for k,v in info.items() if k!='drug_ids'}),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',type=pathlib.Path,required=True);a.add_argument('--out',type=pathlib.Path,required=True);a.add_argument('--download',action='store_true');args=a.parse_args()
    if args.download:acquire(args.input)
    extract(args.input,args.out)
