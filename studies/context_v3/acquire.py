"""Bounded public Tahoe pseudobulk acquisition; no model fitting or test scoring."""
from pathlib import Path
import urllib.request, hashlib, json, time
import h5py
import numpy as np
import pandas as pd
REV='6c54c8eb0321cceff4f888c54e199077d055e20b'
SHA='09ac7c8f63a77dc730156c78e361273c3369959640120e5c54241f1a6f82d7b2'
SIZE=4473155146
URL=f'https://huggingface.co/datasets/theislab/chem-perturbridge/resolve/{REV}/tahoe/tahoe_processed.h5ad'

def readcol(v):
    if isinstance(v,h5py.Dataset):
        return v.asstr()[:] if h5py.check_string_dtype(v.dtype) else v[:]
    if set(['codes','categories']) <= set(v):
        cat=readcol(v['categories']); ix=v['codes'][:]
        return np.array([cat[i] if i>=0 else None for i in ix],dtype=object)
    raise ValueError('Unsupported metadata encoding')

def order(xs,salt):
    return sorted(xs,key=lambda s:hashlib.sha256((salt+str(s)).encode()).hexdigest())

def subset_x(f,rows,cols):
    x=f['X'];cols=np.asarray(cols,int)
    if isinstance(x,h5py.Dataset):
        return np.stack([np.asarray(x[i,:],dtype=np.float32)[cols] for i in rows])
    if x.attrs.get('encoding-type')!='csr_matrix':raise ValueError('Expected CSR or dense X')
    ptr=x['indptr'][:];g=int(x.attrs['shape'][1]);m=np.full(g,-1,dtype=np.int32);m[cols]=np.arange(len(cols))
    out=np.zeros((len(rows),len(cols)),np.float32)
    for j,i in enumerate(rows):
        lo,hi=int(ptr[i]),int(ptr[i+1]);ids=x['indices'][lo:hi];keep=m[ids]>=0
        out[j,m[ids[keep]]]=x['data'][lo:hi][keep]
    return out

def main():
    start=time.time(); out=Path('acquired_v3');out.mkdir(exist_ok=False);p=Path('tahoe_source.h5ad')
    h=hashlib.sha256();total=0
    req=urllib.request.Request(URL,headers={'User-Agent':'predict-research-public-asset/3'})
    with urllib.request.urlopen(req,timeout=180) as r,p.open('wb') as f:
        while True:
            b=r.read(8*1024*1024)
            if not b:break
            total+=len(b)
            if total>SIZE:raise RuntimeError('Expected size exceeded')
            h.update(b);f.write(b)
    if total!=SIZE or h.hexdigest()!=SHA:raise RuntimeError('Public file hash/size mismatch')
    print('SOURCE_VERIFIED',total,SHA,flush=True)
    with h5py.File(p,'r',rdcc_nbytes=64*1024*1024) as f:
        obs=pd.DataFrame({k:readcol(v) for k,v in f['obs'].items()})
        var=pd.DataFrame({k:readcol(v) for k,v in f['var'].items()})
        gene_key=f['var'].attrs.get('_index','_index');genes=np.asarray(readcol(f['var'][gene_key]),dtype=str)
        if len(set(genes))!=len(genes):raise ValueError('Nonunique genes')
        obs.insert(0,'source_row',np.arange(len(obs)))
        obs.to_csv(out/'full_obs.csv.gz',index=False)
        ctrl=obs.is_control.astype(str).str.lower().isin(['true','1'])
        usable=~ctrl & (obs.psbulk_cells.astype(float)>=20)
        times=obs.loc[usable,'pert_time_h'].value_counts()
        tm=24. if 24. in times.index else float(times.index[0])
        usable &= np.isclose(obs.pert_time_h.astype(float),tm)
        controls=ctrl & np.isclose(obs.pert_time_h.astype(float),tm)
        keys=['cell_type','plate']
        ctrl_n=obs[controls].groupby(keys,observed=True).size()
        pairable=set(ctrl_n[ctrl_n>=2].index)
        paired=np.array([tuple(v) in pairable for v in obs[keys].itertuples(index=False,name=None)])
        usable &= paired; controls &= paired
        cell_n=obs[usable].groupby('cell_type',observed=True).perturbagen.nunique()
        cells=order(cell_n[cell_n>=100].index,'context-v3-split-20260908:')
        if len(cells)<24:raise ValueError(f'Only {len(cells)} eligible contexts')
        split={'test':cells[:8],'validation':cells[8:16],'train':cells[16:]}
        train_rows=usable & obs.cell_type.isin(split['train'])
        dose_counts=obs[train_rows].groupby('pert_dose_uM',observed=True).apply(lambda q:len(q[['cell_type','perturbagen']].drop_duplicates()),include_groups=False)
        dose=float(dose_counts.sort_values(ascending=False,kind='stable').index[0])
        usable &= np.isclose(obs.pert_dose_uM.astype(float),dose)
        coverage=obs[usable & obs.cell_type.isin(split['train'])].groupby('perturbagen',observed=True).cell_type.nunique()
        drugs=order(coverage[coverage>=max(4,int(np.ceil(.8*len(split['train']))))].index,'context-v3-drugs-20260908:')[:256]
        if len(drugs)<50:raise ValueError(f'Only {len(drugs)} broad-support compounds at mode dose')
        cr=np.flatnonzero(controls & obs.cell_type.isin(split['train']))
        counts=subset_x(f,cr,np.arange(len(genes)))
        libs=obs.iloc[cr].psbulk_counts.to_numpy(float)
        if np.any(libs<=0) or np.any(counts<0) or not np.isfinite(counts).all():raise ValueError('Invalid source control matrix')
        means=(counts/libs[:,None]*1e6).mean(0)
        candidates=np.flatnonzero(means>=1.)
        panel=candidates[np.argsort(-means[candidates],kind='stable')[:2000]]
        del counts
        selected=(controls | (usable & obs.perturbagen.isin(drugs))) & obs.cell_type.isin(cells)
        rows=np.flatnonzero(selected)
        counts=subset_x(f,rows,panel)
        np.savez_compressed(out/'panel_counts.npz',counts=counts,genes=genes[panel],source_rows=rows)
        obs.iloc[rows].to_csv(out/'obs.csv',index=False)
        var.iloc[panel].to_csv(out/'var.csv',index=False)
        design={'schema':1,'input_url':URL,'revision':REV,'input_sha256':SHA,'input_bytes':total,
                'source_shape':[len(obs),len(genes)],'export_shape':list(counts.shape),
                'split':split,'compounds':drugs,'dose_uM':dose,'time_h':tm,'min_cells':20,
                'gene_selection':'top 2000 mean CPM>=1 using only training-context controls',
                'gene_source_controls':obs.iloc[cr].source_row.tolist(),
                'dose_rule':'largest number of training-context compound pairs; metadata only',
                'cohort_rule':'at least100 compounds per context before dose restriction; at least2 time-matched controls per context/plate; compound in >=80% training contexts; first256 by fixed hash',
                'split_rule':'sha256(context-v3-split-20260908:+cell_id), first8 test, next8 validation, remaining train',
                'study_license':'CC0 Tahoe; CC-BY-4.0 Chem-PerturBridge packaging; retain citations',
                'counts_are_pseudobulk_not_single_cells':True,'test_outcomes_used_for_selection':False,
                'elapsed_seconds':time.time()-start}
        design['outputs_sha256']={a.name:hashlib.file_digest(a.open('rb'),'sha256').hexdigest() for a in out.iterdir() if a.is_file()}
        (out/'design.json').write_text(json.dumps(design,indent=2))
        print(json.dumps({k:v for k,v in design.items() if k not in ['compounds','gene_source_controls']},indent=2),flush=True)
    p.unlink()
if __name__=='__main__':main()
