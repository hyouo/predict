"""Row-allowlisted access to the pinned OP3 pseudobulk counts.

The numerical endpoint is explicitly descriptive log2(1+CPM), not limma logFC.
The store has no generic public 'read all expression' method. Fit construction
cannot access evaluation controls or target held-out treatments.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import h5py
import numpy as np
import pandas as pd
from tools.assets import verify_file
from tools.op3_inventory import read_column

SOURCE_CONTEXTS = ('CL_0000084', 'CL_0000623')
TARGET_CONTEXTS = ('CL_0000236', 'CL_0000763')
NAMES = {'CL_0000084':'T', 'CL_0000623':'NK', 'CL_0000236':'B', 'CL_0000763':'Myeloid'}
OP3_SHA = '8eaf7e63adc68029e88184726b8545bc3836c3ef2a25b9185dc4e8b819537e59'

@dataclass
class SourceData:
    drugs: list[str]
    means: np.ndarray               # source context, drug, gene
    donor_effects: np.ndarray       # source context, donor, drug, gene; NaN is missing
    baselines: np.ndarray           # source context, gene
    donors: tuple[str, ...]
    omitted_source_compounds: tuple[str, ...]=()

@dataclass
class AnchorData:
    indices: np.ndarray
    reference_c: np.ndarray
    reference_d: np.ndarray
    mean: np.ndarray
    donor_effects: np.ndarray
    baseline: np.ndarray
    observation_count: int

class OP3Store:
    def __init__(self, path: Path, donor_map: dict[str,str], *, verify: bool=True,
                 denominator: str='psbulk_counts', pseudocount: float=1.0):
        self.path=Path(path)
        if verify: verify_file(self.path, OP3_SHA, 32*1024*1024)
        if denominator not in ('psbulk_counts','retained_sum'):
            raise ValueError('Unsupported denominator')
        if not np.isfinite(pseudocount) or pseudocount<=0:
            raise ValueError('Positive finite CPM pseudocount required')
        self.denominator=denominator; self.pseudocount=float(pseudocount)
        with h5py.File(self.path) as f:
            self.obs=pd.DataFrame({k:read_column(v) for k,v in f['obs'].items()})
            self.var=pd.DataFrame({k:read_column(v) for k,v in f['var'].items()})
            self.shape=tuple(int(v) for v in f['X'].attrs['shape'])
        o=self.obs
        required={'sample_id','cell_type','is_control','plate','well','perturbagen','split','psbulk_counts','psbulk_cells','pubchem_cid'}
        if not required <= set(o) or len(o)!=self.shape[0]: raise ValueError('Unexpected metadata')
        if o.sample_id.isna().any() or not o.sample_id.is_unique: raise ValueError('Invalid sample IDs')
        if o.is_control.isna().any(): raise ValueError('Missing controls')
        o['is_control']=o.is_control.astype(bool)
        o['donor']=o.plate.map(donor_map)
        if o.donor.isna().any(): raise ValueError('Unverified/missing donor mapping')
        o['well_row']=o.well.str[0]
        self.access: list[dict] = []
        self._cache={}

    def _read(self, rows, role: str, *, evaluation_split: str | None=None):
        rows=np.asarray(rows, dtype=int)
        if rows.ndim!=1 or np.any(rows<0) or np.any(rows>=len(self.obs)): raise ValueError('Bad row IDs')
        o=self.obs.iloc[rows]; ctrl=o.is_control.to_numpy(); src=o.cell_type.isin(SOURCE_CONTEXTS).to_numpy()
        allowed = {
            'source_treatment': (~ctrl)&src&(o['split'].to_numpy()=='train'),
            'source_control': ctrl&src,
            'baseline': ctrl&o.well_row.isin(list('AB')).to_numpy(),
            'anchor_control': ctrl&o.well_row.isin(list('CD')).to_numpy(),
            'anchor_treatment': (~ctrl)&(~src)&(o['split'].to_numpy()=='train'),
            'evaluation_control': ctrl&o.well_row.isin(list('EFGH')).to_numpy(),
            'evaluation_treatment': (~ctrl)&(~src)&(o['split'].to_numpy()==evaluation_split)
        }
        if role not in allowed or not allowed[role].all(): raise PermissionError(f'Rows not permitted for {role}')
        if role.startswith('evaluation_') and evaluation_split not in ('public_test','private_test'):
            raise PermissionError('Evaluation phase required')
        ans=np.empty((len(rows),self.shape[1]),dtype=float)
        with h5py.File(self.path) as f:
            x=f['X']; ptr=x['indptr']
            for j,i in enumerate(rows):
                key=(int(i),self.denominator,self.pseudocount)
                if key in self._cache: ans[j]=self._cache[key]; continue
                a,b=int(ptr[i]),int(ptr[i+1]); indices=x['indices'][a:b]; vals=x['data'][a:b]
                if not np.isfinite(vals).all() or np.any(vals<0): raise ValueError('Invalid counts')
                total=float(self.obs.iloc[i].psbulk_counts) if self.denominator=='psbulk_counts' else float(vals.sum())
                if not np.isfinite(total) or total<=0: raise ValueError('Invalid library total')
                row=np.full(self.shape[1],np.log2(self.pseudocount),dtype=float)
                row[indices]=np.log2(self.pseudocount+1e6*vals/total)
                self._cache[key]=row; ans[j]=row
        self.access.append({'role':role,'sample_ids':o.sample_id.tolist()})
        return ans

    def _controls(self, context: str, donors, letters: str, role: str, evaluation_split=None):
        o=self.obs
        mask=(o.cell_type==context)&o.donor.isin(donors)&o.is_control&o.well_row.isin(list(letters))
        rows=np.flatnonzero(mask)
        x=self._read(rows,role,evaluation_split=evaluation_split)
        meta=o.iloc[rows]
        result={p:x[(meta.plate==p).to_numpy()].mean(axis=0) for p in meta.plate.unique()}
        return result, x, meta

    def source(self, donors=('Donor 1','Donor 2','Donor 3')) -> SourceData:
        donors=tuple(donors);o=self.obs
        available=[set(o.loc[(o.cell_type==c)&o.donor.isin(donors)&(~o.is_control)&(o['split']=='train'),'perturbagen']) for c in SOURCE_CONTEXTS]
        drugs=sorted(set.intersection(*available))
        omitted=tuple(sorted(set.union(*available)-set(drugs)))
        drug_index={d:i for i,d in enumerate(drugs)}
        means=[]; all_reps=[]; baselines=[]
        for context in SOURCE_CONTEXTS:
            controls,b,_=self._controls(context,donors,'ABCDEFGH','source_control')
            mask=(o.cell_type==context)&o.donor.isin(donors)&(~o.is_control)&(o['split']=='train')&o.perturbagen.isin(drugs)
            rows=np.flatnonzero(mask); z=self._read(rows,'source_treatment')
            tensor=np.full((len(donors),len(drugs),self.shape[1]),np.nan)
            for j,i in enumerate(rows):
                r=o.iloc[i]; ix=(donors.index(r.donor),drug_index[r.perturbagen])
                if np.isfinite(tensor[ix]).any(): raise ValueError('Duplicate donor/compound observation')
                tensor[ix]=z[j]-controls[r.plate]
            counts=np.isfinite(tensor[...,0]).sum(axis=0)
            # Do not silently impute a completely unavailable condition.
            mean=np.divide(np.nansum(tensor,axis=0),counts[:,None],out=np.full(tensor.shape[1:],np.nan),where=counts[:,None]>0)
            means.append(mean);all_reps.append(tensor);baselines.append(b.mean(axis=0))
        return SourceData(drugs,np.asarray(means),np.asarray(all_reps),np.asarray(baselines),donors,omitted)

    def baseline(self, context, donors):
        return self._controls(context,donors,'AB','baseline')[1].mean(axis=0)

    def anchors(self, context, source: SourceData) -> AnchorData:
        if context not in TARGET_CONTEXTS: raise ValueError('Unknown target context')
        o=self.obs; donors=source.donors;lookup={d:i for i,d in enumerate(source.drugs)}
        controls_c,_,_=self._controls(context,donors,'C','anchor_control')
        controls_d,_,_=self._controls(context,donors,'D','anchor_control')
        mask=(o.cell_type==context)&o.donor.isin(donors)&(~o.is_control)&(o['split']=='train')&o.perturbagen.isin(source.drugs)
        rows=np.flatnonzero(mask);x=self._read(rows,'anchor_treatment')
        drugs=sorted(o.loc[mask,'perturbagen'].unique());dix={d:i for i,d in enumerate(drugs)}
        xc=np.full((len(donors),len(drugs),self.shape[1]),np.nan);xd=xc.copy()
        for j,i in enumerate(rows):
            r=o.iloc[i];ix=(donors.index(r.donor),dix[r.perturbagen])
            if np.isfinite(xc[ix]).any(): raise ValueError('Duplicate anchor unit')
            xc[ix]=x[j]-controls_c[r.plate];xd[ix]=x[j]-controls_d[r.plate]
        count=np.isfinite(xc[...,0]).sum(axis=0)
        yc=np.nansum(xc,axis=0)/count[:,None]; yd=np.nansum(xd,axis=0)/count[:,None]
        ids=np.array([lookup[d] for d in drugs]);keep=np.isfinite(source.means[:,ids,:]).all(axis=(0,2))
        if not keep.all(): raise ValueError('Anchor source coverage missing; explicit protocol amendment needed')
        return AnchorData(ids,yc,yd,(yc+yd)/2,(xc+xd)/2,self.baseline(context,donors),len(rows))

    def query_indices(self, context, split, source: SourceData, evaluation_donors=None):
        o=self.obs; donors=source.donors if evaluation_donors is None else tuple(evaluation_donors)
        names=set(o.loc[(o.cell_type==context)&o.donor.isin(donors)&(~o.is_control)&(o['split']==split),'perturbagen'])
        available=np.isfinite(source.means).all(axis=(0,2))
        keep=np.array([i for i,d in enumerate(source.drugs) if d in names and available[i]],dtype=int)
        omitted=sorted(names-set(source.drugs[i] for i in keep))
        return keep,omitted

    def truth(self, context, split, source: SourceData, *, prediction_lock: Path, evaluation_donors=None):
        """Call only after the prediction lock file has been written and hashed."""
        if not Path(prediction_lock).is_file(): raise PermissionError('Predictions must be sealed first')
        lock=json.loads(Path(prediction_lock).read_text())
        if not lock.get('predictions_sealed_before_evaluation'): raise PermissionError('Invalid prediction lock')
        o=self.obs;donors=source.donors if evaluation_donors is None else tuple(evaluation_donors)
        ids,omitted=self.query_indices(context,split,source,donors);drugs=[source.drugs[i] for i in ids]
        controls,_,_=self._controls(context,donors,'EFGH','evaluation_control',split)
        mask=(o.cell_type==context)&o.donor.isin(donors)&(~o.is_control)&(o['split']==split)&o.perturbagen.isin(drugs)
        rows=np.flatnonzero(mask);x=self._read(rows,'evaluation_treatment',evaluation_split=split)
        tensor=np.full((len(donors),len(drugs),self.shape[1]),np.nan);di={d:i for i,d in enumerate(drugs)}
        for j,i in enumerate(rows):
            r=o.iloc[i];tensor[donors.index(r.donor),di[r.perturbagen]]=x[j]-controls[r.plate]
        n=np.isfinite(tensor[...,0]).sum(axis=0)
        return ids,np.nansum(tensor,axis=0)/n[:,None],tensor,omitted
