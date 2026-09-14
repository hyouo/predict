"""Pinned OP3 count reader with explicit split and physical-control access gates.

This reader supports the prefiltered published panel only. It never opens private
expression. Donor IDs are not reconstructed. Integer counts are not individual cells.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import h5py
import numpy as np
import pandas as pd
from tools.op3_inventory import read_column

SHA256 = '8eaf7e63adc68029e88184726b8545bc3836c3ef2a25b9185dc4e8b819537e59'
SOURCES = ('CL_0000084', 'CL_0000623')
TARGETS = ('CL_0000236', 'CL_0000763')
NAMES = {'CL_0000084':'T', 'CL_0000623':'NK', 'CL_0000236':'B', 'CL_0000763':'Myeloid'}

class Reader:
    def __init__(self, path: Path, verify: bool = True):
        self.path = Path(path)
        if verify and hashlib.sha256(self.path.read_bytes()).hexdigest() != SHA256:
            raise ValueError('OP3 asset hash mismatch')
        with h5py.File(self.path, 'r') as f:
            self.obs = pd.DataFrame({k: read_column(v) for k,v in f['obs'].items()})
            self.var = pd.DataFrame({k: read_column(v) for k,v in f['var'].items()})
            self.shape = tuple(f['X'].attrs['shape'])
        if not self.obs.sample_id.is_unique or len(self.obs) != self.shape[0]:
            raise ValueError('Unexpected profile identity/shape')
        self.access = []

    def partition(self, rotation: int = 0) -> np.ndarray:
        if rotation not in (0,1,2):
            raise ValueError('Only registered control rotations are allowed')
        part = np.full(len(self.obs), -1, dtype=int)
        ctr = self.obs[self.obs.is_control.astype(bool)]
        # Physical wells receive identical assignments across all cell types.
        for plate, rows in ctr.groupby('plate', observed=True):
            wells = sorted(rows.well.unique())
            if len(wells) != 8:
                raise ValueError('Expected eight control wells per plate')
            assign = {w: (i % 3 - rotation) % 3 for i,w in enumerate(wells)}
            part[rows.index] = rows.well.map(assign).to_numpy()
        return part

    def read(self, ids, allowed, stage: str, denominator: str = 'psbulk_counts'):
        ids = np.asarray(ids, dtype=int)
        if ids.ndim != 1 or len(set(ids)) != len(ids) or not set(ids) <= set(allowed):
            raise PermissionError('Rows outside the explicit access grant')
        if np.any(self.obs.iloc[ids].split == 'private_test'):
            raise PermissionError('Private expression is sealed')
        if stage != 'evaluation' and np.any(self.obs.iloc[ids].split == 'public_test'):
            raise PermissionError('Public outcomes cannot enter training/prediction')
        if stage not in ('fit', 'evaluation') or denominator not in ('psbulk_counts','row_sum'):
            raise ValueError('Unregistered stage or denominator')
        x = np.zeros((len(ids), self.shape[1]), dtype=np.float64)
        with h5py.File(self.path, 'r') as f:
            ptr = f['X/indptr'][:]
            for j,i in enumerate(ids):
                lo,hi = int(ptr[i]),int(ptr[i+1])
                values = f['X/data'][lo:hi]
                if np.any(values < 0) or not np.isfinite(values).all():
                    raise ValueError('Invalid counts')
                x[j, f['X/indices'][lo:hi]] = values
        d = (self.obs.iloc[ids].psbulk_counts.to_numpy(dtype=float)
             if denominator == 'psbulk_counts' else x.sum(1))
        if np.any(d <= 0) or not np.isfinite(d).all():
            raise ValueError('Invalid library exposure')
        self.access.append({'stage':stage, 'row_ids':ids.tolist(),
                            'sample_ids':self.obs.iloc[ids].sample_id.tolist()})
        return np.log2(1 + 1e6*x/d[:,None])

@dataclass
class Bundle:
    reader: Reader
    compounds: list[str]
    source: np.ndarray
    source_validation: np.ndarray
    baseline: dict[str,np.ndarray]
    anchors: dict[str,np.ndarray]
    y_fit: dict[str,np.ndarray]
    y_validation: dict[str,np.ndarray]
    queries: dict[str,np.ndarray]
    evaluation_ids: dict[str,np.ndarray]
    partition: np.ndarray
    missing_source: list[dict]
    denominator: str


def _means(obs, x, ids, keys):
    table = obs.iloc[ids].copy()
    out = {}
    for key, ix in table.groupby(keys, observed=True, sort=True).indices.items():
        out[key] = x[np.asarray(ix)].mean(0)
    return out


def prepare(path, rotation=0, denominator='psbulk_counts') -> Bundle:
    r = Reader(path); o=r.obs; part=r.partition(rotation)
    compounds = sorted(o.loc[~o.is_control.astype(bool), 'perturbagen'].unique())
    lookup = {p:i for i,p in enumerate(compounds)}
    fit_ids = np.flatnonzero((o.split=='train') | (part==0) | (part==1))
    # Never grant final control partition C, public or private responses.
    log = r.read(fit_ids, fit_ids, 'fit', denominator)
    byid = {int(i):log[j] for j,i in enumerate(fit_ids)}
    ctrs = []
    for group in (0,1):
        ids = np.flatnonzero(part==group)
        ctrs.append(_means(o,np.stack([byid[int(i)] for i in ids]),ids,['cell_type','plate']))
    baseline={ct:np.stack([v for (c,p),v in ctrs[0].items() if c==ct]).mean(0)
              for ct in SOURCES+TARGETS}
    def effects(ct, ref):
        ids = np.flatnonzero((o.cell_type==ct)&(o.split=='train'))
        delta=np.stack([byid[int(i)]-ctrs[ref][(ct,o.iloc[i].plate)] for i in ids])
        return _means(o,delta,ids,'perturbagen')
    sm = [effects(ct,0) for ct in SOURCES]
    source=np.full((2,len(compounds),r.shape[1]),np.nan)
    missing=[]
    for s in range(2):
        for p,v in sm[s].items(): source[s,lookup[p]]=v
    for j,p in enumerate(compounds):
        good=np.isfinite(source[:,j,0])
        if not good.any(): raise ValueError('Query lacks all source measurements')
        for s in np.flatnonzero(~good):
            source[s,j]=source[good,j].mean(0)
            missing.append({'source':SOURCES[s],'compound':p,'rule':'other_observed_source'})
    source_validation=source.copy()
    for s,ct in enumerate(SOURCES):
        for p,v in effects(ct,1).items(): source_validation[s,lookup[p]]=v
    anchors={};yf={};yv={};queries={};eval_ids={}
    for ct in TARGETS:
        a,b=effects(ct,0),effects(ct,1)
        names=sorted(a)
        anchors[ct]=np.array([lookup[p] for p in names])
        yf[ct]=np.stack([a[p] for p in names]);yv[ct]=np.stack([b[p] for p in names])
        eval_ids[ct]=np.flatnonzero((o.cell_type==ct)&(o.split=='public_test'))
        queries[ct]=np.array([lookup[p] for p in sorted(o.iloc[eval_ids[ct]].perturbagen.unique())])
        if set(anchors[ct])&set(queries[ct]):raise ValueError('Anchor/query drug overlap')
    return Bundle(r,compounds,source,source_validation,baseline,anchors,yf,yv,queries,eval_ids,part,missing,denominator)


def evaluate_truth(bundle: Bundle, seal_path: Path):
    """Open public outcomes and scoring-only controls AFTER a prediction seal exists."""
    import json
    seal=json.loads(Path(seal_path).read_text())
    if not seal.get('predictions_complete') or not seal.get('arrays_sha256'):
        raise PermissionError('Prediction seal missing')
    arrays=Path(seal_path).parent / seal['arrays_file']
    if hashlib.sha256(arrays.read_bytes()).hexdigest()!=seal['arrays_sha256']:
        raise PermissionError('Predictions changed after seal')
    r=bundle.reader;o=r.obs;ctr=np.flatnonzero(bundle.partition==2)
    ids=np.concatenate([ctr]+list(bundle.evaluation_ids.values()))
    log=r.read(ids,ids,'evaluation',bundle.denominator)
    byid={int(i):log[j] for j,i in enumerate(ids)}
    cb=_means(o,np.stack([byid[int(i)] for i in ctr]),ctr,['cell_type','plate'])
    truth={};replicates={}
    for ct,rows in bundle.evaluation_ids.items():
        delta=np.stack([byid[int(i)]-cb[(ct,o.iloc[i].plate)] for i in rows])
        agg=_means(o,delta,rows,'perturbagen')
        truth[ct]=np.stack([agg[bundle.compounds[j]] for j in bundle.queries[ct]])
        replicates[ct]={'rows':rows,'effects':delta}
    return truth,replicates
