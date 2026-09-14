"""Verified auxiliary phospho-signaling data; never silently substitutes for RNA."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

EXPECTED_BLOB = '89861234c3300d46b2e5573416ab910f8fa32030'
EXPECTED_SHA256 = 'da843e47b33b1a4485244efd59d61176fd7d0554089c3fe4785ebac12cdacc6c'
CUE_COLS = ['TR:IL1a','TR:TGFa','TR:IL6','TR:INS','TR:TNFa',
            'TR:PI3Ki','TR:MEK12i','TR:IKKi']
CONTEXTS = ['HepG2','PriHu','Hep3B','Huh7','Focus']

@dataclass(frozen=True)
class Panel:
    effects: np.ndarray  # context x condition x measured signal
    cues: np.ndarray    # condition x eight experimentally administered cues
    raw: np.ndarray
    controls: np.ndarray # context x two distinct untreated records x signal
    contexts: tuple[str, ...]
    conditions: tuple[str, ...]
    signals: tuple[str, ...]
    audit: dict

def load_panel(path: str | Path, control_mode: str = 'mean') -> Panel:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Measured dataset missing: {path}; no synthetic fallback.')
    b = path.read_bytes()
    blob = hashlib.sha1(b'blob ' + str(len(b)).encode() + b'\0' + b).hexdigest()
    if blob != EXPECTED_BLOB:
        raise ValueError('Input does not match the pinned original Git blob.')
    df = pd.read_csv(path)
    ccols = [f'TR:{c}:CellLine' for c in CONTEXTS]
    times = df.filter(regex='^DA:').to_numpy()
    signals = df.filter(regex='^DV:').columns.tolist()
    if df.shape != (480,52) or not np.all(times == times[:,[0]]):
        raise ValueError('Unexpected assay shape/time schema.')
    flags = df[ccols].fillna(0).to_numpy()
    if not np.all(flags.sum(1)==1): raise ValueError('Ambiguous context flags.')
    cues_all = df[CUE_COLS].fillna(0).to_numpy(dtype=int)
    if not np.isin(cues_all, [0,1]).all(): raise ValueError('Cues must be binary.')
    ctx_id = flags.argmax(1)
    all_signals = df[signals].to_numpy(dtype=float)
    if not np.isfinite(all_signals).all() or np.any(all_signals<0):
        raise ValueError('Invalid measured signal values.')
    # Primary auxiliary experiment uses ONLY ligand-present, time-code-30 rows.
    keep = (times[:,0]==30) & (cues_all[:,:5].sum(1)>0)
    keys_all = np.array([''.join(map(str,r)) for r in cues_all])
    condition_keys = tuple(sorted(set(keys_all[keep])))
    if len(condition_keys)!=72: raise ValueError('Expected 72 ligand-present conditions.')
    raw = np.empty((5,72,18)); controls = np.empty((5,2,18)); x = None
    duplicated_pairs = []
    for ci, context in enumerate(CONTEXTS):
        idx = np.where(keep & (ctx_id==ci))[0]
        if len(idx)!=72 or len(set(keys_all[idx]))!=72: raise ValueError('Condition duplication.')
        order = np.array([idx[np.where(keys_all[idx]==key)[0][0]] for key in condition_keys])
        raw[ci] = all_signals[order]
        xc = cues_all[order]
        if x is not None and not np.array_equal(x,xc): raise ValueError('Unaligned cues.')
        x = xc
        # NO inhibitor-treated time-zero sample is counted as an untreated control.
        ctrl_idx = np.where((ctx_id==ci) & (times[:,0]==0) & (cues_all.sum(1)==0))[0]
        if len(ctrl_idx)!=2: raise ValueError('Expected two distinct untreated records.')
        controls[ci] = all_signals[ctrl_idx]
        for row_id, block in df.loc[ctx_id==ci].groupby(df.columns[0]):
            if len(block)==2:
                if not np.array_equal(block[signals].iloc[0].to_numpy(),block[signals].iloc[1].to_numpy()):
                    raise ValueError('Duplicated ID has inconsistent measurements.')
                duplicated_pairs.append(str(row_id))
    logctrl = np.log2(controls+1)
    if control_mode=='mean': baseline=logctrl.mean(1)
    elif control_mode=='first': baseline=logctrl[:,0]
    elif control_mode=='second': baseline=logctrl[:,1]
    else: raise ValueError('control_mode must be mean, first or second.')
    effects = np.log2(raw+1)-baseline[:,None,:]
    audit = {
      'biological_modality':'auxiliary phosphoprotein signaling; NOT RNA',
      'raw_shape':list(df.shape), 'git_blob_sha1':blob,
      'sha256':hashlib.sha256(b).hexdigest(), 'bytes':len(b),
      'retained_ligand_rows':int(keep.sum()), 'conditions_per_context':72,
      'contexts':CONTEXTS, 'signal_count':18,
      'duplicate_id_pairs':duplicated_pairs, 'duplicate_id_pair_count':len(duplicated_pairs),
      'control_records_per_context':2, 'control_mode':control_mode,
      'effect_scale':'log2(signal+1) minus mean log2(untreated signal+1)',
      'technical_caveat':'No independent stimulated biological replicate is available per retained condition in this file. Time-0/30 duplicated rows are not independent repeats.',
      'analysis_status':'retrospective masked prediction, not prospective blinded validation'
    }
    return Panel(effects,x,raw,controls,tuple(CONTEXTS),condition_keys,
                 tuple(c.removeprefix('DV:') for c in signals),audit)
