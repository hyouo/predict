"""Bounded first acquisition of public upstream OP3 metadata; no model training."""
from __future__ import annotations
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen, Request
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from tools.assets import ROOT, sha256_file
from tools.op3_inventory import read_column

URL = 'https://openproblems-data.s3.amazonaws.com/resources/task_perturbation_prediction/datasets/neurips-2023-data/pseudobulk_filtered_with_uns.h5ad'
CAP = 512 * 1024 * 1024
EXPECTED_SHA256 = "f9c33a0f06fae8c53b299d66e6b8926cb5d5117b96737b6dfce7ef718ceb79e2"

def main():
    out = ROOT / 'runs/op3-upstream-audit'
    out.mkdir(parents=True, exist_ok=False)
    record = {'source_url': URL, 'max_bytes': CAP, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'prior_checksum_available': True, 'expected_sha256': EXPECTED_SHA256, 'purpose': 'metadata_and_count_semantics_audit_only', 'trained_model': False}
    path = None
    try:
        with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as f:
            path = Path(f.name)
            req = Request(URL, headers={'User-Agent': 'predict-research-audit/2'})
            with urlopen(req, timeout=60) as r:
                size = 0
                while chunk := r.read(1024*1024):
                    size += len(chunk)
                    if size > CAP:
                        raise ValueError('Upstream file exceeds 512 MiB cap')
                    f.write(chunk)
        record.update(bytes=size, sha256=sha256_file(path))
        if record['sha256'] != EXPECTED_SHA256:
            raise ValueError('Upstream OP3 SHA256 changed; refusing silent dataset drift')
        with h5py.File(path) as f:
            obs = pd.DataFrame({k: read_column(v) for k,v in f['obs'].items()})
            var = pd.DataFrame({k: read_column(v) for k,v in f['var'].items()})
            obs.to_csv(out / 'upstream_obs.csv.gz', index=False)
            var.to_csv(out / 'upstream_var.csv.gz', index=False)
            x = f['X']
            if x.attrs.get('encoding-type') != 'csr_matrix':
                raise ValueError('Expected CSR upstream count matrix')
            matrix = csr_matrix((x['data'][:], x['indices'][:], x['indptr'][:]), shape=tuple(x.attrs['shape']))
            record.update(shape=list(matrix.shape), obs_columns=list(obs), var_columns=list(var),
                          count_dtype=str(matrix.dtype), min_count=float(matrix.data.min()), max_count=float(matrix.data.max()),
                          integer_valued=bool(np.equal(matrix.data, np.floor(matrix.data)).all()),
                          uns_keys=list(f.get('uns',{})))
            np.savez_compressed(out / 'upstream_counts.npz', data=matrix.data, indices=matrix.indices, indptr=matrix.indptr, shape=matrix.shape)
            if 'donor_id' in obs:
                cols=[k for k in ['plate_name','donor_id','library_id'] if k in obs]
                record['donor_mapping']=obs[cols].drop_duplicates().to_dict('records')
        record['component_sha256']={name:sha256_file(out/name) for name in ('upstream_obs.csv.gz','upstream_var.csv.gz','upstream_counts.npz')}
        record['completed']=True
    except Exception as e:
        record.update(completed=False, error=f'{type(e).__name__}: {e}')
        raise
    finally:
        record['finished_utc']=datetime.now(timezone.utc).isoformat()
        (out/'audit.json').write_text(json.dumps(record, indent=2, default=str)+'\n')
        if path is not None: path.unlink(missing_ok=True)
        print(json.dumps(record, indent=2, default=str))

if __name__ == '__main__': main()
