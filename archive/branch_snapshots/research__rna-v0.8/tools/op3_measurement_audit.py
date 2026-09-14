"""Reconcile count and donor identities. No response statistic is fitted here."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from tools.assets import ROOT,sha256_file,verify_file
from tools.op3_inventory import read_column
from tools.upstream_metadata import EXPECTED_SHA256
from src.rna_data import OP3_SHA


def reconcile(processed:Path,upstream:Path):
    verify_file(processed,OP3_SHA,32*1024**2)
    urec=json.loads((upstream/'audit.json').read_text())
    if urec.get('sha256')!=EXPECTED_SHA256 or not urec.get('completed'):raise ValueError('Unexpected upstream provenance')
    component_hashes={name:sha256_file(upstream/name) for name in ('upstream_obs.csv.gz','upstream_var.csv.gz','upstream_counts.npz')}
    if 'component_sha256' in urec and component_hashes!=urec['component_sha256']:raise ValueError('Upstream component changed')
    uo=pd.read_csv(upstream/'upstream_obs.csv.gz');uv=pd.read_csv(upstream/'upstream_var.csv.gz')
    with np.load(upstream/'upstream_counts.npz',allow_pickle=False) as n:
        ux=csr_matrix((n['data'],n['indices'],n['indptr']),shape=tuple(n['shape']))
    mp={'T cells':'CL_0000084','NK cells':'CL_0000623','B cells':'CL_0000236','Myeloid cells':'CL_0000763'}
    uo['standard_plate']=uo.plate_name.str.replace('-','_',regex=False);uo['standard_cell']=uo.cell_type.map(mp)
    uo['key']=uo.standard_plate+'|'+uo.well+'|'+uo.standard_cell
    if not uo.key.is_unique or uo.key.isna().any():raise ValueError('Upstream unit keys invalid')
    with h5py.File(processed) as f:
        o=pd.DataFrame({k:read_column(v) for k,v in f['obs'].items()});v=pd.DataFrame({k:read_column(v) for k,v in f['var'].items()})
        a=f['X'];x=csr_matrix((a['data'][:],a['indices'][:],a['indptr'][:]),shape=tuple(a.attrs['shape']))
    o['key']=o.plate+'|'+o.well+'|'+o.cell_type
    if not o.key.is_unique:raise ValueError('Processed unit keys invalid')
    rowids=uo.reset_index().set_index('key').loc[o.key,'index'].to_numpy();match=uo.iloc[rowids].reset_index(drop=True)
    if not np.array_equal(match.sm_name.to_numpy(),o.perturbagen.to_numpy()) or not np.array_equal(match.split.to_numpy(),o.split.to_numpy()):raise ValueError('Drug/split identity changed')
    if not uv['_index'].is_unique:raise ValueError('Original gene symbols are ambiguous')
    gids=uv.reset_index().set_index('_index')['index'];exact=v.symbol.isin(gids.index)&(~v.is_merged.astype(bool))
    if not exact.all():raise ValueError('Merged or unmatched genes need a separate reconciliation')
    residual=x-ux[rowids][:,gids.loc[v.symbol].to_numpy()]
    rawsum=np.asarray(ux[rowids].sum(axis=1)).ravel();retsum=np.asarray(x.sum(axis=1)).ravel()
    if not np.array_equal(rawsum,o.psbulk_counts.to_numpy()) or residual.nnz:raise ValueError('Counts or denominator mismatch')
    if np.any(x.data<0) or not np.isfinite(x.data).all() or not np.equal(x.data,np.floor(x.data)).all():raise ValueError('Not nonnegative integer counts')
    if match.groupby('standard_plate').donor_id.nunique().ne(1).any():raise ValueError('Plate has multiple donors')
    donors=match.set_index('standard_plate').donor_id.to_dict()
    return {'upstream_file_sha256':EXPECTED_SHA256,'processed_sha256':OP3_SHA,'component_sha256':component_hashes,
            'component_hashes_recorded_at_download':'component_sha256' in urec,
            'upstream_shape':list(ux.shape),'processed_shape':list(x.shape),'mapped_rows':len(o),
            'exactly_matched_genes':len(v),'nonzero_count_differences':int(residual.nnz),
            'psbulk_counts_equals_original_gene_sum':True,'nonnegative_integer_counts':True,
            'donor_mapping':donors,'retained_count_fraction_quantiles':np.quantile(retsum/rawsum,[0,.25,.5,.75,1]).tolist(),
            'gene_universe_limitation':'Upstream filter used all treatments and contexts. Conditional-on-release analysis only.',
            'unit':'Plate/well/cell-type raw-count pseudobulk according to original pipeline, confirmed by full entrywise reconciliation',
            'approved_scope':'descriptive log2CPM prediction under frozen sample/control roles',
            'approved_for_confirmatory_sota_claim':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--upstream',type=Path,default=ROOT/'runs/op3-upstream-audit')
    p.add_argument('--processed',type=Path,default=ROOT/'data/rna/op3_standardized_processed.h5ad');a=p.parse_args()
    if a.output.exists():raise FileExistsError('Do not overwrite an audit')
    audit=reconcile(a.processed,a.upstream);a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
