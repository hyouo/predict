"""Metadata-only inventory of the pinned OP3 asset; no expression outcomes are read."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from tools.assets import ROOT, inspect_h5ad, verify_file


def read_column(node) -> np.ndarray:
    """Decode dense or AnnData categorical columns, preserving missing categories."""
    if isinstance(node, h5py.Dataset):
        if node.ndim != 1:
            raise ValueError('Observation columns must be one-dimensional')
        return node.asstr()[:] if h5py.check_string_dtype(node.dtype) else node[:]
    if not isinstance(node, h5py.Group) or not {'categories', 'codes'} <= set(node):
        raise ValueError('Unsupported observation column encoding; refusing to guess')
    categories = read_column(node['categories'])
    codes = read_column(node['codes'])
    if not np.issubdtype(codes.dtype, np.integer) or np.any(codes < -1) or np.any(codes >= len(categories)):
        raise ValueError('Invalid categorical indices')
    out = np.empty(len(codes), dtype=object)
    out[:] = None
    known = codes >= 0
    out[known] = categories[codes[known]]
    return out


def inventory(path: Path) -> dict:
    asset = json.loads((ROOT / 'assets/manifest.json').read_text())['assets']['op3']
    digest = verify_file(path, asset['sha256'], asset['max_bytes'])
    structure = inspect_h5ad(path)
    with h5py.File(path, 'r') as f:
        obs = pd.DataFrame({name: read_column(node) for name, node in f['obs'].items()})
    required = {'sample_id', 'cell_type', 'is_control', 'split', 'plate', 'well',
                'perturbagen', 'pert_time_h', 'pert_dose_uM', 'psbulk_cells', 'psbulk_counts'}
    if not required <= set(obs) or len(obs) != structure['X']['shape'][0]:
        raise ValueError('Unexpected OP3 metadata schema or row count')
    if obs['is_control'].isna().any() or not obs['is_control'].map(lambda x: isinstance(x, (bool, np.bool_))).all():
        raise ValueError('Control status must be explicitly Boolean and complete')
    if obs['sample_id'].isna().any() or not obs['sample_id'].is_unique:
        raise ValueError('Sample IDs must be complete and unique')
    controls = obs['is_control'].astype(bool)
    columns = ['cell_type', 'is_control', 'split', 'plate', 'pert_time_h', 'pert_dose_uM']
    counts = {name: {str(k): int(v) for k, v in obs[name].value_counts(dropna=False).items()}
              for name in columns}
    return {'asset': 'op3', 'source_url': asset['url'], **digest,
            'metadata_only': True, 'expression_values_read': False,
            'approved_for_training': False, 'rna_training_performed': False,
            'structure': structure, 'observation_count': len(obs),
            'unique_sample_ids': int(obs['sample_id'].nunique()),
            'counts': counts,
            'control_perturbagens': sorted(obs.loc[controls, 'perturbagen'].dropna().unique().tolist()),
            'noncontrol_perturbagen_labels': int(obs.loc[~controls, 'perturbagen'].nunique()),
            'noncontrol_cell_type_perturbagen_pairs': int(len(obs.loc[~controls, ['cell_type', 'perturbagen']].drop_duplicates())),
            'unique_plate_well_pairs': int(len(obs[['plate', 'well']].drop_duplicates())),
            'unique_plate_well_cell_type_triples': int(len(obs[['plate', 'well', 'cell_type']].drop_duplicates())),
            'explicit_donor_column_present': 'donor' in obs,
            'all_missing_obs_columns': sorted(obs.columns[obs.isna().all()].tolist()),
            'measurement_contract_status': 'pending_upstream_preprocessing_and_biological_unit_audit',
            'caveats': [
                'Pseudobulk-related metadata are present; do not call these rows individual cells or verified raw counts.',
                'Distinct perturbagen labels are not an independently verified count of chemically unique entities.',
                'No donor mapping is inferred from demographics, plate, well or library labels.',
                'Stored train/public_test/private_test labels must be preserved; metadata inspection is not permission to tune on test outcomes.',
                'This processed asset is not asserted to contain the full original experiment.'
            ]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--path', type=Path, default=ROOT / 'data/rna/op3_standardized_processed.h5ad')
    parser.add_argument('--output', type=Path, required=True, help='New JSON file; existing audit is not overwritten')
    args = parser.parse_args()
    result = inventory(args.path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2, allow_nan=False))
