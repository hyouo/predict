"""Create the narrowly scoped OP3 measurement contract; never read test expression.

Upstream semantics are documented using fixed source-code references. Numeric checks
are limited to train and allowed A/B control rows. A readable file is not by itself
approval for causal, donor-independent or full-transcriptome claims.
"""
from __future__ import annotations
import argparse, json, hashlib
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from src.op3_data import Reader, SHA256, SOURCES, TARGETS, NAMES
from tools.op3_pilot import ROOT, dump


def contract(path: Path) -> dict:
    reader=Reader(path);o=reader.obs;p=reader.partition(0)
    ids=np.flatnonzero((o.split=='train')|(p==0)|(p==1));ratios=[]
    nonnegative=True;integral=True
    with h5py.File(path,'r') as f:
        ptr=f['X/indptr'][:]
        dtype=str(f['X/data'].dtype)
        for i in ids:
            v=f['X/data'][int(ptr[i]):int(ptr[i+1])]
            nonnegative &= bool(np.all(v>=0) and np.isfinite(v).all())
            integral &= bool(np.all(v==np.floor(v)))
            ratios.append(float(v.sum()/o.iloc[i].psbulk_counts))
    rows=[]
    for ct in SOURCES+TARGETS:
        for split,part in o[o.cell_type==ct].groupby('split',observed=True):
            rows.append({'context':NAMES[ct],'cell_type':ct,'split':split,
                         'profiles':len(part),'compounds':int(part.perturbagen.nunique())})
    treat=o[~o.is_control.astype(bool)]
    profiles_per_pair=treat.groupby(['cell_type','perturbagen'],observed=True).size()
    result={'version':'op3-exploratory-002','data_sha256':SHA256,'bytes':path.stat().st_size,
      'shape':[int(n) for n in reader.shape],'units':'nonnegative integer raw-count-derived pseudobulk on an upstream filtered gene panel',
      'aggregation_unit':'plate x well x reannotated cell type; NOT individual cells',
      'training_approval':{'approved_for_exploratory_pilot':True,'approved_for_pristine_inductive_benchmark':False,
          'approved_for_donor_independent_validation':False,'approved_for_causal_mechanism_claims':False,
          'permitted_splits':['train','control_A','control_B'], 'evaluation_only':['public_test','control_C'],
          'forbidden_expression_split':'private_test'},
      'audit_expression_scope':'train plus A/B controls only; test metadata are known',
      'expression_row_ids_checked':ids.tolist(),'expression_dtype':dtype,'nonnegative_finite':nonnegative,
      'all_checked_entries_integral':integral,
      'matrix_row_sum_over_psbulk_counts_range':[min(ratios),max(ratios)],
      'normalization':'log2(1 + 1e6 * counts / psbulk_counts); denominator was computed after initial upstream gene filtering',
      'normalization_not':'not absolute RNA/cell, not verified full-transcriptome CPM, not official signed-p-value outcome',
      'source_contexts':{c:NAMES[c] for c in SOURCES},'target_contexts':{c:NAMES[c] for c in TARGETS},
      'profile_coverage':rows,'profiles_per_cell_compound':{str(k):int(v) for k,v in profiles_per_pair.value_counts().items()},
      'cell_counts_are_not_replicate_n':True,'explicit_donor_id_present':'donor' in o or 'donor_id' in o,
      'donor_reconstruction_performed':False,'plate_count':int(o.plate.nunique()),
      'unique_plate_wells':int(o[['plate','well']].drop_duplicates().shape[0]),
      'paired_cell_types_share_treatment_wells':True,'nominal_dose_uM':sorted(treat.pert_dose_uM.unique().tolist()),
      'nominal_time_hours':sorted(treat.pert_time_h.unique().tolist()),
      'control_design':{'physical_wells_per_plate':8,'same_assignment_across_cell_types':True,
        'primary_A_indices':[0,3,6],'primary_B_indices':[1,4,7],'primary_C_indices':[2,5],
        'control_profiles_per_target_A_B_C':[18,18,12],
        'physical_disjointness_not_donor_independence':True},
      'global_gene_selection_uses_all_conditions':True,
      'upstream_provenance':[
        {'repository':'openproblems-bio/task_perturbation_prediction','commit':'626f94521fd236e6b79cc60696159172dd07d2bd',
         'path':'src/process_dataset/compute_pseudobulk/script.py','blob':'71c8ddb6a9dccc8abc0fc021d3d7c53adadf3884',
         'evidence':'X assigned from raw.X and summed by plate_well_celltype_reannotated'},
        {'repository':'openproblems-bio/task_perturbation_prediction','commit':'626f94521fd236e6b79cc60696159172dd07d2bd',
         'path':'src/process_dataset/filter_vars/script.R','blob':'7566f1cfac94034991d5c7ed23aab7b5d549f364',
         'evidence':'filterByExpr per cell type, intersection, before split-specific analysis'},
        {'repository':'theislab/Chem-PerturBridge','commit':'994399c79048c4bb92118363839ec6a2d36afb8f',
         'path':'src/pseudobulking/datasets/op3/standardization.py','blob':'ff55571ba01c8437a42325880c0e86df7a799e30',
         'evidence':'DMSO control filter; psbulk_counts calculation; strict metadata schema; X int64 copy'}],
      'limitations':['Gene universe was globally filtered before our split and cannot support target-only absent-gene discovery.',
                     'No donor IDs: plate/well replicates are not asserted independent donors.',
                     'Source/target profiles can share the same treated well and assay/culture variation.',
                     'Only two source and two evaluated target contexts; no external study.',
                     'No ATAC or target mechanism measurements in this asset.']}
    if not (nonnegative and integral) or reader.shape!=(1813,5288):raise ValueError('Unexpected measured asset')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--path',type=Path,default=ROOT/'data/rna/op3_standardized_processed.h5ad')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.parent.mkdir(parents=True,exist_ok=True);dump(a.output,contract(a.path))
