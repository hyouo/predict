"""Verify frozen model bytes, numerical reference scores, and profile access logs."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from src.op3_data import Reader
from tools.op3_pilot import ROOT, hashes


def verify(run: Path | None = None):
    frozen=json.loads((ROOT/'docs/experiments/002-op3-frozen-manifest.json').read_text())
    if hashes()!=frozen['code_sha256']:raise ValueError('Experiment-002 core differs from preregistered hashes')
    result={'core_hashes_match':True}
    if run is not None:
        expected=pd.read_csv(ROOT/'reports/op3/experiment-002/primary_macro.csv').set_index('method').sort_index()
        cols=['mse','centered_mse','bias_mse','row_rmse','mae','cosine','plate_orthogonal_mse','plate_orthogonal_cosine']
        current=pd.read_csv(run/'summary.csv').groupby('method')[cols].mean().sort_index()
        if list(expected.index)!=list(current.index):raise ValueError('Model identities differ')
        np.testing.assert_allclose(current[cols],expected[cols],rtol=1e-6,atol=1e-8)
        r=Reader(ROOT/'data/rna/op3_standardized_processed.h5ad');part=r.partition()
        fit=json.loads((run/'fit_access.json').read_text());ev=json.loads((run/'evaluation_access.json').read_text())
        fitids={i for x in fit for i in x['row_ids']};allids=fitids|{i for x in ev for i in x['row_ids']}
        if any(r.obs.iloc[i].split=='private_test' for i in allids):raise ValueError('Private data access')
        if any(part[i]==2 or r.obs.iloc[i].split=='public_test' for i in fitids):raise ValueError('Test data in fitting')
        result.update({'reference_scores_match':True,'max_score_abs_difference':float(np.max(abs(current[cols].to_numpy()-expected[cols].to_numpy()))),
                       'private_expression_used':False,'fit_rows':len(fitids),'evaluation_plus_fit_rows':len(allids)})
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path);a=p.parse_args();verify(a.run)
