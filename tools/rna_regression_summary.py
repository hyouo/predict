"""Compact numerical regression record; not an additional scientific experiment."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from tools.assets import ROOT,sha256_file

METHODS=['zero_shot__zero','zero_shot__source_consensus_shrink','few_shot__mixed_stack_cv','few_shot__mixed_stack_adjusted','few_shot__simple_stack_adjusted']
EXPECTED={
 'op3':dict(zip(METHODS,[0.1951018811187269, 0.1720730599495276, 0.1734205461441141, 0.1661479929976906, 0.1662719517703888])),
 'sciplex':dict(zip(METHODS,[2.4217097849183173,2.4114579153889215,2.7435998218058906,2.4171921725618115,2.4171749077305047]))}


def main():
    out=ROOT/'runs/ci-rna-summary';out.mkdir(exist_ok=False)
    matrices={};record={}
    for dataset,d in [('op3',ROOT/'runs/ci-op3/evaluation_private_test'),('sciplex',ROOT/'runs/ci-sciplex/evaluation')]:
        frame=pd.read_csv(d/'macro_summary.csv').set_index('method');frame.to_csv(out/(dataset+'_macro.csv'))
        pd.read_csv(d/'by_context.csv').to_csv(out/(dataset+'_contexts.csv'),index=False)
        rec={}
        for method in METHODS:
            actual=float(frame.loc[method,'mse']);expected=EXPECTED[dataset][method]
            if not np.isclose(actual,expected,rtol=2e-7,atol=2e-8):raise AssertionError(f'{dataset}/{method} unexpected regression: {actual} != {expected}')
            rec[method]={'mse':actual,'expected':expected,'absolute_difference':abs(actual-expected)}
        record[dataset]=rec
        locations=[(c,ROOT/f'runs/ci-op3/{c}_predictions.npz') for c in ('B','Myeloid')] if dataset=='op3' else [(c,ROOT/f'runs/ci-sciplex/{c}/predictions.npz') for c in ('K562','A549','MCF7')]
        for context,path in locations:
            with np.load(path,allow_pickle=False) as p:
                for method in METHODS:
                    a=p[method];rows=np.arange(0,len(a),11);cols=np.unique(np.r_[0,127,128,255,256,511,1023,a.shape[1]-1])
                    matrices[f'{dataset}__{context}__{method}']=a[np.ix_(rows,cols)]
    np.savez_compressed(out/'sampled_predictions.npz',**matrices)
    record['sample_array_sha256']=sha256_file(out/'sampled_predictions.npz')
    record['purpose']='cross-machine regression only; source runs and outcomes are already consumed'
    (out/'comparison.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))

if __name__=='__main__':main()
