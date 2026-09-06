"""Fit, seal and evaluate the registered quantitative OP3 pilot in distinct phases.

Private outcomes are never read. Final-reference wells and public treatment counts
are only opened by `evaluate`, after saved predictions have passed hash verification.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
from src.op3_data import prepare, evaluate_truth, TARGETS, NAMES, SHA256
from src.op3_models import fewshot, zero_shot, LAMBDAS, GATE_LAMBDAS
from src.op3_metrics import reference_design, contrast_metrics

ROOT=Path(__file__).resolve().parents[1]

def dump(path,obj):
    with Path(path).open('x') as f:json.dump(obj,f,indent=2,allow_nan=False)

def hashes():
    paths=['src/op3_data.py','src/op3_models.py','src/op3_metrics.py','tools/op3_pilot.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def fit(output: Path, data: Path, rotation=0, denominator='psbulk_counts', permutation=None,
        frozen_commit='uncommitted-development'):
    if output.exists():raise FileExistsError('Run directories are immutable')
    output.mkdir(parents=True)
    bundle=prepare(data,rotation,denominator)
    if permutation is not None:
        order=np.random.default_rng(permutation).permutation(len(bundle.compounds))
        bundle.source=bundle.source[:,order]
        bundle.source_validation=bundle.source_validation[:,order]
        for missing in bundle.missing_source:
            old=bundle.compounds.index(missing['compound'])
            missing['compound']=bundle.compounds[int(np.flatnonzero(order==old)[0])]
    else:order=np.arange(len(bundle.compounds))
    zs,zparams=zero_shot(bundle)
    predictions={};parameters={'zero_shot':zparams};info={}
    for ct in TARGETS:
        q=bundle.queries[ct];a=bundle.anchors[ct]
        fs,pars=fewshot(bundle.source,a,bundle.y_fit[ct],bundle.y_validation[ct])
        for (cell,method),p in zs.items():
            if cell==ct:predictions[ct+'__'+method]=p[q]
        for method,p in fs.items():predictions[ct+'__'+method]=p[q]
        parameters[ct]=pars
        info[ct]={'name':NAMES[ct],'anchors':[bundle.compounds[i] for i in a],
                  'queries':[bundle.compounds[i] for i in q],
                  'evaluation_rows':bundle.evaluation_ids[ct].tolist(),
                  'anchor_fit_sha256':hashlib.sha256(bundle.y_fit[ct].tobytes()).hexdigest()}
    np.savez_compressed(output/'predictions.npz',**predictions)
    dump(output/'parameters.json',parameters)
    dump(output/'fit_access.json',bundle.reader.access)
    seal={'predictions_complete':True,'private_expression_read':False,
          'public_expression_read':False,'final_controls_read':False,
          'protocol_commit':'f89de0a0bfbf9c6358b0b0ecb99e45eb479124f4',
          'frozen_code_commit':frozen_commit,'code_sha256':hashes(),
          'data_sha256':SHA256,'rotation':rotation,'denominator':denominator,
          'source_permutation_seed':permutation,'source_permutation':order.tolist(),
          'missing_source':bundle.missing_source,'targets':info,
          'arrays_file':'predictions.npz',
          'arrays_sha256':hashlib.sha256((output/'predictions.npz').read_bytes()).hexdigest(),
          'python':platform.python_version(),'numpy':np.__version__,
          'models':sorted(predictions)}
    dump(output/'seal.json',seal)
    print(json.dumps({'phase':'fit','output':str(output),'models':len(predictions),
                      'arrays_sha256':seal['arrays_sha256'],'test_outcomes_read':False}),flush=True)
    return bundle


def metrics(y,p):
    err=p-y;bias=err.mean(0)
    norms=np.linalg.norm(y,axis=1)*np.linalg.norm(p,axis=1)
    cosine=np.divide(np.sum(y*p,axis=1),norms,out=np.zeros(len(y)),where=norms>0)
    return {'mse':float(np.mean(err**2)),
            'centered_mse':float(np.mean((err-bias)**2)),
            'bias_mse':float(np.mean(bias**2)),
            'row_rmse':float(np.sqrt(np.mean(err**2,axis=1)).mean()),
            'mae':float(np.mean(abs(err))),'cosine':float(cosine.mean())}


def score(output: Path,data: Path):
    seal=json.loads((output/'seal.json').read_text())
    if seal['code_sha256']!=hashes():raise RuntimeError('Code changed after prediction')
    if (output/'summary.csv').exists():raise FileExistsError('Results already exist')
    bundle=prepare(data,seal['rotation'],seal['denominator'])
    truth,reps=evaluate_truth(bundle,output/'seal.json')
    saved=np.load(output/'predictions.npz',allow_pickle=False)
    rows=[];per=[]
    for key in saved.files:
        ct,method=key.split('__');y=truth[ct];p=saved[key]
        # Retrieval removes only TRAIN-estimated target background, not public mean.
        bg=bundle.y_fit[ct].mean(0)
        yp=y-bg;pp=p-bg
        yn=np.maximum(np.linalg.norm(yp,axis=1),1e-12)
        pn=np.maximum(np.linalg.norm(pp,axis=1),1e-12)
        sim=(pp/pn[:,None])@(yp/yn[:,None]).T
        rank=1+np.sum(sim>np.diag(sim)[:,None]+1e-12,axis=1)
        tied=np.sum(abs(sim-np.diag(sim)[:,None])<=1e-12,axis=1)
        rank=rank+(tied-1)/2
        row={'context':NAMES[ct],'method':method,'rotation':seal['rotation'],
             'denominator':seal['denominator'],'n_compounds':len(y),
             'n_genes':y.shape[1],**metrics(y,p),
             **contrast_metrics(y,p,reference_design(bundle.reader.obs,bundle.evaluation_ids[ct],seal['targets'][ct]['queries'])),
             'mean_reciprocal_rank':float(np.mean(1/rank)),
             'retrieval_top1':float(np.mean((rank==1).astype(float)))}
        rows.append(row)
        for j,drug in enumerate(seal['targets'][ct]['queries']):
            per.append({'context':NAMES[ct],'compound':drug,'method':method,
                        'mse':float(np.mean((p[j]-y[j])**2)),
                        'retrieval_rank':float(rank[j])})
    pd.DataFrame(rows).to_csv(output/'summary.csv',index=False)
    pd.DataFrame(per).to_csv(output/'per_compound.csv',index=False)
    np.savez_compressed(output/'evaluation_truth.npz',**truth)
    dump(output/'evaluation_access.json',bundle.reader.access)
    print(pd.DataFrame(rows).groupby('method')[['mse','centered_mse','bias_mse','cosine']].mean().sort_values('mse').to_string(),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['fit','evaluate'])
    p.add_argument('--data',type=Path,default=ROOT/'data/rna/op3_standardized_processed.h5ad')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rotation',type=int,default=0,choices=[0,1,2])
    p.add_argument('--denominator',choices=['psbulk_counts','row_sum'],default='psbulk_counts')
    p.add_argument('--permutation',type=int)
    p.add_argument('--frozen-commit',default='uncommitted-development')
    a=p.parse_args()
    if a.phase=='fit':fit(a.output,a.data,a.rotation,a.denominator,a.permutation,a.frozen_commit)
    else:score(a.output,a.data)

if __name__=='__main__':main()
