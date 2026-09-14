"""Protocol 006 source-only outputs on existing, explicitly exploratory folds."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from tools.assets import ROOT,sha256_file
from src.spectral_consensus import spectral_consensus
from src.rna_metrics import effect_metrics,row_correlation


def run(a):
    out=a.output;out.mkdir(parents=True,exist_ok=False);sources={};labels={};input_hashes={}
    if a.dataset=='op3':
        lock=json.loads((a.frozen/'prediction_lock.json').read_text())
        inp=a.frozen/'source_inputs.npz';input_hashes[str(inp)]=sha256_file(inp)
        source=np.load(inp,allow_pickle=False)['means']
        for name in ('B','Myeloid'):sources[name]=source;labels[name]=lock['drugs']
        truth_dir=a.frozen/a.evaluation
    else:
        for name in ('K562','A549','MCF7'):
            d=a.frozen/name;lock=json.loads((d/'prediction_lock.json').read_text())
            inp=d/'inputs.npz';input_hashes[str(inp)]=sha256_file(inp)
            sources[name]=np.load(inp,allow_pickle=False)['source'];labels[name]=lock['cids']
        truth_dir=a.frozen/'evaluation'
    predictions={};fits={};hashes={}
    for name,source in sources.items():
        pred={};record={}
        for ridge in (.1,1.,10.):
            for weighted in (False,True):
                method=f"rho{ridge:g}_"+('gene' if weighted else 'plain')
                pred[method],record[method]=spectral_consensus(source,ridge=ridge,gene_weighted=weighted)
        predictions[name]=pred;fits[name]=record
        file=out/f'{name}_predictions.npz';np.savez_compressed(file,**pred);hashes[file.name]=sha256_file(file)
    seal={'scientific_status':'exploratory after previous test analyses','nominated_default':'rho1_gene',
          'zero_target_treatment_inputs':True,'zero_target_baseline_inputs':True,
          'input_hashes':input_hashes,'prediction_hashes':hashes,
          'code_sha256':{p:sha256_file(ROOT/p) for p in ('src/spectral_consensus.py','tools/spectral_experiment.py')},
          'sealed_utc':datetime.now(timezone.utc).isoformat()}
    (out/'prediction_lock.json').write_text(json.dumps(seal,indent=2)+'\n')
    (out/'fit.json').write_text(json.dumps(fits,indent=2)+'\n')
    summary=[];rows=[]
    for name,pred in predictions.items():
        t=np.load(truth_dir/f'{name}_truth.npz',allow_pickle=False)
        truth=t['mean'] if a.dataset=='op3' else t['truth']
        ix=t['source_indices'] if a.dataset=='op3' else t['query_indices']
        for method,p in pred.items():
            p=p[ix];summary.append({'context':name,'method':method,**effect_metrics(p,truth)})
            cor=row_correlation(p,truth)
            for j,q in enumerate(ix):rows.append({'context':name,'method':method,'compound':labels[name][q],
                                                  'mse':float(np.mean((p[j]-truth[j])**2)),'pearson':float(cor[j])})
    df=pd.DataFrame(summary);df.to_csv(out/'by_context.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'by_compound.csv',index=False)
    df.groupby('method').mean(numeric_only=True).to_csv(out/'macro_summary.csv')
    print(df[['context','method','mse','centered_mse','bias_mse']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',choices=['op3','sciplex'],required=True)
    p.add_argument('--frozen',type=Path,required=True);p.add_argument('--evaluation',default='evaluation_private_test')
    p.add_argument('--output',type=Path,required=True);run(p.parse_args())
