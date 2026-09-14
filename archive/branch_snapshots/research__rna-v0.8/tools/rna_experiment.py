"""Fit first, seal numerical predictions, then separately evaluate OP3 holdouts."""
from __future__ import annotations
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import platform
import sys
import numpy as np
import pandas as pd
from tools.assets import ROOT, sha256_file
from src.rna_data import OP3Store, TARGET_CONTEXTS, NAMES
from src.rna_models import fit_rna_bank, zero_shot
from src.rna_metrics import effect_metrics, row_correlation

CODE_FILES=['src/rna_data.py','src/rna_models.py','src/rna_metrics.py','tools/rna_experiment.py']

def store(args):
    mapping=json.loads((ROOT/'reports/data/op3_donors.json').read_text())
    return OP3Store(args.path,mapping,denominator=args.denominator,pseudocount=args.pseudocount)

def fit(args):
    args.output.mkdir(parents=True,exist_ok=False)
    s=store(args);donors=tuple(args.donors.split(','));source=s.source(donors)
    record={'started_utc':datetime.now(timezone.utc).isoformat(),'purpose':'retrospective_quantitative_RNA',
            'data_sha256':sha256_file(args.path),'code_sha256':{p:sha256_file(ROOT/p) for p in CODE_FILES},
            'protocol_sha256':sha256_file(ROOT/'docs/experiments/002-rna-protocol.md'),
            'donors_for_fitting':list(donors),'source_contexts':['T','NK'],'target_contexts':['B','Myeloid'],
            'target_query_outcomes_received_by_fit':False,'evaluation_controls_used_by_fit':False,
            'denominator':args.denominator,'pseudocount_CPM':args.pseudocount,
            'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'files':{},'contexts':{}}
    np.savez_compressed(args.output/'source_inputs.npz',means=source.means,baselines=source.baselines,donor_effects=source.donor_effects)
    for context in TARGET_CONTEXTS:
        baseline=s.baseline(context,donors)
        zp,zrec=zero_shot(source.means,source.baselines,baseline)
        # This access occurs only AFTER zero-shot predictions have been computed.
        anchor=s.anchors(context,source)
        fp,frec=fit_rna_bank(source.means,anchor.indices,anchor.reference_c,anchor.reference_d,block_size=args.block_size)
        allp={**{'zero_shot__'+k:v for k,v in zp.items()},**{'few_shot__'+k:v for k,v in fp.items()}}
        name=NAMES[context];filename=f'{name}_predictions.npz'
        np.savez_compressed(args.output/filename,**allp)
        (args.output/f'{name}_fit.json').write_text(json.dumps({'zero_shot':zrec,'few_shot':frec},indent=2)+'\n')
        record['files'][filename]=sha256_file(args.output/filename)
        record['contexts'][context]={'anchor_compounds':[source.drugs[i] for i in anchor.indices],
                                     'anchor_count':len(anchor.indices),'target_treatment_observations':anchor.observation_count}
    forbidden=[a for a in s.access if a['role'].startswith('evaluation')]
    if forbidden:raise RuntimeError('Fit accessed evaluation data')
    (args.output/'fit_access.json').write_text(json.dumps(s.access,indent=2)+'\n')
    record['omitted_source_compounds']=source.omitted_source_compounds
    record['drugs']=source.drugs;record['genes']=s.var.ensembl_id.tolist()
    record['predictions_sealed_before_evaluation']=True
    record['finished_utc']=datetime.now(timezone.utc).isoformat()
    (args.output/'prediction_lock.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'stage':'fit','directory':str(args.output),'files':record['files'],'context_budgets':record['contexts']},indent=2))

def evaluate(args):
    root=args.output;lockpath=root/'prediction_lock.json';lock=json.loads(lockpath.read_text())
    for name,digest in lock['files'].items():
        if sha256_file(root/name)!=digest:raise ValueError('A sealed prediction file changed')
    for name,digest in lock['code_sha256'].items():
        if sha256_file(ROOT/name)!=digest:raise ValueError('Code changed since prediction sealing; use a new fit directory')
    if args.denominator!=lock['denominator'] or args.pseudocount!=lock['pseudocount_CPM']:raise ValueError('Endpoint mismatch')
    out=root/('evaluation_'+args.split+('_'+args.eval_donors.replace(',','_').replace(' ','') if args.eval_donors else ''))
    out.mkdir(exist_ok=False)
    s=store(args);source=s.source(tuple(lock['donors_for_fitting']))
    eval_donors=tuple(args.eval_donors.split(',')) if args.eval_donors else tuple(lock['donors_for_fitting'])
    summary=[];rows=[];truthfiles={};repstats=[]
    for context in TARGET_CONTEXTS:
        ix,t,reps,omitted=s.truth(context,args.split,source,prediction_lock=lockpath,evaluation_donors=eval_donors)
        p=np.load(root/f'{NAMES[context]}_predictions.npz',allow_pickle=False)
        np.savez_compressed(out/f'{NAMES[context]}_truth.npz',mean=t,donor_effects=reps,source_indices=ix)
        truthfiles[context]={'queries':[source.drugs[i] for i in ix],'omitted_no_source_coverage':omitted,
                             'independent_biological_donors':len(eval_donors)}
        for method in p.files:
            pred=p[method][ix];metrics=effect_metrics(pred,t)
            summary.append({'context':NAMES[context],'method':method,**metrics})
            pear=row_correlation(pred,t);cos=row_correlation(pred,t,False)
            for j,i in enumerate(ix):rows.append({'context':NAMES[context],'method':method,'compound':source.drugs[i],
                                                'mse':float(np.mean((pred[j]-t[j])**2)),'pearson':float(pear[j]),'cosine':float(cos[j])})
        for a in range(len(reps)):
            for b in range(a):
                ok=np.isfinite(reps[a]).all(axis=1)&np.isfinite(reps[b]).all(axis=1)
                repstats.append({'context':NAMES[context],'donor_a':eval_donors[a],'donor_b':eval_donors[b],
                                 'compound_count':int(ok.sum()),'mean_cross_donor_pearson':float(row_correlation(reps[a,ok],reps[b,ok]).mean()),
                                 'mean_cross_product':float(np.mean(reps[a,ok]*reps[b,ok])),
                                 'note':'Independent donors, not technical replicates; heterogeneous effects and finite controls remain.'})
    pd.DataFrame(summary).to_csv(out/'by_context.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'by_compound.csv',index=False)
    d=pd.DataFrame(summary);numeric=[c for c in d.columns if c not in ('context','method')]
    d.groupby('method')[numeric].mean().sort_values('mse').to_csv(out/'macro_summary.csv')
    (out/'evaluation_manifest.json').write_text(json.dumps({'prediction_lock_sha256':sha256_file(lockpath),'split':args.split,
                              'evaluation_donors':eval_donors,'contexts':truthfiles,'cross_donor_diagnostics':repstats,
                              'sample_access':s.access,'finished_utc':datetime.now(timezone.utc).isoformat()},indent=2)+'\n')
    print(d[['context','method','mse','centered_mse','bias_mse','effect_pearson','centered_retrieval_top1']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['fit','evaluate'])
    p.add_argument('--path',type=Path,default=ROOT/'data/rna/op3_standardized_processed.h5ad')
    p.add_argument('--output',type=Path,required=True);p.add_argument('--donors',default='Donor 1,Donor 2,Donor 3')
    p.add_argument('--eval-donors');p.add_argument('--split',choices=['public_test','private_test'],default='public_test')
    p.add_argument('--denominator',choices=['psbulk_counts','retained_sum'],default='psbulk_counts')
    p.add_argument('--pseudocount',type=float,default=1.0);p.add_argument('--block-size',type=int,default=256)
    a=p.parse_args();fit(a) if a.stage=='fit' else evaluate(a)
