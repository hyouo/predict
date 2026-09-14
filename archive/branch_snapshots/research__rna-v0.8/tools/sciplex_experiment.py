"""External cross-assay-panel test; estimator frozen before external outcomes.

Panel names refer to published assay IDs, not asserted biological replicates.
All three folds must be sealed before evaluation of any panel B target.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from tools.assets import ROOT,verify_file,sha256_file
from tools.op3_inventory import read_column
from src.rna_models import fit_rna_bank,zero_shot
from src.rna_metrics import effect_metrics,row_correlation

SHA='8ad35a20e44bf41cb030a0f063eb6a5727965b8f27d0eda4380bda0e294ccc4c'
CONTEXTS=('CVCL_0004','CVCL_0023','CVCL_0031')
NAMES={'CVCL_0004':'K562','CVCL_0023':'A549','CVCL_0031':'MCF7'}
PROTOCOL_COMMIT='13f3f4a4dded7a023c4fa3c3505f58dff42110e7'
MODEL_SHA='059ca5276f80307c566b5b735f3077d4ff4b3474358ecaa8871f782d3d00f46c'
CODE=['tools/sciplex_experiment.py','src/rna_models.py','src/rna_metrics.py']

class Store:
    def __init__(self,path,target):
        verify_file(path,SHA,160*1024**2)
        self.path=path;self.target=target;self.sources=[c for c in CONTEXTS if c!=target]
        with h5py.File(path) as f:
            self.o=pd.DataFrame({k:read_column(v) for k,v in f['obs'].items()})
            self.v=pd.DataFrame({k:read_column(v) for k,v in f['var'].items()})
            self.ng=int(f['X'].attrs['shape'][1])
        o=self.o
        if not o.sample_id.is_unique or not self.v['_index'].is_unique:raise ValueError('Nonunique IDs')
        o['is_control']=o.is_control.astype(bool)
        o['plate_number']=o.plate.str.replace('plate','',regex=False).astype(int)
        o['panel']=np.where(o.plate_number.between(1,24),'A',np.where(o.plate_number.between(25,48),'B','excluded'))
        o['column']=o.well.str.extract(r'([0-9]+)$').astype(int)
        o['cid']=pd.to_numeric(o.pubchem_cid,errors='coerce').astype('Int64').astype('string')
        o['qualifying']=(~o.is_control)&np.isclose(o.pert_dose_uM,10)&np.isclose(o.pert_time_h,24)&(o.psbulk_cells>=10)&(o.psbulk_counts>0)&o.cid.notna()
        pools=[set(o.loc[o.qualifying&(o.cell_type==c)&o.panel.eq(p),'cid']) for c in CONTEXTS for p in ('A','B')]
        self.cids=sorted(set.intersection(*pools));self.omitted=sorted(set.union(*pools)-set(self.cids))
        ranked=sorted(self.cids,key=lambda c:(hashlib.sha256(('predict-sciplex-v08|'+c).encode()).hexdigest(),c))
        self.anchors=ranked[:16];self.queries=sorted(set(self.cids)-set(self.anchors));self.access=[]
        if len(self.queries)<20:raise ValueError('Insufficient common eligible chemicals')

    def read(self,rows,role,genes=None):
        rows=np.asarray(rows,dtype=int);o=self.o.iloc[rows];ctrl=o.is_control
        valid={
          'source_control':ctrl&o.panel.eq('A')&o.cell_type.isin(self.sources),
          'source_treatment':o.qualifying&o.panel.eq('A')&o.cell_type.isin(self.sources)&o.cid.isin(self.cids),
          'target_baseline':ctrl&o.panel.eq('A')&o.cell_type.eq(self.target),
          'anchor_control':ctrl&o.panel.eq('A')&o.cell_type.eq(self.target),
          'anchor_treatment':o.qualifying&o.panel.eq('A')&o.cell_type.eq(self.target)&o.cid.isin(self.anchors),
          'evaluation_control':ctrl&o.panel.eq('B')&o.cell_type.eq(self.target),
          'evaluation_treatment':o.qualifying&o.panel.eq('B')&o.cell_type.eq(self.target)&o.cid.isin(self.queries)}
        if role not in valid or not valid[role].all():raise PermissionError('Unpermitted sample access: '+role)
        if genes is None:genes=np.arange(self.ng)
        lookup=np.full(self.ng,-1,dtype=int);lookup[genes]=np.arange(len(genes))
        ans=np.zeros((len(rows),len(genes)))
        with h5py.File(self.path) as f:
            x=f['X'];ptr=x['indptr']
            for j,i in enumerate(rows):
                a,b=int(ptr[i]),int(ptr[i+1]);v=x['data'][a:b];ix=x['indices'][a:b]
                if np.any(v<0) or not np.isfinite(v).all() or not np.equal(v,np.floor(v)).all():raise ValueError('Input not raw integer counts')
                total=float(self.o.iloc[i].psbulk_counts)
                if total<=0 or not np.isfinite(total) or float(v.sum())>total*1.00001:raise ValueError('Invalid count denominator')
                pos=lookup[ix];keep=pos>=0
                ans[j,pos[keep]]=np.log2(1+1e6*v[keep]/total)
        self.access.append({'role':role,'sample_ids':o.sample_id.tolist()})
        return ans

    def controls(self,context,panel,role,genes=None):
        o=self.o;mask=o.is_control&o.cell_type.eq(context)&o.panel.eq(panel)&np.isclose(o.pert_time_h,24)
        ids=np.flatnonzero(mask);x=self.read(ids,role,genes);result={}
        for plate in o.iloc[ids].plate.unique():
            positions=np.flatnonzero(o.iloc[ids].plate.eq(plate));positions=positions[np.argsort(o.iloc[ids].iloc[positions]['column'])]
            if len(positions)!=2:raise ValueError('Expected exactly two controls per assay plate')
            result[plate]=x[positions]
        return result,x,ids

    def effects(self,context,panel,cids,control,role,genes,which=None):
        o=self.o;mask=o.qualifying&o.cell_type.eq(context)&o.panel.eq(panel)&o.cid.isin(cids)
        ids=np.flatnonzero(mask);z=self.read(ids,role,genes);groups={c:[] for c in cids}
        for j,i in enumerate(ids):
            row=o.iloc[i];refs=control[row.plate];ref=refs.mean(axis=0) if which is None else refs[which]
            groups[row.cid].append(z[j]-ref)
        if any(not v for v in groups.values()):raise ValueError('A chemical has no allowed observation')
        return np.asarray([np.mean(groups[c],axis=0) for c in cids]),len(ids)

def fit(args):
    if sha256_file(ROOT/'src/rna_models.py')!=MODEL_SHA:raise ValueError('Estimator changed after external protocol registration')
    out=args.output/NAMES[args.target];out.mkdir(parents=True,exist_ok=False)
    s=Store(args.path,args.target);control_all=[]
    for c in s.sources:control_all.append(s.controls(c,'A','source_control')[1])
    score=np.vstack(control_all).mean(axis=0)
    genes=np.sort(np.argsort(-score,kind='stable')[:4000])
    means=[];bases=[];source_rows=0
    for c in s.sources:
        ctrl,z,_=s.controls(c,'A','source_control',genes)
        e,n=s.effects(c,'A',s.cids,ctrl,'source_treatment',genes);means.append(e);bases.append(z.mean(axis=0));source_rows+=n
    source=np.asarray(means);baselines=np.asarray(bases)
    _,z,ctrlids=s.controls(args.target,'A','target_baseline',genes)
    zp,zrec=zero_shot(source,baselines,z.mean(axis=0));zrec['fit_contexts']=s.sources
    # No target treatment was accessed before the zero-shot outputs above.
    ctrl,_,_=s.controls(args.target,'A','anchor_control',genes)
    ac,n=s.effects(args.target,'A',s.anchors,ctrl,'anchor_treatment',genes,0)
    ad,n2=s.effects(args.target,'A',s.anchors,ctrl,'anchor_treatment',genes,1)
    ix=np.array([s.cids.index(c) for c in s.anchors])
    fp,frec=fit_rna_bank(source,ix,ac,ad,block_size=128)
    pred={**{'zero_shot__'+k:v for k,v in zp.items()},**{'few_shot__'+k:v for k,v in fp.items()}}
    # Uncompressed NPZ avoids costly repeated compression during development.
    np.savez(out/'predictions.npz',**pred)
    np.savez_compressed(out/'inputs.npz',source=source,source_baselines=baselines,target_baseline=z.mean(axis=0),anchor_c=ac,anchor_d=ad,anchor_indices=ix)
    (out/'fit.json').write_text(json.dumps({'zero_shot':zrec,'few_shot':frec},indent=2)+'\n')
    if any(a['role'].startswith('evaluation') for a in s.access):raise AssertionError('Premature test access')
    record={'data_sha256':SHA,'protocol_commit':PROTOCOL_COMMIT,'estimator_freeze':'4779f0ed03c44a1e16757945ec2e81567420fb62',
       'code_sha256':{p:sha256_file(ROOT/p) for p in CODE},'prediction_sha256':sha256_file(out/'predictions.npz'),
       'target':args.target,'sources':s.sources,'cids':s.cids,'anchor_cids':s.anchors,'query_cids':s.queries,'omitted_cids':s.omitted,
       'gene_indices':genes.tolist(),'gene_ids':s.v['_index'].iloc[genes].tolist(),'gene_selection':'4000 highest mean source-panel-A CONTROL logCPM; stable tie break',
       'source_treatment_rows':source_rows,'target_anchor_chemicals':16,'target_anchor_treatment_rows':n,'target_control_rows':len(ctrlids),
       'target_low_cell_controls':s.o.iloc[ctrlids].loc[s.o.iloc[ctrlids].psbulk_cells<10,['sample_id','psbulk_cells']].to_dict('records'),
       'predictions_sealed_before_evaluation':True,'access':s.access,'python':sys.version,'finished_utc':datetime.now(timezone.utc).isoformat()}
    (out/'prediction_lock.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ('target','source_treatment_rows','target_anchor_chemicals','target_anchor_treatment_rows','prediction_sha256')}))
    print('ELIGIBLE',len(s.cids),'QUERIES',len(s.queries))

def evaluate(args):
    locks={c:json.loads((args.output/NAMES[c]/'prediction_lock.json').read_text()) for c in CONTEXTS}
    for c,l in locks.items():
        if not l.get('predictions_sealed_before_evaluation'):raise ValueError('Unsealed fold')
        if sha256_file(args.output/NAMES[c]/'predictions.npz')!=l['prediction_sha256']:raise ValueError('Changed prediction file')
        if any(sha256_file(ROOT/p)!=v for p,v in l['code_sha256'].items()):raise ValueError('Changed code')
    out=args.output/'evaluation';out.mkdir(exist_ok=False);summary=[];details=[];records={}
    for c,lock in locks.items():
        s=Store(args.path,c);genes=np.asarray(lock['gene_indices']);ctrl,_,ids=s.controls(c,'B','evaluation_control',genes)
        truth,n=s.effects(c,'B',s.queries,ctrl,'evaluation_treatment',genes)
        pred=np.load(args.output/NAMES[c]/'predictions.npz',allow_pickle=False);ix=np.array([s.cids.index(q) for q in s.queries])
        np.savez_compressed(out/(NAMES[c]+'_truth.npz'),truth=truth,query_indices=ix)
        for method in pred.files:
            p=pred[method][ix];summary.append({'context':NAMES[c],'method':method,**effect_metrics(p,truth)})
            correlations=row_correlation(p,truth)
            for j,cid in enumerate(s.queries):details.append({'context':NAMES[c],'method':method,'cid':cid,'mse':float(np.mean((p[j]-truth[j])**2)),'pearson':float(correlations[j])})
        records[c]={'test_treatment_rows':n,'test_control_rows':len(ids),'query_chemicals':len(s.queries),'access':s.access}
    d=pd.DataFrame(summary);d.to_csv(out/'by_context.csv',index=False);pd.DataFrame(details).to_csv(out/'by_compound.csv',index=False)
    d.groupby('method').mean(numeric_only=True).sort_values('mse').to_csv(out/'macro_summary.csv')
    (out/'evaluation_manifest.json').write_text(json.dumps({'folds':records,'finished_utc':datetime.now(timezone.utc).isoformat()},indent=2)+'\n')
    print(d[['context','method','mse','centered_mse','bias_mse','effect_pearson','centered_retrieval_top1']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['fit','evaluate']);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--target',choices=CONTEXTS,default=CONTEXTS[0]);p.add_argument('--path',type=Path,default=ROOT/'data/rna/srivatsan20_sciplex3_processed.h5ad')
    a=p.parse_args()
    try:fit(a) if a.stage=='fit' else evaluate(a)
    except Exception as e:
        if a.output.exists():
            failure=a.output/f'failure_{a.stage}_{a.target}_{datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")}.json'
            failure.write_text(json.dumps({'error':f'{type(e).__name__}: {e}','stage':a.stage,'target':a.target},indent=2)+'\n')
        raise
