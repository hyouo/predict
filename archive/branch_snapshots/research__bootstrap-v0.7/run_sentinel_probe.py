"""A targeted missing-factor observation vs an equal-cost covered observation.

All candidate probe rows are excluded from BOTH arms' query sets before outcomes
are examined. The comparison changes information, not only the prediction model.
"""
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np,pandas as pd
from src.data import load_panel
from src.covariance import build_kernels,fit_all
from src.aggregation import fit_stage2
from src.kriging import acquire
from run_covariance import metrics
ROOT=Path(__file__).resolve().parent
FACTORS=['IL1a','TGFa','IL6','INS','TNFa','PI3Ki','MEK12i','IKKi']

def probe_sets(x,j):
    held=x[:,j]==1
    if j<5:
        probe=np.flatnonzero(held & (x[:,:5].sum(1)==1) & (x[:,5:].sum(1)==0))
    else:
        probe=np.flatnonzero(held & (x[:,:5].sum(1)==1) & (x[:,5:].sum(1)==1))
    q=np.setdiff1d(np.flatnonzero(held),probe)
    pool=np.flatnonzero(~held)
    assert len(probe)==(1 if j<5 else 5) and len(q)>0
    return pool,q,probe

def choose_probe(K,a,q,pool,noise=.1):
    c=K+1.
    for i in a:
        v=c[:,i].copy();c=c-np.outer(v,v)/(c[i,i]+noise)
    gains=np.mean(c[np.ix_(q,pool)]**2,axis=0)/(np.diag(c)[pool]+noise)
    return int(pool[np.argmax(gains)]),gains

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--context',type=int,required=True);ci=ap.parse_args().context
    p=load_panel(ROOT/'data/hepatocyte_signaling_raw.csv');Y=p.effects;x=p.cues
    src=np.delete(Y,ci,axis=0);k=build_kernels(src,x);rows=[];splits=[];preds={};t0=time.monotonic()
    names=['source_offset','cue','output_disagreement_select','output_marginal_select','expanded_bootstrap','expanded_stacking']
    for j,factor in enumerate(FACTORS):
        pool,q,probes=probe_sets(x,j)
        order=acquire(k['cue_original'],pool,q,25,criterion='proper')[0]
        for b in [8,16,24]:
            a=order[:b];chosen,gains=choose_probe(k['cue_original'],a,q,probes)
            arms={'base':a,'covered_extra':order[:b+1],'missing_factor_probe':np.r_[a,chosen]}
            for arm,aa in arms.items():
                sid=f'{ci}_{factor}_{b}_{arm}'
                assert set(aa).isdisjoint(q)
                meta={'split_id':sid,'context':p.contexts[ci],'held_factor':factor,'base_budget':b,'total_budget':len(aa),'arm':arm,'n_queries':len(q)}
                splits.append({**meta,'anchors':aa.tolist(),'queries':q.tolist(),'reserved_probe_candidates':probes.tolist(),
                    'chosen_probe':chosen,'probe_gains_working_model':gains.tolist()})
                fits=fit_all(src,x,aa,Y[ci,aa].copy(),k);fits.update(fit_stage2(src,x,aa,Y[ci,aa].copy(),k))
                for name in names:
                    pp=fits[name].prediction[q]
                    rows.append({**meta,'model':name,**metrics(Y[ci,q],pp),'prediction_sha256':hashlib.sha256(pp.astype('<f8').tobytes()).hexdigest()})
                    if b==16:preds[sid+'_'+name]=pp
        print(ci,factor,len(rows),round(time.monotonic()-t0,2),flush=True)
    out=ROOT/f'results/sentinel_probe/context{ci}';out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'scores.csv',index=False);(out/'splits.json').write_text(json.dumps(splits,indent=2));np.savez_compressed(out/'predictions.npz',**preds)
    (out/'execution.json').write_text(json.dumps({'completed':True,'n_scores':len(rows),'elapsed_seconds':time.monotonic()-t0},indent=2))
    print(pd.DataFrame(rows).query('base_budget==16').groupby(['model','arm']).mse.mean().unstack().to_string())
if __name__=='__main__':main()
