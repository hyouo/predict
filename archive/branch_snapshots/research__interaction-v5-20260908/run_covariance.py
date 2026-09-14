from __future__ import annotations
import json,hashlib,argparse,time
from pathlib import Path
import numpy as np,pandas as pd
from src.data import load_panel
from src.covariance import build_kernels,fit_all
from src.kriging import acquire
ROOT=Path(__file__).resolve().parent

def metrics(y,p):
    e=p-y;ec=e-e.mean(0);yc=y-y.mean(0);pc=p-p.mean(0)
    den=np.linalg.norm(yc,axis=0)*np.linalg.norm(pc,axis=0)
    corr=np.divide(np.sum(yc*pc,0),den,out=np.zeros(y.shape[1]),where=den>1e-12)
    return {'mse':float(np.mean(e*e)),'centered_mse':float(np.mean(ec*ec)),
            'bias_mse':float(np.mean(e.mean(0)**2)), 'pattern_corr':float(np.mean(corr)),
            'sign_accuracy':float(np.mean(np.sign(y)==np.sign(p)))}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seeds',type=int,default=20)
    ap.add_argument('--start-seed',type=int,default=71000);ap.add_argument('--out',default='results/stage1')
    ap.add_argument('--designs',nargs='+',default=['random','cue_aopt','disagreement_aopt','universal_cue','universal_disagreement'])
    ap.add_argument('--budgets',type=int,nargs='+',default=[8,16,24]);args=ap.parse_args()
    out=ROOT/args.out;out.mkdir(parents=True,exist_ok=True)
    p=load_panel(ROOT/'data/hepatocyte_signaling_raw.csv');Y=p.effects;x=p.cues
    (out/'data_audit.json').write_text(json.dumps(p.audit,indent=2))
    rows=[];tuning=[];splits=[];cond=[];preds={};t0=time.monotonic()
    for ci,cname in enumerate(p.contexts):
        src=np.delete(Y,ci,axis=0);k=build_kernels(src,x)
        for seed in range(args.start_seed,args.start_seed+args.seeds):
            rng=np.random.default_rng(np.random.SeedSequence([seed,ci]))
            perm=rng.permutation(72);q=perm[:24];pool=np.sort(perm[24:]);ra=rng.permutation(pool)
            allorders={'random':ra}
            configs={'cue_aopt':(k['cue_original'],'proper'),
                'disagreement_aopt':((k['cue']+k['disagreement'])/2,'proper'),
                'universal_cue':(k['cue_original'],'universal'),
                'universal_disagreement':((k['cue']+k['disagreement'])/2,'universal')}
            for d in args.designs:
                if d!='random':allorders[d]=acquire(*configs[d][:1],pool,q,max(args.budgets),criterion=configs[d][1])[0]
            for b in args.budgets:
                for d in args.designs:
                    a=allorders[d][:b]
                    sid=f'{ci}_{seed}_{b}_{d}'
                    meta={'split_id':sid,'context':cname,'context_index':ci,'seed':seed,'budget':b,'design':d}
                    splits.append({**meta,'anchors':a.tolist(),'queries':q.tolist()})
                    fits=fit_all(src,x,a,Y[ci,a].copy(),k)
                    for name,z in fits.items():
                        yhat=z.prediction[q]
                        rows.append({**meta,'model':name,**metrics(Y[ci,q],yhat),
                                     'prediction_sha256':hashlib.sha256(yhat.astype('<f8').tobytes()).hexdigest()})
                        tuning.append({**meta,'model':name,'parameters':json.dumps(z.parameters)})
                        for j,qi in enumerate(q):
                            cond.append({**meta,'model':name,'query':int(qi),'mse':float(np.mean((Y[ci,qi]-yhat[j])**2))})
                        if seed==args.start_seed:preds[sid+'_'+name]=yhat
        print(cname,len(rows),round(time.monotonic()-t0,1),flush=True)
    pd.DataFrame(rows).to_csv(out/'scores.csv',index=False)
    pd.DataFrame(tuning).to_csv(out/'tuning.csv',index=False)
    pd.DataFrame(cond).to_csv(out/'per_condition.csv',index=False)
    (out/'splits.json').write_text(json.dumps(splits,indent=2))
    np.savez_compressed(out/'sample_predictions.npz',**preds)
    summary=pd.DataFrame(rows).groupby(['budget','design','model'])[['mse','centered_mse','pattern_corr']].mean()
    summary.to_csv(out/'summary.csv')
    (out/'execution.json').write_text(json.dumps({'elapsed_seconds':time.monotonic()-t0,'n_scores':len(rows),'args':vars(args),'completed':True},indent=2))
    print(summary.to_string(),flush=True)
if __name__=='__main__':main()
