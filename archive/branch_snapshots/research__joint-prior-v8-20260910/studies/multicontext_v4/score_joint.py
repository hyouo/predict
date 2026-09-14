"""Single post-freeze joint scoring of all original and full-feature comparators."""
from pathlib import Path
import argparse,json,time,traceback
import numpy as np,pandas as pd
from src.data import Study,sha
from src.metrics import score_cells,score_pairs

def main():
 p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--plan',required=True);p.add_argument('--primary',required=True);p.add_argument('--comparators',required=True);p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True)
 try:
  roots=[Path(a.primary),Path(a.comparators)];freezes=[]
  for root in roots:
   f=json.loads((root/'FREEZE.json').read_text());assert f['complete'] and not f['test_scored']
   if sha(root/'predictions.npz')!=f['prediction_sha256'] or sha(a.plan)!=f['plan_sha256']:raise ValueError('Frozen predictions or plan modified')
   freezes.append({'file':str(root/'FREEZE.json'),'sha256':sha(root/'FREEZE.json')})
  (out/'JOINT_FREEZE.json').write_text(json.dumps({'time_unix':time.time(),'frozen_inputs':freezes,'code_sha256':sha(__file__)},indent=2))
  d=Study(a.data,a.plan);d.allow_score(out/'JOINT_FREEZE.json');truth,records=d.effects(d.plan['test'],role='score')
  count=np.load(roots[0]/'source_inputs.npz',allow_pickle=False)['source_counts'];rows=[];drugrows=[];pairs=[]
  for root in roots:
   with np.load(root/'predictions.npz',allow_pickle=False) as arrays:
    for name in arrays.files:
     pred=arrays[name];r,b=score_cells(pred,truth,d.plan['test'],d.drugs,3,count);ps=score_pairs(pred,truth,d.plan['test'],count)
     rows.extend(dict(model=name,**v) for v in r);drugrows.extend(dict(model=name,**v) for v in b);pairs.extend(dict(model=name,**v) for v in ps)
  frame=pd.DataFrame(rows);frame.to_csv(out/'by_context.csv',index=False);pd.DataFrame(drugrows).to_csv(out/'by_drug.csv',index=False);pd.DataFrame(pairs).to_csv(out/'paired_contrasts.csv',index=False)
  frame.groupby('model').mean(numeric_only=True).reset_index().to_csv(out/'macro.csv',index=False);pd.DataFrame(pairs).groupby('model').mean(numeric_only=True).reset_index().to_csv(out/'contrast_macro.csv',index=False)
  np.savez_compressed(out/'truth.npz',truth=truth,contexts=np.asarray(d.plan['test']),drugs=np.asarray(d.drugs),genes=np.asarray(d.genes))
  (out/'score_access.json').write_text(json.dumps(d.log));(out/'target_records.json').write_text(json.dumps(records));(out/'COMPLETE.json').write_text(json.dumps({'time_unix':time.time(),'complete':True,'truth_sha256':sha(out/'truth.npz'),'sources':freezes},indent=2))
  print(frame.groupby('model').mse.mean().sort_values().to_string())
 except Exception as e:(out/'FAILED.json').write_text(json.dumps(dict(error=str(e),traceback=traceback.format_exc())));raise
if __name__=='__main__':main()
