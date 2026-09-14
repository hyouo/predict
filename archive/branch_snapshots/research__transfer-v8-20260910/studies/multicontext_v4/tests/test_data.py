import unittest,tempfile,json
from pathlib import Path
import numpy as np,pandas as pd
from src.data import metadata_plan,Study,sha,SOURCE_SHA
from src.metrics import score_cells,score_pairs
class DataTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);records=[];x=[];r=np.random.default_rng(0)
  for c in range(10):
   for well in range(4):
    records.append(dict(sample=f'{c}_C{well}',original_row=len(records),cell_type=f'C{c}',plate=f'P{c}',well=f'A{well}',is_control=True,perturbagen='ctl'));x.append(r.normal(size=5))
   for p in range(51):
    for k in range(2):
     records.append(dict(sample=f'{c}_D{p}_{k}',original_row=len(records),cell_type=f'C{c}',plate=f'P{c}',well=f'D{p}_{k}',is_control=False,perturbagen=str(1000+p)));x.append(r.normal(size=5))
  self.obs=pd.DataFrame(records).set_index('sample');self.obs.to_csv(self.root/'slice_obs.csv.gz')
  pd.DataFrame({'symbol':[f'g{i}' for i in range(5)]},index=[f'g{i}' for i in range(5)]).to_csv(self.root/'slice_var.csv.gz')
  np.savez_compressed(self.root/'part.npz',X=np.array(x),original_rows=np.arange(len(x)))
  (self.root/'acquisition.json').write_text(json.dumps({'sha256':SOURCE_SHA,'parts':[{'file':'part.npz','sha256':sha(self.root/'part.npz')}]}))
  self.plan=metadata_plan(self.root,self.root/'plan');self.d=Study(self.root,self.plan)
 def tearDown(self):self.temp.cleanup()
 def test_disjoint_contexts(self):
  a,b,c=[set(self.plan[k]) for k in ['train','validation','test']];self.assertFalse(a&b or a&c or b&c)
 def test_test_values_blocked(self):
  rows=self.obs[self.obs.cell_type.isin(self.plan['test'])&~self.obs.is_control].original_row
  with self.assertRaises(PermissionError):self.d.values(rows,'fit')
 def test_test_score_refs_blocked(self):
  c=self.plan['test'][0];v=next(v for k,v in self.plan['refs'].items() if k.startswith(c+'|'))
  with self.assertRaises(PermissionError):self.d.values(v['B'],'reference_B')
  self.assertEqual(self.d.values(v['A'],'baseline_A').shape,(2,5))
 def test_direct_fit_guard(self):
  with self.assertRaises(PermissionError):self.d.effects(self.plan['test'])
 def test_scale_preserved(self):
  y,_=self.d.effects(self.plan['train'][:1]);self.assertEqual(y.shape,(1,51,5));self.assertTrue(np.isfinite(y).all())
 def test_missing_freeze_rejected(self):
  with self.assertRaises(ValueError):self.d.allow_score(self.root/'absent')
 def test_corrupt_part_rejected(self):
  with (self.root/'part.npz').open('ab') as f:f.write(b'x')
  with self.assertRaises(ValueError):Study(self.root,self.plan)
 def test_plan_overwrite_rejected(self):
  with self.assertRaises(FileExistsError):metadata_plan(self.root,self.root/'plan')
 def test_metrics_zero_ties_and_decomposition(self):
  r=np.random.default_rng(2);y=r.normal(size=(2,20,5));p=np.zeros_like(y)
  rows,_=score_cells(p,y,['a','b'],np.arange(20),3,np.repeat(4,20))
  for v in rows:self.assertAlmostEqual(v['retrieval_top1'],.05);self.assertAlmostEqual(v['mse'],v['centered_mse']+v['bias_mse'])
  q=score_pairs(p,y,['a','b'],np.repeat(4,20));self.assertEqual(q[0]['contrast_skill'],0.)
if __name__=='__main__':unittest.main()
