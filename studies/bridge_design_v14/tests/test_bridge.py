import unittest,sys,itertools
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bridge import enumerate_worlds,conditional_interval,one_link_design,greedy_one_per_condition

class BridgeTests(unittest.TestCase):
 def test_equal_lengths_required(self):
  with self.assertRaises(ValueError):enumerate_worlds([1,2],[1,2,3])
 def test_nonfinite_rejected(self):
  with self.assertRaises(ValueError):enumerate_worlds([1,np.nan],[1,2])
 def test_dimension_rejected(self):
  with self.assertRaises(ValueError):enumerate_worlds([[1,2]],[[1,2]])
 def test_exhaustive_world_count(self):
  w,v=enumerate_worlds([0,1,3],[2,3,9]);self.assertEqual(w.shape,(6,3));self.assertEqual(len(v),6)
 def test_centered_covariance(self):
  x=np.array([0.,1.,3.]);y=np.array([2.,3.,9.]);w,v=enumerate_worlds(x,y)
  for p,c in zip(w,v):self.assertAlmostEqual(c,np.mean((x-x.mean())*(y[p]-y.mean())))
 def test_additive_shift_invariance(self):
  a=one_link_design([0,1,3],[2,3,9]);b=one_link_design([70,71,73],[-8,-7,-1]);self.assertAlmostEqual(a['worst_width'],b['worst_width'])
 def test_constant_margin(self):
  a=one_link_design([1,1,1],[0,3,4]);self.assertEqual(a['initial_width'],0);self.assertEqual(a['guaranteed_reduction'],0)
 def test_all_links_reduce_range(self):
  w,v=enumerate_worlds([0,1,3],[2,3,9]);base=np.ptp(v)
  for i in range(3):
   for j in range(3):
    l,h,n=conditional_interval(w,v,[(i,j)]);self.assertLessEqual(h-l,base+1e-12);self.assertEqual(n,2)
 def test_two_links_fix_third(self):
  w,v=enumerate_worlds([0,1,3],[2,3,9])
  for p in w:
   l,h,n=conditional_interval(w,v,[(0,p[0]),(2,p[2])]);self.assertEqual(l,h);self.assertEqual(n,1)
 def test_conflicting_links_rejected(self):
  w,v=enumerate_worlds([0,1,3],[2,3,9])
  with self.assertRaises(ValueError):conditional_interval(w,v,[(0,1),(2,1)])
 def test_duplicate_record_ids_rejected(self):
  with self.assertRaises(ValueError):one_link_design([0,1,3],[2,3,9],['a','a','c'])
 def test_triplicate_half_width_guarantee(self):
  rng=np.random.default_rng(33)
  for _ in range(300):
   x=rng.normal(size=3);y=rng.normal(size=3);d=one_link_design(x,y)
   self.assertLessEqual(d['worst_width'],.5*d['initial_width']+1e-12)
 def test_closed_form_triplet(self):
  rng=np.random.default_rng(28)
  for _ in range(100):
   x=np.sort(rng.normal(size=3));y=np.sort(rng.normal(size=3));d=one_link_design(x,y)
   exact=min(x[1]-x[0],x[2]-x[1])*(y[2]-y[0])/3
   self.assertAlmostEqual(d['worst_width'],exact)
 def test_target_order_does_not_change_policy(self):
  d=one_link_design([0,1,4],[2,7,9]);e=one_link_design([0,1,4],[9,2,7]);self.assertEqual(d['chosen_index'],e['chosen_index'])
 def test_not_a_recovered_pairing(self):
  d=one_link_design([0,1,4],[2,7,9]);self.assertTrue(d['never_an_inferred_pairing']);self.assertNotIn('actual_pair',d)
 def test_cost_budget_optimizer(self):
  gs=[.2,.1,.7,.4];sel=greedy_one_per_condition(gs,2,['a','b','c','d'])
  self.assertEqual(sum(gs[i] for i in sel),max(sum(gs[i] for i in q) for q in itertools.combinations(range(4),2)))
 def test_budget_failure(self):
  with self.assertRaises(ValueError):greedy_one_per_condition([.1],2,['a'])
 def test_transport_face_extrema(self):
  x=np.array([0,1,4.]);y=np.array([2,7,9.]);w,v=enumerate_worlds(x,y);c=((x-x.mean())[:,None]*(y-y.mean())).ravel()
  a=np.zeros((6,9))
  for i in range(3):a[i,i*3:(i+1)*3]=1;a[3+i,i::3]=1
  for i,j in itertools.product(range(3),repeat=2):
   b=[(0,None)]*9;b[i*3+j]=(1/3,1/3)
   lo=linprog(c,A_eq=a,b_eq=np.ones(6)/3,bounds=b,method='highs');hi=linprog(-c,A_eq=a,b_eq=np.ones(6)/3,bounds=b,method='highs')
   self.assertTrue(lo.success and hi.success);l,h,_=conditional_interval(w,v,[(i,j)])
   self.assertAlmostEqual(lo.fun,l);self.assertAlmostEqual(-hi.fun,h)
if __name__=='__main__':unittest.main()
