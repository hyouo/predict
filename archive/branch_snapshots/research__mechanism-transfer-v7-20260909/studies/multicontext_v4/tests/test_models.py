import unittest
import numpy as np
from src.models import Predictor, consensus, ridge, project,pca_fit
class Models(unittest.TestCase):
 def setUp(self):
  r=np.random.default_rng(9);self.b=r.normal(size=(8,12));self.y=r.normal(size=(8,20,12))
 def test_missing_not_zero(self):
  y=self.y.copy();y[0,0]=np.nan;f,n,m=consensus(y);np.testing.assert_allclose(f[0],y[1:,0].mean(0));self.assertEqual(n[0],7)
 def test_partial_missing_rejected(self):
  y=self.y.copy();y[0,0,0]=np.nan
  with self.assertRaises(ValueError):Predictor().fit(y,self.b)
 def test_no_target_input_to_fit(self):
  m=Predictor().fit(self.y,self.b);a=m.predict(self.b[:2]);b=m.predict(self.b[:2]);np.testing.assert_array_equal(a,b)
 def test_shared_model_identical_targets(self):
  for kind in ['zero','source_mean','scalar','shared_ridge']:
   a=Predictor(kind=kind).fit(self.y,self.b).predict(self.b[:2]);np.testing.assert_array_equal(a[0],a[1])
 def test_all_candidates_finite(self):
  for kind in ['source_mean','context_rbf','context_linear','scalar','shared_ridge','interaction_centered','interaction_uncentered','context_additive']:
   a=Predictor(kind=kind).fit(self.y,self.b).predict(self.b[:2]);self.assertTrue(np.isfinite(a).all())
 def test_context_row_permutation(self):
  ix=np.array([4,3,7,1,6,2,5,0]);m=Predictor().fit(self.y,self.b);n=Predictor().fit(self.y[ix],self.b[ix]);np.testing.assert_allclose(m.predict(self.b[:2]),n.predict(self.b[:2]),atol=1e-10)
 def test_zero_strength_equals_trunk(self):
  a=Predictor(strength=0).fit(self.y,self.b).predict(self.b[:2]);b=Predictor(kind='shared_ridge').fit(self.y,self.b).predict(self.b[:2]);np.testing.assert_allclose(a,b,atol=1e-12)
 def test_support_center_identity(self):
  r=np.random.default_rng(10);z=r.normal(size=(8,3));a=r.normal(size=(20,12));w=r.normal(size=(20,3,12));y=a[None]+np.einsum('ck,pkg->cpg',z,w)
  ix=np.array([0,2,4,7]);f=y[ix].mean(0);zt=r.normal(size=3);truth=a+np.einsum('k,pkg->pg',zt,w);pred=f+np.einsum('k,pkg->pg',zt-z[ix].mean(0),w);np.testing.assert_allclose(pred,truth,atol=1e-12)
 def test_pairwise_loss_identity(self):
  e=self.y[:,0];pair=np.mean([(e[i]-e[j])**2 for i in range(8) for j in range(i)],axis=0);center=np.mean((e-e.mean(0))**2,axis=0);np.testing.assert_allclose(pair,2*8/7*center)
if __name__=='__main__':unittest.main()
