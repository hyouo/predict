import unittest
import numpy as np
from src.kernels import affine_weights,Geometry,SupportTransport,shuffled_membership

class TestSupport(unittest.TestCase):
 def setUp(self):self.rng=np.random.default_rng(5)
 def test_primal_dual(self):
  x=self.rng.normal(size=(8,4));q=self.rng.normal(size=(3,4));y=self.rng.normal(size=(8,6));lam=.7
  w=affine_weights(x@x.T,q@x.T,lam)
  xc=x-x.mean(0);yc=y-y.mean(0)
  beta=np.linalg.solve(xc.T@xc+lam*np.eye(4),xc.T@yc)
  np.testing.assert_allclose(w@y,y.mean(0)+(q-x.mean(0))@beta,atol=1e-12)
 def test_affine_weights(self):
  x=self.rng.normal(size=(6,4));q=self.rng.normal(size=(3,4))
  w=affine_weights(x@x.T,q@x.T,.1)
  np.testing.assert_allclose(w.sum(1),1,atol=1e-12)
 def test_feature_translation(self):
  x=self.rng.normal(size=(6,4));q=self.rng.normal(size=(3,4));shift=self.rng.normal(size=4)
  np.testing.assert_allclose(affine_weights(x@x.T,q@x.T,.3),affine_weights((x+shift)@(x+shift).T,(q+shift)@(x+shift).T,.3),atol=1e-12)
 def test_outcome_shift(self):
  x=self.rng.normal(size=(6,4));q=self.rng.normal(size=(3,4));y=self.rng.normal(size=(6,2));w=affine_weights(x@x.T,q@x.T,.3)
  np.testing.assert_allclose(w@(y+4),w@y+4,atol=1e-12)
 def test_unseen_outcome_mask(self):
  b=self.rng.normal(size=(5,4));y=self.rng.normal(size=(5,6,4));y[2,3]=np.nan
  m=SupportTransport().fit(y,b);p=m.predict(b[:2],.5)
  self.assertTrue(np.isfinite(p).all())
 def test_missing_whole_drug(self):
  b=self.rng.normal(size=(5,4));y=self.rng.normal(size=(5,6,4));y[:,3]=np.nan
  p=SupportTransport().fit(y,b).predict(b[:2],.5);self.assertTrue(np.isnan(p[:,3]).all())
 def test_partial_missing_rejected(self):
  b=self.rng.normal(size=(5,4));y=self.rng.normal(size=(5,6,4));y[1,1,1]=np.nan
  with self.assertRaises(ValueError):SupportTransport().fit(y,b)
 def test_homoskedastic_equivalence(self):
  b=self.rng.normal(size=(5,4));y=self.rng.normal(size=(5,6,4));v=np.ones((5,6));q=b[:2]
  np.testing.assert_allclose(SupportTransport().fit(y,b).predict(q,.5),SupportTransport(noise_gamma=1).fit(y,b,v).predict(q,.5),atol=1e-12)
 def test_permutation_preserves_structure(self):
  a=self.rng.uniform(size=(100,10))<.3;z,p=shuffled_membership(a,9101)
  np.testing.assert_array_equal(z.sum(1),a.sum(1));np.testing.assert_array_equal(z.sum(0),a.sum(0))
  np.testing.assert_array_equal(z.astype(int).T@z,a.astype(int).T@a)
 def test_kernel_psd(self):
  b=self.rng.normal(size=(6,10));a=self.rng.uniform(size=(10,5))<.5
  k=Geometry(.5,a).fit(b).k;self.assertGreater(np.linalg.eigvalsh(k).min(),-1e-10)
 def test_bad_lambda(self):
  with self.assertRaises(ValueError):affine_weights(np.eye(3),np.ones((2,3)),0)
 def test_bad_baseline(self):
  with self.assertRaises(ValueError):Geometry().fit(np.array([[1,np.nan],[2,3]]))
 def test_single_support(self):
  np.testing.assert_allclose(affine_weights(np.array([[1.]]),np.array([[3.],[2.]]),.01),np.ones((2,1)))
if __name__=='__main__':unittest.main()
