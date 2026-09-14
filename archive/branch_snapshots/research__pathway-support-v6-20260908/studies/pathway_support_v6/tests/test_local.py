import unittest
import numpy as np
from src.kernels import SupportTransport,shuffled_membership
from src.local import LocalTransport,neighbourhood
class TestLocal(unittest.TestCase):
 def setUp(self):self.r=np.random.default_rng(14)
 def test_uniform_equals_global(self):
  b=self.r.normal(size=(7,8));y=self.r.normal(size=(7,9,8));y[2,:3]=np.nan;q=self.r.normal(size=(3,8))
  x=LocalTransport(np.ones((8,8))/8).fit(y,b).predict(q,.4,.75,.5)
  z=SupportTransport().fit(y,b).predict(q,.4,.75,.5)
  np.testing.assert_allclose(x,z,atol=2e-13)
 def test_prior_rowstochastic(self):
  a=self.r.uniform(size=(10,6))<.4;a[3]=0
  w=neighbourhood(a);np.testing.assert_allclose(w.sum(1),1);self.assertTrue((w>=0).all());np.testing.assert_allclose(w[3],.1)
 def test_prior_equivariance(self):
  a=self.r.uniform(size=(40,6))<.4;s,p=shuffled_membership(a,3)
  np.testing.assert_allclose(neighbourhood(s),neighbourhood(a)[np.ix_(p,p)])
 def test_local_single_gene(self):
  b=self.r.normal(size=(7,1));y=self.r.normal(size=(7,9,1));q=self.r.normal(size=(3,1))
  np.testing.assert_allclose(LocalTransport(np.ones((1,1))).fit(y,b).predict(q,.4),SupportTransport().fit(y,b).predict(q,.4),atol=1e-12)
 def test_local_outcome_shift(self):
  b=self.r.normal(size=(7,8));y=self.r.normal(size=(7,9,8));q=self.r.normal(size=(3,8));w=np.eye(8)
  np.testing.assert_allclose(LocalTransport(w).fit(y+2,b).predict(q,.4),LocalTransport(w).fit(y,b).predict(q,.4)+2,atol=1e-12)
if __name__=='__main__':unittest.main()
