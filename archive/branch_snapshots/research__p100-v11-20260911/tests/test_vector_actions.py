import unittest
import numpy as np
from src.contrast_design import Action,VectorAction,posterior_update,vector_posterior_update,vector_action_gain,multiplex_actions,risk

class VectorActionTests(unittest.TestCase):
    def test_correlated_readout_gain_equals_actual_variance_drop(self):
        r=np.random.default_rng(41);z=r.normal(size=(6,6));C=z@z.T;H=r.normal(size=(3,6));L=r.normal(size=(4,6))
        noise=np.eye(3)+.3*np.ones((3,3));a=VectorAction('single_well',H,noise,2.)
        gain=vector_action_gain(C,L,a);_,C1=vector_posterior_update(np.zeros(6),C,a,np.ones(3))
        self.assertAlmostEqual(gain['absolute_risk_reduction'],risk(C,L)-risk(C1,L),places=10)
        self.assertAlmostEqual(gain['risk_reduction_per_cost'],gain['absolute_risk_reduction']/2.)
    def test_independent_readouts_equal_sequential_updates(self):
        C=np.eye(5);m=np.zeros(5);H=np.array([[1,1,0,0,0],[0,0,1,1,0]],float);R=np.diag([.2,.8]);y=np.array([.3,-1.])
        mm,CC=vector_posterior_update(m,C,VectorAction('well',H,R),y)
        for i in range(2):m,C=posterior_update(m,C,Action('component',H[i],R[i,i]),y[i])
        np.testing.assert_allclose(mm,m,atol=1e-12);np.testing.assert_allclose(CC,C,atol=1e-12)
    def test_well_cost_not_multiplied_by_readout_count(self):
        acts=multiplex_actions(3,2,np.eye(2),np.eye(2),[2,3,4,5])
        self.assertEqual(len(acts),4);self.assertEqual(acts[1].cost,3)
        self.assertEqual(acts[1].design.shape,(2,8))
        for i,g in enumerate([0,4]):
            self.assertEqual(acts[1].design[i,g],1);self.assertEqual(acts[1].design[i,g+1],1)
    def test_singular_noise_is_rejected(self):
        with self.assertRaises(ValueError):vector_action_gain(np.eye(3),np.eye(3),VectorAction('bad',np.ones((2,3)),np.ones((2,2))))
if __name__=='__main__':unittest.main()
