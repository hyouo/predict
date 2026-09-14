import unittest
import numpy as np
from src.reference_risk import reference_correction, fit_reference_ensembles
from src.rna_models import fit_rna_bank

class ReferenceRiskTests(unittest.TestCase):
    def test_shape_nonfinite_and_bad_reconstruction(self):
        h = np.zeros((2,3,3)); f = np.zeros((4,3,2))
        self.assertEqual(reference_correction(h,f),0.)
        with self.assertRaises(ValueError): reference_correction(h,f[:,:,:1])
        f[0,0,0]=np.nan
        with self.assertRaises(ValueError): reference_correction(h,f)
        s=np.ones((2,5,2));a=np.ones((3,2))
        with self.assertRaises(ValueError):
            fit_reference_ensembles(s,np.arange(3),a,a,np.ones((2,3,2)),np.eye(71))

    def test_unit_order_and_sign_invariance(self):
        r=np.random.default_rng(21);h=r.normal(size=(5,4,4));f=r.normal(size=(3,4,5))
        v=reference_correction(h,f)
        self.assertAlmostEqual(v,reference_correction(h,f[::-1]),places=11)
        f[1]*=-1
        self.assertAlmostEqual(v,reference_correction(h,f),places=11)

    def test_single_unit_recovers_original_ensemble(self):
        r=np.random.default_rng(3);s=r.normal(size=(2,8,6));a=r.normal(size=(4,6));b=a+r.normal(size=(4,6))*.2
        p,d=fit_rna_bank(s,np.arange(4),a,b,block_size=3)
        q,e=fit_reference_ensembles(s,np.arange(4),a,b,(a-b)[None],d['gram'],block_size=3)
        np.testing.assert_allclose(d['control_correction'],e['correction_before_half_factor'],atol=1e-12)
        for k in q:np.testing.assert_allclose(q[k],p[k+'_adjusted'],atol=1e-8)

    def test_design_symmetrization_mean_and_variance(self):
        r=np.random.default_rng(2718);u,n,trials=4,8,50000
        z=r.normal(size=(n,u));h=.7*(np.ones((n,n))-np.eye(n))/(n-1)
        delta=r.normal(size=(u,trials));f=z.T[:,:,None]*delta[:,None,:]
        pooled=f.sum(axis=0)
        a=.5*np.sum(pooled*(h@pooled),axis=0)/n
        b=.5*np.einsum('uig,ij,ujg->g',f,h,f,optimize=True)/n
        expected=.5*np.trace(h@(z@z.T))/n
        self.assertLess(abs(a.mean()-expected),.01)
        self.assertLess(abs(b.mean()-expected),.01)
        self.assertLess(b.var(),a.var()*.6)

if __name__=='__main__': unittest.main()
