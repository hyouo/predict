import unittest
import numpy as np
from src.correlated_anchors import correlated_smoother,fit_correlated_bank
from src.rna_models import smoother

class CorrelatedAnchorTests(unittest.TestCase):
    def test_correlated_loo_equals_direct_refits(self):
        rng=np.random.default_rng(731);z=rng.normal(size=(9,4));k=(z@z.T)[None]
        ix=np.array([0,2,3,5,8]);q=rng.normal(size=(5,3));v=(q@q.T)[None]*.3;y=rng.normal(size=5)
        for mode in ('none','shrink','free'):
            h,final=correlated_smoother(k,ix,.4,mode,v)
            for i in range(5):
                keep=np.delete(np.arange(5),i)
                _,pr=correlated_smoother(k,ix[keep],.4,mode,v[:,keep][:,:,keep])
                self.assertAlmostEqual(float(h[0,i]@y),float(pr[0,ix[i]]@y[keep]),places=10)
            self.assertEqual(h.shape,(1,5,5));self.assertEqual(final.shape,(1,9,5))

    def test_zero_covariance_recovers_frozen_smoother(self):
        rng=np.random.default_rng(2);z=rng.normal(size=(3,8,4));k=z@z.transpose(0,2,1);a=np.array([1,3,5,6])
        for mode in ('none','shrink','free'):
            h,p=correlated_smoother(k,a,.3,mode,np.zeros((3,4,4)));h0,p0=smoother(k,a,.3,mode)
            np.testing.assert_allclose(h,h0,atol=1e-11);np.testing.assert_allclose(p,p0,atol=1e-11)

    def test_free_intercept_constant_equivariance(self):
        rng=np.random.default_rng(42);z=rng.normal(size=(10,3));k=(z@z.T)[None];a=np.arange(5)
        q=rng.normal(size=(5,2));v=(q@q.T)[None]
        h,p=correlated_smoother(k,a,.1,'free',v)
        np.testing.assert_allclose(h.sum(axis=-1),1,atol=1e-10);np.testing.assert_allclose(p.sum(axis=-1),1,atol=1e-10)

    def test_invalid_noise_rejected(self):
        k=np.eye(4)[None];a=np.arange(3)
        with self.assertRaises(ValueError):correlated_smoother(k,a,.1,'none',-np.eye(3)[None])
        bad=np.eye(3)[None];bad[0,0,1]=1
        with self.assertRaises(ValueError):correlated_smoother(k,a,.1,'none',bad)
        with self.assertRaises(ValueError):correlated_smoother(k,a,0,'none',np.eye(3)[None])

    def test_factor_reconstruction_and_model_weights(self):
        rng=np.random.default_rng(9);s=rng.normal(size=(2,8,4));a=rng.normal(size=(4,4));f=rng.normal(size=(2,4,4))*.1;b=a-f.sum(0)
        for mode in ('full','diagonal'):
            p,d=fit_correlated_bank(s,np.arange(4),a,b,f,noise_mode=mode,block_size=2)
            for name,w in d['weights'].items():
                self.assertAlmostEqual(sum(w),1,places=8);self.assertTrue(np.isfinite(p[name]).all())
        with self.assertRaises(ValueError):fit_correlated_bank(s,np.arange(4),a,b,f+1)

if __name__=='__main__':unittest.main()
