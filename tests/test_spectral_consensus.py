import unittest
import numpy as np
from src.spectral_consensus import spectral_consensus

class SpectralConsensusTests(unittest.TestCase):
    def test_identical_views_recover_common_response(self):
        rng=np.random.default_rng(18);x=rng.normal(size=(8,30));p,r=spectral_consensus(np.stack([x,x]))
        np.testing.assert_allclose(p,x,atol=1e-8)

    def test_opposite_views_cancel(self):
        x=np.random.default_rng(1).normal(size=(7,11));p,_=spectral_consensus(np.stack([x,-x]))
        np.testing.assert_allclose(p,0,atol=1e-12)

    def test_swaps_and_condition_gene_permutations(self):
        rng=np.random.default_rng(32);x=rng.normal(size=(2,9,14));p,_=spectral_consensus(x)
        q,_=spectral_consensus(x[::-1]);np.testing.assert_allclose(q,p,atol=1e-9)
        i=rng.permutation(9);j=rng.permutation(14);q,_=spectral_consensus(x[:,i][:,:,j])
        np.testing.assert_allclose(q,p[i][:,j],atol=1e-9)

    def test_zero_input_and_invalid_inputs(self):
        p,_=spectral_consensus(np.zeros((2,5,8)));np.testing.assert_equal(p,0)
        with self.assertRaises(ValueError):spectral_consensus(np.zeros((3,5,8)))
        with self.assertRaises(ValueError):spectral_consensus(np.full((2,5,8),np.nan))
        with self.assertRaises(ValueError):spectral_consensus(np.zeros((2,5,8)),ridge=0)

if __name__=='__main__':unittest.main()
