import unittest
import numpy as np
from src.covariance import build_kernels
from src.aggregation import bank,aggregate
from src.streaming import fit_streaming

class StreamingTests(unittest.TestCase):
    def test_same_model_with_multiple_readout_blocks(self):
        rng=np.random.default_rng(991);s=rng.normal(size=(4,14,11));x=rng.integers(0,2,size=(14,8));a=np.arange(7);y=rng.normal(size=(7,11));k=build_kernels(s,x)
        p,e,ps=bank(s,x,k,a,y)
        for mode in ['bootstrap','stacking']:
            full=aggregate(p,e,ps,mode);out=np.zeros((14,11),float)
            z=fit_streaming(s,x,a,y,block_size=4,mode=mode,kernels=k,out=out)
            np.testing.assert_allclose(z.prediction,full.prediction,atol=1e-6,rtol=1e-6)
            self.assertIs(z.prediction,out)
    def test_invalid_block_size(self):
        with self.assertRaises(ValueError):fit_streaming(np.ones((4,10,3)),np.ones((10,8)),np.arange(4),np.ones((4,3)),block_size=0)
if __name__=='__main__':unittest.main()
