import unittest
import numpy as np
from pathlib import Path
from src.data import load_panel
from src.contrast_design import Action,action_gain,posterior_update,risk,continuous_replication_allocation
from run_sentinel_probe import probe_sets,choose_probe
from src.covariance import cue_kernel
from src.kriging import acquire

class PracticalDesignTests(unittest.TestCase):
    def test_probe_reserves_all_candidate_test_rows(self):
        p=load_panel(Path(__file__).resolve().parents[1]/'data/hepatocyte_signaling_raw.csv')
        for j in range(8):
            pool,q,probes=probe_sets(p.cues,j)
            self.assertTrue(set(pool).isdisjoint(q));self.assertTrue(set(probes).isdisjoint(q))
            self.assertEqual(set(q)|set(probes),set(np.flatnonzero(p.cues[:,j])))
            k=cue_kernel(p.cues);a,_=acquire(k,pool,q,17,criterion='proper')
            chosen,gains=choose_probe(k,a[:16],q,probes)
            self.assertIn(chosen,probes);self.assertTrue(np.isfinite(gains).all())
            self.assertTrue(set(np.r_[a[:16],chosen]).isdisjoint(q))
    def test_nonfinite_noise_cost_rejected(self):
        for nv,cost in [(np.nan,1.),(1.,np.nan),(np.inf,1.)]:
            a=Action('invalid',np.ones(3),nv,cost)
            with self.assertRaises(ValueError):action_gain(np.eye(3),np.eye(3),a)
    def test_nonfinite_contrasts_mean_and_allocation_rejected(self):
        with self.assertRaises(ValueError):risk(np.eye(2),np.array([[np.nan,1.]]))
        with self.assertRaises(ValueError):risk(np.eye(2),np.eye(2),[1.,np.nan])
        with self.assertRaises(ValueError):posterior_update([np.nan,0.],np.eye(2),Action('well',np.ones(2),1.),0.)
        with self.assertRaises(ValueError):continuous_replication_allocation([1.],[1.],1.,np.nan)
        with self.assertRaises(ValueError):continuous_replication_allocation([],[],1.,10.)
    def test_empty_or_incompatible_contrasts_rejected(self):
        for L in [np.ones((0,2)),np.ones((2,3)),np.ones(2)]:
            with self.assertRaises(ValueError):action_gain(np.eye(2),L,Action('well',np.ones(2),1.))

if __name__=='__main__':unittest.main()
