import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
import pandas as pd
from src.op3_data import Reader, SOURCES, TARGETS, evaluate_truth
from src.op3_models import krr_path, slope_fit, stack_weights, source_kernels, fewshot, zero_shot
from src.op3_metrics import orthogonal_projector, reference_design, contrast_metrics
from tools.op3_pilot import metrics

class OP3PilotTests(unittest.TestCase):
    def setUp(self):self.rng=np.random.default_rng(97)

    def fixture(self,path):
        obs=[]
        for ct in SOURCES+TARGETS:
            for i in range(8):
                obs.append(dict(sample_id=f'{ct}-C{i}',cell_type=ct,plate='P1',well=f'{i+1}C',
                                split='control',is_control=True,psbulk_counts=10))
            for split in ['train','public_test','private_test']:
                obs.append(dict(sample_id=f'{ct}-{split}',cell_type=ct,plate='P1',well=split,
                                split=split,is_control=False,psbulk_counts=10))
        obs=pd.DataFrame(obs);n=len(obs)
        with h5py.File(path,'w') as f:
            og=f.create_group('obs')
            for key in obs:
                values=obs[key].to_numpy()
                if values.dtype.kind=='O':og.create_dataset(key,data=values,dtype=h5py.string_dtype())
                else:og.create_dataset(key,data=values)
            f.create_group('var').create_dataset('symbol',data=np.array(['a','b'],object),dtype=h5py.string_dtype())
            x=f.create_group('X');x.attrs['shape']=(n,2);x.attrs['encoding-type']='csr_matrix'
            x.create_dataset('data',data=np.tile([3,7],n));x.create_dataset('indices',data=np.tile([0,1],n))
            x.create_dataset('indptr',data=np.arange(0,2*n+1,2))
        return Reader(path,verify=False)

    def test_private_hard_rejection_even_with_grant(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.fixture(Path(d)/'a.h5ad');ids=np.flatnonzero(r.obs.split=='private_test')
            with self.assertRaises(PermissionError):r.read(ids,ids,'evaluation')
    def test_public_not_fit(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.fixture(Path(d)/'a.h5ad');ids=np.flatnonzero(r.obs.split=='public_test')
            with self.assertRaises(PermissionError):r.read(ids,ids,'fit')
    def test_final_control_not_granted(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.fixture(Path(d)/'a.h5ad');part=r.partition();c=np.flatnonzero(part==2)
            with self.assertRaises(PermissionError):r.read(c,np.flatnonzero(part!=2),'fit')
    def test_physical_control_partition_shared_across_cells(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.fixture(Path(d)/'a.h5ad')
            for rotation in range(3):
                o=r.obs.copy();o['partition']=r.partition(rotation)
                self.assertEqual(o[o.is_control].groupby(['plate','well']).partition.nunique().max(),1)
                self.assertEqual(sum(o.partition>=0),32)
    def test_reader_normalizes_expected_counts(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.fixture(Path(d)/'a.h5ad');x=r.read([0],[0],'fit')
            np.testing.assert_allclose(x,np.log2(1+np.array([[300000.,700000.]])))
    def test_unknown_library_scale_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            r=self.fixture(Path(d)/'a.h5ad')
            with self.assertRaises(ValueError):r.read([0],[0],'fit','raw')
    def test_no_evaluation_without_seal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'seal.json';p.write_text('{}')
            with self.assertRaises(PermissionError):evaluate_truth(None,p)
    def test_nuisance_projector_annihilates_controls(self):
        h=self.rng.normal(size=(20,4));m,rank=orthogonal_projector(h)
        np.testing.assert_allclose(m@h,0,atol=1e-13)
        np.testing.assert_allclose(m@m,m,atol=1e-13);self.assertEqual(rank,4)
    def test_contrast_score_reference_invariance(self):
        h=self.rng.normal(size=(20,4));y=self.rng.normal(size=(20,9));p=self.rng.normal(size=y.shape)
        c=self.rng.normal(size=(4,9));a=contrast_metrics(y,p,h);b=contrast_metrics(y+h@c,p,h)
        self.assertAlmostEqual(a['plate_orthogonal_mse'],b['plate_orthogonal_mse'],12)
    def test_no_residual_contrast_rejected(self):
        with self.assertRaises(ValueError):orthogonal_projector(np.eye(4))
    def test_metric_bias_decomposition(self):
        y=self.rng.normal(size=(20,9));p=self.rng.normal(size=y.shape);s=metrics(y,p)
        self.assertAlmostEqual(s['mse'],s['bias_mse']+s['centered_mse'],13)
    def test_analytic_loo_matches_refits(self):
        x=self.rng.normal(size=(9,4));k=x@x.T/4;y=self.rng.normal(size=(5,7));a=np.arange(5)
        prior=self.rng.normal(size=(9,7));p,loo=krr_path(k,a,y,prior,[.1,10])
        for i in range(5):
            ix=np.arange(5)!=i;q,_=krr_path(k,a[ix],y[ix],prior,[.1,10])
            np.testing.assert_allclose(q[:,i],loo[:,i],atol=1e-10)
    def test_krr_response_translation(self):
        x=self.rng.normal(size=(9,4));k=x@x.T/4;y=self.rng.normal(size=(5,7));a=np.arange(5)
        prior=np.zeros((9,7));p,loo=krr_path(k,a,y,prior);q,lq=krr_path(k,a,y+2,prior)
        np.testing.assert_allclose(p+2,q,atol=1e-10);np.testing.assert_allclose(loo+2,lq,atol=1e-10)
    def test_gene_slope_approaches_shared(self):
        x=self.rng.normal(size=(7,9,2));y=self.rng.normal(size=(7,9));q=self.rng.normal(size=(4,9,2))
        a=slope_fit(x,y,q,1e10);b=slope_fit(x,y,q,kind='shared')
        np.testing.assert_allclose(a,b,atol=1e-8)
    def test_slope_nonfinite_rejected(self):
        with self.assertRaises(ValueError):slope_fit(np.full((5,8,2),np.nan),np.zeros((5,8)),np.zeros((2,8,2)))
    def test_stacking_weights_form_simplex(self):
        w=stack_weights(self.rng.normal(size=(5,7,8)))
        self.assertAlmostEqual(w.sum(),1,12);self.assertTrue(np.all(w>=0))
    def test_source_kernels_psd(self):
        for k in source_kernels(self.rng.normal(size=(2,12,9))).values():
            self.assertGreater(np.linalg.eigvalsh(k).min(),-1e-10)
    def test_fewshot_no_query_outcome_input(self):
        s=self.rng.normal(size=(2,10,8));y=self.rng.normal(size=(5,8));a=np.arange(5)
        p,pars=fewshot(s,a,y,y+.1)
        self.assertTrue(all(x.shape==(10,8) for x in p.values()))
    def test_zero_shot_needs_no_target_treatments(self):
        b=SimpleNamespace(source=self.rng.normal(size=(2,20,8)),source_validation=self.rng.normal(size=(2,20,8)),
                          compounds=[f'd{i}' for i in range(20)],missing_source=[],
                          baseline={ct:self.rng.normal(size=8) for ct in SOURCES+TARGETS})
        p,pars=zero_shot(b);self.assertEqual(len(p),10)
    def test_reference_design_uses_metadata_counts(self):
        o=pd.DataFrame({'plate':['p1','p2','p2'],'perturbagen':['a','a','b']})
        h=reference_design(o,[0,1,2],['a','b'])
        np.testing.assert_allclose(h,[[.5,.5],[0,1]])

if __name__=='__main__':unittest.main()
