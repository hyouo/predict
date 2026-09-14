import unittest
import numpy as np
from src.rna_models import smoother,fit_weights,fit_rna_bank,zero_shot,shared_kernels
from src.rna_metrics import effect_metrics

class RNAStatisticsTests(unittest.TestCase):
    def test_loo_refit_all_intercept_modes(self):
        r=np.random.default_rng(172);x=r.normal(size=(7,3));k=x@x.T;y=r.normal(size=(5,2));ix=np.arange(5)
        for mode in ('none','shrink','free'):
            h,f=smoother(np.repeat(k[None],2,axis=0),ix,.3,mode)
            pred=np.einsum('gij,jg->ig',h,y)
            for leave in range(5):
                keep=np.delete(ix,leave);_,test=smoother(np.repeat(k[None],2,axis=0),keep,.3,mode)
                explicit=np.einsum('gij,jg->ig',test,y[np.arange(5)!=leave])
                np.testing.assert_allclose(pred[leave],explicit[leave],atol=1e-10)

    def test_paired_control_correction_unbiased_in_working_model(self):
        r=np.random.default_rng(640);n=8;b=40000
        treatment=.3*r.normal(size=(n,b));c=r.normal(size=(1,b));d=r.normal(size=(1,b));testctrl=r.normal(size=(1,b))/np.sqrt(2)
        yc=treatment-c;yd=treatment-d;y=(yc+yd)/2;diff=yc-yd
        h=.7*(np.ones((n,n))-np.eye(n))/(n-1)
        cv=np.mean((h@y-y)**2);correction=.5*np.mean((h@diff)*diff)
        true_noisy=np.mean((h@y-(treatment-testctrl))**2)
        self.assertLess(abs(cv+correction-true_noisy),.015)
        self.assertGreater(true_noisy-cv,.60)

    def test_constant_reference_correction_cancels_between_free_intercepts(self):
        n=6;r=np.random.default_rng(66);z=r.normal(size=(n,2));k=z@z.T
        ds=np.full((n,1),2.0);values=[]
        for lam in (.01,.5,10):
            h,_=smoother(k[None],np.arange(n),lam,'free')
            values.append(float(np.mean((h[0]@ds)*ds)))
        np.testing.assert_allclose(values,4.,atol=1e-10)

    def test_adjusted_stacking_uses_linear_correction(self):
        g=np.diag([.2,.5]);c=np.array([.8,-.2]);w,rec=fit_weights(g,c,adjusted=True)
        grid=np.linspace(0,1,10001);risk=.2*grid**2+.5*(1-grid)**2+.5*(.8*grid-.2*(1-grid))
        self.assertAlmostEqual(w[0],float(grid[np.argmin(risk)]),places=3)
        self.assertTrue(rec['success']);self.assertAlmostEqual(w.sum(),1)

    def test_zero_shot_is_indeed_anchor_free(self):
        r=np.random.default_rng(5);s=r.normal(size=(2,10,12));b=r.normal(size=(2,12))+7;bt=np.ones(12)*7
        p,rec=zero_shot(s,b,bt)
        self.assertEqual(rec['target_treatment_rows_used'],0)
        np.testing.assert_allclose(p['source'],s.mean(axis=0))

    def test_centered_metric_exposes_drug_agnostic_prediction(self):
        r=np.random.default_rng(3);t=r.normal(size=(6,10));p=np.repeat(t.mean(axis=0)[None],6,axis=0)
        m=effect_metrics(p,t)
        self.assertAlmostEqual(m['centered_retrieval_top1'],1/6)
        self.assertAlmostEqual(m['mse'],m['centered_mse']+m['bias_mse'])

    def test_bank_block_equivalence_and_finite_output(self):
        r=np.random.default_rng(71);s=r.normal(size=(2,9,7));a=r.normal(size=(4,7));b=a+.2*r.normal(size=a.shape)
        p1,q1=fit_rna_bank(s,np.arange(4),a,b,block_size=3)
        p2,q2=fit_rna_bank(s,np.arange(4),a,b,block_size=7)
        for k in p1:np.testing.assert_allclose(p1[k],p2[k],atol=2e-6)
        np.testing.assert_allclose(q1['control_correction'],q2['control_correction'],atol=1e-12)
        self.assertFalse(q1['target_test_outcomes_received'])

    def test_invalid_shapes_and_missingness_rejected(self):
        r=np.random.default_rng(3);s=r.normal(size=(2,8,4));a=r.normal(size=(4,4));s[0,0,0]=np.nan
        with self.assertRaises(ValueError):fit_rna_bank(s,np.arange(4),a,a)
        with self.assertRaises(ValueError):effect_metrics(a,np.full_like(a,np.nan))

if __name__=='__main__':unittest.main()
