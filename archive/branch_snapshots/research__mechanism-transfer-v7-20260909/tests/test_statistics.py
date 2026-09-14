import unittest
import numpy as np
from src.covariance import build_kernels,normalize_kernel,kernel_path,cue_kernel,tune
from src.aggregation import path,output_kernels,aggregate,bank
from src.kriging import uk_covariance,uk_predict,acquire
from src.contrast_design import Action,posterior_update,risk,action_gain,control_and_treatment_actions,continuous_replication_allocation
from src.weighted_risk import source_weights
from src.deployment_risk import fit_deployment

class StatisticsTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(707)
        self.source=self.rng.normal(size=(4,12,3));self.x=self.rng.integers(0,2,size=(12,8)).astype(float)
        self.a=np.array([0,1,3,5,8,10]);self.y=self.rng.normal(size=(6,3))
        self.k=build_kernels(self.source,self.x)
    def test_kernel_psd(self):
        for k in self.k.values(): self.assertGreater(np.linalg.eigvalsh(k).min(),-1e-8)
    def test_disagreement_constant_offsets_invariant(self):
        b=self.source+self.rng.normal(size=(4,1,3))*10
        np.testing.assert_allclose(self.k['disagreement'],build_kernels(b,self.x)['disagreement'],atol=1e-12)
    def test_zero_disagreement_finite(self):
        s=np.repeat(self.source[:1],4,axis=0)
        np.testing.assert_array_equal(build_kernels(s,self.x)['disagreement'],np.zeros((12,12)))
    def test_normalized_mean_diagonal(self):
        k=normalize_kernel(self.k['cue_original']);self.assertAlmostEqual(np.trace(k)/len(k),1.)
        np.testing.assert_allclose(k.sum(0),0,atol=1e-12)
    def test_batched_equals_shared(self):
        p,e=path(self.k['cue_original'],self.a,self.y)
        pp,ee=kernel_path(self.k['cue_original'],self.a,self.y)
        np.testing.assert_allclose(p,pp,atol=1e-10);np.testing.assert_allclose(e*e,ee,atol=1e-9)
    def test_loo_matches_refitting(self):
        K=self.k['cue_original'];p,e=path(K,self.a,self.y,lambdas=[.1])
        for j in range(len(self.a)):
            ids=np.delete(np.arange(len(self.a)),j)
            pr=uk_predict(K,self.a[ids],self.y[ids],.1)
            np.testing.assert_allclose(self.y[j]-e[0,j],pr[self.a[j]],atol=1e-10)
    def test_output_specific_loo(self):
        K=output_kernels(self.source);p,e=path(K,self.a,self.y,lambdas=[.1])
        for g in range(3):
            for j in range(len(self.a)):
                ids=np.delete(np.arange(len(self.a)),j)
                pr=uk_predict(K[g],self.a[ids],self.y[ids],.1)
                np.testing.assert_allclose(self.y[j,g]-e[0,j,g],pr[self.a[j],g],atol=1e-10)
    def test_output_marginal_equals_primal_ridge(self):
        K=output_kernels(self.source,'marginal');pp,_=path(K,self.a,self.y,lambdas=[.1])
        for g in range(3):
            z=self.source[:,:,g].T;z=z-z.mean(0)
            z=z/np.sqrt(np.mean(np.sum(z*z,1)))
            za=z[self.a];xc=za-za.mean(0);yc=self.y[:,g]-self.y[:,g].mean()
            coef=np.linalg.solve(xc.T@xc+.1*np.eye(4),xc.T@yc)
            p=self.y[:,g].mean()+(z-za.mean(0))@coef
            np.testing.assert_allclose(pp[0,:,g],p,atol=1e-10)
    def test_uk_equal_kernel_prediction(self):
        pr,_=path(self.k['cue_original'],self.a,self.y,lambdas=[.01,.1,1.])
        for j,lam in enumerate([.01,.1,1.]):np.testing.assert_allclose(pr[j],uk_predict(self.k['cue_original'],self.a,self.y,lam),atol=1e-10)
    def test_uk_psd(self):
        c=uk_covariance(self.k['cue_original'],self.a);self.assertGreater(np.linalg.eigvalsh(c).min(),-1e-8)
    def test_uk_sequential_covariance_identity(self):
        K=self.k['cue_original'];a=self.a[:3];i=7;c=uk_covariance(K,a,.1)
        nxt=c-np.outer(c[:,i],c[:,i])/(c[i,i]+.1)
        np.testing.assert_allclose(nxt,uk_covariance(K,np.r_[a,i],.1),atol=1e-10)
    def test_uk_constant_kernel_invariant(self):
        np.testing.assert_allclose(uk_covariance(self.k['cue_original'],self.a),uk_covariance(self.k['cue_original']+9,self.a),atol=1e-10)
    def test_design_disjoint(self):
        q=np.array([9,10,11]);pool=np.arange(9)
        for typ in ['proper','universal']:
            a,_=acquire(self.k['cue_original'],pool,q,5,criterion=typ)
            self.assertEqual(len(set(a)),5);self.assertEqual(len(np.intersect1d(a,q)),0)
    def test_bad_design_raises(self):
        with self.assertRaises(ValueError):acquire(self.k['cue_original'],np.arange(8),np.array([7,8]),3)
    def test_bank_prediction_reference_equivariance(self):
        p,e,pars=bank(self.source,self.x,self.k,self.a,self.y)
        shift=np.array([3.,-2.,4.]);p1,e1,_=bank(self.source,self.x,self.k,self.a,self.y+shift)
        np.testing.assert_allclose(p1,p+shift,atol=1e-10);np.testing.assert_allclose(e1,e,atol=1e-9)
    def test_aggregation_simplex(self):
        p,e,pars=bank(self.source,self.x,self.k,self.a,self.y,families=('cue',))
        for mode in ['bootstrap','stacking','uniform']:
            fit=aggregate(p,e,pars,mode);w=np.array(fit.parameters['weights'])
            self.assertAlmostEqual(w.sum(),1);self.assertTrue((w>=0).all());self.assertTrue(np.isfinite(fit.prediction).all())
    def test_bootstrap_structural_tie_prior(self):
        # Identical validation losses cannot rank the two penalties. Use declared
        # strong-penalty tie rule rather than numerical round-off.
        e=np.ones((2,5,1));e[0]-=1e-13
        pp=np.stack([np.zeros((8,1)),np.ones((8,1))]);pars=[{'lambda':.001},{'lambda':1000.}]
        z=aggregate(pp,e,pars,'bootstrap')
        np.testing.assert_allclose(z.prediction,1.,atol=0.)
    def test_weight_reference_invariance(self):
        w,_=source_weights(self.source,self.a,self.y)
        ss=self.source+self.rng.normal(size=(4,1,3))*5
        w1,_=source_weights(ss,self.a,self.y+np.array([2,-1,3]))
        np.testing.assert_allclose(w,w1,atol=1e-12)
    def test_target_query_not_passed_meta(self):
        # Only target anchors are in the public fitting signature. A separate
        # target query array can change arbitrarily with no effect on the fit.
        q=np.array([2,4]);z=fit_deployment(self.source,self.x,self.a,self.y,q,self.k,families=('cue',),weights=(.5,))
        target_query=self.rng.normal(size=(2,3))*1e6
        z1=fit_deployment(self.source,self.x,self.a,self.y.copy(),q,self.k,families=('cue',),weights=(.5,))
        np.testing.assert_array_equal(z['meta_0.5'].prediction,z1['meta_0.5'].prediction)
    def test_contrast_action_exact_gain(self):
        Z=self.rng.normal(size=(6,6));C=Z@Z.T;L=self.rng.normal(size=(3,6));h=self.rng.normal(size=6)
        act=Action('well',h,.2,2.)
        gain=action_gain(C,L,act);_,C1=posterior_update(np.zeros(6),C,act,.4)
        self.assertAlmostEqual(gain['absolute_risk_reduction'],risk(C,L)-risk(C1,L),places=10)
        self.assertAlmostEqual(gain['risk_reduction_per_cost'],gain['absolute_risk_reduction']/2)
    def test_posterior_psd(self):
        C=np.eye(4);a=Action('well',np.array([1,0,1,0]),.1)
        _,C1=posterior_update(np.zeros(4),C,a,1.)
        self.assertGreaterEqual(np.linalg.eigvalsh(C1).min(),-1e-12)
    def test_common_control_covariance(self):
        nt=4;v=np.diag([.5]+[1.]*nt);D=np.zeros((nt,nt+1));D[:,0]=-1;D[:,1:]=np.eye(nt)
        np.testing.assert_allclose(D@v@D.T,np.eye(nt)+.5*np.ones((nt,nt)))
    def test_control_cancels_direct_pairwise_contrast(self):
        L=np.array([[1.,-1.]]);V=np.eye(2)+7*np.ones((2,2))
        np.testing.assert_allclose(L@V@L.T,L@np.eye(2)@L.T)
    def test_symmetric_balanced_control_no_pairwise_gain(self):
        C=np.diag([10.,1.,1.]);m=np.zeros(3);actions=control_and_treatment_actions(2,1.,[1.,1.])
        for a in actions[1:]:m,C=posterior_update(m,C,a,0.)
        L=np.array([[0.,1.,-1.]])
        self.assertLess(abs(action_gain(C,L,actions[0])['absolute_risk_reduction']),1e-15)
    def test_replication_allocation(self):
        z=continuous_replication_allocation(np.ones(16)/16,np.ones(16),1.,96)
        self.assertAlmostEqual(z['control'],19.2)
        np.testing.assert_allclose(z['treatments'],4.8)
        self.assertAlmostEqual(z['control']/z['treatments'][0],4.)
    def test_allocation_budget_cost(self):
        c=np.array([2.,1.,3.]);z=continuous_replication_allocation([.2,.8],[1.,2.],.5,30,c)
        self.assertAlmostEqual(c@np.r_[z['control'],z['treatments']],30.)
    def test_invalid_covariance_raises(self):
        with self.assertRaises(ValueError):risk(np.diag([1.,-2.]),np.eye(2))

if __name__=='__main__':unittest.main(verbosity=2)
