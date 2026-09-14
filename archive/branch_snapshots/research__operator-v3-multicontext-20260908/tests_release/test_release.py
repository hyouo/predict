"""v0.1.0 unit/CLI tests. Synthetic arrays here are test fixtures, not evidence."""
from __future__ import annotations
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from perturb_predict import TransferModel
from perturb_predict.model import coefficients, PENALTIES, RHOS
from perturb_predict.io import (write_bundle, load_bundle, save_archive, load_archive,
                               run_directory, alignment, digest, identifiers)
from perturb_predict.op3 import guard_rows, prepare, fetch, DATA_SHA256
from perturb_predict.cli import train, predict, evaluate, read_predictions, main, benchmark, score_op3
from perturb_predict.evaluation import metrics

CONTRACT = {'effect_space': 'test_fixture_effect', 'source_context': 'source', 'target_context': 'target',
            'organism': 'test', 'dose': 'fixed', 'time': 'fixed', 'gene_namespace': 'test_gene'}


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        rng = np.random.default_rng(20260907)
        self.x = rng.normal(size=(6, 12))
        self.y = .7*self.x + rng.normal(size=(6, 12))*.1
        self.q = rng.normal(size=(3, 12))
        self.ids = np.array([f'd{i}' for i in range(6)])
        self.qids = np.array(['q0', 'q1', 'q2'])
        self.genes = np.array([f'g{i}' for i in range(12)])

    def fitted(self):
        return TransferModel().fit(self.x, self.y, ids=self.ids, genes=self.genes, contract=CONTRACT)

    def bundles(self):
        write_bundle(self.root/'train.npz', 'train', self.ids, self.genes, CONTRACT, source=self.x, target=self.y)
        write_bundle(self.root/'query.npz', 'query', self.qids, self.genes, CONTRACT, source=self.q)
        return self.root/'train.npz', self.root/'query.npz'

    def test_coefficients_independent_normal_equation(self):
        b, a = coefficients(self.x, self.y, 10.)
        xc, yc = self.x-self.x.mean(0), self.y-self.y.mean(0)
        total = float(np.sum(xc**2))
        prior = np.sum(xc*yc)/(total+max(total*1e-8, 1e-10))
        lam = 10*max(float(np.median(np.sum(xc**2, 0))), 1e-8)
        for j in range(12):
            direct = np.linalg.solve(np.array([[np.dot(xc[:,j], xc[:,j])+lam]]),
                                     np.array([np.dot(xc[:,j], yc[:,j])+lam*prior]))[0]
            self.assertAlmostEqual(b[j], direct, places=13)
        np.testing.assert_allclose(a, self.y.mean(0)-b*self.x.mean(0))

    def test_full_compound_loo_grid(self):
        m = self.fitted()
        self.assertEqual(len(m.report['candidates']), 25)
        for c in m.report['candidates']:
            preds = []
            for i in range(6):
                keep = np.arange(6) != i
                b, a = coefficients(self.x[keep], self.y[keep], c['penalty'])
                preds.append(self.x[i]*b+c['rho']*a)
            self.assertAlmostEqual(c['loo_mse'], float(np.mean((np.array(preds)-self.y)**2)), places=14)

    def test_roundtrip_bitwise(self):
        m = self.fitted(); m.save(self.root/'m.npz'); n = TransferModel.load(self.root/'m.npz')
        kwargs = dict(genes=self.genes, ids=self.qids, contract=CONTRACT)
        np.testing.assert_array_equal(m.predict(self.q, **kwargs)[0], n.predict(self.q, **kwargs)[0])

    def test_reordered_genes(self):
        m = self.fitted(); order = np.arange(12)[::-1]
        p, _ = m.predict(self.q, genes=self.genes, ids=self.qids, contract=CONTRACT)
        q, _ = m.predict(self.q[:,order], genes=self.genes[order], ids=self.qids, contract=CONTRACT)
        np.testing.assert_array_equal(p, q)

    def test_reordered_query_rows(self):
        m = self.fitted(); order = [2,0,1]
        p, _ = m.predict(self.q, genes=self.genes, ids=self.qids, contract=CONTRACT)
        q, _ = m.predict(self.q[order], genes=self.genes, ids=self.qids[order], contract=CONTRACT)
        np.testing.assert_array_equal(p[order], q)

    def test_missing_gene_rejected(self):
        with self.assertRaises(ValueError):
            self.fitted().predict(self.q[:,:-1], genes=self.genes[:-1], ids=self.qids, contract=CONTRACT)

    def test_extra_gene_rejected(self):
        with self.assertRaises(ValueError): alignment(np.append(self.genes,'extra'), self.genes,'genes')

    def test_duplicate_gene_rejected(self):
        g = self.genes.copy(); g[-1] = g[0]
        with self.assertRaises(ValueError): self.fitted().predict(self.q, genes=g, ids=self.qids, contract=CONTRACT)

    def test_duplicate_compound_rejected(self):
        with self.assertRaises(ValueError): identifiers(['d','d'], 'ids')

    def test_empty_identifier_rejected(self):
        with self.assertRaises(ValueError): identifiers(['d',' '], 'ids')

    def test_train_query_overlap_rejected(self):
        with self.assertRaises(ValueError):
            self.fitted().predict(self.q, genes=self.genes, ids=['d0','q1','q2'], contract=CONTRACT)

    def test_context_contract_mismatch(self):
        c = dict(CONTRACT, target_context='other')
        with self.assertRaises(ValueError): self.fitted().predict(self.q, genes=self.genes, ids=self.qids, contract=c)

    def test_units_contract_mismatch(self):
        c = dict(CONTRACT, effect_space='raw_counts')
        with self.assertRaises(ValueError): self.fitted().predict(self.q, genes=self.genes, ids=self.qids, contract=c)

    def test_minimum_calibration(self):
        with self.assertRaises(ValueError):
            TransferModel().fit(self.x[:3], self.y[:3], ids=self.ids[:3], genes=self.genes, contract=CONTRACT)

    def test_nonfinite_train_rejected(self):
        for bad in (np.nan, np.inf):
            x = self.x.copy(); x[0,0] = bad
            with self.assertRaises(ValueError): TransferModel().fit(x,self.y,ids=self.ids,genes=self.genes,contract=CONTRACT)

    def test_nonfinite_query_rejected(self):
        q = self.q.copy(); q[0,0] = np.nan
        with self.assertRaises(ValueError): self.fitted().predict(q,genes=self.genes,ids=self.qids,contract=CONTRACT)

    def test_constant_source_finite(self):
        m = TransferModel().fit(np.ones_like(self.x),self.y,ids=self.ids,genes=self.genes,contract=CONTRACT)
        self.assertTrue(np.isfinite(m.slopes).all())

    def test_structural_tie_prefers_shrinkage(self):
        m = TransferModel().fit(np.zeros_like(self.x),np.zeros_like(self.y),ids=self.ids,genes=self.genes,contract=CONTRACT)
        self.assertEqual(m.report['selected']['rho'],0.)
        self.assertEqual(m.report['selected']['penalty'],1000.)

    def test_unfitted_rejected(self):
        with self.assertRaises(ValueError): TransferModel().save(self.root/'m.npz')

    def test_missing_file_rejected(self):
        with self.assertRaises(FileNotFoundError): load_bundle(self.root/'missing.npz','train')

    def test_object_arrays_prohibited(self):
        with self.assertRaises(ValueError): save_archive(self.root/'a.npz',{},x=np.array([{}],dtype=object))

    def test_loading_object_array_prohibited(self):
        np.savez(self.root/'a.npz', metadata_json=np.array('{}'), x=np.array([{}],dtype=object))
        with self.assertRaises(ValueError): load_archive(self.root/'a.npz')

    def test_unknown_schema_rejected(self):
        save_archive(self.root/'a.npz',{'schema_version':99},x=np.zeros(2))
        with self.assertRaises(ValueError): load_archive(self.root/'a.npz')

    def test_query_cannot_contain_labels(self):
        with self.assertRaises(ValueError):
            write_bundle(self.root/'q.npz','query',self.qids,self.genes,CONTRACT,source=self.q,target=self.q)

    def test_archive_written_once(self):
        self.fitted().save(self.root/'a.npz'); before=digest(self.root/'a.npz')
        with self.assertRaises(FileExistsError): self.fitted().save(self.root/'a.npz')
        self.assertEqual(before,digest(self.root/'a.npz'))

    def test_run_written_once(self):
        with run_directory(self.root/'run'): pass
        with self.assertRaises(FileExistsError):
            with run_directory(self.root/'run'): pass

    def test_failed_run_not_complete(self):
        with self.assertRaises(ValueError):
            with run_directory(self.root/'run'): raise ValueError('test_failure')
        self.assertTrue((self.root/'run/FAILED.json').exists())
        self.assertFalse((self.root/'run/COMPLETE.json').exists())

    def test_private_locked_even_for_scoring(self):
        o={'split':np.array(['private_test']), 'is_control':np.array([False]), 'well':np.array(['A1'])}
        for purpose in ('prepare','score'):
            with self.assertRaises(PermissionError): guard_rows(o, np.array([0]), purpose)

    def test_public_locked_during_prepare(self):
        o={'split':np.array(['public_test']), 'is_control':np.array([False]), 'well':np.array(['A1'])}
        with self.assertRaises(PermissionError): guard_rows(o,np.array([0]),'prepare')
        np.testing.assert_array_equal(guard_rows(o,np.array([0]),'score'),[0])

    def test_reserved_controls_locked(self):
        o={'split':np.array(['control']), 'is_control':np.array([True]), 'well':np.array(['E1'])}
        with self.assertRaises(PermissionError): guard_rows(o,np.array([0]),'prepare')

    def test_invalid_row_index(self):
        o={'split':np.array(['train']), 'is_control':np.array([False]), 'well':np.array(['A1'])}
        for ids in ([-1],[1],[0,0],[.5]):
            with self.assertRaises(ValueError): guard_rows(o,np.array(ids),'prepare')

    def test_conditional_requires_explicit_ack(self):
        with self.assertRaises(ValueError): prepare('missing',self.root/'out')
        with self.assertRaises(ValueError): benchmark('missing',self.root/'out')
        with self.assertRaises(ValueError): score_op3('missing','missing',self.root/'out')

    def test_invalid_cache_not_overwritten(self):
        p=self.root/'raw.h5ad';p.write_bytes(b'not_real_data')
        with self.assertRaises(ValueError): fetch(p)
        self.assertEqual(p.read_bytes(),b'not_real_data')

    def test_failed_download_cleans_partial(self):
        with patch('perturb_predict.op3.urlopen',side_effect=OSError('offline')):
            with self.assertRaises(OSError): fetch(self.root/'raw.h5ad')
        self.assertFalse((self.root/'raw.h5ad').exists())
        self.assertFalse((self.root/'raw.h5ad.partial').exists())

    def test_download_wrong_hash(self):
        fake=io.BytesIO(b'wrong');fake.geturl=lambda:'https://example.test'
        with patch('perturb_predict.op3.urlopen',return_value=fake):
            with self.assertRaises(ValueError): fetch(self.root/'raw.h5ad')
        self.assertFalse((self.root/'raw.h5ad').exists())

    def test_metric_error_decomposition(self):
        m=metrics(self.x,self.y)
        self.assertAlmostEqual(m['mse'],m['centered_mse']+m['bias_mse'],places=14)

    def test_metric_ties_not_first_index(self):
        m=metrics(np.zeros_like(self.x),self.y)
        self.assertAlmostEqual(m['retrieval_top1'],1/len(self.x))

    def test_standalone_pipeline_and_id_alignment(self):
        t,q=self.bundles();train(t,self.root/'m');predict(self.root/'m/model.npz',q,self.root/'p')
        # Evaluation labels deliberately permute both rows and genes.
        order=[2,0,1];g=np.arange(12)[::-1]
        write_bundle(self.root/'truth.npz','truth',self.qids[order],self.genes[g],CONTRACT,target=(.7*self.q)[order][:,g])
        scores=evaluate(self.root/'p/predictions.npz',self.root/'truth.npz',self.root/'s')
        _,a=read_predictions(self.root/'p/predictions.npz')
        self.assertAlmostEqual(scores['model']['mse'],np.mean((a['prediction']-.7*self.q)**2),places=14)
        self.assertTrue((self.root/'p/predictions.csv').exists())

    def test_incomplete_predictions_fail(self):
        t,q=self.bundles();train(t,self.root/'m');predict(self.root/'m/model.npz',q,self.root/'p')
        (self.root/'p/COMPLETE.json').unlink()
        with self.assertRaises(ValueError): read_predictions(self.root/'p/predictions.npz')

    def test_tampered_predictions_fail(self):
        t,q=self.bundles();train(t,self.root/'m');predict(self.root/'m/model.npz',q,self.root/'p')
        with (self.root/'p/predictions.npz').open('ab') as f:f.write(b'changed')
        with self.assertRaises(ValueError): read_predictions(self.root/'p/predictions.npz')

    def test_cli_actionable_missing_file(self):
        with contextlib.redirect_stderr(io.StringIO()) as err:
            code=main(['train','--data',str(self.root/'none.npz'),'--out',str(self.root/'m')])
        self.assertEqual(code,2);self.assertIn('FileNotFoundError',err.getvalue())

    def test_training_shape_rejected(self):
        with self.assertRaises(ValueError):
            TransferModel().fit(self.x,self.y[:,:-1],ids=self.ids,genes=self.genes,contract=CONTRACT)

    def test_model_coefficient_corruption(self):
        m=self.fitted();m.save(self.root/'m.npz');meta,a=load_archive(self.root/'m.npz')
        a['slopes']=np.array([np.nan]*12);save_archive(self.root/'bad.npz',meta,**a)
        with self.assertRaises(ValueError): TransferModel.load(self.root/'bad.npz')

if __name__=='__main__': unittest.main()
