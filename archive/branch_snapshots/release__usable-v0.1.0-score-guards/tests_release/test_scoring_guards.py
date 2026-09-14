"""Review regressions: scoring completeness and audited truth-access lifecycle."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch
import numpy as np
from perturb_predict.cli import score_op3, evaluate
from perturb_predict.op3 import Reader, SOURCE, NK, TARGETS, DATA_SHA256, contract, expected_public_queries


class ScoringGuardTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ct = TARGETS['B_cells']
        self.ids = np.array(['q0', 'q1', 'q2'])
        self.genes = np.array(['g0', 'g1'])
        self.meta = {'contract': contract(self.ct),
                     'input_provenance': {'dataset_sha256': DATA_SHA256}}
        self.arrays = {'ids': self.ids, 'genes': self.genes,
                       'prediction': np.ones((3,2)), 'source_copy': np.zeros((3,2))}
        self.reader = MagicMock()
        self.reader.genes = self.genes
        self.reader.log = []
        self.reader.ids.return_value = np.array([0,1,2])
        self.reader.effects.return_value = ({d: np.ones(2) for d in self.ids}, np.arange(3))

    def invoke(self, **kwargs):
        arrays = kwargs.pop('arrays', self.arrays)
        with patch('perturb_predict.cli.read_predictions',return_value=(self.meta, arrays)), \
             patch('perturb_predict.cli.digest',return_value='unit_test_hash'), \
             patch('perturb_predict.op3.Reader',return_value=self.reader), \
             patch('perturb_predict.op3.expected_public_queries',return_value=self.ids.tolist()):
            return score_op3('fixture.h5ad','fixture.npz',self.root/'out',accept_public_reuse=True)

    def test_expected_queries_use_all_metadata_intersections(self):
        r = Reader.__new__(Reader)
        r.obs = {'perturbagen':np.array(['q0','q1','q0','q1','q2','q0','q1','q2']),
                 'cell_type':np.array([SOURCE]*2+[NK]*3+[self.ct]*3),
                 'split':np.array(['train']*5+['public_test']*3),
                 'is_control':np.zeros(8,dtype=bool), 'well':np.array(['A1']*8)}
        self.assertEqual(expected_public_queries(r,self.ct),['q0','q1'])

    def test_op3_subset_rejected_before_truth(self):
        a = dict(self.arrays,ids=self.ids[:1],prediction=np.ones((1,2)),source_copy=np.zeros((1,2)))
        with self.assertRaisesRegex(ValueError,'complete frozen query set'):self.invoke(arrays=a)
        self.reader.effects.assert_not_called()
        self.assertTrue((self.root/'out/FAILED.json').exists())
        self.assertFalse((self.root/'out/COMPLETE.json').exists())
        self.assertFalse((self.root/'out/access_started.json').exists())

    def test_op3_extra_query_rejected_before_truth(self):
        a = dict(self.arrays,ids=np.append(self.ids,'extra'),prediction=np.ones((4,2)),source_copy=np.zeros((4,2)))
        with self.assertRaisesRegex(ValueError,'extra=1'):self.invoke(arrays=a)
        self.reader.effects.assert_not_called()

    def test_existing_op3_run_rejected_before_reader(self):
        (self.root/'out').mkdir()
        with patch('perturb_predict.cli.read_predictions',return_value=(self.meta,self.arrays)), \
             patch('perturb_predict.op3.Reader') as reader:
            with self.assertRaises(FileExistsError):
                score_op3('missing','missing',self.root/'out',accept_public_reuse=True)
            reader.assert_not_called()

    def test_failed_op3_read_retains_started_and_final_audit(self):
        def fail(*a,**k):
            self.reader.log.append({'ids':[0], 'purpose':'score', 'status':'failed'})
            raise ValueError('simulated_decoder_failure')
        self.reader.effects.side_effect=fail
        with self.assertRaisesRegex(ValueError,'decoder_failure'):self.invoke()
        out=self.root/'out'
        self.assertTrue((out/'access_started.json').exists())
        self.assertTrue((out/'FAILED.json').exists())
        self.assertFalse((out/'COMPLETE.json').exists())
        self.assertEqual(json.loads((out/'access_log.json').read_text())['reader_log'][0]['status'],'failed')

    def test_existing_generic_run_rejected_before_truth(self):
        (self.root/'out').mkdir()
        with patch('perturb_predict.cli.read_predictions',return_value=(self.meta,self.arrays)), \
             patch('perturb_predict.cli.load_bundle') as truth:
            with self.assertRaises(FileExistsError):evaluate('unused','unused',self.root/'out')
            truth.assert_not_called()

    def test_reader_records_failed_decode_attempt(self):
        r=Reader.__new__(Reader);r.path='fixture';r.shape=(1,2);r.log=[]
        r.obs={'split':np.array(['public_test']),'is_control':np.array([False]),'well':np.array(['A1'])}
        with patch('perturb_predict.op3.h5py.File',side_effect=OSError('decode_failure')):
            with self.assertRaises(OSError):r.logcpm(np.array([0]),'score')
        self.assertEqual(r.log[0]['status'],'failed')
        self.assertEqual(r.log[0]['ids'],[0])

    def test_complete_op3_query_scoring_succeeds(self):
        scores=self.invoke()
        self.assertEqual(scores['model']['n_compounds'],3)
        self.assertEqual(scores['model']['mse'],0.)
        self.assertTrue((self.root/'out/COMPLETE.json').exists())
        self.assertTrue((self.root/'out/access_started.json').exists())
        self.assertTrue((self.root/'out/access_log.json').exists())

if __name__=='__main__':unittest.main()
