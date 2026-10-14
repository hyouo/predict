"""The offline example is a software check, never a biological benchmark."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from perturb_predict import cli
from perturb_predict.demo import DEMO_SEED, EVIDENCE, run_demo
from perturb_predict.io import digest, load_archive, load_bundle


class DemoTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)

    def test_complete_offline_pipeline_preserves_provenance_and_hashes(self):
        out = self.root / 'demo'
        summary = run_demo(out)
        self.assertEqual(summary, json.loads((out / 'summary.json').read_text()))
        self.assertEqual(summary['evidence'], EVIDENCE)
        self.assertIn('do not measure biological performance', summary['disclaimer'])
        self.assertFalse(summary['evaluation_is_new_independent_evidence'])
        self.assertFalse(summary['query_target_labels_read'])
        self.assertTrue(summary['all_predictions_frozen_before_scoring'])
        self.assertEqual(summary['n_calibration_compounds'], 6)
        self.assertEqual(summary['n_predictions'], 3)
        self.assertEqual(summary['n_genes'], 8)
        self.assertEqual(set(summary['scores']), {'model', 'source_copy', 'zero'})
        metrics = json.loads((out / 'evaluation/metrics.json').read_text())
        self.assertEqual(summary['scores'], metrics['scores'])
        for stage in ('', 'inputs', 'model', 'prediction', 'evaluation'):
            completion = json.loads((out / stage / 'COMPLETE.json').read_text())
            self.assertEqual(completion['status'], 'complete')
            self.assertIn('demo.py', completion['code_sha256'])
            self.assertFalse((out / stage / 'FAILED.json').exists())
        for path in ('protocol.json', 'model/fit.json', 'prediction/freeze.json',
                     'evaluation/access_started.json', 'evaluation/metrics.json'):
            record = json.loads((out / path).read_text())
            self.assertIn('package_version', record)
            self.assertIn('python', record)
            self.assertIn('numpy', record)
            self.assertIn('utc', record)
            self.assertIn('demo.py', record['code_sha256'])
        self.assertTrue((out / 'prediction/predictions.csv').is_file())
        self.assertTrue((out / 'evaluation/per_compound.csv').is_file())
        artifact_paths = {p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file()}
        self.assertEqual(set(summary['artifact_sha256']),
                         artifact_paths - {'summary.json', 'COMPLETE.json'})
        for path, expected in summary['artifact_sha256'].items():
            self.assertEqual(digest(out / path), expected, path)
        for kind, expected in summary['input_sha256'].items():
            self.assertEqual(digest(out / f'inputs/{kind}.npz'), expected)

    def test_existing_run_is_refused_without_touching_outputs(self):
        out = self.root / 'demo'
        run_demo(out)
        before = {p.relative_to(out): digest(p) for p in out.rglob('*') if p.is_file()}
        with patch('perturb_predict.demo.write_bundle') as writer, \
             patch('perturb_predict.cli.train') as train:
            with self.assertRaises(FileExistsError):
                run_demo(out)
            writer.assert_not_called()
            train.assert_not_called()
        after = {p.relative_to(out): digest(p) for p in out.rglob('*') if p.is_file()}
        self.assertEqual(before, after)

    def test_query_labels_are_isolated_until_frozen_prediction_is_evaluated(self):
        out = self.root / 'demo'
        reads = []

        def read(path, kind):
            reads.append(kind)
            if kind == 'truth':
                self.assertEqual(reads, ['train', 'query', 'truth'])
                freeze = json.loads((out / 'prediction/freeze.json').read_text())
                self.assertFalse(freeze['query_target_labels_read'])
                self.assertEqual(freeze['prediction_sha256'],
                                 digest(out / 'prediction/predictions.npz'))
                self.assertTrue((out / 'prediction/COMPLETE.json').is_file())
                self.assertTrue((out / 'evaluation/access_started.json').is_file())
            return load_bundle(path, kind)

        with patch('perturb_predict.cli.load_bundle', side_effect=read):
            run_demo(out)
        self.assertEqual(reads, ['train', 'query', 'truth'])
        _, train = load_bundle(out / 'inputs/train.npz', 'train')
        query_meta, query = load_bundle(out / 'inputs/query.npz', 'query')
        _, truth = load_bundle(out / 'inputs/truth.npz', 'truth')
        self.assertEqual(set(query), {'source', 'ids', 'genes'})
        self.assertEqual(set(truth), {'target', 'ids', 'genes'})
        self.assertFalse(set(train['ids']) & set(query['ids']))
        np.testing.assert_array_equal(query['ids'], truth['ids'])
        self.assertEqual(query_meta['provenance']['evidence'], EVIDENCE)
        self.assertEqual(query_meta['provenance']['target_information_budget']
                         ['target_query_compounds_available_to_fit'], 0)
        model_meta, _ = load_archive(out / 'model/model.npz')
        self.assertFalse(model_meta['fit_report']['private_or_query_target_labels_used_by_fit'])

    def test_fixed_seed_reproduces_inputs_predictions_and_scores(self):
        first, second = self.root / 'first', self.root / 'second'
        a, b = run_demo(first), run_demo(second)
        self.assertEqual(a['seed'], DEMO_SEED)
        self.assertEqual(b['seed'], DEMO_SEED)
        self.assertEqual(a['input_sha256'], b['input_sha256'])
        self.assertEqual(a['selected_parameters'], b['selected_parameters'])
        self.assertEqual(a['scores'], b['scores'])
        for path in ('model/model.npz', 'prediction/predictions.npz'):
            _, left = load_archive(first / path)
            _, right = load_archive(second / path)
            for key in left:
                np.testing.assert_array_equal(left[key], right[key], err_msg=f'{path}: {key}')

    def test_evaluation_failure_keeps_frozen_prediction_and_failure_marker(self):
        out = self.root / 'failed'
        with patch.object(cli, 'evaluate', side_effect=OSError('demo evaluation failure')):
            with self.assertRaisesRegex(OSError, 'demo evaluation failure'):
                run_demo(out)
        self.assertFalse((out / 'COMPLETE.json').exists())
        self.assertFalse((out / 'summary.json').exists())
        failure = json.loads((out / 'FAILED.json').read_text())
        self.assertEqual(failure['status'], 'failed')
        self.assertIn('demo evaluation failure', failure['error'])
        self.assertTrue((out / 'prediction/COMPLETE.json').is_file())
        cli.read_predictions(out / 'prediction/predictions.npz')


if __name__ == '__main__':
    unittest.main()
