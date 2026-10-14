"""Deterministic offline walkthrough of the stable software pipeline.

These invented arrays check installation, file contracts and the separation of
training, prediction and evaluation. Their scores are not biological evidence.
"""
from __future__ import annotations

import numpy as np

from . import cli
from .io import digest, json_write, provenance, run_directory, write_bundle


DEMO_SEED = 20261014
EVIDENCE = 'synthetic_software_demo_only'
DISCLAIMER = ('All compounds, genes and effects are invented. Scores check the '
              'software workflow only and do not measure biological performance.')


def run_demo(out):
    """Run a fixed synthetic example in a new directory, without network access.

    Input arrays and resulting scores are reproducible for the recorded NumPy
    version. Runtime provenance includes timestamps, so complete run directories
    are intentionally not byte-identical. Existing output directories are refused.
    """
    with run_directory(out) as root:
        calibration_ids = np.array([f'synthetic_calibration_{i:02d}' for i in range(1, 7)])
        query_ids = np.array([f'synthetic_query_{i:02d}' for i in range(1, 4)])
        genes = np.array([f'synthetic_gene_{i:02d}' for i in range(1, 9)])
        contract = {
            'effect_space': 'synthetic_unitless_effect',
            'source_context': 'synthetic_source',
            'target_context': 'synthetic_target',
            'organism': 'synthetic_no_organism',
            'dose': 'synthetic_not_applicable',
            'time': 'synthetic_not_applicable',
            'gene_namespace': 'synthetic_gene',
        }
        budget = {
            'target_calibration_compounds': len(calibration_ids),
            'target_query_compounds_available_to_fit': 0,
            'target_query_compounds_available_to_predict': 0,
            'target_query_compounds_for_evaluation': len(query_ids),
            'biological_samples': 0,
        }
        recipe = {
            'evidence': EVIDENCE,
            'disclaimer': DISCLAIMER,
            'seed': DEMO_SEED,
            'random_generator': 'numpy.random.Generator(PCG64)',
            'generator_version': 1,
            'numpy': np.__version__,
            'n_calibration_compounds': len(calibration_ids),
            'n_query_compounds': len(query_ids),
            'n_genes': len(genes),
            'source_distribution': 'standard_normal; calibration rows then query rows',
            'target_formula': 'source * linspace(0.6, 1.2, 8) + '
                              'linspace(-0.15, 0.15, 8) + 0.03 * standard_normal',
            'split_ids': {'calibration': calibration_ids.tolist(), 'query': query_ids.tolist()},
            'target_information_budget': budget,
        }
        # Persist the fixed software-check recipe before generating or scoring data.
        json_write(root / 'protocol.json', {
            **recipe, 'contract': contract,
            'algorithm': 'single_source_partial_pooling_v1',
            'evaluation_is_new_independent_evidence': False,
            **provenance(),
        })
        rng = np.random.Generator(np.random.PCG64(DEMO_SEED))
        source = rng.normal(size=(len(calibration_ids) + len(query_ids), len(genes)))
        target = (source * np.linspace(0.6, 1.2, len(genes))
                  + np.linspace(-0.15, 0.15, len(genes))
                  + 0.03 * rng.normal(size=source.shape))
        n = len(calibration_ids)
        with run_directory(root / 'inputs') as inputs:
            write_bundle(inputs / 'train.npz', 'train', calibration_ids, genes, contract,
                         source=source[:n], target=target[:n], provenance_info=recipe)
            # Query labels live exclusively in a separate truth bundle. The train
            # and predict entrypoints receive no truth path or query target array.
            write_bundle(inputs / 'query.npz', 'query', query_ids, genes, contract,
                         source=source[n:], provenance_info=recipe)
            write_bundle(inputs / 'truth.npz', 'truth', query_ids, genes, contract,
                         target=target[n:], provenance_info=recipe)

        fit = cli.train(inputs / 'train.npz', root / 'model')
        prediction = cli.predict(root / 'model/model.npz', inputs / 'query.npz',
                                 root / 'prediction')
        scores = cli.evaluate(root / 'prediction/predictions.npz', inputs / 'truth.npz',
                              root / 'evaluation')
        # Hash every artifact already produced, including stage completion records.
        # summary.json and the outer COMPLETE.json are written after this manifest.
        artifact_hashes = {
            path.relative_to(root).as_posix(): digest(path)
            for path in sorted(root.rglob('*')) if path.is_file()
        }
        summary = {
            'evidence': EVIDENCE,
            'disclaimer': DISCLAIMER,
            'seed': DEMO_SEED,
            'contract': contract,
            'n_calibration_compounds': len(calibration_ids),
            **prediction,
            'target_information_budget': budget,
            'selected_parameters': fit['selected'],
            'scores': scores,
            'evaluation_is_new_independent_evidence': False,
            'all_predictions_frozen_before_scoring': True,
            'input_sha256': {name: digest(inputs / f'{name}.npz')
                             for name in ('train', 'query', 'truth')},
            'artifact_sha256': artifact_hashes,
            **provenance(),
        }
        json_write(root / 'summary.json', summary)
    return summary
