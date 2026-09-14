"""Installable, separable prepare/train/predict/evaluate workflow."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from . import __version__
from .model import TransferModel
from .io import (load_bundle, load_archive, save_archive, check_contract, identifiers, matrix,
                 alignment, run_directory, json_write, provenance, digest, SCHEMA)
from .evaluation import metrics


def train(data, out):
    meta, a = load_bundle(data, 'train')
    with run_directory(out) as root:
        model = TransferModel().fit(a['source'], a['target'], ids=a['ids'], genes=a['genes'],
                                   contract=meta['contract'])
        model.save(root / 'model.npz')
        json_write(root / 'fit.json', {**model.report, 'training_input_sha256': digest(data),
                      'model_sha256': digest(root / 'model.npz'), 'contract': model.contract, **provenance()})
    return model.report


def predict(model_path, data, out):
    model = TransferModel.load(model_path)
    meta, a = load_bundle(data, 'query')
    with run_directory(out) as root:
        p, info = model.predict(a['source'], genes=a['genes'], ids=a['ids'], contract=meta['contract'])
        ix = alignment(a['genes'], model.genes, 'genes')
        source = a['source'][:, ix]
        metadata = {'schema_version': SCHEMA, 'kind': 'prediction', 'contract': model.contract,
                    'model_sha256': digest(model_path), 'query_sha256': digest(data),
                    'model_calibration_ids': model.ids.tolist(), 'query_target_labels_read': False,
                    'package_version': __version__, 'input_provenance': meta.get('provenance', {})}
        save_archive(root / 'predictions.npz', metadata, prediction=p, source_copy=source,
                     ids=a['ids'], genes=model.genes)
        with (root / 'predictions.csv').open('x', encoding='utf-8', newline='') as f:
            w = csv.writer(f); w.writerow(['compound_id', *model.genes])
            for name, row in zip(a['ids'], p): w.writerow([name, *row])
        diagnostics = [{'compound_id': str(d), 'relative_source_distance': float(v),
                         'warning': 'far_from_calibration' if v > 3 else ''}
                        for d, v in zip(a['ids'], info['distance_to_calibration_center'])]
        # Freeze manifest produced before a separate evaluator is allowed to read truth.
        json_write(root / 'freeze.json', {**metadata, 'prediction_sha256': digest(root / 'predictions.npz'),
                  'diagnostics': diagnostics, 'confidence_intervals_available': False,
                  'diagnostic_threshold_is_heuristic': True, **provenance()})
    return {'n_predictions': len(p), 'n_genes': p.shape[1], 'query_target_labels_read': False}


def read_predictions(path):
    p = Path(path)
    completion = p.parent / 'COMPLETE.json'
    if not completion.exists() or json.loads(completion.read_text()).get('status') != 'complete':
        raise ValueError('Prediction run is incomplete or failed')
    freeze_path = p.parent / 'freeze.json'
    freeze = json.loads(freeze_path.read_text(encoding='utf-8'))
    if freeze.get('prediction_sha256') != digest(p): raise ValueError('Frozen prediction SHA256 mismatch')
    meta, a = load_archive(p)
    if meta.get('kind') != 'prediction' or set(a) != {'prediction', 'source_copy', 'ids', 'genes'}:
        raise ValueError('Malformed prediction archive')
    check_contract(meta)
    a['ids'], a['genes'] = identifiers(a['ids'], 'ids'), identifiers(a['genes'], 'genes')
    for k in ('prediction', 'source_copy'):
        a[k] = matrix(a[k], k)
        if a[k].shape != (len(a['ids']), len(a['genes'])): raise ValueError('Prediction shape mismatch')
    if set(a['ids']) & set(meta.get('model_calibration_ids', [])):
        raise ValueError('Calibration/query overlap')
    if freeze.get('contract') != meta['contract']:
        raise ValueError('Freeze/array metadata contract mismatch')
    return meta, a


def score_arrays(out, pred, truth, ids, *, record):
    with run_directory(out) as root:
        scores = {name: metrics(p, truth) for name, p in pred.items()}
        json_write(root / 'metrics.json', {'scores': scores, **record, **provenance()})
        with (root / 'per_compound.csv').open('x', encoding='utf-8', newline='') as f:
            w = csv.writer(f); w.writerow(['model', 'compound_id', 'mse', 'rmse'])
            for name, p in pred.items():
                for d, mse in zip(ids, np.mean((p-truth)**2, axis=1)):
                    w.writerow([name, d, float(mse), float(np.sqrt(mse))])
    return scores


def evaluate(prediction, truth_path, out):
    pm, p = read_predictions(prediction)
    tm, t = load_bundle(truth_path, 'truth')
    if pm['contract'] != tm['contract']: raise ValueError('Truth and prediction contract mismatch')
    rows, genes = alignment(t['ids'], p['ids'], 'ids'), alignment(t['genes'], p['genes'], 'genes')
    truth = t['target'][rows][:, genes]
    return score_arrays(out, {'model': p['prediction'], 'source_copy': p['source_copy'],
                    'zero': np.zeros_like(truth)}, truth, p['ids'],
                    record={'prediction_sha256': digest(prediction), 'truth_sha256': digest(truth_path),
                            'evidence': 'user_supplied_truth; independence_not_certified'})


def score_op3(data, prediction, out, *, accept_public_reuse=False):
    if not accept_public_reuse:
        raise ValueError('Pass --accept-public-reuse: public labels were previously used in development')
    from .op3 import Reader, TARGETS, DATA_SHA256, contract
    meta, p = read_predictions(prediction)  # check freeze before any public expression access
    ct = meta['contract']['target_context']
    if ct not in TARGETS.values() or meta['contract'] != contract(ct):
        raise ValueError('Prediction is not for the supported OP3 contract')
    if meta.get('input_provenance', {}).get('dataset_sha256') != DATA_SHA256:
        raise ValueError('Prediction lacks pinned OP3 query provenance')
    r = Reader(data)
    if not np.array_equal(r.genes, p['genes']): raise ValueError('Pinned OP3 gene panel/order mismatch')
    allowed = set(r.obs['perturbagen'][r.ids(ct, split='public_test')])
    if not set(p['ids']) <= allowed: raise ValueError('Predictions contain non-public query compounds')
    effects, _ = r.effects(ct, 'public_test', 'EFGH', purpose='score')
    truth = np.stack([effects[d] for d in p['ids']])
    return score_arrays(out, {'model': p['prediction'], 'source_copy': p['source_copy'],
                            'zero': np.zeros_like(truth)}, truth, p['ids'],
                  record={'prediction_sha256': digest(prediction), 'dataset_sha256': DATA_SHA256,
                          'evidence': 'engineering_acceptance_on_reused_public_OP3; NOT_new_validation',
                          'private_expression_read': False, 'reader_log': r.log})


def benchmark(data, out, *, accept_conditional=False, accept_public_reuse=False):
    if not (accept_conditional and accept_public_reuse):
        raise ValueError('Both --accept-conditional and --accept-public-reuse are required')
    from .op3 import prepare, TARGETS
    with run_directory(out) as root:
        audit = prepare(data, root / 'prepared', accept_conditional=True)
        for name in TARGETS:
            train(root/'prepared'/name/'train.npz', root/name/'model')
            predict(root/name/'model/model.npz', root/'prepared'/name/'query.npz', root/name/'prediction')
        # Both target predictions are frozen before either target public label is read.
        scores = {name: score_op3(data, root/name/'prediction/predictions.npz', root/name/'evaluation',
                                 accept_public_reuse=True) for name in TARGETS}
        macro = {m: {k: float(np.mean([scores[n][m][k] for n in TARGETS]))
                     for k in ('mse', 'centered_mse', 'bias_mse', 'retrieval_top1')}
                 for m in ('model', 'source_copy', 'zero')}
        json_write(root/'benchmark.json', {'by_context': scores, 'macro': macro,
                  'targets': audit['targets'], 'evaluation_is_new_independent_evidence': False,
                  'all_predictions_frozen_before_scoring': True, **provenance()})
    return macro


def main(argv=None):
    ap = argparse.ArgumentParser(description='Predict v0.1.0: auditable few-shot RNA effect baseline, not SOTA.')
    ap.add_argument('--version', action='version', version=f'perturb-predict {__version__}')
    sub = ap.add_subparsers(dest='command', required=True)
    fetch_p = sub.add_parser('fetch-op3', help='Download and verify fixed public OP3 data, <=32 MiB')
    fetch_p.add_argument('--out', type=Path, required=True)
    for command in ('prepare-op3', 'benchmark-op3'):
        p = sub.add_parser(command)
        p.add_argument('--data', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
        p.add_argument('--accept-conditional', action='store_true')
        if command == 'benchmark-op3': p.add_argument('--accept-public-reuse', action='store_true')
    p = sub.add_parser('train'); p.add_argument('--data', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('predict'); p.add_argument('--model', type=Path, required=True)
    p.add_argument('--data', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('evaluate'); p.add_argument('--prediction', type=Path, required=True)
    p.add_argument('--truth', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('score-op3'); p.add_argument('--data', type=Path, required=True)
    p.add_argument('--prediction', type=Path, required=True); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--accept-public-reuse', action='store_true')
    a = ap.parse_args(argv)
    try:
        if a.command == 'fetch-op3':
            from .op3 import fetch
            result = fetch(a.out)
        elif a.command == 'prepare-op3':
            from .op3 import prepare
            result = prepare(a.data, a.out, accept_conditional=a.accept_conditional)
            result = {'targets': result['targets'], 'private_expression_read': False}
        elif a.command == 'benchmark-op3':
            result = benchmark(a.data, a.out, accept_conditional=a.accept_conditional,
                               accept_public_reuse=a.accept_public_reuse)
        elif a.command == 'train': result = train(a.data, a.out)
        elif a.command == 'predict': result = predict(a.model, a.data, a.out)
        elif a.command == 'evaluate': result = evaluate(a.prediction, a.truth, a.out)
        else: result = score_op3(a.data, a.prediction, a.out, accept_public_reuse=a.accept_public_reuse)
        print(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2))
        return 0
    except (OSError, ValueError, KeyError, PermissionError) as exc:
        print(json.dumps({'status': 'error', 'type': type(exc).__name__, 'message': str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
