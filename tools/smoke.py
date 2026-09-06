"""One measured-data integration check, not a new benchmark or scientific validation."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
import pandas as pd
from src.data import load_panel
from src.covariance import build_kernels
from src.aggregation import bank, aggregate
from src.streaming import fit_streaming
from src.kriging import acquire
from run_covariance import metrics
from tools.assets import sha256_file

ROOT = Path(__file__).resolve().parents[1]


def run(output: Path) -> dict:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    panel = load_panel(ROOT / 'data/hepatocyte_signaling_raw.csv')
    source = np.delete(panel.effects, 0, axis=0)
    kernels = build_kernels(source, panel.cues)
    perm = np.random.default_rng(20260906).permutation(72)
    queries = perm[:24]
    anchors, _ = acquire(kernels['cue_original'], np.sort(perm[24:]), queries, 16, criterion='proper')
    if not set(queries).isdisjoint(anchors):
        raise AssertionError('Anchor/query overlap')
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                         text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    record = {'purpose': 'integration_check_only', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'git_commit': commit, 'git_dirty': dirty,
              'code_sha256': {str(p.relative_to(ROOT)): sha256_file(p)
                              for folder in ['src', 'tools'] for p in sorted((ROOT / folder).glob('*.py'))},
              'environment': {'python': sys.version, 'platform': platform.platform(),
                              'numpy': np.__version__, 'scipy': scipy.__version__, 'pandas': pd.__version__,
                              'threads': {k: os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']}},
              'data_sha256': panel.audit['sha256'], 'context': panel.contexts[0], 'seed': 20260906,
              'anchors': anchors.tolist(), 'queries': queries.tolist(),
              'source_query_outcomes_allowed': True, 'target_query_outcomes_passed_to_fit': False,
              'rna_training_performed': False, 'new_independent_biological_studies': 0, 'completed': False}
    path = output / 'run.json'
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    try:
        pp, errors, params = bank(source, panel.cues, kernels, anchors, panel.effects[0, anchors].copy())
        rows = []
        for mode in ['bootstrap', 'stacking']:
            full = aggregate(pp, errors, params, mode)
            blocked = fit_streaming(source, panel.cues, anchors, panel.effects[0, anchors].copy(),
                                    block_size=5, mode=mode, kernels=kernels)
            np.testing.assert_allclose(full.prediction, blocked.prediction, rtol=1e-6, atol=1e-6)
            prediction = full.prediction[queries]
            np.save(output / f'{mode}.npy', prediction, allow_pickle=False)
            rows.append({'model': mode, **metrics(panel.effects[0, queries], prediction),
                         'prediction_sha256': hashlib.sha256(prediction.astype('<f8').tobytes()).hexdigest(),
                         'streaming_max_abs_difference': float(np.max(np.abs(full.prediction - blocked.prediction)))})
        record.update({'completed': True, 'scores': rows})
    except Exception as exc:
        record['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New directory; existing results are never overwritten')
    print(json.dumps(run(parser.parse_args().output), indent=2, allow_nan=False))
