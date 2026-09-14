"""Descriptive metrics. No biological independence or calibrated uncertainty claims."""
from __future__ import annotations
import numpy as np
from .io import matrix


def metrics(pred, truth):
    p, y = matrix(pred, 'prediction'), matrix(truth, 'truth')
    if p.shape != y.shape: raise ValueError('Prediction/truth shape mismatch')
    e = p - y
    bias = e.mean(0)
    pc, yc = p - p.mean(0), y - y.mean(0)
    norm = float(np.linalg.norm(pc) * np.linalg.norm(yc))
    ps = pc / np.maximum(np.linalg.norm(pc, axis=1, keepdims=True), 1e-12)
    ys = yc / np.maximum(np.linalg.norm(yc, axis=1, keepdims=True), 1e-12)
    similarity = ps @ ys.T
    top, reciprocal_rank = [], []
    for i, row in enumerate(similarity):
        same = abs(row - row[i]) <= 1e-12
        rank = 1 + np.sum(row > row[i] + 1e-12) + (same.sum() - 1) / 2
        top.append(float(abs(row[i] - row.max()) <= 1e-12) / max(int(np.sum(abs(row-row.max()) <= 1e-12)), 1))
        reciprocal_rank.append(1 / rank)
    return {'mse': float(np.mean(e**2)), 'mean_row_rmse': float(np.sqrt(np.mean(e**2, axis=1)).mean()),
            'bias_mse': float(np.mean(bias**2)), 'centered_mse': float(np.mean((e-bias)**2)),
            'centered_pattern_cosine': float((pc*yc).sum()/norm) if norm > 1e-15 else 0.,
            'retrieval_top1': float(np.mean(top)), 'retrieval_mrr': float(np.mean(reciprocal_rank)),
            'n_compounds': len(p), 'n_genes': p.shape[1]}
