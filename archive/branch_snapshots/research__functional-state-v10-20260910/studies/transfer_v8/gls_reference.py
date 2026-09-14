"""Small independent reference for profiled, ridge-regularized GLS.

Each block is (design, observed_source_response, working_kernel). This mathematical
core is not the complete data pipeline and does not establish biological noise
calibration, novelty, or superiority. No target response is needed for prediction.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


def profiled_precision(kernel, residual_penalty):
    k = np.asarray(kernel, dtype=float)
    if k.ndim != 2 or k.shape[0] != k.shape[1] or len(k) == 0:
        raise ValueError('A nonempty square kernel is required')
    if not np.isfinite(k).all() or not np.isfinite(residual_penalty) or residual_penalty <= 0:
        raise ValueError('Finite kernel and positive residual penalty required')
    c = (k + k.T) / 2 + residual_penalty * np.eye(len(k))
    inverse = cho_solve(cho_factor(c, lower=True), np.eye(len(k)))
    a = inverse.sum(axis=1)
    q = inverse - np.outer(a, a) / a.sum()
    return (q + q.T) / 2


def fit_profiled_gls(blocks, residual_penalty=1., prior_penalty=1.):
    if not np.isfinite(prior_penalty) or prior_penalty <= 0:
        raise ValueError('Positive prior penalty required')
    gram = rhs = None
    for design, response, kernel in blocks:
        x, y = np.asarray(design, float), np.asarray(response, float)
        if x.ndim != 2 or y.ndim != 2 or len(x) != len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError('Aligned finite design/response matrices required')
        q = profiled_precision(kernel, residual_penalty)
        if len(q) != len(x):
            raise ValueError('Kernel and response dimensions differ')
        if gram is None:
            gram = np.zeros((x.shape[1], x.shape[1]))
            rhs = np.zeros((x.shape[1], y.shape[1]))
        if gram.shape[0] != x.shape[1] or rhs.shape[1] != y.shape[1]:
            raise ValueError('Block feature/output dimensions differ')
        gram += x.T @ q @ x
        rhs += x.T @ q @ y
    if gram is None:
        raise ValueError('At least one block required')
    positive = np.diag(gram)
    positive = positive[positive > 1e-12]
    scale = float(positive.mean()) if len(positive) else 1.
    a = (gram + gram.T) / 2 + prior_penalty * scale * np.eye(len(gram))
    coefficient = cho_solve(cho_factor(a, lower=True), rhs)
    return coefficient, {'penalty_scale': scale, 'gram': gram, 'rhs': rhs}
