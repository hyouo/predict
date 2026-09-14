"""Reference-unit-aware risk refinement; existing frozen RNA models are unchanged.

Independent, centered reference errors are a WORKING covariance assumption.
A donor and an assay plate are not interchangeable biological replication units.
"""
from __future__ import annotations
import numpy as np
from src.rna_models import (candidate_bank, shared_kernels, kernel_blocks,
                           maps_for, prior_block, fit_weights, _finite_shape)


def reference_correction(loo: np.ndarray, factors: np.ndarray) -> float:
    """Unscaled sum_u,g F[u,:,g]^T H[g] F[u,:,g], without cross-unit terms."""
    h = np.asarray(loo, dtype=float)
    f = np.asarray(factors, dtype=float)
    if h.ndim != 3 or f.ndim != 3 or h.shape != (f.shape[2], f.shape[1], f.shape[1]):
        raise ValueError('Reference factors and smoothers have incompatible shapes')
    if not np.isfinite(h).all() or not np.isfinite(f).all():
        raise ValueError('Nonfinite reference factors or smoother')
    return float(np.einsum('uig,gij,ujg->', f, h, f, optimize=True))


def fit_reference_ensembles(source, anchors, reference_c, reference_d, factors,
                            gram, *, block_size: int = 128):
    """Change only covariance correction; all 71 candidates and labels are fixed.

    factors: reference_unit x anchor_compound x gene, whose sum equals C-D
    anchor effects. `gram` is the ordinary LOO error Gram from the frozen bank.
    """
    source, yc, yd, ix = _finite_shape(source, reference_c, reference_d, anchors)
    f = np.asarray(factors, dtype=float)
    if f.ndim != 3 or f.shape[1:] != yc.shape or not len(f) or not np.isfinite(f).all():
        raise ValueError('Invalid reference factors')
    if not np.allclose(f.sum(axis=0), yc-yd, rtol=1e-8, atol=1e-9):
        raise ValueError('Factors do not reconstruct the paired reference difference')
    if isinstance(block_size, bool) or not isinstance(block_size, int) or block_size < 1:
        raise ValueError('Invalid block size')
    bank = candidate_bank(); m = len(bank); n, ng = yc.shape; nc = source.shape[1]
    g = np.asarray(gram, dtype=float)
    if g.shape != (m, m) or not np.isfinite(g).all():
        raise ValueError('Ordinary Gram has the wrong candidate order or shape')
    shared = shared_kernels(source); mean = source.mean(axis=0); y = (yc+yd)/2
    correction = np.zeros(m)
    for start in range(0, ng, block_size):
        stop = min(start+block_size, ng); b = stop-start
        kernels = kernel_blocks(source, start, stop, shared)
        for j, can in enumerate(bank):
            loo, _ = maps_for(can, kernels, ix, b, nc)
            correction[j] += reference_correction(loo, f[:, :, start:stop])
    correction /= n*ng
    choices = {}; diagnostics = {}
    for label, mask in (
        ('simple_stack', np.array([not c.kernel.startswith('mixed') for c in bank])),
        ('mixed_stack', np.ones(m, dtype=bool)),
    ):
        ids = np.flatnonzero(mask)
        w, rec = fit_weights(g[np.ix_(ids, ids)], correction[ids], adjusted=True)
        full = np.zeros(m); full[ids] = w
        choices[label] = full; diagnostics[label] = rec
    predictions = {k: np.zeros((nc, ng)) for k in choices}
    for start in range(0, ng, block_size):
        stop = min(start+block_size, ng); b = stop-start
        kernels = kernel_blocks(source, start, stop, shared)
        for j, can in enumerate(bank):
            needed = [k for k, w in choices.items() if w[j] > 1e-12]
            if not needed:
                continue
            _, final = maps_for(can, kernels, ix, b, nc)
            prior = prior_block(can, mean, start, stop)
            pred = prior + np.einsum('gij,jg->ig', final, y[:, start:stop]-prior[ix])
            for k in needed:
                predictions[k][:, start:stop] += choices[k][j]*pred
    return predictions, {
        'reference_unit_count': len(f), 'factor_sum_verified': True,
        'correction_before_half_factor': correction.tolist(),
        'weights': {k: w.tolist() for k, w in choices.items()},
        'optimization': diagnostics, 'target_test_outcomes_received': False,
        'scientific_status': 'exploratory_followup_after_original_test_results',
        'assumption': 'Independent centered reference units conditional on fixed source inputs; no selected-model unbiasedness guarantee.',
    }
