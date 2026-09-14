"""Exact link-audit design under a CONDITIONAL one-to-one matching model.

This module never infers an actual pairing. Each enumerated permutation is a
possible world used to quantify design risk, not a biological record.
"""
from __future__ import annotations
from itertools import permutations
import numpy as np


def enumerate_worlds(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError('Equal one-dimensional margins required; identity remains unverified')
    if not 2 <= len(x) <= 6 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Finite margins of length 2..6 required')
    worlds = np.asarray(list(permutations(range(len(y)))), dtype=int)
    xc, yc = x-x.mean(), y-y.mean()
    values = np.mean(xc[None, :] * yc[worlds], axis=1)
    return worlds, values


def conditional_interval(worlds, values, links=()):
    worlds, values = np.asarray(worlds), np.asarray(values)
    keep = np.ones(len(worlds), bool)
    seen_i, seen_j = set(), set()
    for i, j in links:
        if i in seen_i or j in seen_j:
            raise ValueError('Duplicate or contradictory authenticated links')
        if not (0 <= i < worlds.shape[1] and 0 <= j < worlds.shape[1]):
            raise ValueError('Link index out of bounds')
        keep &= worlds[:, i] == j
        seen_i.add(i); seen_j.add(j)
    if not keep.any():
        raise ValueError('No compatible one-to-one correspondence')
    return float(values[keep].min()), float(values[keep].max()), int(keep.sum())


def one_link_design(x, y, native_ids=None):
    worlds, values = enumerate_worlds(x, y)
    n = len(x)
    ids = list(map(str, native_ids)) if native_ids is not None else list(map(str, range(n)))
    if len(ids) != n or len(set(ids)) != n:
        raise ValueError('Unique native IDs required')
    initial = float(np.ptp(values))
    widths = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            low, high, _ = conditional_interval(worlds, values, [(i,j)])
            widths[i,j] = max(0., high-low)
    worst = widths.max(1)
    chosen = min(range(n), key=lambda i:(worst[i], ids[i]))
    # For each possible true world: a deterministic query versus a randomized
    # record query. The latter averages query choice, not presumed world truth.
    achieved = widths[chosen, worlds[:,chosen]]
    random_record = np.mean(np.stack([widths[i,worlds[:,i]] for i in range(n)]),axis=0)
    return dict(chosen_index=chosen, chosen_native_id=ids[chosen],
                initial_low=float(values.min()), initial_high=float(values.max()),
                initial_width=initial, worst_width=float(worst[chosen]),
                guaranteed_reduction=initial-float(worst[chosen]),
                random_record_adversarial_mean_width=float(worst.mean()),
                candidate_worst_widths=worst.tolist(),
                widths_by_possible_partner=widths.tolist(),
                possible_world_count=len(worlds),
                achieved_width_each_world=achieved.tolist(),
                random_record_width_each_world=random_record.tolist(),
                never_an_inferred_pairing=True)


def greedy_one_per_condition(gains, budget, tie_ids):
    gains=np.asarray(gains,float)
    if gains.ndim != 1 or not np.isfinite(gains).all() or (gains < -1e-12).any():
        raise ValueError('Finite nonnegative gains required')
    if len(tie_ids)!=len(gains) or len(set(tie_ids))!=len(gains):
        raise ValueError('Unique condition IDs required')
    if not isinstance(budget,int) or budget<0 or budget>len(gains):
        raise ValueError('Budget out of range')
    return sorted(range(len(gains)),key=lambda i:(-gains[i],tie_ids[i]))[:budget]
