"""Single-source, partially pooled gene slopes with intercept shrinkage.

Classical statistical baseline extracted from v0.8, without spectral truncation.
Tuning uses complete calibration compounds only; no query labels enter this API.
"""
from __future__ import annotations
import numpy as np
from .io import (matrix, identifiers, save_archive, load_archive, check_contract,
                 alignment, SCHEMA, provenance)

PENALTIES = (0.1, 1., 10., 100., 1000.)
RHOS = (0., .25, .5, .75, 1.)


def coefficients(x, y, penalty):
    """Match the one-source v0.8 hierarchical model's pooled prior and scale."""
    x, y = matrix(x, 'source'), matrix(y, 'target')
    if x.shape != y.shape or len(x) < 2 or not np.isfinite(penalty) or penalty <= 0:
        raise ValueError('Matched source/target with >=2 compounds and positive penalty required')
    mx, my = x.mean(0), y.mean(0)
    xc, yc = x - mx, y - my
    xx, xy = (xc * xc).sum(0), (xc * yc).sum(0)
    total = float(xx.sum())
    prior = float(xy.sum()) / (total + max(total * 1e-8, 1e-10))
    ridge = penalty * max(float(np.median(xx)), 1e-8)
    slopes = (xy + ridge * prior) / (xx + ridge)
    intercept = my - mx * slopes
    if not np.isfinite(slopes).all() or not np.isfinite(intercept).all():
        raise ValueError('Ill-scaled input produced invalid coefficients')
    return slopes, intercept


class TransferModel:
    """Fit one target context from aligned, unique calibration-compound effects."""
    def fit(self, source, target, *, ids, genes, contract):
        x, y = matrix(source, 'source'), matrix(target, 'target')
        ids, genes = identifiers(ids, 'ids'), identifiers(genes, 'genes')
        check_contract({'contract': contract})
        if x.shape != y.shape or x.shape != (len(ids), len(genes)) or len(ids) < 4:
            raise ValueError('Need >=4 distinct calibration compounds, aligned source/target/genes')
        risks = []
        for penalty in PENALTIES:
            terms, intercepts = np.empty_like(x), np.empty_like(x)
            for i in range(len(x)):
                keep = np.arange(len(x)) != i
                b, a = coefficients(x[keep], y[keep], penalty)
                terms[i], intercepts[i] = x[i] * b, a
            for rho in RHOS:
                loss = float(np.mean((terms + rho * intercepts - y)**2))
                if not np.isfinite(loss): raise ValueError('Nonfinite validation loss')
                risks.append({'penalty': penalty, 'rho': rho, 'loo_mse': loss})
        best = min(v['loo_mse'] for v in risks)
        selected = sorted((v for v in risks if v['loo_mse'] <= best + 1e-12),
                          key=lambda v: (v['rho'], -v['penalty']))[0]
        b, a = coefficients(x, y, selected['penalty'])
        self.slopes = b
        self.intercept = selected['rho'] * a
        self.genes, self.ids, self.training_source = genes, ids, x.copy()
        self.contract = dict(contract)
        self.report = {'selected': selected, 'candidates': risks,
                       'n_calibration_compounds': len(x), 'n_genes': len(genes),
                       'validation': 'leave_complete_compound_out; same pooled reference',
                       'confidence_intervals_available': False,
                       'private_or_query_target_labels_used_by_fit': False}
        return self

    def predict(self, source, *, genes, ids, contract):
        if not hasattr(self, 'slopes'): raise ValueError('Model is not fitted')
        check_contract({'contract': contract})
        if contract != self.contract:
            raise ValueError('Effect-space/source/target/dose/time contract differs from model')
        q = matrix(source, 'source')
        ix = alignment(genes, self.genes, 'genes')
        ids = identifiers(ids, 'ids')
        if q.shape != (len(ids), len(ix)): raise ValueError('Query shape mismatch')
        if set(ids) & set(self.ids):
            raise ValueError('Query contains calibration compound IDs; in-sample predictions blocked')
        q = q[:, ix]
        pred = q * self.slopes + self.intercept
        if not np.isfinite(pred).all(): raise ValueError('Nonfinite predictions')
        # Descriptive support only: not a calibrated biological confidence score.
        train = self.training_source
        center = train.mean(0)
        scale = max(float(np.median(np.sqrt(np.mean((train - center)**2, axis=1)))), 1e-12)
        relative_distance = np.sqrt(np.mean((q - center)**2, axis=1)) / scale
        return pred, {'distance_to_calibration_center': relative_distance,
                      'distance_is_calibrated_confidence': False}

    def save(self, path):
        if not hasattr(self, 'slopes'): raise ValueError('Model is not fitted')
        save_archive(path, {'schema_version': SCHEMA, 'kind': 'model',
                     'algorithm': 'single_source_partial_pooling_v1', 'contract': self.contract,
                     'fit_report': self.report, 'provenance': provenance()},
                     slopes=self.slopes, intercept=self.intercept, genes=self.genes,
                     calibration_ids=self.ids, training_source=self.training_source)

    @classmethod
    def load(cls, path):
        meta, a = load_archive(path)
        if meta.get('kind') != 'model' or meta.get('algorithm') != 'single_source_partial_pooling_v1':
            raise ValueError('Unsupported model kind/algorithm')
        if set(a) != {'slopes', 'intercept', 'genes', 'calibration_ids', 'training_source'}:
            raise ValueError('Model arrays missing or unexpected')
        obj = cls()
        obj.contract = check_contract(meta)
        obj.genes = identifiers(a['genes'], 'genes')
        obj.ids = identifiers(a['calibration_ids'], 'calibration_ids')
        obj.training_source = matrix(a['training_source'], 'training_source')
        if len(obj.ids) < 4 or obj.training_source.shape != (len(obj.ids), len(obj.genes)):
            raise ValueError('Invalid model training array')
        for name in ('slopes', 'intercept'):
            v = a[name]
            if v.dtype.kind not in 'fiu' or v.shape != (len(obj.genes),) or not np.isfinite(v).all():
                raise ValueError('Invalid model coefficient dimensions/values')
            setattr(obj, name, np.asarray(v, float))
        obj.report = meta.get('fit_report')
        if not isinstance(obj.report, dict): raise ValueError('Missing fit report')
        return obj
