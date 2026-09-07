"""Strict adapter for one pinned OP3 asset; private outcomes always denied.

Only this conditional, already-filtered panel is supported here. For other data,
prepare explicitly aligned effect bundles with perturb_predict.io.write_bundle.
"""
from __future__ import annotations
from pathlib import Path
from urllib.request import urlopen
import hashlib
import h5py
import numpy as np
from .io import digest, json_write, provenance, write_bundle, run_directory

DATA_SHA256 = '8eaf7e63adc68029e88184726b8545bc3836c3ef2a25b9185dc4e8b819537e59'
DATA_URL = ('https://huggingface.co/datasets/theislab/chem-perturbridge/resolve/'
            '6c54c8eb0321cceff4f888c54e199077d055e20b/op3/'
            'op3_standardized_processed.h5ad?download=true')
MAX_DATA_BYTES = 32 * 1024**2
SOURCE = 'CL_0000084'
NK = 'CL_0000623'
TARGETS = {'B_cells': 'CL_0000236', 'Myeloid_cells': 'CL_0000763'}


def fetch(path):
    """Fixed HTTPS download with bounds/hash, no overwriting bad cached files."""
    p = Path(path)
    if p.exists():
        if p.stat().st_size > MAX_DATA_BYTES or digest(p) != DATA_SHA256:
            raise ValueError('Existing OP3 cache is invalid; it was not overwritten')
        return {'cached': True, 'sha256': DATA_SHA256, 'bytes': p.stat().st_size}
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name + '.partial')
    created = False
    try:
        with temp.open('xb') as f:
            created = True
            with urlopen(DATA_URL, timeout=60) as response:
                if not response.geturl().startswith('https://'):
                    raise ValueError('Non-HTTPS redirect rejected')
                total = 0
                while chunk := response.read(1024**2):
                    total += len(chunk)
                    if total > MAX_DATA_BYTES: raise ValueError('Download exceeds 32 MiB')
                    f.write(chunk)
        if digest(temp) != DATA_SHA256: raise ValueError('Downloaded OP3 hash mismatch')
        # Exclusive creation of final path, even under concurrent calls.
        with p.open('xb') as out, temp.open('rb') as inp:
            while chunk := inp.read(1024**2): out.write(chunk)
    finally:
        if created: temp.unlink(missing_ok=True)
    return {'cached': False, 'sha256': DATA_SHA256, 'bytes': p.stat().st_size}


def column(node):
    if isinstance(node, h5py.Dataset):
        if node.ndim != 1: raise ValueError('Metadata columns must be vectors')
        return node.asstr()[:] if h5py.check_string_dtype(node.dtype) else node[:]
    if not isinstance(node, h5py.Group) or not {'categories', 'codes'} <= set(node):
        raise ValueError('Unsupported categorical encoding')
    values, codes = column(node['categories']), column(node['codes'])
    if codes.dtype.kind not in 'iu' or np.any(codes < -1) or np.any(codes >= len(values)):
        raise ValueError('Invalid category codes')
    out = np.empty(len(codes), object)
    out[:] = None
    out[codes >= 0] = values[codes[codes >= 0]]
    return out


def guard_rows(obs, ids, purpose):
    """The scoring capability never enables private-test expression."""
    ids = np.asarray(ids)
    n = len(obs['split'])
    if ids.ndim != 1 or ids.dtype.kind not in 'iu' or np.any(ids < 0) or np.any(ids >= n):
        raise ValueError('Invalid row indices')
    if len(set(ids.tolist())) != len(ids): raise ValueError('Duplicate row indices')
    if purpose not in ('prepare', 'score'): raise ValueError('Unknown access purpose')
    splits = obs['split'][ids]
    if np.any(splits == 'private_test'): raise PermissionError('Private outcomes permanently locked')
    if purpose != 'score' and np.any(splits == 'public_test'):
        raise PermissionError('Public treatment values may only be read after prediction')
    if purpose != 'score':
        bad = [int(i) for i in ids if obs['is_control'][i] and str(obs['well'][i])[0] in 'EFGH']
        if bad: raise PermissionError('Reserved scoring controls cannot enter preparation')
    return ids


class Reader:
    def __init__(self, path):
        self.path = Path(path)
        if self.path.stat().st_size > MAX_DATA_BYTES or digest(self.path) != DATA_SHA256:
            raise ValueError('This adapter only accepts the pinned, SHA256-verified OP3 file')
        with h5py.File(self.path) as f:
            self.obs = {k: column(v) for k, v in f['obs'].items()}
            self.genes = column(f['var']['ensembl_id']).astype(str)
            self.shape = tuple(map(int, f['X'].attrs['shape']))
            if f['X'].attrs.get('encoding-type') != 'csr_matrix' or self.shape != (1813, 5288):
                raise ValueError('Unexpected pinned OP3 shape/encoding')
        self.log = []
        required = ('is_control', 'split', 'cell_type', 'plate', 'well', 'perturbagen',
                    'psbulk_counts', 'sample_id', 'pert_dose_uM', 'pert_time_h')
        if any(k not in self.obs or len(self.obs[k]) != self.shape[0] for k in required):
            raise ValueError('Missing/invalid OP3 metadata')
        if any(not isinstance(x, (bool, np.bool_)) for x in self.obs['is_control']):
            raise ValueError('Control status must be explicit Boolean')
        if len(set(self.obs['sample_id'])) != self.shape[0] or len(set(self.genes)) != self.shape[1]:
            raise ValueError('Duplicate sample or gene identifiers')

    def ids(self, ct, split=None, control_letters=None):
        o = self.obs
        match = o['cell_type'] == ct
        if control_letters is not None:
            match &= o['is_control'].astype(bool)
            match &= np.array([str(w)[0] in control_letters for w in o['well']])
        else:
            match &= ~o['is_control'].astype(bool)
        if split is not None: match &= o['split'] == split
        return np.flatnonzero(match)

    def logcpm(self, ids, purpose='prepare'):
        ids = guard_rows(self.obs, ids, purpose)
        event = {'ids': ids.tolist(), 'purpose': purpose,
                 'splits': sorted(set(map(str, self.obs['split'][ids]))),
                 'status': 'started'}
        self.log.append(event)
        try:
            out = np.zeros((len(ids), self.shape[1]), float)
            with h5py.File(self.path) as f:
                x = f['X']; pointers = x['indptr'][:]
                for j, i in enumerate(ids):
                    a, b = int(pointers[i]), int(pointers[i+1])
                    np.add.at(out[j], x['indices'][a:b], x['data'][a:b])
            totals = np.asarray(self.obs['psbulk_counts'][ids], float)
            if (not np.isfinite(out).all() or np.any(out < 0) or np.any(out != np.floor(out))
                    or not np.isfinite(totals).all() or np.any(totals <= 0)
                    or np.any(out.sum(1) > totals + 1e-6)):
                raise ValueError('Invalid count values or supplied library denominators')
            result = np.log2(1. + 1e6 * out / totals[:, None])
        except Exception as exc:
            event.update(status='failed', error_type=type(exc).__name__)
            raise
        else:
            event['status'] = 'completed'
            return result

    def effects(self, ct, split, letters, purpose='prepare'):
        controls = self.ids(ct, control_letters=letters)
        c = self.logcpm(controls, purpose)
        references = {}
        for plate in sorted(set(self.obs['plate'])):
            ix = self.obs['plate'][controls] == plate
            if int(ix.sum()) != len(letters): raise ValueError('Missing/extra plate-matched controls')
            references[plate] = c[ix].mean(0)
        ids = self.ids(ct, split=split)
        doses, times = self.obs['pert_dose_uM'][ids], self.obs['pert_time_h'][ids]
        if np.any(doses != 1.) or np.any(times != 24.): raise ValueError('Unexpected dose/time')
        values = self.logcpm(ids, purpose) - np.stack([references[p] for p in self.obs['plate'][ids]])
        names = self.obs['perturbagen'][ids]
        effects = {str(d): values[names == d].mean(0) for d in sorted(set(names))}
        return effects, ids


def contract(ct):
    return {'effect_space': 'log2(1+1e6*counts/psbulk_counts)-plate_matched_control_mean.v1',
            'source_context': SOURCE, 'target_context': ct, 'organism': 'human',
            'dose': '1 uM', 'time': '24 h', 'gene_namespace': 'OP3 pinned Ensembl IDs'}


def expected_public_queries(reader, ct):
    """One metadata-only query rule shared by preparation and OP3 evaluation."""
    if ct not in TARGETS.values(): raise ValueError('Unsupported OP3 target context')
    names = reader.obs['perturbagen']
    source = set(names[reader.ids(SOURCE, split='train')])
    nk = set(names[reader.ids(NK, split='train')])
    public = set(names[reader.ids(ct, split='public_test')])
    expected = sorted(source & nk & public)
    if not expected: raise ValueError('Empty expected OP3 query set')
    return expected


def prepare(data, out, *, accept_conditional=False):
    if not accept_conditional:
        raise ValueError('Pass --accept-conditional: upstream-selected panel, not strict zero-shot validation')
    with run_directory(out) as root:
        r = Reader(data)
        source, source_ids = r.effects(SOURCE, 'train', 'CD')
        # Reproduce v0.8's fixed common source panel using NK metadata ONLY.
        nk_names = set(r.obs['perturbagen'][r.ids(NK, split='train')])
        allowed = sorted(set(source) & nk_names)
        audit = {'dataset_sha256': DATA_SHA256, 'shape': list(r.shape),
                 'task': 'few_shot_unseen_target_compound_pair', 'source_context': SOURCE,
                 'source_reference': 'C/D controls per plate',
                 'target_fit_reference': 'A-D controls per plate',
                 'score_reference': 'E-H controls per plate',
                 'conditional_panel': True, 'private_expression_read': False,
                 'source_query_outcomes_allowed': True, 'targets': {}, **provenance(),
                 'limits': ['Upstream gene selection used the published combined panel.',
                            'Shared source/target wells remain; donor independence is not established.',
                            'Public evaluation has already been used for research development.',
                            'Not a new-drug, unseen-context, ATAC or single-cell distribution model.']}
        for name, ct in TARGETS.items():
            target, target_ids = r.effects(ct, 'train', 'ABCD')
            train_ids = sorted(set(target) & set(allowed))
            queries = expected_public_queries(r, ct)  # metadata only
            if set(train_ids) & set(queries): raise ValueError('Train/query compound overlap')
            dest = root / name
            for kind, ids in [('train', train_ids), ('query', queries)]:
                write_bundle(dest / f'{kind}.npz', kind, ids, r.genes, contract(ct),
                             source=np.stack([source[d] for d in ids]),
                             target=np.stack([target[d] for d in ids]) if kind == 'train' else None,
                             provenance_info={'dataset_sha256': DATA_SHA256, 'adapter': 'op3_v1',
                                              'target_labels_present': kind == 'train'})
            audit['targets'][name] = {'cell_type': ct, 'calibration_ids': train_ids, 'query_ids': queries,
                  'n_calibration_compounds': len(train_ids), 'n_query_compounds': len(queries),
                  'calibration_pseudobulk_rows_used': int(np.isin(r.obs['perturbagen'][target_ids], train_ids).sum()),
                  'target_control_rows_used': len(r.ids(ct, control_letters='ABCD')),
                  'source_control_rows_used': len(r.ids(SOURCE, control_letters='CD')),
                  'excluded_target_train_compounds': sorted(set(target) - set(allowed))}
        audit['reader_log'] = r.log
        if any(set(a['splits']) & {'public_test', 'private_test'} for a in r.log):
            raise RuntimeError('Leakage guard failure during prepare')
        json_write(root / 'audit.json', audit)
    return audit
