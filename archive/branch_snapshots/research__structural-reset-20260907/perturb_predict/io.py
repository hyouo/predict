"""Versioned, non-pickle arrays with explicit identifiers and bounded loading."""
from __future__ import annotations
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import zipfile
import numpy as np

MAX_BYTES = 512 * 1024**2
SCHEMA = 1


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def json_write(path, value):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def provenance():
    from . import __version__
    return {'package_version': __version__, 'python': platform.python_version(),
            'numpy': np.__version__, 'utc': datetime.now(timezone.utc).isoformat(),
            'code_sha256': {p.name: digest(p) for p in sorted(Path(__file__).parent.glob('*.py'))}}


@contextmanager
def run_directory(path):
    """Reserve an output once. Failed runs stay marked; never overwrite/reuse."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=False)
    try:
        yield p
    except Exception as exc:
        json_write(p / 'FAILED.json', {'status': 'failed', 'error_type': type(exc).__name__,
                                     'error': str(exc), **provenance()})
        raise
    else:
        json_write(p / 'COMPLETE.json', {'status': 'complete', **provenance()})


def identifiers(value, name):
    a = np.asarray(value)
    if a.ndim != 1 or a.dtype.kind not in 'US' or not len(a):
        raise ValueError(f'{name} must be a nonempty vector of Unicode strings')
    a = a.astype(str)
    if any(not v.strip() for v in a) or len(set(a.tolist())) != len(a):
        raise ValueError(f'{name} contains empty/duplicate identifiers')
    return a


def matrix(value, name):
    a = np.asarray(value)
    if a.ndim != 2 or min(a.shape) == 0 or a.dtype.kind not in 'fiu':
        raise ValueError(f'{name} must be a nonempty numeric matrix')
    if a.size * 8 > MAX_BYTES:
        raise ValueError(f'{name} exceeds the 512 MiB float64 matrix limit')
    a = np.asarray(a, dtype=np.float64)
    if not np.isfinite(a).all():
        raise ValueError(f'{name} contains NaN/infinity; missing is not zero')
    return a


def save_archive(path, metadata, **arrays):
    """Exclusive file creation; object arrays/pickle are never written."""
    for key, value in arrays.items():
        if key == 'metadata_json' or np.asarray(value).dtype.kind == 'O':
            raise ValueError('Reserved name or object array')
    text = json.dumps(metadata, ensure_ascii=False, allow_nan=False, sort_keys=True)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('xb') as f:
        np.savez_compressed(f, metadata_json=np.array(text), **arrays)


def load_archive(path):
    p = Path(path)
    if p.stat().st_size > MAX_BYTES:
        raise ValueError('Archive exceeds 512 MiB limit')
    with zipfile.ZipFile(p) as z:
        entries = z.infolist()
        if len(entries) > 32 or sum(i.file_size for i in entries) > MAX_BYTES:
            raise ValueError('Expanded archive too large')
        names = [i.filename for i in entries]
        if len(set(names)) != len(names) or any('/' in n or not n.endswith('.npy') for n in names):
            raise ValueError('Invalid archive entries')
    with np.load(p, allow_pickle=False, max_header_size=10000) as z:
        arrays = {k: z[k] for k in z.files}
    raw = arrays.pop('metadata_json', None)
    if raw is None or raw.shape != () or raw.dtype.kind not in 'US':
        raise ValueError('Missing scalar metadata_json')
    meta = json.loads(str(raw))
    if not isinstance(meta, dict) or meta.get('schema_version') != SCHEMA:
        raise ValueError('Unsupported schema_version')
    if any(a.dtype.kind == 'O' for a in arrays.values()):
        raise ValueError('Object arrays are prohibited')
    return meta, arrays


CONTRACT_FIELDS = ('effect_space', 'source_context', 'target_context', 'organism',
                   'dose', 'time', 'gene_namespace')


def check_contract(meta):
    contract = meta.get('contract')
    if not isinstance(contract, dict) or any(not isinstance(contract.get(k), str) or
                                           not contract[k].strip() for k in CONTRACT_FIELDS):
        raise ValueError('Explicit effect/context/organism/dose/time/gene contract required')
    return contract


def load_bundle(path, kind):
    meta, a = load_archive(path)
    if meta.get('kind') != kind:
        raise ValueError(f'Expected {kind} bundle, received {meta.get("kind")}')
    expected = {'train': {'source', 'target', 'ids', 'genes'},
                'query': {'source', 'ids', 'genes'},
                'truth': {'target', 'ids', 'genes'}}
    if kind not in expected or set(a) != expected[kind]:
        raise ValueError('Unexpected/missing bundle arrays; query bundles cannot contain target labels')
    check_contract(meta)
    a['ids'] = identifiers(a['ids'], 'ids')
    a['genes'] = identifiers(a['genes'], 'genes')
    for k in ('source', 'target'):
        if k in a:
            a[k] = matrix(a[k], k)
            if a[k].shape != (len(a['ids']), len(a['genes'])):
                raise ValueError('Matrix does not match compound and gene IDs')
    return meta, a


def write_bundle(path, kind, ids, genes, contract, *, source=None, target=None, provenance_info=None):
    a = {'ids': identifiers(ids, 'ids'), 'genes': identifiers(genes, 'genes')}
    if source is not None: a['source'] = matrix(source, 'source')
    if target is not None: a['target'] = matrix(target, 'target')
    meta = {'schema_version': SCHEMA, 'kind': kind, 'contract': contract,
            'provenance': provenance_info or {}}
    check_contract(meta)
    expected = {'train': {'source', 'target', 'ids', 'genes'},
                'query': {'source', 'ids', 'genes'}, 'truth': {'target', 'ids', 'genes'}}
    if kind not in expected or set(a) != expected[kind]: raise ValueError('Bundle kind/array mismatch')
    for key in ('source', 'target'):
        if key in a and a[key].shape != (len(a['ids']), len(a['genes'])):
            raise ValueError('Bundle shape/identifier mismatch')
    save_archive(path, meta, **a)


def alignment(actual, expected, name):
    actual = identifiers(actual, name)
    expected = identifiers(expected, name)
    if set(actual) != set(expected):
        missing, extra = set(expected) - set(actual), set(actual) - set(expected)
        raise ValueError(f'{name} mismatch: missing={len(missing)}, extra={len(extra)}; no imputation')
    lookup = {s: i for i, s in enumerate(actual)}
    return np.array([lookup[s] for s in expected], dtype=int)
