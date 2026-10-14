"""Versioned, non-pickle arrays with explicit identifiers and bounded loading."""
from __future__ import annotations
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import os
import tempfile
import zipfile
import zlib
import numpy as np

MAX_BYTES = 512 * 1024**2
SCHEMA = 1


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


@contextmanager
def atomic_output(path, *, prefix='.output-'):
    """Write in the destination directory, then publish once without replacing it."""
    destination = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w+b', delete=False,
                                         dir=destination.parent, prefix=prefix, suffix='.tmp') as f:
            temporary = Path(f.name)
            yield f
            f.flush()
            os.fsync(f.fileno())
        os.link(temporary, destination)  # Atomic publication, never overwrite an existing result.
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def json_write(path, value):
    """Publish complete JSON exclusively; serialization/I/O errors leave no partial file."""
    text = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
    with atomic_output(path, prefix='.json-') as f:
        f.write(text.encode('utf-8'))


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
        json_write(p / 'COMPLETE.json', {'status': 'complete', **provenance()})
    except (Exception, KeyboardInterrupt) as exc:
        try:
            json_write(p / 'FAILED.json', {'status': 'failed', 'error_type': type(exc).__name__,
                                         'error': str(exc), **provenance()})
        except Exception as marker_error:
            # Preserve the operation's cause even when the disk cannot store its marker.
            exc.add_note(f'Could not write FAILED.json: {marker_error}')
        raise


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
    """Atomically publish a bounded archive; object arrays/pickle are never written."""
    arrays = {key: np.asarray(value) for key, value in arrays.items()}
    for key, value in arrays.items():
        if (not key.isidentifier() or key in {'metadata_json', 'file', 'allow_pickle'}
                or value.dtype.hasobject or value.dtype.kind not in 'biufcUS'):
            raise ValueError('Invalid/reserved name or unsupported array dtype')
    text = json.dumps(metadata, ensure_ascii=False, allow_nan=False, sort_keys=True)
    arrays = {'metadata_json': np.array(text), **arrays}
    if len(arrays) > 32 or sum(a.nbytes for a in arrays.values()) > MAX_BYTES:
        raise ValueError('Expanded archive too large')
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with atomic_output(p) as f:
        # Dtypes above exclude pickle even on older supported NumPy releases.
        np.savez_compressed(f, **arrays)
        f.flush()
        if f.tell() > MAX_BYTES:
            raise ValueError('Archive exceeds 512 MiB limit')
        with zipfile.ZipFile(f) as z:
            if sum(i.file_size for i in z.infolist()) > MAX_BYTES:
                raise ValueError('Expanded archive too large')


def _validate_array_header(stream, size):
    version = np.lib.format.read_magic(stream)
    readers = {(1, 0): np.lib.format.read_array_header_1_0,
               (2, 0): np.lib.format.read_array_header_2_0}
    if version not in readers:
        raise ValueError('Unsupported NPY format; use plain numeric or string arrays')
    shape, _, dtype = readers[version](stream, max_header_size=10000)
    if dtype.hasobject or dtype.kind not in 'biufcUS':
        raise ValueError('Object/structured arrays are prohibited')
    if any(n < 0 for n in shape):
        raise ValueError('Invalid array shape')
    count = math.prod(shape)
    nbytes = count * dtype.itemsize
    if count > MAX_BYTES or nbytes > MAX_BYTES or stream.tell() + nbytes != size:
        raise ValueError('Array shape/size mismatch or expanded array too large')


def _reject_json_constant(value):
    raise ValueError(f'Nonfinite JSON value: {value}')


def load_archive(path):
    p = Path(path)
    if p.stat().st_size > MAX_BYTES:
        raise ValueError('Archive exceeds 512 MiB limit')
    try:
        with zipfile.ZipFile(p) as z:
            entries = z.infolist()
            if len(entries) > 32 or sum(i.file_size for i in entries) > MAX_BYTES:
                raise ValueError('Expanded archive too large')
            names = [i.filename for i in entries]
            if len(set(names)) != len(names) or any(
                    not n.endswith('.npy') or not n[:-4].isidentifier() for n in names):
                raise ValueError('Invalid archive entries')
            # A tiny ZIP member can claim a huge NPY shape. Check all headers
            # against actual payload sizes BEFORE NumPy allocates any array.
            for entry in entries:
                with z.open(entry) as stream:
                    _validate_array_header(stream, entry.file_size)
            arrays = {}
            for entry in entries:
                with z.open(entry) as stream:
                    arrays[entry.filename[:-4]] = np.lib.format.read_array(
                        stream, allow_pickle=False, max_header_size=10000)
    except (zipfile.BadZipFile, EOFError, RuntimeError, NotImplementedError, zlib.error) as exc:
        raise ValueError(f'Invalid or corrupted NPZ archive: {exc}') from exc
    raw = arrays.pop('metadata_json', None)
    if raw is None or raw.shape != () or raw.dtype.kind not in 'US':
        raise ValueError('Missing scalar metadata_json')
    meta = json.loads(str(raw), parse_constant=_reject_json_constant)
    if (not isinstance(meta, dict) or type(meta.get('schema_version')) is not int
            or meta['schema_version'] != SCHEMA):
        raise ValueError('Unsupported schema_version')
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
