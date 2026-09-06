"""Acquire public pinned assets; never substitute synthetic data or overwrite bad caches."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify_file(path: Path, expected_sha256: str, max_bytes: int) -> dict:
    size = path.stat().st_size
    if size > max_bytes:
        raise ValueError('Asset exceeds explicit size cap')
    digest = sha256_file(path)
    if digest != expected_sha256:
        raise ValueError('SHA256 mismatch; refusing to use or overwrite this file')
    return {'bytes': size, 'sha256': digest}


def materialize(url: str, destination: Path, expected_sha256: str,
                max_bytes: int, source: Path | None = None) -> dict:
    """Bounded atomic transfer. An optional local source supports offline reproduction."""
    if urlsplit(url).scheme != 'https':
        raise ValueError('Only HTTPS public asset URLs are accepted')
    if not re.fullmatch(r'[0-9a-f]{64}', expected_sha256):
        raise ValueError('A pinned lowercase SHA256 is required')
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes <= 0:
        raise ValueError('A positive integer byte cap is required')
    destination = Path(destination)
    if destination.exists():
        return {**verify_file(destination, expected_sha256, max_bytes), 'cached': True}
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix='.download-',
                                         suffix='.partial', delete=False) as output:
            temporary = Path(output.name)
            request = Request(url, headers={'User-Agent': 'predict-research-data-audit/1'})
            stream = Path(source).open('rb') if source is not None else urlopen(request, timeout=60)
            with stream:
                if hasattr(stream, 'geturl') and urlsplit(stream.geturl()).scheme != 'https':
                    raise ValueError('Refusing an HTTPS downgrade redirect')
                size = 0
                while chunk := stream.read(1024 * 1024):
                    size += len(chunk)
                    if size > max_bytes:
                        raise ValueError('Asset exceeds explicit size cap')
                    output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        result = verify_file(temporary, expected_sha256, max_bytes)
        os.replace(temporary, destination)
        return {**result, 'cached': False, 'retrieval': 'local_verified_copy' if source else 'https'}
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def inspect_h5ad(path: Path) -> dict:
    """Read structure only. This is not approval of assay or control semantics."""
    import h5py
    def matrix(node):
        if isinstance(node, h5py.Dataset):
            return {'encoding': 'dense', 'shape': list(node.shape), 'dtype': str(node.dtype)}
        if isinstance(node, h5py.Group):
            shape = node.attrs.get('shape')
            enc = node.attrs.get('encoding-type', 'unknown')
            if isinstance(enc, bytes):
                enc = enc.decode()
            return {'encoding': str(enc), 'shape': None if shape is None else [int(n) for n in shape]}
        return None
    with h5py.File(path, 'r') as f:
        return {'top_level_keys': list(f), 'X': matrix(f.get('X')),
                'layers': {k: matrix(v) for k, v in f.get('layers', {}).items()},
                'obs_columns': list(f['obs']) if 'obs' in f else [],
                'var_columns': list(f['var']) if 'var' in f else [],
                'uns_keys': list(f['uns']) if 'uns' in f else [],
                'approved_for_training': False,
                'warning': 'Verify units, normalization, donor/well/replicate identity, controls and splits before training.'}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--asset', choices=['auxiliary', 'op3', 'sciplex'], required=True)
    parser.add_argument('--source', type=Path, help='Optional offline copy, still checked against the pinned hash')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'assets/manifest.json').read_text())
    asset = manifest['assets'][args.asset]
    output = (ROOT / asset['destination']).resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError('Destination must remain inside this checkout')
    report = materialize(asset['url'], output, asset['sha256'], asset['max_bytes'], args.source)
    report.update({'asset': args.asset, 'source_url': asset['url'], 'path': asset['destination'],
                   'modality': asset['modality'], 'trained_model': False})
    if args.asset == 'auxiliary':
        from src.data import load_panel
        report['metadata'] = load_panel(output).audit
    else:
        report['metadata'] = inspect_h5ad(output)
    audit = output.with_suffix('.audit.json')
    audit.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
