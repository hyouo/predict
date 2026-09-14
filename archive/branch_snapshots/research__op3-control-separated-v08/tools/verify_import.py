"""Verify the original v0.7 import. Later algorithm changes should intentionally fail this audit."""
import json
from pathlib import Path
from tools.assets import sha256_file

ROOT = Path(__file__).resolve().parents[1]

if __name__ == '__main__':
    manifest = json.loads((ROOT / 'docs/history/v0.7_import.json').read_text())
    changed = [p for p, digest in manifest['files'].items()
               if not (ROOT / p).is_file() or sha256_file(ROOT / p) != digest]
    print(json.dumps({'matched': not changed, 'files_checked': len(manifest['files']), 'changed': changed}, indent=2))
    raise SystemExit(bool(changed))
