"""Verify preserved branch tips and exact snapshots without remote branch names."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def main():
    manifest = json.loads((ROOT / 'archive/branch-heads-20261014.json').read_text())
    rows = manifest['branches']
    if len(rows) != 36 or len({r['branch'] for r in rows}) != 36:
        raise ValueError('Expected the 36 distinct branch tips recorded before consolidation')
    snapshots = 0
    for row in rows:
        if git('rev-parse', row['commit'] + '^{tree}') != row['tree']:
            raise ValueError('Historical commit tree mismatch: ' + row['branch'])
        subprocess.run(['git', 'merge-base', '--is-ancestor', row['commit'], 'HEAD'],
                       cwd=ROOT, check=True)
        if 'snapshot' in row:
            if git('rev-parse', 'HEAD:' + row['snapshot']) != row['tree']:
                raise ValueError('Historical snapshot changed: ' + row['branch'])
            snapshots += 1
    if snapshots != 29:
        raise ValueError('Expected 29 unchanged original snapshots')
    print(json.dumps({'preserved_branch_tips': len(rows), 'exact_snapshots': snapshots,
                      'all_tips_are_ancestors': True, 'head': git('rev-parse', 'HEAD')}))


if __name__ == '__main__':
    main()
