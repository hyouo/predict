"""Preview or atomically remove preserved branches from hyouo/predict.

Run ``python -m tools.prune_branches`` to inspect the plan. Only ``--apply``
deletes remote branches. Git authentication must already be available locally.
No branch contents are rewritten, and a rejected atomic push is never retried.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path('archive/branch-heads-20261014.json')
INTEGRATION_BRANCH = 'fix/public-repair-20260930'
ALLOWED_URLS = {
    'https://github.com/hyouo/predict',
    'https://github.com/hyouo/predict.git',
    'git@github.com:hyouo/predict.git',
}
SHA1 = re.compile(r'[0-9a-f]{40}')


class CleanupError(ValueError):
    """A safety precondition failed; no fallback deletion is permitted."""


def _git(repo, *args, acceptable=(0,)):
    result = subprocess.run(
        ['git', '--no-replace-objects', '-C', str(repo), *args],
        capture_output=True, text=True, check=False)
    if result.returncode not in acceptable:
        detail = result.stderr.strip() or result.stdout.strip() or f'exit {result.returncode}'
        raise CleanupError(f'git {args[0]} failed: {detail}')
    return result


def _validate_origin(repo):
    urls = {}
    for kind, args in (
        ('fetch', ('remote', 'get-url', '--all', 'origin')),
        ('push', ('remote', 'get-url', '--push', '--all', 'origin')),
    ):
        values = _git(repo, *args).stdout.splitlines()
        if len(values) != 1 or values[0] not in ALLOWED_URLS:
            raise CleanupError(f'origin {kind} URL must point only to hyouo/predict on GitHub')
        urls[kind] = values[0]
    return urls


def _manifest(repo):
    data = json.loads((Path(repo) / MANIFEST).read_text(encoding='utf-8'))
    if (not isinstance(data, dict) or type(data.get('schema_version')) is not int
            or data['schema_version'] != 1 or data.get('canonical_branch') != 'main'
            or not isinstance(data.get('branches'), list)):
        raise CleanupError('Unsupported branch preservation manifest')
    recorded = {}
    for entry in data['branches']:
        if not isinstance(entry, dict):
            raise CleanupError('Malformed branch manifest entry')
        branch, commit = entry.get('branch'), entry.get('commit')
        if (not isinstance(branch, str) or not branch or branch in recorded
                or not isinstance(commit, str) or SHA1.fullmatch(commit) is None):
            raise CleanupError('Invalid or duplicate branch/commit in manifest')
        _git(repo, 'check-ref-format', f'refs/heads/{branch}')
        recorded[branch] = commit
    if 'main' not in recorded:
        raise CleanupError('Manifest must preserve main')
    return recorded


def _remote_heads(repo):
    result = {}
    for line in _git(repo, 'ls-remote', '--heads', 'origin').stdout.splitlines():
        fields = line.split()
        if (len(fields) != 2 or SHA1.fullmatch(fields[0]) is None
                or not fields[1].startswith('refs/heads/')):
            raise CleanupError('Malformed remote branch response')
        branch = fields[1][len('refs/heads/'):]
        if not branch or branch in result:
            raise CleanupError('Empty or duplicate remote branch')
        result[branch] = fields[0]
    if 'main' not in result:
        raise CleanupError('Remote main is missing; no branch may be deleted')
    return result


def _require_default_main(repo, main_commit):
    lines = _git(repo, 'ls-remote', '--symref', 'origin', 'HEAD').stdout.splitlines()
    symbolic = [line.split() for line in lines if line.startswith('ref:')]
    heads = [line.split() for line in lines if not line.startswith('ref:')]
    if symbolic != [['ref:', 'refs/heads/main', 'HEAD']]:
        raise CleanupError('Remote default branch must be main; refusing to delete any default branch')
    if heads != [[main_commit, 'HEAD']]:
        raise CleanupError('Remote main changed while checking its default branch; rerun preview')


def _require_ancestor(repo, commit, main_commit):
    result = _git(repo, 'merge-base', '--is-ancestor', commit, main_commit, acceptable=(0, 1))
    if result.returncode != 0:
        raise CleanupError(f'Commit {commit} is not preserved in remote main; deletion refused')


def prune_branches(repo=ROOT, *, apply=False, integration_branch=INTEGRATION_BRANCH):
    """Validate the fixed manifest and return the preview or applied cleanup record.

    Only the named integration branch may have advanced since the manifest was
    written. Both its recorded and current heads must be ancestors of fetched
    main. A full Git history is required for all preservation checks.
    """
    if integration_branch != INTEGRATION_BRANCH:
        raise CleanupError(f'The only permitted integration branch is {INTEGRATION_BRANCH}')
    urls = _validate_origin(repo)
    recorded = _manifest(repo)
    heads = _remote_heads(repo)
    unknown = set(heads) - set(recorded)
    if unknown:
        raise CleanupError(f'Unrecorded remote branches; review before deleting: {sorted(unknown)}')
    for branch, commit in heads.items():
        if branch not in ('main', integration_branch) and commit != recorded[branch]:
            raise CleanupError(f'Remote head changed for {branch}; deletion refused')
    main_commit = heads['main']
    _require_default_main(repo, main_commit)
    _git(repo, 'fetch', '--no-tags', 'origin', 'main')
    fetched = _git(repo, 'rev-parse', '--verify', 'FETCH_HEAD^{commit}').stdout.strip()
    if fetched != main_commit:
        raise CleanupError('Remote main changed during fetch; rerun preview')
    # Check every recorded commit, including branches already removed elsewhere.
    preserved = set(recorded.values()) | set(heads.values())
    for commit in sorted(preserved):
        _require_ancestor(repo, commit, main_commit)
    deletions = [{
        'branch': branch,
        'recorded_commit': recorded[branch],
        'expected_remote_commit': heads[branch],
        'integration_head_advanced': heads[branch] != recorded[branch],
    } for branch in sorted(heads) if branch != 'main']
    command = ['push', '--atomic', '--no-follow-tags']
    command += [f'--force-with-lease=refs/heads/{item["branch"]}:{item["expected_remote_commit"]}'
                for item in deletions]
    command += ['origin', *[f':refs/heads/{item["branch"]}' for item in deletions]]
    report = {
        'status': 'preview', 'repository': 'hyouo/predict', 'remote': 'origin',
        'remote_urls': urls, 'default_branch': 'main', 'main_commit': main_commit,
        'manifest': MANIFEST.as_posix(), 'preserved_commit_count': len(preserved),
        'delete_count': len(deletions), 'delete': deletions,
        'push_command': ['git', *command] if deletions else [],
        'apply_requested': apply,
    }
    if not apply:
        return report
    # Freeze the checked plan. A newly created/advanced branch needs a new preview.
    if _remote_heads(repo) != heads:
        raise CleanupError('Remote branch set changed after validation; deletion refused')
    _require_default_main(repo, main_commit)
    if deletions:
        _git(repo, *command)  # Exactly one atomic, leased deletion; no retry/fallback.
    remaining = _remote_heads(repo)
    if set(remaining) != {'main'}:
        raise CleanupError('Cleanup verification failed: remote branches other than main remain')
    if remaining['main'] != main_commit:
        raise CleanupError('Remote main changed during cleanup; inspect the remote before continuing')
    _require_default_main(repo, main_commit)
    report['status'] = 'applied' if deletions else 'already_clean'
    report['remaining_branches'] = ['main']
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true',
                        help='Apply the verified atomic deletion; default is preview only')
    parser.add_argument('--integration-branch', choices=[INTEGRATION_BRANCH],
                        default=INTEGRATION_BRANCH,
                        help='Only this integration branch may advance if fully preserved in main')
    args = parser.parse_args(argv)
    try:
        report = prune_branches(apply=args.apply, integration_branch=args.integration_branch)
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0
    except (OSError, ValueError) as exc:
        print(json.dumps({'status': 'error', 'message': str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print(json.dumps({'status': 'error', 'message': 'Interrupted; inspect remote state before retrying.'}),
              file=sys.stderr)
        return 130


if __name__ == '__main__':
    sys.exit(main())
