"""Branch cleanup tests use mocks or disposable local repositories, never GitHub."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from tools import prune_branches as cleanup


OLD_MAIN = '1' * 40
OLD_BRANCH = '2' * 40
OLD_INTEGRATION = '3' * 40
NEW_MAIN = '4' * 40
NEW_INTEGRATION = '5' * 40
OTHER = '6' * 40
ARCHIVED_BRANCH = 'research/archived'


def write_manifest(repo, heads):
    path = repo / cleanup.MANIFEST
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        'schema_version': 1, 'canonical_branch': 'main',
        'branches': [{'branch': branch, 'commit': commit} for branch, commit in heads.items()],
    }), encoding='utf-8')


class FakeGit:
    def __init__(self):
        self.heads = {'main': NEW_MAIN, ARCHIVED_BRANCH: OLD_BRANCH,
                      cleanup.INTEGRATION_BRANCH: OLD_INTEGRATION}
        self.default = 'main'
        self.fetch_url = self.push_url = 'https://github.com/hyouo/predict.git'
        self.ancestors = {OLD_MAIN, OLD_BRANCH, OLD_INTEGRATION, NEW_MAIN, NEW_INTEGRATION}
        self.calls = []
        self.fail_push = False
        self.after_push = None

    def __call__(self, repo, *args, acceptable=(0,)):
        self.calls.append(args)
        output, code = '', 0
        if args[0] == 'remote':
            output = self.push_url if '--push' in args else self.fetch_url
        elif args[:2] == ('ls-remote', '--heads'):
            output = '\n'.join(f'{commit}\trefs/heads/{branch}'
                               for branch, commit in sorted(self.heads.items()))
        elif args[:2] == ('ls-remote', '--symref'):
            output = f'ref: refs/heads/{self.default}\tHEAD\n{self.heads[self.default]}\tHEAD\n'
        elif args[0] == 'rev-parse':
            output = self.heads['main']
        elif args[0] == 'merge-base':
            code = 0 if args[2] in self.ancestors else 1
        elif args[0] == 'push':
            if self.fail_push:
                raise cleanup.CleanupError('atomic push rejected')
            self.heads = {'main': self.heads['main']}
            if self.after_push:
                self.heads.update(self.after_push)
        elif args[0] not in ('check-ref-format', 'fetch'):
            raise AssertionError(f'Unexpected git command: {args}')
        return subprocess.CompletedProcess(args, code, output.rstrip('\n') + '\n', '')

    @property
    def pushes(self):
        return [call for call in self.calls if call[0] == 'push']


class BranchCleanupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name)
        write_manifest(self.repo, {'main': OLD_MAIN, ARCHIVED_BRANCH: OLD_BRANCH,
                                   cleanup.INTEGRATION_BRANCH: OLD_INTEGRATION})
        self.git = FakeGit()

    def invoke(self, **kwargs):
        with patch.object(cleanup, '_git', side_effect=self.git):
            return cleanup.prune_branches(self.repo, **kwargs)

    def test_preview_checks_history_without_pushing(self):
        report = self.invoke()
        self.assertEqual(report['status'], 'preview')
        self.assertFalse(report['apply_requested'])
        self.assertEqual(report['delete_count'], 2)
        self.assertEqual(self.git.pushes, [])
        self.assertIn(('fetch', '--no-tags', 'origin', 'main'), self.git.calls)
        checked = {call[2] for call in self.git.calls if call[0] == 'merge-base'}
        self.assertEqual(checked, {OLD_MAIN, OLD_BRANCH, OLD_INTEGRATION, NEW_MAIN})

    def test_cli_defaults_to_preview_and_prints_json(self):
        with patch.object(cleanup, 'prune_branches', return_value={'status': 'preview'}) as prune, \
             contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(cleanup.main([]), 0)
        self.assertEqual(json.loads(output.getvalue()), {'status': 'preview'})
        prune.assert_called_once_with(apply=False, integration_branch=cleanup.INTEGRATION_BRANCH)

    def test_changed_nonintegration_head_and_unknown_branch_refuse_deletion(self):
        for changed in ({ARCHIVED_BRANCH: OTHER}, {'unrecorded/new': OTHER}):
            with self.subTest(changed=changed):
                self.git = FakeGit()
                self.git.heads.update(changed)
                with self.assertRaises(cleanup.CleanupError):
                    self.invoke(apply=True)
                self.assertEqual(self.git.pushes, [])
                self.assertFalse(any(call[0] == 'fetch' for call in self.git.calls))

    def test_original_and_updated_integration_heads_must_be_preserved(self):
        self.git.heads[cleanup.INTEGRATION_BRANCH] = NEW_INTEGRATION
        report = self.invoke()
        item = next(item for item in report['delete'] if item['branch'] == cleanup.INTEGRATION_BRANCH)
        self.assertTrue(item['integration_head_advanced'])
        self.assertEqual(item['expected_remote_commit'], NEW_INTEGRATION)
        for unpreserved in (OLD_INTEGRATION, NEW_INTEGRATION):
            with self.subTest(commit=unpreserved):
                self.git = FakeGit()
                self.git.heads[cleanup.INTEGRATION_BRANCH] = NEW_INTEGRATION
                self.git.ancestors.remove(unpreserved)
                with self.assertRaisesRegex(cleanup.CleanupError, 'not preserved'):
                    self.invoke(apply=True)
                self.assertEqual(self.git.pushes, [])

    def test_no_unpreserved_manifest_commit_may_be_ignored(self):
        del self.git.heads[ARCHIVED_BRANCH]
        self.git.ancestors.remove(OLD_BRANCH)
        with self.assertRaisesRegex(cleanup.CleanupError, 'not preserved'):
            self.invoke(apply=True)
        self.assertEqual(self.git.pushes, [])

    def test_apply_uses_exact_atomic_deletions_and_expected_head_leases(self):
        report = self.invoke(apply=True)
        self.assertEqual(report['status'], 'applied')
        self.assertEqual(report['remaining_branches'], ['main'])
        self.assertEqual(len(self.git.pushes), 1)
        command = self.git.pushes[0]
        self.assertIn('--atomic', command)
        self.assertIn('--no-follow-tags', command)
        for branch, commit in ((ARCHIVED_BRANCH, OLD_BRANCH),
                               (cleanup.INTEGRATION_BRANCH, OLD_INTEGRATION)):
            self.assertIn(f'--force-with-lease=refs/heads/{branch}:{commit}', command)
            self.assertIn(f':refs/heads/{branch}', command)
        self.assertNotIn(':refs/heads/main', command)
        self.assertNotIn('--force', command)
        refs = [arg for arg in command if arg.startswith(':')]
        self.assertEqual(set(refs), {f':refs/heads/{ARCHIVED_BRANCH}',
                                    f':refs/heads/{cleanup.INTEGRATION_BRANCH}'})

    def test_default_branch_and_different_push_destination_are_rejected(self):
        self.git.default = ARCHIVED_BRANCH
        with self.assertRaisesRegex(cleanup.CleanupError, 'default branch'):
            self.invoke(apply=True)
        self.assertEqual(self.git.pushes, [])
        self.git = FakeGit()
        self.git.push_url = 'https://github.com/another/repository.git'
        with self.assertRaisesRegex(cleanup.CleanupError, 'push URL'):
            self.invoke(apply=True)
        self.assertEqual(self.git.pushes, [])

    def test_push_failure_has_no_retry_or_unprotected_fallback(self):
        self.git.fail_push = True
        with self.assertRaisesRegex(cleanup.CleanupError, 'atomic push rejected'):
            self.invoke(apply=True)
        self.assertEqual(len(self.git.pushes), 1)
        self.assertEqual(len(self.git.heads), 3)

    def test_post_apply_verification_rejects_extra_remote_branch(self):
        self.git.after_push = {'new/concurrent': OTHER}
        with self.assertRaisesRegex(cleanup.CleanupError, 'other than main remain'):
            self.invoke(apply=True)
        self.assertEqual(len(self.git.pushes), 1)

    def test_already_clean_is_noop_and_cannot_select_arbitrary_exception_branch(self):
        self.git.heads = {'main': NEW_MAIN}
        report = self.invoke(apply=True)
        self.assertEqual(report['status'], 'already_clean')
        self.assertEqual(self.git.pushes, [])
        with self.assertRaises(cleanup.CleanupError):
            self.invoke(apply=True, integration_branch=ARCHIVED_BRANCH)
        self.assertEqual(self.git.pushes, [])


class LocalGitCleanupTests(unittest.TestCase):
    """Run real Git leases/atomic pushes only against a disposable bare directory."""
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo, self.remote = self.root / 'checkout', self.root / 'remote.git'
        self.git('init', '--initial-branch=main', str(self.repo), cwd=self.root)
        self.git('init', '--bare', '--initial-branch=main', str(self.remote), cwd=self.root)
        self.git('config', 'user.email', 'fixture@example.invalid')
        self.git('config', 'user.name', 'Branch Cleanup Test')
        (self.repo / 'fixture.txt').write_text('first\n')
        self.git('add', 'fixture.txt')
        self.git('-c', 'commit.gpgsign=false', 'commit', '-m', 'first fixture')
        self.old = self.git('rev-parse', 'HEAD').stdout.strip()
        self.git('branch', ARCHIVED_BRANCH)
        self.git('branch', cleanup.INTEGRATION_BRANCH)
        (self.repo / 'fixture.txt').write_text('second\n')
        self.git('-c', 'commit.gpgsign=false', 'commit', '-am', 'second fixture')
        self.current = self.git('rev-parse', 'HEAD').stdout.strip()
        self.git('remote', 'add', 'origin', str(self.remote))
        self.git('push', '--atomic', 'origin', 'main', ARCHIVED_BRANCH, cleanup.INTEGRATION_BRANCH)
        write_manifest(self.repo, {'main': self.old, ARCHIVED_BRANCH: self.old,
                                   cleanup.INTEGRATION_BRANCH: self.old})

    def git(self, *args, cwd=None):
        return subprocess.run(['git', '-C', str(cwd or self.repo), *args],
                              capture_output=True, text=True, check=True)

    def invoke(self, **kwargs):
        # Only this test bypasses the GitHub URL allowlist for its temporary bare repo.
        with patch.object(cleanup, '_validate_origin', return_value={'fixture': 'local bare repository'}):
            return cleanup.prune_branches(self.repo, **kwargs)

    def test_real_git_preview_then_atomic_cleanup_leaves_main_untouched(self):
        before = self.git('ls-remote', '--heads', 'origin').stdout
        self.assertEqual(self.invoke()['status'], 'preview')
        self.assertEqual(self.git('ls-remote', '--heads', 'origin').stdout, before)
        self.assertEqual(self.invoke(apply=True)['status'], 'applied')
        self.assertEqual(self.git('ls-remote', '--heads', 'origin').stdout.strip(),
                         f'{self.current}\trefs/heads/main')
        self.assertEqual(self.invoke(apply=True)['status'], 'already_clean')

    def test_real_git_lease_race_rejects_all_deletions_without_retry(self):
        original = cleanup._git
        pushes = []

        def raced(repo, *args, **kwargs):
            if args[0] == 'push':
                pushes.append(args)
                self.git('update-ref', f'refs/heads/{ARCHIVED_BRANCH}', self.current, cwd=self.remote)
            return original(repo, *args, **kwargs)

        with patch.object(cleanup, '_git', side_effect=raced):
            with self.assertRaisesRegex(cleanup.CleanupError, 'push failed'):
                self.invoke(apply=True)
        self.assertEqual(len(pushes), 1)
        remaining = cleanup._remote_heads(self.repo)
        self.assertEqual(remaining, {'main': self.current, ARCHIVED_BRANCH: self.current,
                                     cleanup.INTEGRATION_BRANCH: self.old})


if __name__ == '__main__':
    unittest.main()
