import hashlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import h5py
import numpy as np
from tools.assets import materialize, inspect_h5ad
from tools.smoke import run


class RepositoryToolsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = b'verified test fixture, not biological data\n'
        self.digest = hashlib.sha256(self.data).hexdigest()
        self.url = 'https://example.org/public-fixture'

    def test_offline_verified_copy(self):
        source = self.root / 'input'
        source.write_bytes(self.data)
        output = self.root / 'cache/output'
        result = materialize(self.url, output, self.digest, 1024, source)
        self.assertEqual(output.read_bytes(), self.data)
        self.assertFalse(result['cached'])
        self.assertTrue(materialize(self.url, output, self.digest, 1024)['cached'])

    def test_wrong_cached_file_never_overwritten(self):
        output = self.root / 'output'
        output.write_bytes(b'wrong')
        with self.assertRaises(ValueError):
            materialize(self.url, output, self.digest, 1024)
        self.assertEqual(output.read_bytes(), b'wrong')

    def test_wrong_download_removed(self):
        output = self.root / 'output'
        with patch('tools.assets.urlopen', return_value=io.BytesIO(b'wrong')):
            with self.assertRaises(ValueError):
                materialize(self.url, output, self.digest, 1024)
        self.assertFalse(output.exists())
        self.assertFalse(list(self.root.glob('*.partial')))

    def test_transfer_and_cache_size_caps(self):
        output = self.root / 'output'
        with patch('tools.assets.urlopen', return_value=io.BytesIO(self.data)):
            with self.assertRaises(ValueError):
                materialize(self.url, output, self.digest, 4)
        output.write_bytes(self.data)
        with self.assertRaises(ValueError):
            materialize(self.url, output, self.digest, 4)

    def test_non_https_and_bad_hash_rejected(self):
        for url, digest, cap in [('http://example.org/x', self.digest, 1024),
                                  (self.url, 'unverified', 1024),
                                  (self.url, self.digest, -1),
                                  (self.url, self.digest, True)]:
            with self.assertRaises(ValueError):
                materialize(url, self.root / 'output', digest, cap)

    def test_network_failure_cleans_temporary_file(self):
        with patch('tools.assets.urlopen', side_effect=OSError('offline test')):
            with self.assertRaises(OSError):
                materialize(self.url, self.root / 'output', self.digest, 1024)
        self.assertFalse(list(self.root.glob('*.partial')))
        self.assertFalse((self.root / 'output').exists())

    def test_h5ad_structure_does_not_approve_training(self):
        path = self.root / 'fixture.h5ad'
        with h5py.File(path, 'w') as f:
            f.create_dataset('X', data=np.ones((3, 2)))
            f.create_group('obs').create_dataset('_index', data=np.arange(3))
            f.create_group('var').create_dataset('_index', data=np.arange(2))
        report = inspect_h5ad(path)
        self.assertEqual(report['X']['shape'], [3, 2])
        self.assertFalse(report['approved_for_training'])

    def test_run_refuses_existing_results(self):
        output = self.root / 'existing-run'
        output.mkdir()
        (output / 'keep.txt').write_text('old result')
        with self.assertRaises(FileExistsError):
            run(output)
        self.assertEqual((output / 'keep.txt').read_text(), 'old result')


if __name__ == '__main__':
    unittest.main()
