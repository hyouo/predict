"""Offline regression checks for verified asset publication and interrupted writes."""
import hashlib
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tools.assets import materialize


class InterruptedResponse(io.BytesIO):
    def geturl(self):
        return 'https://example.org/fixture'

    def read(self, size=-1):
        if self.tell():
            raise OSError('connection lost after first chunk')
        return super().read(4)


class AssetPublicationTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)
        self.payload = b'small offline asset fixture'
        self.digest = hashlib.sha256(self.payload).hexdigest()
        self.source = self.root / 'source'
        self.source.write_bytes(self.payload)
        self.destination = self.root / 'cache/asset'
        self.url = 'https://example.org/fixture'

    def copy(self):
        return materialize(self.url, self.destination, self.digest, 1024, self.source)

    def test_flush_and_publication_errors_leave_no_cached_or_partial_asset(self):
        for operation in ('fsync', 'link'):
            with self.subTest(operation=operation):
                with patch(f'tools.assets.os.{operation}', side_effect=OSError('storage failed')):
                    with self.assertRaisesRegex(OSError, 'storage failed'):
                        self.copy()
                self.assertEqual(list(self.destination.parent.iterdir()), [])
                self.assertEqual(self.source.read_bytes(), self.payload)

    def test_competing_publication_cannot_replace_existing_destination(self):
        original_link = os.link
        winner = b'previously published by another process'

        def competing_publish(source, destination):
            Path(destination).write_bytes(winner)
            return original_link(source, destination)

        with patch('tools.assets.os.link', side_effect=competing_publish):
            with self.assertRaises(FileExistsError):
                self.copy()
        self.assertEqual(self.destination.read_bytes(), winner)
        self.assertEqual(list(self.destination.parent.iterdir()), [self.destination])
        self.assertEqual(self.source.read_bytes(), self.payload)

    def test_partial_download_failure_cleans_cache_and_temporary_file(self):
        with patch('tools.assets.urlopen', return_value=InterruptedResponse(self.payload)):
            with self.assertRaisesRegex(OSError, 'connection lost after first chunk'):
                materialize(self.url, self.destination, self.digest, 1024)
        self.assertEqual(list(self.destination.parent.iterdir()), [])

    def test_interrupted_download_is_not_published(self):
        class CancelledResponse(InterruptedResponse):
            def read(self, size=-1):
                if self.tell():
                    raise KeyboardInterrupt('cancelled during download')
                return super().read(size)

        with patch('tools.assets.urlopen', return_value=CancelledResponse(self.payload)):
            with self.assertRaises(KeyboardInterrupt):
                materialize(self.url, self.destination, self.digest, 1024)
        self.assertEqual(list(self.destination.parent.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
