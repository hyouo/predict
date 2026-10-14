"""Regression checks for untrusted arrays and exclusive artifact publication."""
import hashlib
import io as stdio
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings
import zipfile

import numpy as np

from perturb_predict import cli, io, op3
from perturb_predict.demo import run_demo


def npy_bytes(array):
    stream = stdio.BytesIO()
    np.lib.format.write_array(stream, np.asarray(array), allow_pickle=True)
    return stream.getvalue()


def declared_array(shape, *, dtype='<f8', payload=b''):
    stream = stdio.BytesIO()
    np.lib.format.write_array_header_1_0(
        stream, {'descr': dtype, 'fortran_order': False, 'shape': shape})
    return stream.getvalue() + payload


class ArchiveSafetyTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)

    def archive(self, name, entries, *, metadata='{"schema_version": 1}'):
        path = self.root / name
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                archive.writestr('metadata_json.npy', npy_bytes(np.array(metadata)))
                for key, payload in entries:
                    archive.writestr(key, payload)
        return path

    def assert_rejected_before_array_read(self, path):
        with patch.object(np.lib.format, 'read_array') as read:
            with self.assertRaises(ValueError):
                io.load_archive(path)
            read.assert_not_called()

    def test_partial_npz_serialization_does_not_publish_or_leave_temporary_file(self):
        def fail(stream, **arrays):
            stream.write(b'partial archive bytes')
            raise OSError('compression write failure')

        with patch.object(np, 'savez_compressed', side_effect=fail):
            with self.assertRaisesRegex(OSError, 'compression write failure'):
                io.save_archive(self.root / 'data.npz', {'schema_version': 1}, x=np.ones(2))
        self.assertEqual(list(self.root.iterdir()), [])

    def test_npz_flush_and_publication_failures_leave_no_final_or_temporary_file(self):
        for operation in ('fsync', 'link'):
            with self.subTest(operation=operation):
                with patch.object(io.os, operation, side_effect=OSError('storage failure')):
                    with self.assertRaisesRegex(OSError, 'storage failure'):
                        io.save_archive(self.root / 'data.npz', {'schema_version': 1}, x=np.ones(2))
                self.assertEqual(list(self.root.iterdir()), [])

    def test_npz_publication_race_preserves_existing_bytes(self):
        path = self.root / 'data.npz'
        original_link = os.link
        winner = b'created by the other writer'

        def competing_publish(source, destination):
            Path(destination).write_bytes(winner)
            return original_link(source, destination)

        with patch.object(io.os, 'link', side_effect=competing_publish):
            with self.assertRaises(FileExistsError):
                io.save_archive(path, {'schema_version': 1}, x=np.ones(2))
        self.assertEqual(path.read_bytes(), winner)
        self.assertEqual(list(self.root.iterdir()), [path])

    def test_object_fields_are_rejected_on_both_write_and_read(self):
        dtypes = [object, [('value', object)], [('nested', [('value', object)])]]
        for index, dtype in enumerate(dtypes):
            with self.subTest(dtype=str(dtype)):
                array = np.zeros(1, dtype=dtype)
                rejected = self.root / f'writer-{index}.npz'
                with self.assertRaises(ValueError):
                    io.save_archive(rejected, {'schema_version': 1}, values=array)
                self.assertFalse(rejected.exists())
                crafted = self.archive(f'object-{index}.npz', [('values.npy', npy_bytes(array))])
                self.assert_rejected_before_array_read(crafted)

    def test_forged_shapes_and_inexact_payloads_are_rejected_before_allocation(self):
        entries = {
            'huge_product': declared_array((2**40, 2**40), payload=b'\0' * 8),
            'large_zero_width': declared_array((io.MAX_BYTES + 1,), dtype='<U0'),
            'negative_dimension': declared_array((-1,)),
            'truncated_payload': declared_array((2,), payload=b'\0' * 8),
            'trailing_payload': declared_array((1,), payload=b'\0' * 16),
        }
        for name, payload in entries.items():
            with self.subTest(name=name):
                path = self.archive(f'{name}.npz', [('values.npy', payload)])
                self.assertLess(path.stat().st_size, 1024)
                self.assert_rejected_before_array_read(path)

    def test_oversized_or_unsupported_npy_headers_are_rejected_before_array_read(self):
        header = b"{'descr': '<f8', 'fortran_order': False, 'shape': (1,), }"
        oversized = header + b' ' * (10001 - len(header) - 1) + b'\n'
        payloads = [np.lib.format.magic(1, 0) + struct.pack('<H', len(oversized))
                    + oversized + b'\0' * 8,
                    np.lib.format.magic(3, 0) + struct.pack('<I', len(header))
                    + header + b'\0' * 8]
        for index, payload in enumerate(payloads):
            with self.subTest(index=index):
                path = self.archive(f'header-{index}.npz', [('values.npy', payload)])
                self.assert_rejected_before_array_read(path)

    def test_invalid_member_names_and_duplicates_are_rejected_before_array_read(self):
        payload = npy_bytes(np.array([1.]))
        cases = [[('../values.npy', payload)], [('values.npy/extra.npy', payload)],
                 [('values.txt', payload)], [('values.npy', payload), ('values.npy', payload)]]
        for index, entries in enumerate(cases):
            with self.subTest(index=index):
                self.assert_rejected_before_array_read(self.archive(f'name-{index}.npz', entries))

    def test_file_expanded_size_and_member_count_limits_precede_array_reads(self):
        path = self.archive('limits.npz', [('values.npy', npy_bytes(np.zeros(4096, dtype=np.uint8)))])
        compressed_size = path.stat().st_size
        with patch.object(io, 'MAX_BYTES', compressed_size - 1):
            self.assert_rejected_before_array_read(path)
        with zipfile.ZipFile(path) as archive:
            expanded_size = sum(member.file_size for member in archive.infolist())
        self.assertLess(compressed_size + 1, expanded_size)
        with patch.object(io, 'MAX_BYTES', compressed_size + 1):
            self.assert_rejected_before_array_read(path)
        too_many = self.archive('members.npz', [
            (f'value_{i}.npy', npy_bytes(np.array([i]))) for i in range(32)])
        self.assert_rejected_before_array_read(too_many)

    def test_nonfinite_metadata_and_noninteger_schema_are_rejected(self):
        cases = ['{"schema_version": 1, "value": NaN}',
                 '{"schema_version": 1, "value": Infinity}',
                 '{"schema_version": true}', '{"schema_version": 1.0}']
        for index, metadata in enumerate(cases):
            with self.subTest(metadata=metadata):
                path = self.archive(f'metadata-{index}.npz', [], metadata=metadata)
                with self.assertRaises(ValueError):
                    io.load_archive(path)

    def test_corrupted_zip_cli_returns_json_exit_two_without_traceback(self):
        non_zip = self.root / 'not-zip.npz'
        non_zip.write_bytes(b'not a zip file')
        truncated = self.archive('truncated.npz', [])
        truncated.write_bytes(truncated.read_bytes()[:-22])
        short_npy = self.archive('short-npy.npz', [('values.npy', b'\x93NUMPY\x01\x00\x7f\x00')])
        for path in (non_zip, truncated, short_npy):
            with self.subTest(path=path.name):
                destination = self.root / path.stem
                completed = subprocess.run(
                    [sys.executable, '-m', 'perturb_predict', 'train', '--data', str(path),
                     '--out', str(destination)],
                    capture_output=True, text=True, check=False)
                self.assertEqual(completed.returncode, 2, completed.stderr)
                self.assertEqual(completed.stdout, '')
                self.assertNotIn('Traceback', completed.stderr)
                error = json.loads(completed.stderr)
                self.assertEqual(error['status'], 'error')
                self.assertEqual(error['type'], 'ValueError')
                self.assertFalse(destination.exists())

    def test_failed_and_complete_prediction_cannot_open_truth_for_scoring(self):
        out = self.root / 'demo'
        run_demo(out)
        io.json_write(out / 'prediction/FAILED.json', {'status': 'failed'})
        with patch.object(cli, 'load_bundle') as truth:
            with self.assertRaisesRegex(ValueError, 'incomplete or failed'):
                cli.evaluate(out / 'prediction/predictions.npz', out / 'inputs/truth.npz',
                             self.root / 'evaluation')
            truth.assert_not_called()
        self.assertFalse((self.root / 'evaluation').exists())

    def test_keyboard_interrupt_marks_run_failed(self):
        out = self.root / 'interrupted'
        with self.assertRaises(KeyboardInterrupt):
            with io.run_directory(out):
                raise KeyboardInterrupt('user cancelled')
        failure = json.loads((out / 'FAILED.json').read_text())
        self.assertEqual(failure['status'], 'failed')
        self.assertEqual(failure['error_type'], 'KeyboardInterrupt')
        self.assertFalse((out / 'COMPLETE.json').exists())

    def test_failure_marker_write_error_preserves_operation_exception(self):
        for index, original in enumerate((ValueError('original cause'),
                                          KeyboardInterrupt('original interruption'))):
            with self.subTest(error=type(original).__name__):
                out = self.root / f'failed-{index}'
                with patch.object(io, 'json_write', side_effect=OSError('disk unavailable')):
                    with self.assertRaises(type(original)) as caught:
                        with io.run_directory(out):
                            raise original
                self.assertIs(caught.exception, original)
                self.assertIn('disk unavailable', '\n'.join(original.__notes__))
                self.assertFalse((out / 'COMPLETE.json').exists())


class Op3DownloadPublicationTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)
        self.payload = b'offline download fixture; not an OP3 dataset'
        self.digest = hashlib.sha256(self.payload).hexdigest()

    def response(self):
        response = stdio.BytesIO(self.payload)
        response.geturl = lambda: 'https://example.org/fixture'
        return response

    def test_successful_download_is_published_once_and_then_cached(self):
        path = self.root / 'op3.h5ad'
        with patch.object(op3, 'DATA_SHA256', self.digest), \
             patch.object(op3, 'urlopen', return_value=self.response()) as download:
            result = op3.fetch(path)
            cached = op3.fetch(path)
        self.assertFalse(result['cached'])
        self.assertTrue(cached['cached'])
        self.assertEqual(path.read_bytes(), self.payload)
        download.assert_called_once()
        self.assertEqual(list(self.root.iterdir()), [path])

    def test_download_publication_failure_cleans_temporary_file(self):
        path = self.root / 'op3.h5ad'
        with patch.object(op3, 'DATA_SHA256', self.digest), \
             patch.object(op3, 'urlopen', return_value=self.response()), \
             patch.object(io.os, 'link', side_effect=OSError('publication failed')):
            with self.assertRaisesRegex(OSError, 'publication failed'):
                op3.fetch(path)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_download_race_preserves_file_created_during_transfer(self):
        path = self.root / 'op3.h5ad'
        original_link = os.link
        winner = b'existing concurrent download'

        def competing_publish(source, destination):
            Path(destination).write_bytes(winner)
            return original_link(source, destination)

        with patch.object(op3, 'DATA_SHA256', self.digest), \
             patch.object(op3, 'urlopen', return_value=self.response()), \
             patch.object(io.os, 'link', side_effect=competing_publish):
            with self.assertRaises(FileExistsError):
                op3.fetch(path)
        self.assertEqual(path.read_bytes(), winner)
        self.assertEqual(list(self.root.iterdir()), [path])

    def test_download_size_limit_leaves_no_cache_or_temporary_file(self):
        with patch.object(op3, 'DATA_SHA256', self.digest), \
             patch.object(op3, 'MAX_DATA_BYTES', len(self.payload) - 1), \
             patch.object(op3, 'urlopen', return_value=self.response()):
            with self.assertRaisesRegex(ValueError, 'exceeds'):
                op3.fetch(self.root / 'op3.h5ad')
        self.assertEqual(list(self.root.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
