import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from perturb_predict import io

class AtomicOutputTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory()
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)

    def test_serialization_failure_creates_no_partial_json(self):
        for value in ({"x": float("nan")}, {"x": object()}):
            path = self.root / 'output.json'
            with self.assertRaises((TypeError, ValueError)):
                io.json_write(path, value)
            self.assertFalse(path.exists())
        self.assertEqual(list(self.root.iterdir()), [])

    def test_existing_output_preserved(self):
        path = self.root / 'output.json'
        path.write_text('original')
        with self.assertRaises(FileExistsError):
            io.json_write(path, {"replacement": True})
        self.assertEqual(path.read_text(), 'original')
        self.assertEqual(list(self.root.glob('.json-*')), [])

    def test_publish_error_does_not_leave_partial_file(self):
        path = self.root / 'output.json'
        with patch.object(io.os, 'link', side_effect=OSError('disk failure')):
            with self.assertRaises(OSError):
                io.json_write(path, {"x": 1})
        self.assertFalse(path.exists())
        self.assertEqual(list(self.root.iterdir()), [])

    def test_complete_run_roundtrip(self):
        with io.run_directory(self.root / 'success') as directory:
            io.json_write(directory / 'result.json', {"value": "中文"})
        self.assertEqual(json.loads((directory / 'COMPLETE.json').read_text())['status'], 'complete')
        self.assertFalse((directory / 'FAILED.json').exists())

    def test_failed_completion_is_marked_failed(self):
        original = io.json_write
        def fail_completion(path, value):
            if Path(path).name == 'COMPLETE.json':
                raise OSError('completion storage failure')
            return original(path, value)
        with patch.object(io, 'json_write', side_effect=fail_completion):
            with self.assertRaisesRegex(OSError, 'completion storage failure'):
                with io.run_directory(self.root / 'failure'):
                    pass
        directory = self.root / 'failure'
        self.assertFalse((directory / 'COMPLETE.json').exists())
        self.assertEqual(json.loads((directory / 'FAILED.json').read_text())['status'], 'failed')

if __name__ == '__main__':
    unittest.main()
