import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np
from tools.op3_inventory import read_column


class MetadataColumnTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.file = h5py.File(Path(self.temp.name) / 'synthetic-format-fixture.h5', 'w')
        self.addCleanup(self.file.close)

    def test_dense_text_and_numeric(self):
        text = self.file.create_dataset('text', data=['a', 'b'], dtype=h5py.string_dtype())
        values = self.file.create_dataset('values', data=np.array([1., 2.]))
        self.assertEqual(read_column(text).tolist(), ['a', 'b'])
        np.testing.assert_array_equal(read_column(values), [1., 2.])

    def test_categories_preserve_boolean_and_missing(self):
        group = self.file.create_group('c')
        group.create_dataset('categories', data=[False, True])
        group.create_dataset('codes', data=np.array([1, -1, 0]))
        self.assertEqual(read_column(group).tolist(), [True, None, False])

    def test_empty_categories_allow_only_missing(self):
        group = self.file.create_group('c')
        group.create_dataset('categories', data=np.array([], dtype=float))
        group.create_dataset('codes', data=np.array([-1, -1]))
        self.assertEqual(read_column(group).tolist(), [None, None])

    def test_invalid_categories_rejected(self):
        group = self.file.create_group('c')
        group.create_dataset('categories', data=[1, 2])
        for values in [[-2], [2], [0.5]]:
            if 'codes' in group:
                del group['codes']
            group.create_dataset('codes', data=values)
            with self.assertRaises(ValueError):
                read_column(group)

    def test_unknown_encoding_or_matrix_rejected(self):
        with self.assertRaises(ValueError):
            read_column(self.file.create_group('unknown'))
        with self.assertRaises(ValueError):
            read_column(self.file.create_dataset('matrix', data=np.ones((2, 2))))


if __name__ == '__main__':
    unittest.main()
