import unittest,json
from pathlib import Path
from unittest.mock import patch
import numpy as np
from tools.assets import ROOT
from src.rna_data import OP3Store

class RNAAccessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (ROOT/'data/rna/op3_standardized_processed.h5ad').is_file():
            raise unittest.SkipTest('Pinned RNA input is an integration-test dependency, not simulated')
    def store(self):
        return OP3Store(ROOT/'data/rna/op3_standardized_processed.h5ad',json.loads((ROOT/'reports/data/op3_donors.json').read_text()))
    def test_target_hidden_expression_rejected_for_training(self):
        s=self.store();o=s.obs
        rows=np.flatnonzero((o.cell_type=='CL_0000236')&o['split'].eq('private_test'))
        for role in ('source_treatment','anchor_treatment','baseline'):
            with self.assertRaises(PermissionError):s._read(rows,role)
    def test_evaluation_control_cannot_be_training_control(self):
        s=self.store();o=s.obs
        rows=np.flatnonzero((o.cell_type=='CL_0000236')&o.is_control&o.well_row.eq('E'))
        with self.assertRaises(PermissionError):s._read(rows,'anchor_control')
        with self.assertRaises(PermissionError):s._read(rows,'baseline')
        with self.assertRaises(PermissionError):s._read(rows,'evaluation_control')
    def test_count_scale_nonnegative_finite(self):
        s=self.store();o=s.obs;rows=np.flatnonzero((o.cell_type=='CL_0000236')&o.is_control&o.well_row.eq('A'))
        x=s._read(rows,'baseline');self.assertTrue(np.isfinite(x).all());self.assertTrue((x>=0).all())
    def test_unverified_donor_mapping_rejected(self):
        with self.assertRaises(ValueError):OP3Store(ROOT/'data/rna/op3_standardized_processed.h5ad',{})
    def test_roles_validate_before_using_cache(self):
        s=self.store();o=s.obs;rows=np.flatnonzero((o.cell_type=='CL_0000236')&o.is_control&o.well_row.eq('A'))
        s._read(rows,'baseline')
        with self.assertRaises(PermissionError):s._read(rows,'anchor_control')
    def test_target_query_unreachable_during_source_and_anchor(self):
        s=self.store();original=s._read
        def guarded(rows,role,**kw):
            if role.startswith('evaluation'):raise AssertionError('Test outcomes were requested')
            return original(rows,role,**kw)
        with patch.object(s,'_read',side_effect=guarded):
            so=s.source(('Donor 1',));a=s.anchors('CL_0000236',so)
            self.assertEqual(a.mean.shape[1],5288)
        self.assertTrue(all(not e['role'].startswith('evaluation') for e in s.access))

if __name__=='__main__':unittest.main()
