import unittest
import numpy as np
from tools.assets import ROOT
from tools.sciplex_experiment import Store,CONTEXTS

class ExternalAccessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (ROOT/'data/rna/srivatsan20_sciplex3_processed.h5ad').is_file():raise unittest.SkipTest('External input not installed')
        cls.s=Store(ROOT/'data/rna/srivatsan20_sciplex3_processed.h5ad',CONTEXTS[0])
    def test_no_chemical_overlap(self):
        self.assertFalse(set(self.s.anchors)&set(self.s.queries));self.assertEqual(len(self.s.anchors),16)
    def test_target_panel_b_unavailable_to_fit(self):
        s=self.s;o=s.o;rows=np.flatnonzero(o.qualifying&o.cell_type.eq(s.target)&o.panel.eq('B')&o.cid.isin(s.queries))
        for role in ('source_treatment','anchor_treatment','target_baseline'):
            with self.assertRaises(PermissionError):s.read(rows,role)
    def test_evaluation_controls_are_not_baseline(self):
        s=self.s;o=s.o;rows=np.flatnonzero(o.is_control&o.cell_type.eq(s.target)&o.panel.eq('B'))
        with self.assertRaises(PermissionError):s.read(rows,'target_baseline')
    def test_single_cid_for_duplicate_drug_alias(self):
        o=self.s.o
        self.assertEqual(o.loc[o.perturbagen.eq('ENMD-2076'),'cid'].nunique(),1)
    def test_source_panel_b_disallowed(self):
        s=self.s;o=s.o;rows=np.flatnonzero(o.qualifying&o.cell_type.isin(s.sources)&o.panel.eq('B'))
        with self.assertRaises(PermissionError):s.read(rows,'source_treatment')
    def test_anchor_hash_order_fixed(self):
        import hashlib
        rank=sorted(self.s.cids,key=lambda c:(hashlib.sha256(('predict-sciplex-v08|'+c).encode()).hexdigest(),c))
        self.assertEqual(self.s.anchors,rank[:16])

if __name__=='__main__':unittest.main()
