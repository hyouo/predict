# Metadata-only lock, 2026-09-08

Acquisition succeeded in run 34179374331 after the first slower CSR-reader run hit its 12-minute limit. Original HDF5 707,799,618 bytes, SHA256 34d198df9eddac5b535a0408794caaf41073fa014c74429d1050e483576fb581, shape 333680 x 978. Fixed 24h / 10uM-plus-controls slice: 65,679 observations, including 18,182 controls; exported without changing expression values. 30 metadata contexts meet the prespecified eligibility; 1,726 treatment labels. No expression or target effect filtering was used for this lock.

Train (18): CVCL_0598, CVCL_0178, CVCL_0027, CVCL_0332, CVCL_2959, CVCL_0031, CVCL_0132, NPC, CVCL_1794, CD34, SKL, CVCL_3383, CVCL_VU89, NEU, CVCL_0065, CVCL_0062, CVCL_0395, CVCL_0035.
Validation (6): NPC.TAK, CVCL_0030, CVCL_B161, CVCL_0033, ASC.C, CVCL_0320.
Test (6): CVCL_0023, CVCL_5136, CVCL_U602, NPC.CAS9, MNEU.E, SKL.C.

All 978 landmark features are retained. Compound IDs must be decoded as strings because the supplied perturbagen field contains numeric and nonnumeric identifiers. There are no duplicate (cell_type, plate, well) observations in the slice. The 860 context/plate control groups each contain 8–28 records, so all pass the >=4-control rule. The file's psbulk_cells and psbulk_counts are sentinel -666, not counts or evidence of biological replication. It is normalized L1000 expression, not single-cell RNA-seq.

Context holdout is not lineage-family holdout: names such as NPC/NPC.TAK/NPC.CAS9 and SKL/SKL.C can represent related backgrounds. The primary experiment holds exact metadata context IDs, not necessarily completely novel ancestral cell lines, tissue families, donors or assay batches. Raw source metadata and the exact per-context drug coverage will be retained. No claim of biological independence for nominal repeats.

Test treatment and test scoring-control values remain inaccessible to scientific fit/selection. All model families will be selected on the six validation contexts, then refit on the 24 train+validation contexts and jointly frozen before any test scoring. The previously registered same-source paired-target contrast evaluation remains mandatory. No test-set model iteration will be presented as confirmation.
