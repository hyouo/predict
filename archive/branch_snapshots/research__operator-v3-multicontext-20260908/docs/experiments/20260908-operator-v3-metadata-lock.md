# LINCS operator-v3 metadata lock

Upstream file: 333,680 profiles x 978 landmark genes, float32, verified SHA256 34d198df9eddac5b535a0408794caaf41073fa014c74429d1050e483576fb581. The 24h/(10uM or control) subset contains 65,679 profiles, 30 metadata context labels, 1,726 noncontrol compounds and 18,182 DMSO controls. All control labels map to DMSO (CID 679); compound labels each map to exactly one name; none of 978 gene columns is marked merged. psbulk_cells/counts are -666 sentinels in this L1000 assay and are NOT count offsets.

Acquisition run 34178506952 exceeded its ten-minute limit without producing a result. The bounded byte-range retry 34179282079 succeeded, reassembled the exact upstream hash and exported the metadata-defined subset. Artifact 10038337172 archive SHA256 78a8fecda2bd19bfb427f3129a969b8d837f27a82fdac0d005dd632a24190219; selected expression archive SHA256 6823d883a781f51acac0beb0de7af30f48dd8a19eff45fc83cad7db9f5bc4507. No target-response summaries or model fits were used for selection.

Applying the previously registered rules gives 26 eligible context labels, 20 development and 6 held out. No threshold revision was needed. Split file SHA256: 96cefff9b45ff650d4885fb10a464f49133602693e4d7110889362e610af2429.

Development labels: ASC.C, CD34, CVCL_0023, CVCL_0030, CVCL_0031, CVCL_0033, CVCL_0062, CVCL_0178, CVCL_0320, CVCL_0332, CVCL_0598, CVCL_1794, CVCL_2959, CVCL_3383, CVCL_U602, CVCL_VU89, NEU, NPC, NPC.CAS9, NPC.TAK.

Final labels and eligible queries: SKL (162), SKL.C (165), CVCL_0027/HepG2 (267), CVCL_5136/HCC515 (268), CVCL_0035/PC3 (1710), CVCL_0132/A375 (1702). These sum to 4,274 context-compound pairs over 1,724 source-covered compound IDs; they are not independent biological experiments. Related SKL/SKL.C labels and NPC sublabels are retained as distinct metadata contexts, not claimed to be independent cell lineages.

Target A-control counts in that order: 104,106,186,184,1004,1023. Target B scoring references: 104,106,184,179,991,1010. Input/control information is explicitly nonzero even though target treatment calibration is zero. Native Level3 normalized values are differenced against within-cell/plate controls, with no CPM or additional log transformation.

Before any development fitting, add an important strong comparator: per-compound baseline-kernel regression with free intercept (linear or RBF, ridge grid as registered, baseline-only scaling, target excluded), with global amplitude selected by the same development-context CV. This is classical kernel regression, not the proposed innovation. It helps distinguish learned response-by-state interaction from straightforward interpolation among available source contexts.

All iterative fitting remains restricted to the twenty development contexts. Final target treatment reads stay blocked until the complete prediction-file manifest is frozen. Upstream quantile normalization, assay noise, selection, potential campaign effects and related context labels limit independence claims.
