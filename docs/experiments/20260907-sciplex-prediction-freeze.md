# Prediction freeze before scoring

All six predeclared quantitative prediction sets have been generated and saved before invoking the separate scoring program: full published panel primary, source-control-expressed panel budget/seed sweep, pooled reference sensitivity, swapped reference sensitivity, pseudocount 0.1 sensitivity, and single-record source/target stress. Primary full panel: six transfers at 16 calibration drugs. Expressed panel: six transfers x four nested budgets x five fixed calibration orders. Four additional expressed-panel sensitivities: six transfers each. Seven same-information methods per transfer. Query drug labels are the fixed 74 labels, with declared record-availability restrictions only for the single-record stress.

Frozen estimator SHA256: 80ddbbc643f03844d589f89c365fb6da00024b38b5cf00d59a36cc590b43b176.
Frozen split SHA256: ddadaae32f834deda8a54193e3a70d73bde58f1129b9f9b0e01cdd54dacb08d5.
Models SHA256: dd082b006fb9d512aed00b74bc37f503a578b14aac965b6f327177bbb9876f39.
Prediction driver SHA256: 3d3e0ca813dd95effd51ddfb73e7d21127496a8ba01654635566e0ca3852d3da.
Scoring program SHA256: 6cc79fe0853a68bd34673f419bcbcf12769f9e1851ea21afcc1a4b139a9d78de.
Data adapter SHA256: 7ec864663939b4a058cfaf5e920681cd4ec413c45b1e66db8bae924a7c2fbff2 (previous adapter preserved separately for the primary full-panel run; extension only adds the predeclared single-record B-reference effect arrays).

Each local prediction manifest includes per-file SHA256, fit hyperparameters, calibration and query IDs, source-specific gene panel size, input SHA256, environment, start and freeze times. The initial timed-out primary attempt is retained and excluded; the complete rerun is `primary_all_complete`. All 23 new implementation/algebra/split checks passed before scoring. Tests are not biological replication.

The frozen estimator can change its effective regularization when zero-variance genes are appended, because its penalty scale uses the median source variance. This behavior has been documented by an algebraic test, not repaired or tuned against external query results. Full-panel and the already predeclared source-control-expressed-panel results will both be reported, including failures. No external prediction score has been used to select an estimator or panel.
