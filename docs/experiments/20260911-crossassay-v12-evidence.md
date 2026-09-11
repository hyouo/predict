# v12 completed evidence checkpoint (overall research incomplete)

This continuation preserves all prior experiments and the interrupted Stage A/B exposure history. No main/default model changed, no wet-lab experiment or modern-method superiority claim. Full executable source and lossless evidence are supplied as verified conversation artifacts; this text does not claim those payloads are all uploaded to Git.

## Data recovery and a new provenance finding
Fifteen original GSE101406 compressed assets were recovered;14 data/metadata files pass publisher SHA512. Actual Level3 matrices: P1001684x96, RNA1667x12328 (only978 directly measured genes used), GCP1712x59. Cross-assay culture/aliquot/preparation identity is unknown despite the original paper reporting biological replicates. det_plate must not be relabelled as culture. Exact dose differences and NPC/NPC.TAK remain explicit.

The bounded pinned Phase-II LPROT6h extraction (run34575308387) completed. All1,102 shared records match exactly at978 measured genes:1,077,756 values, maximum absolute difference0, no fitted transformation.16 same-detection-plate vehicle records absent from the GSE101406 release can be traced in the fixed Phase-II source. The explicit reference registry retains20 original controls and adds only those16, recovering all1,647 RNA treatment references; P100/RNA condition coverage rises301->533 and six-base-cell common conditions30->88. These are shared/recovered measurements, not new independent biology. GCP56 records and cultivation links remain unresolved.

## Actual bounded model comparison
StageB: six base cells x30 exact-dose conditions. Source-only calibrated MSE0.527054367262; magnitude-only early-P1000.526698; full-P100 residual0.527137209915. StageA historical candidate figures are retained but not asserted to be bitwise-reproduced implementations.

Separate pre-registered StageC: five exact-native cancer cells (noNPC/NPC.TAK),88 conditions,978 measured genes. Source-only calibrated MSE0.470554174787; magnitude-only0.470036675555; full-P1000.470721715448. Full P100 does not outperform source RNA. Same-query targetP1003h is explicitly used to predict RNA6h; this is neither query-excluded functional probing nor independent-culture validation. GCP24h is not an input. All32 stage-level predictor arrays include zero/source controls and ten fixed within-processing-block permutations per stage, not32 distinct state-of-the-art methods. No cross-cohort improvement is claimed.

## Verification
27 tests passed. Chosen-parameter replays:96 StageB blocks and80 StageC blocks exactly reproduce. Direct raw-value reconstruction reproduces180 and440 RNA means respectively. The final three archives were independently extracted:106 research/evidence payload files and21 public-data payload files matched SHA256; all27 checks passed again. An independent scorer reproduces176 context/model rows and all reported metric components (maximum difference9.72e-17). All14 publisher hashes and three Phase-II extraction hashes checked again. Full hyperparameter search was not re-run twice; no new biological replication.

Preserved failures include inefficient scoring timeout, object-typed identifier serialization (rebuilt from raw IDs with pickle disabled), and an overstrict overlap-invariance check when newly referenced positive-control observations legitimately changed8 means. No failed-output predictions were used.

## Remaining gate
Resolve original culture/preparation/aliquot links, TableS1 dose discrepancies, NPC provenance and G-0016R meaning before claiming cross-preparation functional state. External phosphorylation signatures also require provenance screening: PTMsigDB includes P100-derived sets, so public knowledge is not automatically an independent prior. Current models did not use PTMsigDB. Do not repeat downloads or unrestricted searches on the now-exposed conditional-mean cohorts.
