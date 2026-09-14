# Next research gates

## P0 — Quantitative RNA data and measurement contract

Acquire the pinned OP3 file and verify SHA256. Inventory matrix layers and units, controls, plate/well/donor/replicate identity, dose/time and missing combinations. Decide whether raw counts or a processed pseudobulk matrix is present; do not guess from the filename. Preserve upstream notices. Deliver a machine-readable data audit and a protocol declaring which observations may be used. A readable HDF5 file is not sufficient to pass this gate.

## P1 — Lock task, splits and same-information baselines

Separate complete unseen-target baseline-only prediction from few-shot calibration. Use grouped biological units and shared-control accounting. Freeze preprocessing and genes using permitted training/baseline information only. Include zero effect, source mean, offset, ridge/low-rank and the v0.7 residual baseline where its information budget is allowed. Later published baselines require actual same-task execution, not borrowed scores. Primary evaluation includes total and centered error, mean bias, intervention discrimination and worst-context performance. Test permutations/near-neighbor exclusion without selecting methods on the final test.

## P2 — Test whether baseline information substitutes for a probe

Compare permitted target baseline RNA/ATAC information with matched-cost target probes under a separately declared few-shot protocol. The claim is incremental information, not merely more flexible fitting. Joint control/treatment/repeat design requires independently supported noise estimates or explicitly labeled sensitivity analysis. A new case cannot be presented as an independent study by changing random seeds.

## Engineering gate

Modernize the remaining historical runners to fail rather than overwrite a run, support explicit paths and record Git/data/protocol hashes. Each PR carries the precise hypothesis, tests performed, data scope, negative results and unresolved limitations. No automated tuning against the final test and no automatic claim of breakthrough.
