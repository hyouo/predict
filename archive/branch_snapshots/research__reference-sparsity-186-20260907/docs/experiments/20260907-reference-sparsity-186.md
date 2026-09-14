# Control-reference and sparse-panel reanalysis (186 labels)

Date: 2026-09-07. This is a reviewable evidence checkpoint, not a software release, journal submission, or claim that the full computational archive has been uploaded to GitHub.

## Provenance and separation from other analyses

Repository baseline: `575c325a4dd6f99eb055a305e789932667835b93`. Public asset acquisition: `652457e1c45cc109db01153109f88177d63a1d18`, Actions run 34082109943. SciPlex file SHA256: `8ad35a20e44bf41cb030a0f063eb6a5727965b8f27d0eda4380bda0e294ccc4c`.

This cohort uses 186 common literal compound labels, 75 fixed query labels, all 56,799 published genes, six directed K562/A549/MCF7 transfers, 16 target calibration compounds, 10uM and 24h. Target fits use reference A and primary scores reference B; source responses pool A/B. The same conditions, model grids and information budgets are used across comparators.

Final synchronization identified prior quantitative SciPlex work at c2dd42122035f31050499020e09e70cc2ca6c9dd and a distinct 157-compound reanalysis in #9. Therefore this is a **quantitative reanalysis of a study external to OP3, not a previously untouched blind validation**. Cohorts and numerical scores must not be silently combined. The local freeze order does not erase earlier project exposure. No estimates, split or model were changed after this provenance correction.

## Primary results (macro MSE, lower is better)

| Method | MSE |
|---|---:|
| Zero effect | 0.900489 |
| Pooled scalar calibration | 0.958991 |
| Frozen partial pooling | 1.013692 |
| Zero-prior gene ridge | 1.014563 |
| Target calibration mean | 1.022538 |
| Multivariate ridge | 1.168610 |
| Source copying | 1.557023 |

Partial pooling is 34.90% better than source copy but 12.57% worse than zero, and loses to zero in every directed transfer. It also loses to zero in all 20 nested calibration-budget/order settings. This is not a SOTA result.

## Statistical diagnostics

- Symmetrized multivariate-ridge MSE changes from 0.647532 with reused target references to 1.180888 with disjoint references; the symmetric zero baseline is 0.922576. Swapping only the scoring reference satisfies an exact fixed-prediction relative-risk identity, verified on all 42 primary arrays. Reference differences cannot automatically be classified as purely technical noise.
- All-gene median calibration variance is zero for K562 and A549. Appending all-zero feature columns to real OP3 inputs changes original-gene predictions at fixed hyperparameters (B-cell RMS difference 0.187018, maximum 6.312158 with twice as many zero columns). Positive-variance-only scaling fixes this particular invariance property.
- Positive scaling plus reference-connected cross-validation lowers 10uM error to 0.898820, only 0.185% below zero. At 1uM it is 0.820483 versus zero 0.821093, only 0.074% lower. These predictions are nearly null, not a recovery of strong treatment-specific mechanisms. The 10uM corrections are post hoc; 1uM is locally predeclared but neither a new study nor wholly unseen in project history.
- OP3 reused-versus-disjoint reference analysis is performed on its already-consumed public data and shows weaker, dataset-dependent reference sensitivity. No private OP3 expression was read.

## Deliverables and checks

A complete 17-page working manuscript in the official Springer Nature sn-jnl / sn-nature template, six vector figures, figure source tables, executed analysis code, fixed protocols, original prediction hashes and an additive repository patch are provided as conversation artifacts. Large primary prediction arrays are separately archived; they are not in this checkpoint commit.

Twenty new numerical/protocol tests pass. All 42 primary prediction matrices were replayed exactly. A fresh primary pipeline from the original H5AD, using uncompressed archive containers, also completed with exact prediction and score-table agreement. Twelve repaired predictions were independently recomputed against frozen hashes. The complete multi-sensitivity wrapper has not been rerun as one fresh monolithic acceptance job after packaging; its constituent analyses were actually executed.

Limits: shared controls/compound identities, unresolved biological-replicate semantics, inherited gene panel, no independently harmonized compound structures, no modern foundation-model comparison, no ATAC, no human-approved authorship/declarations. Compound resampling is conditional sensitivity, not independent biological inference. Human verification and genuinely independent confirmation remain required. Do not merge or release algorithmic superiority claims from this checkpoint.
