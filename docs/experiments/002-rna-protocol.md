# 002 — First quantitative RNA experiments: prospective analysis protocol

Date: 2026-09-06. Base repository commit: a193fc9aa5bc014c63738fb75b7e0d278d44fdd0.
This protocol is registered before inspecting held-out OP3 expression outcomes. It describes retrospective evaluation on a fixed public processed asset, NOT a prospective wet-lab experiment or a reproduction of the competition leaderboard.

## Questions
1. Does gene-specific response shrinkage outperform a shared drug-response kernel on quantitative RNA at the same target-label budget?
2. Does untreated target baseline information improve zero-shot transfer beyond source response shrinkage?
3. Does any improvement survive independent control pools, drug-specific (centered) metrics and leave-donor-out evaluation?

## Measurement gate
Use pinned OP3 SHA256 8eaf7e63adc68029e88184726b8545bc3836c3ef2a25b9185dc4e8b819537e59. Verify nonnegative integer count storage, unique sample/plate/well/context keys, denominator psbulk_counts, upstream count aggregation and gene filtering. Validate donor mapping against original upstream observation metadata where obtainable; otherwise report only a documented demographic reconstruction, not an independent donor verification. Incomplete fields remain explicitly missing.
The released 5,288-gene universe was filtered upstream across cell types; conditional-on-released-universe claims only. No new gene selection or scaling may inspect test expression. Raw counts are never interpreted as already normalized effects. CPM normalization must specify the denominator and pseudocount; use log2(1+10^6 * counts / psbulk_counts), with separately reported retained-gene-sum sensitivity if needed. This is a defined descriptive endpoint, not limma logFC or signed p-values.

## Data use and staged lockbox
T and NK are source contexts. B and myeloid are evaluation contexts. Zero-shot uses no B/myeloid treatment rows (not even their upstream train rows). Few-shot may use only their upstream train conditions, grouped by compound and averaged over allowed donors; every target anchor/control count is reported. Original public_test is development validation; original private_test remains sealed until code and final candidate definitions are committed. Source query interventions are observed, so neither task is a novel-chemical test. Missing conditions remain missing.
Untreated wells are partitioned by their well-row letters: AB baseline input; CD target-anchor reference; EFGH evaluation reference, within plate/context. No evaluation control is passed as a feature or into anchor labels. Source-only transformations may use all source controls. This partition trades library matching for independence and must be compared to the matched-library / pooled-control sensitivity explicitly. Holdout-specific code must enumerate all accessed sample IDs. No normalization fits are estimated across target test samples.
A leave-donor-out sensitivity trains source response summaries and target anchors on two donors and evaluates the third donor; it is not an independent study. Target baseline inputs are restricted to training donors in that sensitivity.

## Preregistered method families and development
Simple baselines: zero response, source mean, drug-agnostic source template, source-only scalar shrinkage, target-anchor mean/offset (few-shot only), shared linear drug-response kernel, and gene-wise affine/ridge source transfer. Candidate family: a mixture of shared and gene-specific kernels with shrinkage, robust selection/averaging across a small fixed regularization grid. Zero-shot baseline-informed multiplicative gates are fitted only on source pseudo-context transfers. Target treatment anchors and original public_test may inform subsequent DEVELOPMENT versions but cannot redefine the first locked result retroactively.
Primary loss: macro mean squared log2-CPM effect error, equal weighting of target contexts and queried compounds. Secondary: rowwise effect correlation/cosine, error split into across-compound mean bias and centered MSE, centered intervention retrieval with ties handled, and per-context/per-donor failures. Test-derived centering is a metric only, never a feature. Source-only responsive-gene metrics must keep the same gene indices for all methods. No target DE-based feature selection.
All fit APIs receive source tensors, allowed untreated summaries and explicit target-anchor arrays, never full target outcome tensors. Baseline-only and few-shot are reported separately. Hyperparameter validation leaves out complete compounds, not individual cells or donor replicas. No external SOTA scores may be compared numerically to this different endpoint.

## Inference and reporting
Preserve every method, failure and subsequent sensitivity; do not promote the best cell type only. Paired resampling by compound is conditional on the observed contexts/donors and is not evidence of population-level biological replication. With only three donors, avoid claims of a calibrated biological confidence bound. Cross-donor signal-energy diagnostics state independence/exchangeability assumptions and do not silently equate donor heterogeneity with technical noise.
External sci-Plex is acquired for a separate untouched experiment. A distinct protocol/code freeze is required before examining its treatment outcomes. No automatic GPU spending, recurring runs or autonomous performance claims.

## Artifact contract
Each scientific run writes to a fresh directory: protocol hash, code hashes, input hashes, sample and compound partitions, predictions before evaluation, selected parameters, metrics, environment and failure status. Counterfactual tamper tests confirm hidden target outcomes cannot change predictions. Tests use synthetic arrays only for correctness; no synthetic sample replaces an unavailable biological row.
