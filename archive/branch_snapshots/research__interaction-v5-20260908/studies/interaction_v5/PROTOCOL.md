# Interaction v5: retrospective, outcome-isolated iteration

This is not a first-use blind test: the same L1000 Phase II dataset and six held contexts were evaluated by multicontext_v4. No new independent study, causal identification, or modern-SOTA comparison is claimed. Main and prior research remain unchanged.

## Source and information contract

Restore the fixed public asset with SHA256 `34d198df9eddac5b535a0408794caaf41073fa014c74429d1050e483576fb581`, using the five successful artifacts of run 34179374331. Independently verify archive and extracted numeric-part hashes. Retain Level-3 normalized landmark values; do not reinterpret them as single-cell or RNA count data. Use the original metadata-only 30-context, 1,726-label, 978-gene, 24h/10uM cohort and 18/6/6 train/validation/test context split. Only baseline A controls from held targets are allowed before the joint prediction freeze. No held-target treatment labels enter fitting or selection. All candidate predictions and stress-test predictions must be frozen before this round's joint scoring.

## Models fixed locally before fitting

Statistical comparators: zero, source mean, scalar shrinkage, full-feature direct/residual ridge, rank-6 and full-baseline RBF, full-baseline linear kernel, v4 support-centred interaction, leave-source-out anchored and direct interactions. The regression grid is penalties 0.01, 0.1, 1, 10, 100 and interaction/amplitude strengths 0, 0.25, 0.5, 0.75, 1. Separate MSE and interaction-selected predictions must not be conflated.

New candidate 1: positive context-distance metric conditioned on source response programs, with an independently trained static-metric ablation. Source pseudo-target profiles are excluded from their own input means and attention support. Baseline rank 12, drug rank 16, seeds 11/37/71, Adam learning rate 0.02, at most 200 steps, validation every 20. Learned positive diagonal metric has 205 scalar parameters including amplitude. PCA representations use permitted training contexts only; inner leave-source-out features do not imply that every fitted representation is cross-fitted.

New candidate 2: the same metric trained on four-way contrasts: (cell a minus cell b) for drug p minus the same difference for drug q. The source pool excludes both pseudo-target contexts. Only drug pairs with identical rational per-plate control-weight signatures within each target are used; the resulting labels cancel shared reference vectors algebraically. Multiply quartet squared loss by (block size minus one)/(2 times block size) to match within-block centred loss. Static/dynamic versions, seeds 11/37/71, 400 steps maximum, Adam 0.01, batches 64, validation every 40. No late choice of favourable seed; report three-seed mean and all individual outcomes.

## Repeatability and interaction evaluation

Partition entire physical treatment/control plates by the recorded _Xn replicate tags: X1/X3 versus X2. These partitions have unequal precision and are not assumed to be independent biological replicates. Within each target pair, centre predictions and both measured difference arrays in groups sharing exact reference-weight signatures in both partitions. Remove singleton groups and disclose coverage. Report raw split cross-products, error against each split, prediction energy, and gain P*(A+B)-P^2 versus zero interaction. The gain has an exact observed-risk interpretation; a latent biological-risk interpretation additionally needs unbiasedness and suitable prediction/error independence. Cross-covariance between partitions is not assumed zero. Never clip negative cross-products or silently turn missing profiles into zeros.

## Additional fixed stresses

1. Remove source codes ASC.C, NPC, NPC.TAK, NEU and SKL for all targets, then reselect parameters on the four remaining validation contexts. This is a code-neighbour exclusion, not verified lineage holdout. Record lost support/coverage.
2. Hold out the first 20% of drug labels sorted by SHA256 of `interaction-v5-drug-holdout:<label>` from parameter learning, PCA fitting and validation selection. Supply these drugs' permitted source measurements only at prediction. This tests new support queries, not unseen chemicals without source experiments. Compare on the identical query subset.

All raw MSE, centred error, target differences, interaction scores, retrieval, coverage and adverse results are retained. Context/compound resampling fixes fitted models and upstream processing and is descriptive, not independent biological inference. No hypothesis is judged from the number of genes, pairs or repeated fits. State, MAP, PrePR-CT and TranSiGen remain unrun under this exact task.

## Execution record

Original v4 source reference: 3dc6849385452ef1baf2ef79b37ad2912466f0ea. Local v5 protocol/implementation commits before scoring include c791912, 1485491, 697454a, 3bfeb55, 9ce6cac, b8d385c and 80b61f8. An initial pandas index ambiguity failed before training and was repaired; the failure record is retained. A prior continuation attachment contained no recovered data/results, so this round rebuilds from the actual verified numerical input rather than treating that attachment as a completed experiment.

This file is a reviewable protocol checkpoint. Complete executable code, frozen arrays and verification records are separate artifacts until explicitly uploaded. No automatic merge, scheduled research or paid compute is enabled.
