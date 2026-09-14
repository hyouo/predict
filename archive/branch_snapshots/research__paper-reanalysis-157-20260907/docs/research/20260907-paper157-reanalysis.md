# Paper evidence record: 157-compound quantitative sci-Plex reanalysis

Date: 2026-09-07. Status: computational manuscript draft; not submitted; no SOTA or Nature-acceptance claim. This record is separate from the earlier 16-calibration/74-query workflow on `research/paper-validation-20260907`.

## Actual analysis

Pinned public sci-Plex3 processed counts: 4,974 pseudobulks, 56,799 genes, 139,848,567 bytes. SHA256: `8ad35a20e44bf41cb030a0f063eb6a5727965b8f27d0eda4380bda0e294ccc4c`. Data revision: `6c54c8eb0321cceff4f888c54e199077d055e20b`, theislab/chem-perturbridge. Acquisition run 34082109943, artifact 10004012749.

Primary: 10uM, 24h, at least 10 contributing cells, 157 compounds after documented prior-signature exclusions; all six ordered A549/K562/MCF7 transfers. Select 3,000 genes using source controls only. Five complete-compound folds, target calibration budgets 8/16/32/64, 32 primary. Calibration and scoring use distinct hash-split target control records; record disjointness does not establish biological independence. Source query responses are allowed. This is few-shot missing-pair prediction, not unseen-drug or target-baseline-only zero-shot prediction.

The fixed partial-pooling core has SHA256 `80ddbbc643f03844d589f89c365fb6da00024b38b5cf00d59a36cc590b43b176`; coefficients were checked against repository commit `575c325a4dd6f99eb055a305e789932667835b93`. Grids and estimator were not adapted to current test scores. Prediction archives and splits were saved before scoring. The offline loader does materialize target arrays; the safeguard is exclusion of query outcomes from fitting/tuning, not physical blindness.

## Main result: do not cherry-pick the comparator

| Method | Equal-direction MSE, budget 32 |
|---|---:|
| No effect | 1.002071 |
| Partial pooling | 1.054363 |
| Global slope plus shrunken offset | 1.054797 |
| Gene-wise ridge | 1.055384 |
| Linear KRR | 1.103987 |
| Target mean | 1.109117 |
| Source copy | 1.551221 |
| Residual KRR | 1.623019 |

Partial pooling reduces MSE 32.03% relative to copying, but is 5.22% WORSE than no effect and only 0.04% better than global-slope calibration. No effect wins every direction; MCF7-to-A549 also favors copying over partial pooling. Drug retrieval is 1.38% for partial pooling versus 6.79% for copying. The 942 direction/compound outcomes are not 942 independent experiments.

Paired 5,000-compound-bootstrap intervals: source-copy reduction 30.94% to 33.17%; increase over no effect 4.31% to 6.12%. These are conditional on fixed fitted models and the current cohort, excluding fitting/plate/study uncertainty.

All eight sensitivity settings and 100 within-fold query-source permutations were completed. Partial pooling loses to no effect in seven sensitivities. Plate-shared references change the ranking, but ordinary linear KRR is better in that setting. Permutation shows a small drug-specific increment, not no-effect superiority. Source-response-cluster blocking evaluates only 686/942 possible queries at the 32-compound budget; missing coverage is explicitly reported.

## Prior-exposure correction

The local current protocol was time-stamped at 05:00:18Z, SHA256 `726c16cf55980933f92fd1d82c865dc53fd98615eb1a0b84d61ba589fdfe6bec`. During final synchronization, repository commit `c2dd42122035f31050499020e09e70cc2ca6c9dd` (04:30:58Z) was found to document an earlier same-day quantitative SciPlex evaluation. Signature exclusions therefore do NOT establish a wholly untouched quantitative cohort. This work is external relative to OP3 but is an external reanalysis, not first-use blinded confirmation or public preregistration. The original protocol, splits, estimator and numerical results are unchanged; the limitation is included in both manuscript and supplement. Do not silently combine the two workflows.

## Reproducibility and deliverables

Sixteen new mathematical/integrity checks passed. A separate process reproduced all eight models for one complete primary direction bitwise; all 48 primary method/direction MSEs were recalculated from arrays. Both OP3 development replay arrays were reproduced bitwise. Interrupted computations and a concurrent-memory failure were retained and excluded from final scores.

A 14-page English main manuscript and 6-page supplement were compiled and visually checked using the unmodified official Springer Nature sn-jnl template, sn-nature option, December 2024 release. Six vector figures and their source tables accompany the manuscript. Human authorship, affiliations and disclosures are pending; AI assistance is disclosed. Full LaTeX, analysis code, data manifests, detailed scores and prediction archives are supplied as conversation deliverables and an additive repository patch. This narrow PR publishes the evidence record, not the complete analysis tree or a new trained-model release.

The scientifically supportable claim is a calibration/evaluation limitation, not a breakthrough response model. Modern same-information model comparisons, genuinely new controlled validation and stronger biological replication remain necessary for a superiority manuscript.
