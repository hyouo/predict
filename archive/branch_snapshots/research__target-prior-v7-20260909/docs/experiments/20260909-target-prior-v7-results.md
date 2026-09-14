# v7 external-target prior: completed retrospective evidence record

Date: 2026-09-09. Continuation of the verified annotation acquisition and protocols on this branch. This is not a SOTA release, independent biological confirmation or journal submission. The existing L1000 cohort and six target contexts have been used in earlier project development.

## Actual new information and implementation

Exact unique PubChem CID joins to the Broad Drug Repurposing Hub 2020-03-24 drug/sample files: 1032/1726 mapped, 904 with target labels, 276 with a measured target among 978 landmarks, six ambiguous CIDs excluded from annotation but not overall scoring. Preserve non-commercial research notices. No Reactome data were acquired; both requests returned HTTP403. Nominal targets/drug-level MOAs are not target engagement, and the target-by-MOA Cartesian-product encoding is not an established per-target action map.

Compared fixed full-feature RBF, target-gene distance mixtures, rank16 and full2852-feature baseline-by-target residual operators, corresponding fixed annotation permutations, and post-score shared-only target/ID/constant corrections. Source-LOO residuals and actual-source-support centring are used. Target treated calibration is zero. Full and rank16 resolution comparisons share the extended penalty grid. All 34 stage1/2 arrays froze jointly before new v7 target scoring; the 16 shared-only ablation arrays were explicitly post-score development and froze separately. All adverse/forced-component results are retained.

## Results (same information budget)

| Method | Original split MSE | Same-tissue-source-excluded MSE |
|---|---:|---:|
| Fixed RBF | 0.304908744 | 0.348013976 |
| Rank16 target-conditioned correction | 0.302789731 | 0.348603513 |
| Rank16 shuffled annotation | 0.304419196 | 0.348043208 |
| Full target-conditioned correction | 0.302344569 | 0.349437179 |
| Full shuffled annotation | 0.303169663 | 0.349075980 |
| Shared target correction (added term has no target-state input) | 0.306354401 | 0.350838266 |
| Same-covered drug-ID shared correction | 0.305432189 | 0.353278879 |

Original split: rank16 improves0.695% against RBF, all six contexts; fixed-model two-axis reweighting range0.108–1.600%. Full-feature point improvement0.841%, four of six contexts; range-0.114–2.187%. These are descriptive current-cohort intervals, not multiplicity-corrected inferential or study-level claims. Harder-context errors increase0.169% and0.409%. Retrieval is not improved. Full-feature interaction risk does not improve despite its lower overall MSE. Direct target-only distance is rejected by validation (eta0); forced use increases error in both scenarios.

Shared-only corrections do not explain away the selected conditional models by matching their scores, but they are not exhaustive alternatives. Their additions provably cancel from same-source target-pair contrasts, verified on real arrays. No target-conditioned causal mechanism has been identified.

## Verification and delivery

26 mathematical/interface checks passed.50 stored prediction files were hash-checked and rescored,40 saved coefficient replays matched exactly, and8 real refits matched bitwise. Independently extracted evidence reproduced52 model/scenario checks (including two derived zero baselines),156 subset summary rows, all978 genes. Fresh code archive and additive patch each passed26 tests. These checks are not biological replicates.

Conversation deliverables:
- predict_target_prior_v7_research.zip:1528862 bytes; SHA25603eab7ca0561fac7dc3e5aa52721ab78595d19cd13209e5d71d5646346293fbb
- predict_target_prior_v7_evidence.zip:285022268 bytes; SHA256cc271236b8f01f4d65cc6e50457a1005d88fe9899b050ee648a2ed9a3718f796
- predict_target_prior_v7.patch.gz: additive studies/target_v7 source/results package; SHA25612be9c5a350f46d0bb7922257575742fcb7124e22bf44a139609d23756110746

This repository record does not claim that complete analysis source, raw annotation tables or prediction arrays have been uploaded into Git. They are supplied through the verified code/evidence archives and additive patch. State, MAP, PrePR-CT and TranSiGen were not executed under this exact protocol. Keep the branch/PR in draft; do not automatically merge or change released predictions or manuscript superiority claims.
