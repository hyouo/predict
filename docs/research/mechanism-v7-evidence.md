# Mechanism-transfer v7 evidence checkpoint — 2026-09-09

Scope: retrospective research on the previously opened six L1000 Phase II test contexts. Zero target-treatment calibration; test baseline controls and measured source responses to query compounds are allowed. Not unseen-compound prediction, a first-use independent test, causal-mechanism recovery or a modern-SOTA claim.

## Actually completed

Acquired Broad Drug Repurposing Hub dated 2020-03-24 annotation tables via the completed bounded Actions run 34317384691. Exact integer PubChem matching: 1038/1726 drugs matched; 1032 have unambiguous target/MOA features, 1038 have unique chemical fingerprints; 276 have at least one annotated target among the 978 measured genes. Preserve the original dated notices; missing annotation is not a no-effect label.

Tested a support-centred bilinear mean prior plus free-intercept source-residual RBF interpolation. This is a classical mean-function/kernel construction. Eight information settings: full support and same-tissue-source exclusion, plus at-most-three-source-context-per-drug masks at three fixed seeds in both scenarios. All methods receive the same retained source response profiles before any learned response embedding. Source baseline controls are not reduced.

## Verified full-support point results

| Method | Primary MSE | Tissue-excluded MSE |
|---|---:|---:|
| RBF | 0.304908744 | 0.348013976 |
| Source-response mean prior + RBF residual | 0.297048311 | 0.344211707 |
| External target/MOA mean prior + residual | 0.302128351 | 0.348678549 |
| Response + correct mechanism | 0.297165823 | 0.344306042 |
| Response + shuffled mechanism, same capacity | 0.296864196 | 0.344154291 |

Source-response prior gain versus RBF is 2.578% primary and 1.093% tissue-excluded. Exploratory fixed-fit context/compound reweighting ranges: [0.954%,4.559%] and [-0.144%,3.274%]. No multiple-testing adjustment, retraining uncertainty, new replicate or new study inference is claimed. Primary drug retrieval decreases from 30.805% to 29.464%. Correct mechanism alone has a small primary benefit but fails the tissue-excluded comparison; joint correct annotations do not reliably beat same-capacity shuffled annotations.

At most three source profiles per drug, response-prior point gains are 2.128/2.239/2.037% primary and 0.822/1.176/1.147% tissue-excluded across seeds 2702/2703/2704. These masks are not independent experiments or prediction ensembles.

Same-selected-kernel-and-amplitude controls show the response innovation, not a new amplitude choice, explains the main improvement. Explicitly post-scoring within-source-count-stratum response-row permutations worsen MSE and interaction risk; they are association diagnostics, not causal tests.

## Verification and failures

25 unit/math/portability checks passed. All 128 original saved arrays hash-verified and context MSE independently recomputed (maximum discrepancy 8.33e-17); these include 32 duplicated baseline controls, not 128 distinct methods. Twelve selected full prediction refits were bitwise identical. Independent Cholesky checks of prior coefficients and kernel computations passed. Six post-scoring identity controls and eight pre-scoring matched-parameter controls completed. One initial posthoc command timed out; its failed directory is retained and not counted as success. A fresh replay completed.

## Publication boundary

This remote branch contains the actual acquisition workflow, protocols and this evidence record. The attempted core-source write was blocked by the connector and was not retried using another mechanism. Full executable source, all tables, numerical evidence and an additive patch are delivered as conversation artifacts, not claimed as uploaded here. Do not merge a superiority claim or replace the default released model. No paid compute or unattended service was started. State, MAP, PrePR-CT and TranSiGen were not executed under this exact protocol.
