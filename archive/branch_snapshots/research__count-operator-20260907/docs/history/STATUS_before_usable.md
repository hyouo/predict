# Research status — repository onboarding

Updated: 2026-09-06. Baseline: v0.7 core import. Stage: reproducible repository setup, not a new algorithm result.

## Evidence carried forward

The auxiliary panel contains five cellular backgrounds, 72 ligand-present conditions per background and 18 phosphoprotein readouts. No independently replicated stimulated response is available per retained condition. The task permits measured source responses of query interventions and target calibration conditions.

Historical expanded-stacking MSE is 0.254994 in the declared pair-mask study, versus 0.263091 for the stronger simple readout-marginal comparator (about 3.08% relative improvement). This does not establish general dominance; the candidate loses to simple alternatives under some factor-blocked tasks.

The equal-cost probe study changed information: 17 covered conditions gave MSE 0.476562, while base16 plus one missing-factor probe gave 0.356214. Most of the gain concerned query-family mean error. It is not strict unseen-factor or RNA zero-shot success. Two historical CSV summaries are copied to reports/v0.7; they are not newly measured during onboarding.

## Checks performed during onboarding

The original 37 tests and 8 new repository utility tests passed locally (45 total). All 15 imported Python files matched original archive SHA256. A fixed real-panel integration case ran with both bootstrap and stacking; the largest full-versus-blocked prediction difference was approximately 1.33e-14. This is implementation validation, not another biological study.

The initial remote Research CI run 34015151595 also passed all 45 checks and the measured-data integration step. Local/remote floating-point predictions are numerically consistent but have different byte hashes; no cross-machine bitwise identity is claimed. Five additional metadata-decoder tests pass locally, giving 50 checks in the updated checkout.

OP3 acquisition SUCCEEDED in run 34015151607. The 23,700,750-byte fixed HDF5 file was retrieved through the repository artifact, and both artifact and matrix SHA256 were rechecked locally. The CSR matrix has shape (1813, 5288); metadata include four cell-type codes, 138 non-control perturbagen labels, 192 control observations and six plates. These are not asserted to be individual cells or raw counts. There is no explicit donor column. Read reports/data/op3_acquisition.json and op3_inventory.json. The artifact expires on 2026-09-13; fixed upstream acquisition remains reproducible. Full measurement semantics and RNA training are still pending; sci-Plex was not downloaded during onboarding.

## Not completed

No quantitative RNA/ATAC model training; no complete target-context zero-shot biological validation; no same-task execution of the latest large baselines; no independent prospective experiment; no calibrated biological confidence claim.

See docs/ROADMAP.md and open issues for the next bounded work units. No autonomous research service or scheduled paid computation is enabled.
