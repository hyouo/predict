# Multi-context v4 — source-support-centred response transfer

Research experiment, not a replacement for the released `perturb-predict` package.

The primary cohort is a fixed public LINCS L1000 Phase II Level-3 slice: 30 exact metadata contexts, 1,726 compound labels, 978 supplied landmark genes, 24 h and 10 uM. These are normalized L1000 measurements, NOT raw RNA-seq counts or individual cells. The split is 18 training contexts, 6 validation contexts, 6 test contexts. Hyperparameters are selected on validation contexts; the final model is refitted on the 24 permitted contexts. Test contexts supply untreated baseline controls but no treatment calibration. Query compounds remain observed in at least 3 source contexts.

## Primary results

The prespecified MSE-selected support-centred interaction has held-context macro MSE 0.319051 versus 0.331512 for full-feature direct ridge, 0.322008 for a baseline-conditioned RBF comparator, and 0.409854 for zero effect. Same-model, same-source target-difference MSE is 0.592252 versus 0.634743 for zero difference (6.69% reduction; 14/15 pairs improve). Centre each target difference over compounds: 4.80% improvement remains. These are finite-cohort prediction results, not causal-mechanism evidence or a SOTA claim.

The RBF comparison remains inconclusive under conditional two-axis resampling. Drug retrieval is worse for the candidate (9.36%) than RBF (27.96%). Related backgrounds exist across exact-context splits; this is not lineage-family holdout. Post-primary diagnostics and a full-spectrum development-only prototype are explicitly separated; the latter was NOT scored on the six opened test contexts.

## Data

The acquisition workflow `.github/workflows/l1000-v4-acquire.yml` pins revision, original 707,799,618-byte file SHA256, extraction rules and per-part hashes. Run 34179374331 succeeded; its five artifacts provide metadata and four numeric parts. Put `acquisition.json`, `slice_obs.csv.gz`, `slice_var.csv.gz` and `part_000.npz` through `part_003.npz` in one directory. The conversation's data-slice archive contains the same files and upstream attribution. Rebuild after artifact expiry using the fixed workflow, not a moving main-branch dataset.

## Reproduce

Python 3.13.5 was used locally. Work inside this directory. Public dependencies require internet for initial installation; numerical fitting itself is local.

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PREDICT_L1000_DATA=/absolute/path/to/l1000_slice
python -m unittest discover -s tests -v
python run.py plan --data "$PREDICT_L1000_DATA" --out runs/plan
python run.py fit --data "$PREDICT_L1000_DATA" --plan runs/plan/plan.json --out runs/primary
python compare_full.py --data "$PREDICT_L1000_DATA" --plan runs/plan/plan.json --out runs/comparators
python score_joint.py --data "$PREDICT_L1000_DATA" --plan runs/plan/plan.json --primary runs/primary --comparators runs/comparators --out runs/final
```

Outputs never overwrite existing directories. Both sets of predictors must be frozen before the joint score command. Scientific access guards reject target treatment and target scoring-reference requests before freezing. Numeric NPZ decompression/integrity checking can touch all bytes; no target outcome is supplied to fitting/selection.

## Audit and scope

The protocol, metadata lock, expanded comparator rule and prediction hashes were committed before the primary joint score. The conditional bootstrap fixes the fitted models and processing and does not include biological-study or retraining uncertainty. The 720 assignment enumeration is a diagnostic, not an experimentally randomized biological p-value. Existing L1000 normalization can induce dependencies among wells. Raw slice reconstruction and all 13 model refits matched; 78 context-by-model matrices replayed bitwise locally. Eighteen unit/guard checks passed. See the complete report and tables for adverse metrics and all variants.

No State, MAP, PrePR-CT, PRnet or TranSiGen result is claimed: none was run under this exact task. The comparison set here is statistical. Cite LINCS L1000/GSE70138 and Chem-PerturBridge, preserve their notices, and do not treat sentinel pseudobulk count fields as actual counts.
