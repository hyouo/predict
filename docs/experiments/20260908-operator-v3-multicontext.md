# Operator v3: test whether more cellular contexts make untreated-state information useful

Registered before local access to the new quantitative LINCS subset or model fitting. Previous operator_v2 used 2-3 source states per outer fold; baseline ablations did not establish incremental information. Preserve those negative results. No algorithmic superiority is assumed.

## Resource and task
Use the fixed public Chem-PerturBridge LINCS phase-II Level-3 landmark file at revision 6c54c8eb0321cceff4f888c54e199077d055e20b, SHA256 34d198df9eddac5b535a0408794caaf41073fa014c74429d1050e483576fb581. Metadata-only acquisition fixes 24h and 10uM compounds, retaining the matching controls. This is an expression assay, NOT RNA-seq counts; never apply CPM/count normalization by guessing the units. Audit the upstream scale and sample mapping before training.

Task: zero treated target-context labels, target untreated observations allowed, query compounds measured in development source contexts. This is not new-compound prediction. The new resource is distinct from OP3/SciPlex, but any connection to earlier L1000 pretraining or external data must be disclosed.

## Metadata-only cohort rules
Keep treatment rows only when their cell/plate/time group has at least four control rows. Require at least two treatment rows per cell/compound, at least 64 qualifying compounds per cell, and a nonmissing context identifier. Split eligible contexts using SHA256('operator-v3-20260908|' + context): reserve the first max(2, floor(C/4)) contexts, capped at six, as final holdouts. At least six development contexts are required; otherwise revise the metadata protocol transparently before fitting, not after seeing results. A query compound must have measured outcomes in at least three development contexts. Missing outcomes remain missing. All measured landmark columns are retained; no target-label gene selection.

Split controls within each cell/plate by sorted sample ID into alternating A and B records. Only A controls form target baseline inputs. Effect estimates subtract plate-matched controls; the final primary score uses B target controls. Same-reference scores are secondary, never substituted for the primary. Different control records do not automatically imply independent biological experiments. Quantile normalization and upstream selection remain limitations.

## Development and candidate families
All iterative design, model selection, scaling, projections and hyperparameter tuning use development contexts only. Validation leaves complete development contexts out, never merely cells or compounds in a known context. Candidate families: zero; measured source mean/median; globally shrunk source mean; baseline-neighbor interpolation; source-response-only multivariate ridge; baseline-additive ridge; baseline-by-response bilinear ridge; nonlinear baseline-kernel interaction ridge. Preserve matched no-baseline and permutation controls. Optional additional prototypes must be documented as development extensions and frozen before accessing final target treatment values.

The interaction predictor anchors on available-source compound effects and adds a low-rank residual depending on both the untreated target state and the source response program. Classical ridge/kernel/tensor methods are not claimed to be invented here. The hypothesis is incremental target-state information that survives a strong same-information response-only comparator.

Initial grids: response rank <=16, baseline rank <=8, ridge penalties 0.001,0.01,0.1,1,10,100 and source-anchor weights 0,0.5,1. Neighbor bandwidth/mixing use a finite recorded grid. Numerical fallbacks and source-coverage counts must be explicit. Nothing may select parameters from final target responses.

## Final gate
After development, save all candidate predictions, parameters, code/data/split hashes and access logs before final target scoring. Score the complete eligible query set: total and centered MSE, bias, response identity retrieval, nonzero amplitude and per-context failures. Compare to zero and response-only strong controls, not just source copy. Conditional compound/context resampling is not a guarantee of population-level generalization. Do not reopen the holdout for another round of model selection. Report all variants and any failures; not-run modern methods remain not-run. No paid compute or unattended research loop.
