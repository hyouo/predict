# Secondary lower-dose sensitivity

Timing: declared after scoring the frozen 10uM primary and prespecified sensitivity analyses. This is a post-primary extension, NOT part of the original primary protocol and NOT a new independent study. No estimator, hyperparameter grid, or preprocessing threshold is changed.

Evaluate the same six directed SciPlex cell-line transfers at 0.01, 0.1 and 1uM, always 24h. Use the same fixed first calibration order (202609070), budget 16, source-control-only expressed gene panels and recorded library totals. Metadata show that YM155 is absent from at least one cell line at 1uM. Exclude this one query label from ALL three lower-dose sensitivity sets, leaving the same 73 query labels; the 16 calibration labels are unaffected. Do not impute the missing condition as zero. The 10uM primary remains unchanged.

For each lower dose, freeze all seven method predictions before reading that dose's target query outcomes for scoring. Source effects pool both control wells. Target fits use reference A, scores use A and B separately; each role has exactly one target control per standardized plate. Same-reference and cross-reference scores are paired sensitivity analyses with identical predicted arrays, equal control precision, and the same treated query measurements. No claim of biological independence between doses is made: controls and compounds are shared.

The main 10uM results, including failure to beat zero in the crossed-reference setting, will not be replaced. The lower-dose analysis tests robustness of this observation, not a new model selected using the primary test answers. Prediction/code/data/split checksums and all failures must be retained.
