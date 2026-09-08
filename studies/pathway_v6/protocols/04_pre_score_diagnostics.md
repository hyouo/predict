# Additional named diagnostics before v6 target scoring

All learned full-gene, real-pathway and randomized-pathway families selected mixture zero on validation; pathway prior weights therefore do not enter their selected predictions. We will not interpret equality of those test predictions as a positive biological annotation result.

Before any v6 target score, also freeze these explicitly secondary diagnostics:
- Force mixture 0.5 for real and degree-preserving randomized pathway aggregate/local kernels; choose penalty and output strength from the already completed validation grids only. Keep the default families unchanged.
- For the learned-metric full/path/null families, also force mixture 0.5 and select penalty/strength by validation MSE. Refit feature weights on allowed final sources at the same 120-step budget.
- Evaluate a fixed permuted-self kernel using exactly the parameters selected for the true self-gene kernel (both MSE and interaction selectors), not parameters selected on test outcomes. Keep the independently tuned null-self model too.

Refit and freeze both 24-source main and predeclared 18-source tissue-label-exclusion predictions before opening the held-target outcomes in v6. Reproduce the same negative/null comparisons in the stress condition. Tissue-label exclusion includes the literal unknown category for continuity with v5; this is explicitly a metadata stress, not evidence that unknown tissues form a biological lineage.

All cases remain development on previously exposed data. Comparing many exploratory families after primary validation is not a new independent confirmatory test. No best-on-test method substitution or journal-superiority claim.
