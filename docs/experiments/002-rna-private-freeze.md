# OP3 private split: first estimator/code freeze

2026-09-06. This amendment is recorded AFTER inspecting public_test development effects, and BEFORE evaluating private_test effects. Initial protocol: de6c0649a6b653b7874a816c0d5be34d4d3d70c2. No claim that the public split is a fresh confirmation.

The current source-complete universe is 137 compounds. R428 is absent in NK and is explicitly omitted, never imputed as zero. This leaves 11 B-cell anchor compounds (33 treatment observations) and 10 myeloid anchors (30 observations), across three donors. All 5,288 released genes are used, conditionally on upstream expression filtering. Untreated features AB, anchor references CD, evaluation references EFGH are disjoint in target wells. Source cells may share physical wells with targets; donor-held-out sensitivity is therefore essential, and control-error exchangeability is only a working assumption.

Current public-test result: adjusted stacking improves ordinary stacking, but the mixed gene/shared bank adds very little to the simpler adjusted bank. The source-only baseline gate does not outperform gene-specific source-consensus shrinkage. All methods are retained, not just the winners.

## Frozen definitions
`src/rna_models.py` includes 71 affine candidate configurations and 25 named zero/few-shot output methods. Lambdas .03,.3,3,30; no, shrunken, and free intercepts; shared linear/RBF, gene-linear, and 50:50 mixed kernels. Paired-control CV correction is one half of mean(d * H d), where d is the difference between effects computed with C and D reference wells. Candidate weights minimize an ordinary error Gram quadratic plus the LINEAR reference correction, on the simplex. The estimator is inspired by existing correlated-data CV and stacking; novelty is not asserted.

Code SHA256 at this freeze:
- src/rna_data.py aed0d64ccd572ff83ad47f7a39a8864b9ea8596b2cfa2aa26662f59a9d6b5b02
- src/rna_models.py 059ca5276f80307c566b5b735f3077d4ff4b3474358ecaa8871f782d3d00f46c
- src/rna_metrics.py 03bb416e3d8b7fcc0ddfeffbd7cdae8140ee08387da04b93396d085b3515ace2
- tools/rna_experiment.py ffbd1f42477bb69378e8693382d1e322179566248774f6bba79ed2de377ea486

Already-sealed all-compound predictions (fit without public/private target outcomes): B NPZ 448e77db22fa2cd4cbec3a19f689812e20a06bac0eff8245f5bba9cd1114a21c; Myeloid NPZ 3b774b5028bab423316f4e3705b36e5e9870afd60c453093a00bfcd9e466c353. Private_test is next evaluated once using these files. Subsequently evaluate the SAME frozen method families in all three leave-donor-out folds; source/baseline/anchor data restricted to other two donors. No donor-level hyperparameter tuning with held-out donor outcomes.

Primary comparison: mixed_stack_adjusted versus mixed_stack_cv (same candidates/data); also report simple_stack_adjusted and source_consensus_shrink to separate architecture from regularization. Macro MSE equally weights B and myeloid and compounds within each. Report centered MSE, bias, retrieval and every failure. A source-only method uses zero target treatment anchors, not the same information budget as few-shot. Paired compound intervals are conditional on the observed donors/contexts, not independent biological validation. Counts were accessed for exact file reconciliation, not for private-response modeling or selection before this freeze.

64 implementation/math/access tests passed locally before the freeze. Models currently have no external sci-Plex performance evaluation; external protocol follows separately.
