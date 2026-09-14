# v6 post-score diagnostic: can a drug-shared baseline slope explain the own-gene gain?

Timing: declared after all 120 primary/stress predictor configurations were frozen and scored. This is a new, explicitly post-score exploratory comparison; it must not replace or relabel the first-stage results. Hyperparameters still use only allowed validation contexts, not target outcomes.

A matched own-gene baseline kernel improved main-context MSE, but its advantage disappeared in the tissue-excluded stress test and its interaction risk became positive. To distinguish additive response drift from drug-specific context effects, fit a single slope per output gene across permitted source profiles, after centering baseline and response within each drug's available source support. Use ridge penalty fixed at 1 times the median positive source centred sum-of-squares, no added slope hyperparameter search. This is a classical regularized drug fixed-effect model.

Evaluate (i) this pooled slope alone plus source-drug means, and (ii) it as the prior mean of the global affine kernel, predicting residuals with the same kernel grid [0.01,0.1,1,10,100] and the same common/correction grids as v6. Neither model uses target treatments. Save all configurations, fit state, and predictions before re-scoring.

The pooled-only model's difference between two targets is constant across drugs and must have exactly zero reference-block-centred interaction. Check that property numerically. Any MSE gain from it is NOT evidence of predicting drug-specific target interactions. Residual-kernel results must retain MSE, centred and bias components, retrieval, and interaction excess in both source settings. No new independent study or SOTA claim is implied.
