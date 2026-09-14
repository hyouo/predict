# 004 — Design-aware reference covariance refinement

Registered 2026-09-06 AFTER OP3 public/private and external sci-Plex results were inspected/recomputed. All subsequent performance values in this experiment are exploratory follow-up, not a fresh confirmatory test. Existing estimator/protocol/results in 002/003 are not overwritten.

## Observation
The pooled paired-reference correction uses 0.5*d.T@H@d/(n*G), where d averages control differences across independent reference units. It contains products between different reference units. Those cross-products have zero expectation under independent, centered reference errors, but add sampling variation. A treatment panel has few reference wells, not G independent biological replicates.

## Fixed proposed change
Keep candidate kernels, priors, regularization grid, target anchor conditions, controls and primary endpoint unchanged. Represent the reference difference as d=sum_u F_u, where F_u is the weighted contribution from a reference unit. Replace the pooled correction by 0.5*sum_u(F_u.T@H@F_u)/(n*G). For OP3 report both donor units (retaining dependence between a donor's two plates) and donor-plate units. For sci-Plex report assay-plate units; plates are NOT called independent biological donors. Same 71 candidates and ordinary error Gram; only the LINEAR covariance correction entering stacking changes. Compute both simple and full candidate-bank outputs, without choosing the best grouping using test scores.

Under independent sign-symmetric zero-mean reference differences, this is the sign-symmetrization of the pooled quadratic correction and has no larger variance. The claim is conditional on fixed affine smoothers and the stated covariance structure; selection of weights, source-target shared technical errors, drift/nonzero means and cross-plate dependence can invalidate a naive unbiased-risk interpretation. This is an application/refinement of correlated-data CV, not a claim to have invented CV or stacking.

## Execution
Use already-sealed OP3 all-donor inputs and three leave-donor-out constructions, and the already-sealed external sci-Plex folds. Reconstruct reference factors only from permitted training controls and anchor metadata. Assert factor sums exactly reproduce the paired anchor-reference difference. Predictions are written and hashed before their exploratory evaluation. Preserve original models, ordinary/pooled-correction models, new donor/plate-correction models and all worsening cases.

Report macro MSE, centered and mean-bias components, centered retrieval and each target context. Conditional compound bootstrap is not biological replication. No SOTA/Nature claim follows from beating the ordinary-CV comparison, especially when a simple source-only baseline is equally good. Tests must verify one-unit equivalence, factor-sum validation, block relabeling invariance and expected correction/variance in a labeled Gaussian simulation.
