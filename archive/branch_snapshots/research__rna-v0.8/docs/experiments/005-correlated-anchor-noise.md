# 005 — Correlated anchor measurement noise, not only a CV correction

2026-09-06. Exploratory follow-up AFTER the original OP3/private/leave-donor-out and external sci-Plex results, and after experiment 004. Experiment 004's reference-unit symmetrization has not materially improved prediction. No previously evaluated split is called fresh confirmation.

## Hypothesis
The existing candidates fit anchor residuals using K_AA + lambda I even though anchors sharing a control have correlated measurement errors. Correcting CV alone does not correct that fitting likelihood. Keep the original 71 candidate definitions, kernel shapes, priors, lambda grid, controls, compounds, gene sets and endpoints. Add a measured-reference covariance N_g = 0.25 sum_u F_ug F_ug^T, where factors F sum to C-D reference differences. Run a full covariance and a diagonal-only covariance variant (same marginal variances) to distinguish correlation from heteroskedasticity. OP3 units are donors, retaining within-donor plate dependencies; sci-Plex units are assay plates, not asserted independent biological replicates.

For each kernel candidate, fit with K_AA + N_g + lambda I; fixed/offset candidates remain unchanged. Compare ordinary and reference-adjusted stacking of the identical candidate bank. Risk adjustment uses the FULL reference-unit covariance for both fitting variants. Do not choose a variant after seeing test scores. Empirical control covariance is a noisy working estimate; source-target shared technical noise, nonzero control drift and selection effects remain limitations.

## Exact LOO derivation to test
Let B=(K_AA+N+lambda I)^-1. With no free intercept, S=K_AA B; with a free intercept replace B by Q=B-B11^T B/(1^T B1) in the latent fit and add the GLS intercept. The latent leave-one-anchor-out map is H=S-diag(diag(S)/diag(Q)) Q (Q=B without a free intercept), not the ordinary I-Q/diag(Q) formula unless noise is scalar. Directly refit each held-out anchor in tests with non-diagonal noise, for all three intercept modes. Cross-covariance between measurement errors is NOT added to the query latent-signal kernel.

## Execution and reporting
Fit all variants on allowed anchors only; seal and hash outputs before exploratory evaluation. Run OP3 all donors and external sci-Plex first, then the SAME fixed variants on all three OP3 leave-donor-out folds. Report MSE, centered/bias components, retrieval, every context, same-information ordinary/corrected controls, and source-only/zero references. Negative results remain. This is a GLS/GP application, not a claim to have invented Gaussian processes or established a new biological mechanism.
