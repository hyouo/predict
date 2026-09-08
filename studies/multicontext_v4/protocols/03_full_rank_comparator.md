# Pre-test comparator expansion

Time: after validation-context model selection for the initially registered families, before any held-test treatment or scoring-reference analysis. This is a transparent expansion of the comparator set, not a change to the original candidates or a reaction to held-test scores.

The source-only comparator must not be limited to 16 principal components. Add a full-978-feature multivariate ridge fitted on exactly the same leave-source-context-out input pairs, with context-balanced weights and a fitted intercept. Use the existing penalty grid {0.01,0.1,1,10,100}; independently select the direct-response and residual-around-source-mean versions on the six validation contexts. Both preserve all 978 features and contain no target-baseline features. Refit on the 24 training+validation contexts, save predictions, then jointly score these models alongside all original frozen families. No target treatment calibration is introduced. These are ordinary classical strong comparators, not novel models.

The already-saved original prediction arrays remain unchanged. New comparator arrays and selection logs must be frozen before the single combined test evaluation. The candidate's claimed advantage, if any, must be stated against these full-feature comparators too. No SOTA claim is supported merely by beating them.
