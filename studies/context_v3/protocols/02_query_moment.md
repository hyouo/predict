# Validation-stage extension: retain query-specific source-state/response associations

Recorded 2026-09-08 after development-context results, before any final8 target treatment scoring. Mean-first bilinear baseline conditioning did not improve validation MSE. A full-gene untreated-control acquisition write was blocked; it was not retried or rerouted. This study continues with the existing2000-gene baseline panel and does not claim full-genome input.

Hypothesis: averaging source responses before conditioning discards the association of source state with response for the particular query drug. Estimate ridge slopes of source response coordinates against source baseline coordinates separately for each compound, using exact fixed-penalty leave-one-out formulas to remove a training recipient's own response. Combine the predicted conditional deviation with source mean response coordinates in a shared cross-gene decoder. All PCA bases and regressions use only34 training contexts; final8 contexts provide controls only.

Development grid: (response rank,baseline rank)=(16,4),(16,8),(32,8); state regression tau=0.1,1,10; decoder ridge=(.001,.01,.1,1,10,100); direct source coefficient=(0,.25,.5,.75,1). Select on8 validation contexts. Refit recipient-baseline-muted control at selected geometry with separate decoder tuning. Compare strong no-baseline residual ridge and baseline-kernel interpolation. All components are classical; the experiment tests information content, not priority of ridge or analytic LOO.

All model choices, parameters and final predictions will be committed/frozen before the single final test scoring pass. Any later analysis is explicitly post-hoc; no replacement of failures. The original task, response panel, dose/time and query cohort remain unchanged.
