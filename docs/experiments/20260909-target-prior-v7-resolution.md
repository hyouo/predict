# v7 full-annotation resolution sensitivity

This extension was implemented and locally committed before the first new v7 target scoring call. Existing L1000 outcomes were exposed in earlier versions, so this remains retrospective development. Stage-1 predictions are retained unchanged.

Metadata-only SVD inspection found rank 16 retains 12.8856375% of squared target-MOA incidence norm; this is not a percentage of biological information. The sparse 2852-feature matrix has 4743 nonzeros and 360 connected components. Fit the same source-LOO, support-centred residual model with all annotation features using sparse normal equations and full 978-gene outputs, no drug rank truncation. Source baseline rank remains 6. Lambda {.1,1,10,100,1000}, correction strength {0,.25,.5,1}. Repeat rank-16 models with the identical grid. At both resolutions include the same seed7313 row-permuted annotation null, retain selected and forced-strength-one predictions. Do not use target outcomes to choose annotation resolution.

Both scenarios and both stages have now completed local fitting and saved individual prediction checksums. All four freezes will be jointly checked and committed locally before unlocking scoring. New external annotations retain Broad non-commercial research notices; no Reactome features, modern-method comparison or mechanistic interpretation is claimed.
