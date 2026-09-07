# Validation-stage response alignment extension

The original global context-kernel validation candidates outperformed the shared low-rank operator. This extension is informed by validation results, not by final test outcomes. All nine test-cell treatment arrays remain unopened by the fitting/evaluation processes. The original seven families and their failures remain reportable.

For training-only program r, center context maps A_r and untreated features Z over the 32 training cells. Define input-gene association s_rg = z_g^T A_r A_r^T z_g / (||z_g||^2 ||A_r||_F^2). Set raw weight max(s_rg - 1/(C-1), 0), normalize and mix with 10% uniform mass. This is a regularizing heuristic, not an unbiased causal/noise estimator. It asks which basal genes are associated with differences in the response map, rather than assuming only strongly responding output genes determine susceptibility.

Compare aligned program kernels, a shared aligned kernel (mean of 16 training program metrics, full response residual output), and matched permuted input-gene metrics. Preserve ranks 4/8/16, RBF bandwidths .3/1/3, ridges .01/.1/1/10, and amplitudes .25/.5/1/2. Validation objective and source/control information unchanged. Keep every candidate score and select parameters on the nine validation contexts only.

Select the overall family by validation BEFORE opening test treatment arrays. Freeze all original and extension family predictions, parameters and code/data hashes first. Test outcomes will be scored once without post-test tuning. Report comparisons with mean source, global kernel, shared-program, original program kernel, and all permutations; do not substitute a test-winning family for the validation-selected family.

These are statistical context-response mappings, not verified regulatory circuits or a drug-target prior. State/MAP/PrePR-CT remain unrun in this specific count-composition experiment and cannot be claimed beaten.
