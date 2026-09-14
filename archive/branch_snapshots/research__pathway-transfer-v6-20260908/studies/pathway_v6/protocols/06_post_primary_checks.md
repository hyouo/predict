# Post-primary checks (not independent confirmation)

The frozen primary scores are now open. MSE-selected self-gene transport has MSE 0.300428633 versus universal-full 0.306387807 and matched permuted-self 0.311538106. All six context point errors improve versus universal-full. The tissue-label exclusion reverses the average self-versus-full MSE comparison and worsens interaction risk. These adverse results will not be replaced.

No model or selector will be retuned. Add only these diagnostic comparisons:
1. Analytically match the global-kernel attenuation implicit in self mixing: remove the self rank-one term but retain (1-gamma) times the full kernel, equivalently lambda/(1-gamma), with the exact chosen gamma=0.25, lambda=1 and output strengths 0.75/0.5. This rules out interpreting a finer effective regularization value as same-gene information.
2. Replay matched self permutations using fixed seeds 2026090802 through 2026090811, unchanged parameters and source/target data. These are input-assignment diagnostics, not biological randomization tests or newly independent data.
3. Swap the two pre-existing control partitions A/B for baseline construction and effect references across all contexts, refit the same rules with fixed parameters, freeze before the swapped-reference score. Treatment rows and total control availability are unchanged. This is a reference sensitivity on the same experiment, not an independent study; odd control counts may make half-sizes unequal.
4. Verify the gene-specific kernel against its independent rank-one-update / residual-regression expression and the source metric Gram loss against explicit predictions. No claim of inventing these linear-algebra identities.

Compute paired conditional resampling summaries on fixed predictions and fixed reference-matched blocks. Make explicit that only six held contexts exist and that intervals omit retraining, preprocessing, reference-selection and study-level uncertainty. Preserve original MSE/interaction/retrieval and all failures.
