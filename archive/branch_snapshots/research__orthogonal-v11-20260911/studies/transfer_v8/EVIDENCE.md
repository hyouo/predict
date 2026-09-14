# v8 evidence checkpoint — 2026-09-10

This is retrospective development on the already exposed LINCS L1000 Phase II cohort. It is not a new independent validation, superiority release, or change to the default model. Target treated calibration remains zero. Source-query measurements and target untreated A baselines are allowed. The primary and source-tissue-excluded tasks use fixed masks:1,080 and1,079 target/compound profiles,978 genes.

## Statistical hypothesis and results

We tested whether jointly learning the mean with the same-drug working residual precision improves over v7's separate context-balanced mean regression. An observation-equal OLS control isolates the change in observation weighting. All methods use the same feature ranks and same validation rules.

| Method | Primary MSE | Tissue-excluded MSE |
|---|---:|---:|
| RBF |0.304908743797|0.348013976357|
| v7 context-balanced response prior |0.297048310624|0.344211706996|
| Observation-equal OLS prior |0.297707592858|0.344081239265|
| Joint GLS prior |0.298310790496|0.344216918662|

The GLS candidate does not replace v7: primary error is0.425% higher, and the difficult setting is effectively tied. A working covariance based on baseline similarity is not identified biological noise. The published gls_reference.py is a small independent algebraic core, numerically checked against the grouped implementation; the complete pipeline is in the additive patch/conversation source package, not claimed to be fully present here.

## Actual published-author-core execution, carefully scoped

Pinned TranSiGen author class e7582e7794f13d4dfddfb15e95e3ba7cd24bc07a was loaded unchanged internally, with its original architecture,forward and six-term loss. ECFP4/random initialization,encoder1200,decoder800,latent100,feature embedding400,5,486,956 parameters,seed364039.64-epoch selection runs chose24/32 epochs; final source refits completed accordingly. This is NOT the original300-epoch LINCS2020 study,MODZ processing,KPGT or shRNA-pretrained configuration; it is a single-seed bounded task adapter. Do not use these scores to claim the original published model has been defeated.

On the exact chemical-coverage subsets (703/702 test profiles), calibrated standalone native MSE is0.423055/0.427396. A validation-informed frozen-weight basal reconstruction correction gives0.382840/0.391784. At the same amplitude, this correction changes only a context-constant component of predictions within the covered set; drug retrieval remains9.528%/8.233%. It is not new drug-specific mechanism recovery.

To match direct source-response access, native means were also combined with the same residualRBF. Source validation selected gamma=0 for original inference in both tasks, returning RBF exactly. A baseline-corrected fusion still did not improve final RBF (full primary MSE0.305220, tissue exactly RBF).

## Verification and failure record

All32 final predictor/scenario arrays hash-checked and rescored; repeated controls are not32 independent methods.26 unit/math/adapter tests passed.6 selected statistical refits and2 safe neural-checkpoint inference replays were bitwise identical. No second complete neural retraining is claimed. Independent augmented-systemRBF matches all target genes within4.89e-15/2.23e-15. MonteCarlo8→128 source-validation checks did not change any frozen predictions. An initial concurrent-memory interruption and one packaging timeout were retained as failures, then recovered with new sequential outputs.

The primary protocol preceded new target scoring. The basal diagnostic was added after source-validation inspection but before target scoring. Post-score fixed-amplitude attribution is explicitly not a new fitted model. Dataset exposure from earlier iterations remains disclosed; conditional resampling fixes models and is not biological-study inference.

## Delivery status

Remote branch contains the source-acquisition workflow, protocols, this evidence record and an independently checked mathematical core. Full research code, exact scalar tables, checkpoints and lossless scoring-row evidence are delivered as conversation artifacts and an additive patch. Keep review in draft; no automatic merge, paid compute, journal submission, or unattended scientific service.
