# Research and software status

Updated 2026-09-07. First usable research software: v0.1.0. Historical research numbering (v0.7/v0.8) remains distinct.

## Implemented

Installable `perturb-predict` package and CLI; pinned OP3 prepare; train/save/load/predict; generic aligned-effect bundles and evaluation; frozen prediction records; controlled private/public/control row access; explicit units/context/dose/time; input/hash checks, non-overwriting outputs and failed-run states. Reference model is single-source partial pooling with pooled target controls and no source-spectrum truncation.

## Actually tested locally

51 release tests plus the unchanged repository's 50 tests passed. Installed wheel run outside source checkout produced full B-cell and myeloid models/predictions, repeated the entire public OP3 benchmark, and ran portable inference using only saved models and unlabeled queries. Replayed predictions were bitwise equal. Comparison to the archived v0.8 reference: max absolute differences 1.78e-15 and 8.88e-16. Local dependency download failed because DNS was unavailable; the isolated wheel reused provisioned numerical dependencies. The remote release workflow tests clean dependency installation separately; actual GitHub checks are authoritative.

## Scientific status

OP3 units are conditional pseudobulk counts with the declared log-CPM effect, not individual cells or competition significance statistics. Release uses 11 B-cell and 10 myeloid calibration compounds and 49/47 public query pairs, not the alternative v0.8 larger-label track. Mean MSE 0.1980665 versus source-copy 0.2197485; drug retrieval 0.438124 versus 0.438993. These are already-used public development data, not new independent validation or SOTA evidence.

Explicit donor identity is unavailable; upstream gene-panel selection and shared source/target wells remain limitations. Strict end-to-end unseen-target and new-drug scientific approval remains withheld. No ATAC or modern-SOTA same-task comparison is claimed. No calibrated biological confidence intervals.

## Next research gate

Freeze an independent dataset/experiment-unit protocol and compare this stable release baseline before adding more architecture. Separate prediction-error reduction, drug-specific response discrimination and incremental baseline-information claims. Issues #2/#3 remain scientific tracking gates, not reasons to prevent clearly labeled conditional research use.

Original onboarding README/STATUS remain in docs/history/*_before_usable.md. No unattended research service or paid compute is enabled.
