# context_v3: larger-context zero-target-treatment validation

Recorded 2026-09-08 before obtaining any numerical Tahoe response in this session.

## Objective
Test whether baseline-conditioned cross-gene response transport gains predictive information when many more cellular backgrounds are available. This is NOT a new-compound task: query compounds are measured in source contexts. No target-context treatment labels may enter training or model selection. Tahoe has not been numerically modeled in the provided prior study archives; it is new to the local development sequence, not necessarily unseen by public foundation models.

## Data and split
Pinned Chem-PerturBridge revision `6c54c8eb0321cceff4f888c54e199077d055e20b`, Tahoe SHA256 `09ac7c8f63a77dc730156c78e361273c3369959640120e5c54241f1a6f82d7b2`. Eligibility, dose/time, compounds, and panel follow acquire.py committed before inspection. Fixed hash order assigns first8 eligible contexts to final test, next8 to development validation, remaining contexts to training. The 2000-gene panel uses TRAINING controls only. Missing labels remain missing. Copying test count arrays during acquisition is not their use for training or selection.

Split controls A/B alternately within context and plate after stable well/sample sorting. Features use log2(1+CPM) treatment-minus-A; training/validation response labels use treatment-minus-B with the same physical treatment measurements. Baseline features use A. Normalize by supplied full library totals. All observations are pseudobulk, not individual cells. A separate scoring gate reads final test responses only after all predictions/parameters/code are frozen. Test control values and condition-availability metadata are allowed. Shared cell-village wells and assay batches remain limitations, not independent biological replications.

## Development models
Use the eight validation contexts, never final test, to develop/select:
- zero and equal-context source mean;
- global source shrinkage;
- baseline-distance source neighbors (1,3,8,all);
- cross-gene residual ridge using source-response PCA without target baseline;
- baseline-additive regression without drug-by-baseline interactions;
- bilinear interaction transport with baseline-coordinate/source-response-coordinate products and ridge regularization;
- baseline-kernel interpolation for each compound, shrunk toward the source mean.
Learn embeddings only on training contexts. Include no-baseline refitting and a fixed shuffled-context baseline control. The building blocks are classical; no novelty claim based only on naming or combining them. Additional architectures after validation inspection are explicitly developmental.

## Freeze/evaluation
Select one setting per family by validation macro MSE. Freeze all final predictions, parameters, code/data/protocol hashes and access logs before scoring eight test contexts once. Designated primary comparison: interaction vs no-baseline residual ridge. Do not turn later refitting into independent confirmation.
Primary metric: average MSE within context, then equal average across test contexts. Secondary: drug-centered error, mean/bias error, signed correlation, cosine identity retrieval with ties handled, predicted/observed RMS, top-magnitude gene errors (scoring-only diagnostic), and all per-context failures. Zero/scale baselines are mandatory. Near-zero predictions are not evidence of recovering nonzero mechanisms. Conditional compound resampling is descriptive; report context jackknife and individual contexts. No claim of study-level biological independence.

Metadata/schema adaptations must be recorded before outcome inspection; post-test analyses are post-hoc. No paid compute, scheduled research, automatic merging, or algorithmic-superiority release.
