# Orthogonal assay v11: continuation, not a restart

Date 2026-09-11. Inherit v10's conclusions and retained v7/RBF/calibration baselines. Do not rerun exposed L1000 model sweeps. The unresolved gate is whether public P100/GCP measurements provide genuinely useful, correctly matched functional information across preparations.

## Acquisition and evidence gate
Retrieve only linked public processed GCT files and official metadata. Retain URLs, acquisition times, byte sizes, content hashes and every failure. Initial scope: the exact MCF7 Plate29 3/6/24-hour PRM series, relevant GCP files, Level2/Level3 definitions and portal inventories. Download bounds do not imply any permission to use private data or paid compute. Scientific use requires checking assay mode, units, missingness, compound identity, dose, time, controls, culture/sample/plate/run/replicate identifiers and possible remeasurement of the same material. Never count PRM and DIA acquisitions of the same material as independent cultures.

## Questions to resolve before modeling
1. Are P100/GCP and existing L1000 matched by physical sample/preparation or only treatment labels? Does full query-drug isolation remain possible?
2. Which time points and cell backgrounds actually exist? Do not infer balanced coverage from a portal overview.
3. Did upstream normalization/filtering use held-out conditions? Prefer permitted lower-level data and frozen source-only filtering when possible; retain published versions for provenance, not as leak-free benchmarks.
4. Can a same-drug early measurement predict a late outcome? This is a distinct diagnostic task and must not substitute for unmeasured-drug probe prediction.
5. Is independent preparation identified, or only different MS runs/assays? Missing identity stays unknown.

No new model or primary predictive result is preregistered by this acquisition note. A separate fixed protocol will precede any new numeric predictive test, after metadata feasibility is established. Any condition already read for exploratory QC will be labelled as such. Archive counterevidence and record why a proposed test cannot be supported. No default-model replacement, automatic merge, journal submission, credentials or unattended research.
