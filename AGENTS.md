# Research working agreement

Begin each session by reading STATUS.md, docs/ROADMAP.md, open issues, and recent commits. The user's primary language is Chinese. Prefer substantive results to routine progress notifications. This repository stores public research decisions and execution records, not private reasoning transcripts or personal credentials.

## Scientific contract

- Original goal: quantitative chemical-response prediction in an unseen target context, using only allowed untreated baseline/source information. Current v0.7 code is an auxiliary phosphoprotein few-shot baseline, not an RNA zero-shot model.
- Distinguish unseen pairs, unseen contexts, unseen chemicals, mechanism/scaffold blocking and few-shot calibration. Count every target treatment/control/probe in the relevant budget. Do not silently change the task.
- Never turn missing labels into biological zeros. Audit units, normalization, gene identifiers, biological replicates, wells, donors, batches and shared controls before modeling.
- Preserve same-information simple baselines. Report absolute error, centered contrasts, mean/bias error, response-specific discrimination and worst-context failure. Do not turn an engineering speedup, extra data or a classical method into a novelty claim.
- Only target anchor outcomes may enter existing fit APIs. Target test outcomes must not determine preprocessing, tuning, probes or model selection. Any label-assisted/oracle diagnostic is explicitly non-deployable.
- Masking seeds and individual readouts do not create independent biological replicates. Shared controls and overlapping samples can induce correlated errors.
- Separate theory under stated assumptions, synthetic implementation checks, retrospective real data and independent/prospective evidence. Record failures; never promise SOTA or a publication tier.

## Execution contract

Use Python 3.13.5 and requirements.txt for baseline reproduction. Set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1. Run `python -m unittest discover -s tests -v` after preparing the auxiliary asset. Use `python -m tools.smoke --output runs/<new-id>` for the integration check.

Every new scientific run needs a new directory and a protocol committed before test evaluation. Record Git commit/dirty state, input SHA256, code hashes, seeds, split IDs, target information budget, parameters, metrics and failure status. Do not overwrite reports/v0.7 or docs/history/v0.7_import.json. Historical entrypoints overwrite outputs: use an isolated worktree until they are modernized.

Use branches/PRs for subsequent changes. Do not force-push, change repository privacy, add credentials, publish private data, or start paid compute without authorization. Retain third-party attribution. No automatic release of performance claims and no scheduled autonomous research. CI executes finite tests; a chat ending does not imply continuing investigation.
