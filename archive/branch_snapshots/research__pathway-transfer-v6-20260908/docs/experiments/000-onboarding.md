# 000 — Repository onboarding, 2026-09-06

Purpose: move the existing research from chat attachments into a versioned, runnable repository. This is not a new scientific experiment.

Input archive: perturbation_research_v0_7.zip, 26,687,094 bytes, SHA256 cb702dabab55317268467f8a1e67e6d252d6ecb3566ff1b8fa36082a5fb03682. Imported core/model/test/entrypoint map: ../history/v0.7_import.json. Full historical arrays and other old runners are not claimed to be in Git.

Local checks: 45 tests passed, comprising all 37 baseline tests plus 8 downloader/audit/non-overwrite tests. Fifteen Python files matched the baseline hashes. The new smoke run used one HepG2 target, four source contexts, 16 anchors, 24 queries and fixed seed 20260906; it tests successful fitting and full/blocked implementation equivalence only. Bootstrap MSE 0.19826231884624448 and stacking MSE 0.1958483616715042 are one integration case, not paper-ready benchmark estimates. Maximum implementation difference: 1.326716514427062e-14.

No raw numerical matrix is committed. Public assets have fixed source revisions, hashes and caps. Bad cached files are not overwritten; failed partial transfers are removed; missing data never triggers a synthetic fallback. Test fixtures are explicitly synthetic software fixtures.

CI and OP3 acquisition workflows, when permitted by the connected GitHub app, are finite CPU jobs with read-only repository permissions. No schedule, remote write-back, paid GPU or autonomous researcher is configured. Execution success is determined from actual logs. Acquisition does not approve assay semantics or constitute RNA training.
