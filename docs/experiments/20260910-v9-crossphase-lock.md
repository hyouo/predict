# v9 cross-phase lock

Both registered Phase II fitting runs have completed. No Phase I treatment value has been supplied to a local fitting, selection or scoring program. The public acquisition worker performed only metadata filtering, byte/hash verification and numeric export. Phase I artifacts are run34443905550: metadata10139085883 and numeric10139087382. Source file2,165,954,369bytes has the registered SHA256. The metadata archive has been read; numeric archive has not yet been downloaded for analysis.

All seven predictor families are fixed: ordinary RBF; identity-nugget v7 response prior; count-only nugget; source-context-average nugget; compound-average nugget; context-by-compound nugget; within-compound shuffled nugget. Every family can select the identity fallback on the original Phase II validation contexts. No new target-guided tuning, filtering or calibration will occur.

Lock hashes:
- PRE_PHASE1_FIT_LOCK.json: bf1e0bdd4b28552e4af395e66c147e155aa3be81e6fafc0f3247b2277d36d3f2
- Phase II primary fitting freeze:11e03fef362075788bc3ad36b9897e390cfee2e00e97502a14f7d940f1641e9b
- Phase II tissue-excluded fitting freeze:6962292e99ec3e902dd2fdc758c05a00c85a1e98a4a7bcecb7af419994deadbe
- frozen transport source:1bfe45f0737f776ecc560b3d58bd74ed89ae518b757f13d61659dd7b2d303b1b
- PHASE1_METADATA_PLAN.json:5a164de63113e82f53c4caa62617b6980a488a968ab66b48b3b9e631d5218431

Metadata eligibility retains3071 profiles in7 targets: CVCL_0023(656),CVCL_2235(678),CVCL_5136(396),CVCL_U602(455),CVCL_UK07(89),PHH(352),SKB(445). Every profile has>=2records and>=3 original Phase II source contexts. CVCL_0336 fails the registered>=50supported-drug threshold and is excluded without inspecting expression. All978 genes match exactly including order. No sample ID or physical plate ID overlaps with the available Phase II input slice. All targets exclude the24 fitting context IDs. Three are historical Phase II held targets; four exact contexts are new to the30-context Phase II cohort. Report both these metadata-defined strata in addition to the whole cohort.

The cross-phase main comparison uses the fixed ordinary Phase II models trained on24 source contexts, not models retuned for the new targets' tissue labels. Source responses remain from Phase II; Phase I contributes untreated A controls only before prediction freeze and treated/B-reference rows only after freeze. This jointly challenges cell-context and collection-phase transport; it is not an independent laboratory, time-forward study or randomized biological validation. Preserve failures, all eligible contexts and all model results.