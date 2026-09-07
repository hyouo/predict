# Paper validation protocol, 2026-09-07

Status: fixed before reading the quantitative SciPlex matrix in this research session. This is a prospective computational protocol, not an externally preregistered study. Primary model is the already delivered full-spectrum single-source partial-pooling estimator (local distribution perturb-predict 0.9.0; repository release numbering is separate). No new architecture or hyperparameter grid will be selected using external query outcomes. Record exact source hashes in the analysis package.

## External-study acquisition
Fixed source: theislab/chem-perturbridge revision 6c54c8eb0321cceff4f888c54e199077d055e20b, sciplex/srivatsan20_sciplex3_processed.h5ad. Expected SHA256 8ad35a20e44bf41cb030a0f063eb6a5727965b8f27d0eda4380bda0e294ccc4c. This is quantitative replication in a study independent of OP3, NOT a claim that no SciPlex-derived gene lists were ever inspected in earlier exploration. Respect CC BY 4.0 and original Sci-Plex/scPerturb/Chem-PerturBridge attribution.

## Decisions allowed after metadata audit
Determine cell-line labels, documented treatment replicates, control wells, common nominal doses and exposure time from metadata and upstream processing code, before outcome-dependent fitting. Select the highest common nonzero dose/time stratum (single dose/time per drug) supported by the metadata. Retain the published gene panel as fixed; disclose any upstream selection. Do not infer missing measurements as zero. If replicate structure or normalization is not supported, report that restriction and do not invent independence.

## Fixed analysis
Use all six directed transfers if three cell lines are present. Conditions/drugs are the split unit. Deterministically reserve approximately 40% of eligible drugs for queries using SHA256 order; remaining drugs form a calibration pool. Primary budget 16 distinct calibration drugs, nested budgets 4/8/32 as sensitivity. Five predeclared calibration orders (seeds 202609070 through 202609074); first seed is primary, repeated seeds are sensitivity, not biological replication. Exact IDs saved before test outcomes enter scoring.

Primary: source and target effect = log2(1 + count/library_total*1e6), treatment minus matched control, averaged over eligible recorded replicates. Preserve known original library totals; if only filtered-panel totals are available, explicitly limit estimand. Disjoint target controls for calibration and scoring when metadata permits; no artificial claim that splitting cells creates biological repeats. Matched-reference analysis is sensitivity, not a replacement result.

Compare frozen partial pooling with source-copy, zero, target-calibration mean, pooled scalar-plus-offset regression, zero-centred gene ridge and multivariate kernel ridge. All fitted baselines receive identical calibration drug outcomes and source information. Fixed penalties 0.1/1/10/100/1000 and offset fractions 0/0.25/0.5/0.75/1 where applicable. Hyperparameters selected only by full-drug leave-one-out within calibration. Queries never used for fit or model selection. Additional residual/cross-control variants are labelled exploratory.

Metrics: per-drug full-panel MSE (primary), centred MSE, bias MSE, effect correlations, drug retrieval and within-experiment replicate reproducibility. Report each directed transfer, equal-weight aggregate, per-query source-copy comparisons and the zero/calibration-mean baselines. Bootstrap drug IDs synchronously across transfers; distinguish conditional intervals from uncertainty over cell lines/studies. No gene-level pseudoreplication. No independent external-study/SOTA superiority claim from OP3 resplits.

Stress tests: where valid metadata exist, cross recorded replicates/wells rather than mixing source and target from the same physical unit; source-drug-label permutation; pseudocount sensitivity. Quantitative data provenance, failed attempts, exact split IDs, frozen predictions and all negative results are retained.

## Manuscript
Use official Springer Nature LaTeX authoring template with Nature reference style, plus compiled PDF, reproducible tables/figures, bibliography, methods, limitations and AI-assistance disclosure. No invented authors, affiliations, ethics approvals, funding or model superiority. No submission or editor contact is authorized. Paper-readiness is assessed from results, not promised from formatting. No scheduled jobs, automatic research loop, paid GPU, private data upload or repository permission changes.
