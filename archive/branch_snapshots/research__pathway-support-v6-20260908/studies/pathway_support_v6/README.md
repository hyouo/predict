# Pathway-support v6: core research estimators and evidence checkpoint

This branch contains the actual support-aware kernel implementations, gene-local priors, additive-drift diagnostic and 18 self-contained core checks. Full data preparation, all 28 checks, retrospective execution records, selection grids, reports, and lossless prediction evidence are supplied in the conversation's v6 research/evidence archives. They are NOT all uploaded to this repository. The existing released model and main branch are unchanged.

## Main conclusion

On the already-used LINCS L1000 Phase II slice (30 contexts, 1,726 compound labels, 978 landmark genes), a same-gene baseline-weighted kernel has held-context MSE 0.301940 versus 0.306660 for the same-calibration full-baseline affine kernel. This is a 1.539% finite-cohort improvement, not SOTA evidence. In the tissue-label-excluded source setting, its centred drug-by-context interaction risk is positive (+0.005463), worse than zero interaction. The global Reactome prior does not beat the generic baseline kernel; gene-local Reactome produces only a small gain and is not better than the same-gene comparator.

Original vs separate mean/correction calibration: 0.319912 vs 0.306930. This calibration effect must not be credited to a biological mechanism or new architecture. Source variability weights are imperfect working proxies. All target treatment calibration budgets are zero; target untreated controls and source observations of the queried drugs remain available.

## Core usage

Run from this directory:

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s tests -v
```

In Python, fit `SupportTransport` or `LocalTransport` on aligned arrays: source effects `(contexts, drugs, genes)`, source baselines `(contexts, genes)`, then predict from target baselines `(targets, genes)`. Missing source profiles must be all-NaN; genes must already be aligned. These research classes do not provide the released software's identifier-level input contract.

```python
import numpy as np
from src.local import LocalTransport
# y_source, b_source and b_target are real, consistently preprocessed arrays.
# Parameters below illustrate the retained primary configuration, not universal defaults.
model = LocalTransport(np.eye(y_source.shape[-1])).fit(y_source, b_source)
prediction = model.predict(b_target, penalty=1.0, common=0.75, correction=0.5)
```

The affine solver is classical ordinary kriging / ridge regression with an unpenalized intercept. Global and local pathway features are unsigned memberships, not signed causal mechanisms or drug targets. `DriftTransport` is a post-score diagnostic, not a replacement selected by target performance.

## Evidence and limitations

Current source-only fit/validation split: 18 training and 6 validation contexts, then refit 24 and evaluate 6 fixed targets. The source stress setting uses 15/3, then 18 sources. The target dataset is previously exposed retrospective development. The 1,080/1,079 response rows differ slightly, but the 1,817 reference-matched interaction rows are identical. Tissue metadata include unknown and do not define validated lineage families.

120 first-stage configurations and 8 explicitly post-score drift configurations were retained. The 28 local mathematical/access checks passed; 24 selected prediction configurations were refit within 1e-12 (8 bitwise identical, max difference 1.33e-15). Independent primal SVD checks differ by at most 1.81e-15. A portability replay reproduced 14 local-model prediction arrays bitwise. These are numerical checks, not biological replications.

No State, MAP, PrePR-CT or other modern SOTA same-budget execution is claimed. Do not present this as universal superiority, independent biological validation, or a journal submission. Keep the branch/PR for review, without automatic merge.

## Sources

LINCS L1000 Phase II / GSE70138; processed theislab/chem-perturbridge revision 6c54c8eb0321cceff4f888c54e199077d055e20b. Original source SHA256 34d198df9eddac5b535a0408794caaf41073fa014c74429d1050e483576fb581. These are normalized landmark measurements, not RNA counts.

Reactome annotation data are CC0; Ragueneau et al., The Reactome Knowledgebase 2026, DOI 10.1093/nar/gkaf1223. Acquired GMT SHA256 89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b, 669 eligible deduplicated pathways over 743 measured genes. The acquisition workflow and four timing-specific protocols are recorded on this branch. No font files, credentials, or private data are included.
