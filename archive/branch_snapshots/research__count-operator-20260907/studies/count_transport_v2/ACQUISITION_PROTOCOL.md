# Fixed new multi-context resource acquisition

This step does not change the completed SciPlex/OP3 results. The small-context baseline-conditioned neural operator has NOT outperformed the stronger scalar baseline. More contexts are needed to test its training assumptions; extra source data must be made available to all matched comparisons.

Download exactly the published Chem-PerturBridge Tahoe pseudobulk file at revision 6c54c8eb0321cceff4f888c54e199077d055e20b, SHA256 09ac7c8f63a77dc730156c78e361273c3369959640120e5c54241f1a6f82d7b2, 4,473,155,146 bytes. This is not the 100-million-cell raw single-cell matrix. Never turn missing pairs into zero.

Before selecting count values, hash-order eligible contexts with count-operator-20260907, reserve up to eight entire contexts for testing and up to four for validation, and use the remainder for training. Time is the metadata-modal time. Fix up to 192 compounds/doses by TRAINING-context metadata coverage (at least 80% of training contexts), not response strength. Select at most two existing treatment records per pair deterministically. Retain baseline and separate audit controls; differing rows/plates do not establish biological independence.

Gene panel: 5,000 highest average baseline-control count proportions across TRAINING contexts only, stable index tie-break. All other recorded counts form one explicit OTHER bin. Export test expression to a separate test_sealed folder. Test treatment arrays must not enter training, tuning, preprocessing, architecture choice or effect inspection. Hashes may be verified. Target baseline controls are allowed at inference; final predictions must be hashed before scoring target treatment counts.

Source data are upstream CC0 1.0; collection packaging/curation is CC BY 4.0. Preserve Tahoe Bio and Chem-PerturBridge notices. A bounded ordinary CPU Action acquires the public file with read-only repository permissions, no secrets, no GPU, no schedule and no autonomous conclusions. Delete the raw cache after export; retain only the subset and provenance as temporary Actions artifacts. Biological or algorithmic superiority is not implied by successful acquisition.
