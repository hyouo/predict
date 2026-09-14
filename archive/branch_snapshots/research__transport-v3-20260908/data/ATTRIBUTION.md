# Data attribution and scope

The numerical data are downloaded into an ignored local working directory; they are not committed to this repository.

## Auxiliary signaling panel

Julio Saez-Rodriguez and colleagues, Comparing Signaling Networks between Normal and Transformed Hepatocytes Using Discrete Logical Models, Cancer Research 71(16):5400–5411 (2011), DOI 10.1158/0008-5472.CAN-10-4453.

Pinned source: https://github.com/cellnopt/cellnopt/blob/319c956e890ed7ae88c8e6e540de7061acfaae00/cno/datasets/HepG2CancerRes2011/MD-HepG2CancerRes2011.csv

Expected unchanged file: 89,549 bytes, Git blob 89861234c3300d46b2e5573416ab910f8fa32030, SHA256 da843e47b33b1a4485244efd59d61176fd7d0554089c3fe4785ebac12cdacc6c. The original repository license notice is retained in CELLNOPT_LICENSE.txt. No separate data-specific license was found in the previously inspected directory. That software notice is not represented as a license grant for the journal article.

The loader selects 360 ligand-present observations and applies a log signal/control transformation without changing the raw file. This is phosphoprotein data, not RNA, pseudobulk RNA or an OP3 subset. Time-code duplicates are not independent biological repeats.

## RNA assets

Chem-PerturBridge: https://huggingface.co/datasets/theislab/chem-perturbridge

The pinned OP3 and sci-Plex paths/hashes are in assets/manifest.json. Credit the distributor and the original Open Problems 2023 (GEO GSE279945) and Srivatsan et al. sci-Plex studies; retain dataset-specific upstream notices. Repository-level CC-BY-4.0 metadata does not eliminate the need to check original dataset terms. Optional short-lived Actions artifacts contain only the specified already-public asset and audit manifest, never personal account data.
