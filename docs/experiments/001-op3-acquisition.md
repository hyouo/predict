# 001 — OP3 acquisition and metadata inventory

Date: 2026-09-06. This is a data-access/structure milestone, not a model experiment.

The repository's finite `Pinned OP3 acquisition audit` workflow succeeded: [run 34015151607](https://github.com/hyouo/predict/actions/runs/34015151607). Artifact 9983664163 was downloaded back through the GitHub connector, and the artifact ZIP and contained HDF5 file were independently rehashed locally. Exact provenance and expiration are in `reports/data/op3_acquisition.json`.

## What was actually read

The fixed Chem-PerturBridge OP3 processed asset is 23,700,750 bytes, SHA256 `8eaf7e63adc68029e88184726b8545bc3836c3ef2a25b9185dc4e8b819537e59`. It contains a CSR X matrix of shape **(1813, 5288)**, no additional layers and no `uns` entries. `tools.op3_inventory` reads structure and observation metadata only, not X expression values.

There are 1,813 unique sample IDs, 1,621 non-control and 192 control observations, four observed cell-type codes, six plate labels, and 138 non-control perturbagen labels. All stored treatment-time fields are 24 h; stored non-control doses are 1 uM, controls 0. There are 545 distinct non-control cell-type/perturbagen pairs. These are descriptive metadata counts, not evidence of chemical uniqueness, independent biological replication or full original-study coverage.

The file contains `psbulk_cells` and `psbulk_counts`. Do not call its rows individual cells or infer raw counts/normalization from the filename. No explicit `donor` column exists: plate/library or demographic fields cannot simply be relabeled as donors without a validated upstream mapping. `guide` and `stimulation` are entirely missing. Existing train/public_test/private_test/control annotations are preserved; no test expression was used for model development in this audit.

## Remaining gate

Issue #2 stays open. Verify the upstream preprocessing and gene-selection procedure, numerical measurement scale, plate/well/biological-replicate mapping, shared-control scheme and split provenance. `approved_for_training` remains false. Metadata-only acquisition does not validate RNA/ATAC performance and is not a new independent biological study.

## Reproduction

```bash
python -m tools.assets --asset op3
python -m tools.op3_inventory --output runs/op3-inventory.json
```

The output path must be new. Raw matrices are ignored by Git. The seven-day artifact is a transfer mechanism, not permanent data hosting; the upstream revision, checksum and capped downloader remain versioned. Five synthetic-format tests validate categorical missingness, Boolean controls, invalid indices and unsupported encodings; no synthetic values replace the public RNA asset.
