# Publication status correction — 2026-09-08

This notice supersedes any implication in this directory's README that the remote estimator files are an exact copy of the executed scientific source.

A post-upload blob check found that two manually transferred files at commit eae4b56ddd5a59585067681b0bc96438622db5e6 are not byte-identical to the actual local files that produced the reported v6 results:

| File | Uploaded blob | Executed local blob |
|---|---|---|
| src/kernels.py | 1762b80da84ae11e600fb57290bbe25d6437a25c | c0f6011c801e975cb9cdb0256efaaeb9c36ebb17 |
| src/local.py | 76c3c2f600cc09aa2022ad8268bae946836956ad | 13e324c8e72df305d55943b026d733ce5c7c0191 |

The subsequent correction write was blocked by the platform. No workaround write has been attempted. Passing the limited remote core tests does not establish that this uploaded variant reproduces the reported scientific outputs.

Use the separately delivered, hash-verified local v6 source archive and numerical evidence for reproduction. The local executed files, frozen predictions, full rescore, and scientific results were not changed by the upload error. The acquisition workflow and timestamped protocol records remain useful provenance, but this branch is not a verified release and must not be merged as one. The repository's main model has not been replaced. This notice intentionally does not replace code or retry the blocked code update.
