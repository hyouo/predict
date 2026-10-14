# v15: primary paired-resource gate after integrated v14

Date: 2026-09-14. Continue Issue #17. Stable main predictor unchanged. The earlier finite metadata-workflow write was blocked and was not retried or routed around. No remote acquisition or model training ran.

## New evidence

GEO GSE297075 (Vivo-seq, Fortmann et al., Cell Reports 2025, doi:10.1016/j.celrep.2025.116006) and native samples GSM8983405/GSM8983406 explicitly describe three cultures A/B/C, split into unstimulated and PMA/ionomycin conditions, then into two antibody panels. Twelve sample groups enter one 10x well. Two GSM entries denote GEX and ADT libraries, not two biological replicates. Independent donor provenance is not established by this design.

The native hash map is explicit: hashtag1=unstimulated, hashtag2=stimulated; hashtags3/4/5=cultures A/B/C in panel1, hashtags6/7/8=A/B/C in panel2. One positive hash in each role is expected; two total positive hashes are not automatically a doublet. Panel1 maps oYo1 to pSTAT3 Y705, oYo2 to p-p65 S536, and oYo3 to rabbit isotype. Panel2 maps oYo1 to pERK1/2 T202/Y204 and oYo3 to pFOS S32. These oligos cannot be renamed globally. Unassayed targets remain missing. The precise panel2 isotype correction remains unresolved.

This is an author-declared, native-decodable design, not a completed local barcode validation. The public H5AD/count files were not acquired because local DNS failed. No new expression matrix or individual biological pairing was reconstructed.

## Additional source checks

QuRIE-seq author QC tables total7449 cells before and6952 after filtering. The aIg and ibrutinib analysis subsets contain4754 and4658 cells and, by the published subset rules, share2460 cells. This is arithmetic on published summaries, not a raw-cell recount. GEO's6976-cell summary requires version reconciliation. The authors did not claim the two analysis subsets were independent tests.

FlexPlex already demonstrates paired-data-driven prediction of pRPS6 from post-perturbation RNA and follow-up experimental validation; it is direct prior art, not our result. VIPerturb-seq's multimodal cell-line pilot is a potential cross-context extension, but its genome-wide bins repeat all non-targeting cells and are not independent datasets. Flobak2019 provides explicit viability replicates, an auxiliary endpoint rather than a replacement for RNA prediction.

## Executed continuation

A local metadata decoder and fixed-prediction within-condition pairing diagnostic were implemented.22 software/mathematical checks passed, including exact finite-permutation verification, panel-dependent feature identity, structural missingness, namespaced barcodes and whole-culture holdout. These checks are not biological validations. New predictive models trained:0. New numerical cell matrices obtained:0. GSE101406 mother-preparation links authenticated:0.

## Next bounded step

Prioritize actual GSE297075 file acquisition and native barcode/hash audit, followed by a culture-held-out within-condition assay-association check. This is a supporting measurement benchmark, not query-excluded drug prediction, baseline-to-future inference, or cross-sequencing-batch evidence. Freeze preprocessing and simple matched comparators before numerical scoring. Retain the original query-excluded, independent-preparation drug-response objective as the separate scientific endpoint. Do not keep reinterpreting old P100 suffixes.

Primary entry points:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297075
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8983406
- https://vanbuggenum.github.io/QuRIE-seq_manuscript/QC.html
- https://zenodo.org/records/20544736
- https://zenodo.org/records/18460279

Full local code, registries and source-access limitations are delivered as conversation artifacts, not claimed to all be present in this branch. No Nature/SOTA claim, automatic merge, email or paid compute.
