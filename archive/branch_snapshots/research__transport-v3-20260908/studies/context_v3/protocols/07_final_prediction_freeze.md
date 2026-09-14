# Final prediction freeze: context_v3

Recorded before any analytical use or scoring of the eight final-context treatment response matrices. Full local freeze created 2026-09-08T01:38:35.316986+00:00; SHA256 `8d448edb5f3bcccd31fd62fd3a3177506cb87c5d241303f58be034c902050aab`. Data design SHA256 `fe295dcd297aaff1d3961074b03c30bc41dc39f003503bfd22d3f1e421538e17`.

## Design and information budget
50 eligible Tahoe contexts:34 fitted,8 development-validation,8 final. 256 compounds at5uM/24h,2000 output/input genes selected using training controls only. Query compounds are measured in source contexts. Each final context has14 A baseline-input control pseudobulks and14 B scoring controls, no treatment calibration labels. Final context order: CVCL_0480,CVCL_1055,CVCL_1635,CVCL_0152,CVCL_1056,CVCL_1724,CVCL_0428,CVCL_1547. This is a within-study cell-context holdout, not new-drug or independent-culture validation. Mechanical acquisition copied the count file, but its final treatment rows did not enter preprocessing, fitting, validation, or model selection.

## Validation-stage additions (not confirmatory results)
Baseline-conditioned gene-program decoder:32 response coordinates,8 baseline coordinates, two64-wide GELU layers and a zero-initialized residual head over the fixed response-ridge model. Source-state features combine local gene baseline, global baseline coordinates, and predefined products with gene program loadings. Three fixed seeds2026090809/810/811, each80epochs maximum, best checkpoint at five-epoch intervals using development contexts. An identically trained muted-input refit is mandatory. Added factorized-affine and additive-only decoder controls with initial seed and identical80epoch maximum. All these additions occurred while final outcomes remained sealed. Full-control-gene acquisition remains unexecuted after its write was blocked; all actual inputs remain the existing2000-gene panel.

Primary candidate: equal-weight conditioned three-seed ensemble. Primary matched comparison: independently fitted muted three-seed ensemble. Classical primary control: source-response rank32 ridge. Development MSEs are0.125173,0.131918,0.133018 respectively; these selected the next hypothesis but are NOT final results. All individual seeds and failed classical candidates remain recorded.

## Scoring contract
Primary: average gene/compound MSE within each final context, then equal context mean, relative to B controls. Secondary: drug-centered/bias errors, amplitude, drug retrieval, top100 observed-effect genes as scoring-only diagnostics, identical-prediction A-reference sensitivity, and two-way context/compound-centered interactions on complete cases. Fixed cyclic reassignment of target-baseline-conditioned predictions is a diagnostic only; it does not replace the muted refit. Report all contexts and failures; no training or selection after final scoring. Conditional compound resampling and context jackknife are descriptive, not independent-study inference. State/MAP/PrePR-CT have not been run under this contract; no SOTA claim.

## Frozen prediction SHA256 (all arrays8x256x2000)
```
zero c3e42a12444f04f76df0625927c77712554457fbc0ecced8295241b4a659089e
source_mean 43d812cd78f14f331767530f21fbef5f73c2bd54fc3c69bae96c38f22eeda315
source_shrink 4a0aa2172019ea803b98af360ff8f0faccb261d16dca9ee09bf69a11d36a3e2c
neighbors 5d5599735c6a428b6257e3c91c7be0c5bbeef8b043ed60dc0374acbaa1a1788c621
baseline_kernel a840bd9606b23b0b1bd8f6ac74c5219ad1a198c5392561b4f3ed87bd6d2364cb
response_ridge 4c779fd35f68f775929c0ad348399c3f3fdb8a39c2ecd38ffc11fd929084c34d
baseline_additive 2db008354daf7336e97d21013f1489ff855949db115e14faca26e6c557960eb5
bilinear 6074f168d4c9a18029eca005ba3b9507e2dc76a1ea47d3fc238b3601222a2b24
query_moment 9772031c10439561c0312570c0fbce1e5f1afaea6ec21586c249c334e2d5e739
moment_no_baseline dfc046cd769aba43af5b3c99056bb57355238003758c78466ac771522767cad1
local_response 96447385f534d57c5a3902a4879e99514fadd88bbc8469c31343f976c67149b2
local_additive 78f4fdd4751d94ca5ee96084ef88e776ca38b9ea9a85e3bc90c0cdd70e3afeb7
local_interaction a81772c25e41ae859eb31d4983b32997196088fb3fc00d3777899462a82b639c
conditioned_seed9 38e3e013199493bcbe42bda323e35208441fec2798746e539eec233411dd5ba1
muted_seed9 65126fbf46779bc791d9ed01ce33eb466b47b7bfc496c78ec2e0a5cc76514555
conditioned_seed10 13c82f401ff9bc6f9a2cfd61dea1a8b74cac36e9358987fc7e965f67f84b1b77
muted_seed10 8881eb35e6e573d8d23c91c7be0c5bbeef8b043ed60dc0374acbaa1a1788c621
conditioned_seed11 f3c6600253fc2b93a682f2266ac5863b1c1a9e750d5a70629f7543fc4c1cda3c
muted_seed11 7eedeeca75afb808d717028d2133289475c9089f9ecc74dfc1760ece16b32e53
affine_decoder de8c0a4344b2f89cdecd68d169f4f7d4517ddcdcd486aa03d61e2d49e80b28ac
additive_only_decoder 1fbaf35ca2dfafce59a7f664cf975b267a09dc0e0e2f2c327da49e8b6b915fed
conditioned_ensemble d582f8e46354f239c0355e3e59b9cbf2a2c092800e8169d3fa9cd22e5b95ba0e
wrong_baseline_matching 34e5a1f6647e798b9049c217ff8be6a7779d669cd66d8f8b8f7f50c7c9887860
muted_ensemble 06dbeb1d6f409199145b933c92b003ad9613764bcaf9f4c1fbd53e00b6533c60
```
The full machine-readable local freeze is authoritative for hashes and exact parameters. Core neural runner SHA256 `08e0b48267d23e21510ee63570f15733df2a779cc4b6ba6e36196af80641a98d`; inference `cc8a73054fe3bc601507617d0e8b44771540260e3aa2961e9ac6e1678d3c64fe`; access gate `756761901c7ce616a5ca277a91e6f0dea53b91010b853225bb729607d6a4edba`.25 local checks passed; all eight neural checkpoints reproduced held-out predictions exactly without treatment labels. One interrupted compression run was resumed only after verifying every existing array; incomplete output was retained separately. No scoring was permitted during that interruption.
