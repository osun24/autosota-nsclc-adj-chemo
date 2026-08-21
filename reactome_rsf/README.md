# Version-pinned MSigDB Reactome RSF

This analysis uses the human MSigDB `C2:CP:REACTOME` collection from release
`2026.1.Hs` (released January 2026). It retrieves the gene-symbol GMT from the
Broad Institute release archive, verifies the pinned SHA-256
`5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341`, and
records the source URL, release documentation, retrieval time, set count, and
unique-gene count in `data/collection_manifest.json`.

The model uses the union of Reactome symbols found in both the training and
validation expression matrices plus all 19 prespecified clinical covariates.
Genes constant in training are dropped; missing predictors are median-imputed
using training values only. The RSF configuration is fixed in `config.py`;
only the validation cohort is used for evaluation, and it is not used for
feature or parameter tuning.
The sealed test CSV is rejected by path and is not part of this program.

Run from the repository root with the environment containing the packages in
`requirements.txt`:

```bash
python -m reactome_rsf.msigdb
python -m reactome_rsf.run
```

The first command is idempotent and validates an existing cache. The second
writes a timestamped directory under `reactome_rsf/runs/` containing the exact
feature list, training imputation medians, collection provenance, model
parameters, software versions, data hashes, validation Harrell C-index, and
60-month delta RMST. Delta RMST is RMST among validation patients whose observed
treatment matches the counterfactual RSF recommendation minus RMST among those
whose treatment does not match. Both metrics have deterministic 1,000-draw
patient-bootstrap 95% percentile intervals. Observed and counterfactual
validation risks are saved for audit; no training/OOB or test metric is reported.
Use `--save-model` only when a potentially large pickle artifact is required.

MSigDB collection files remain subject to the MSigDB license and terms of use.
