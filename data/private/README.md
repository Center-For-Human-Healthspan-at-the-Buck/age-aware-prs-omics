# Local controlled inputs - do not commit participant-level data

This directory is ignored except for this README and the three header-only
`*.placeholder.tsv` files. To rerun the study-specific analyte-wise linear
models, create authorized local files by copying and completing:

- `arivale_analysis_matrix.placeholder.tsv` -> `arivale_analysis_matrix.tsv`
- `arivale_feature_manifest.placeholder.tsv` -> `arivale_feature_manifest.tsv`
- `arivale_prs_manifest.placeholder.tsv` -> `arivale_prs_manifest.tsv`

The analysis matrix should contain one row per participant and columns for age,
covariates, continuous PRSs and omic features. The feature and PRS manifests map
analysis-matrix columns to public display names and omic layers. The generic
workflow is documented in `scripts/01_fit_linear_models.R` and the top-level
README.

The placeholder headers are illustrative rather than the exact locked study
schema. Do not add participant-level Arivale or All of Us data to the public
repository.
