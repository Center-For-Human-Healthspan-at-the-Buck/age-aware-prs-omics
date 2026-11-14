#!/usr/bin/env Rscript

# Validate all public aggregate inputs used by the figure-reproduction scripts.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(file_arg)) else getwd()
repo_guess <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_guess, "R", "common.R"))
paths <- parse_path_args("tables")
source(file.path(paths$repo_root, "R", "model_comparison_plots.R"))
source(file.path(paths$repo_root, "R", "aou_plots.R"))

required <- c(
  "vae/age_residualized.tsv", "vae/age_included.tsv",
  "metadata/metabolomics.tsv", "metadata/proteomics.tsv",
  "metadata/prs_correlation.tsv",
  "manuscript/ols_positive_controls.tsv",
  "manuscript/model_comparison_hits.tsv",
  "aou/replication_summary.tsv", "aou/demographics_summary.csv",
  "aou/qtl_annotations.csv",
  "aou/targets/ols_chemistries.tsv", "aou/targets/ols_proteins.tsv",
  "aou/targets/vae_chemistries.tsv", "aou/targets/vae_proteins.tsv"
)
missing <- required[!file.exists(file.path(paths$data_dir, required))]
if (length(missing) > 0) {
  stop("Missing required repository inputs: ", paste(missing, collapse = ", "), call. = FALSE)
}

static_assets <- c(
  "assets/fig01_workflow_schematic.pdf",
  "assets/fig01_workflow_schematic.png",
  "assets/fig06a_disease_association_workflow.pdf",
  "assets/fig06a_disease_association_workflow.png"
)
missing_assets <- static_assets[!file.exists(file.path(paths$repo_root, static_assets))]
if (length(missing_assets) > 0) {
  stop("Missing required static figure assets: ", paste(missing_assets, collapse = ", "), call. = FALSE)
}

D <- load_discovery_data(paths$data_dir)
comparison_hits <- read_model_comparison_hits(paths$data_dir)
aou <- read_aou_inputs(paths$data_dir)

for (dataset in list(D$ols_positive_controls, comparison_hits)) {
  recorded_n <- suppressWarnings(as.numeric(dataset$N))
  if (any(!is.na(recorded_n) & recorded_n != 1491)) {
    stop("An Arivale discovery-cohort N differs from 1,491.",
         call. = FALSE)
  }
}

if (!isTRUE(all.equal(D$prs_correlation, t(D$prs_correlation), tolerance = 1e-10))) {
  stop("PRS correlation matrix is not symmetric.", call. = FALSE)
}
if (any(abs(diag(D$prs_correlation) - 1) > 1e-10, na.rm = TRUE)) {
  stop("PRS correlation matrix diagonal is not equal to 1.", call. = FALSE)
}

summary <- bind_rows(
  tibble(
    input = c(
      "VAE age-residualized", "VAE age-included", "PRS correlation matrix",
      "OLS positive controls", "Model-comparison hits", "AoU disease-association summary",
      "AoU demographics", "AoU QTL annotations"
    ),
    rows = c(
      nrow(D$vae_residualized), nrow(D$vae_included), nrow(D$prs_correlation),
      nrow(D$ols_positive_controls), nrow(comparison_hits), nrow(aou$results),
      nrow(aou$demographics), nrow(aou$qtl)
    ),
    columns = c(
      ncol(D$vae_residualized), ncol(D$vae_included), ncol(D$prs_correlation),
      ncol(D$ols_positive_controls), ncol(comparison_hits), ncol(aou$results),
      ncol(aou$demographics), ncol(aou$qtl)
    )
  )
)

print(summary)
dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(summary, file.path(paths$output_dir, "input_validation_summary.tsv"))
message("Input validation completed successfully.")
