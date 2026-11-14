#!/usr/bin/env Rscript

# Reproduce Extended Data Figures 1-4 as individual panels.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(file_arg)) else getwd()
repo_guess <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_guess, "R", "common.R"))
paths <- parse_path_args("extended_data")
source(file.path(paths$repo_root, "R", "discovery_plots.R"))
source(file.path(paths$repo_root, "R", "model_comparison_plots.R"))
source(file.path(paths$repo_root, "R", "aou_plots.R"))

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
D <- load_discovery_data(paths$data_dir)
vae_resid <- D$vae_residualized
vae_age <- D$vae_included
ols_resid <- D$ols_positive_controls %>% filter(model == "age_residualized")
ols_age <- D$ols_positive_controls %>% filter(model == "age_included")

# Extended Data Figure 1: domain-level model comparison.
comparison_hits <- read_model_comparison_hits(paths$data_dir)
counts <- make_model_count_table(comparison_hits)
overlap <- make_overlap_table(comparison_hits)
age_context <- make_age_context_table(comparison_hits)

save_plot(plot_model_counts(counts),
          "extended_data_fig01a_domain_model_counts", 105, 82, paths$output_dir)
save_plot(plot_model_overlap(overlap),
          "extended_data_fig01b_pairwise_overlap", 105, 95, paths$output_dir)
save_plot(plot_age_context_counts(age_context),
          "extended_data_fig01c_age_context", 88, 78, paths$output_dir)
write_tsv(counts %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "extended_data_fig01a_data.tsv"))
write_tsv(overlap %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "extended_data_fig01b_data.tsv"))
write_tsv(age_context %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "extended_data_fig01c_data.tsv"))

# Extended Data Figure 2: canonical PRS-biomarker positive controls.
p_positive_resid <- plot_canonical_benchmark(ols_resid, vae_resid, canonical_pairs)
p_positive_age <- plot_canonical_benchmark_age(ols_age, vae_age, canonical_pairs)
save_plot(p_positive_resid, "extended_data_fig02a_positive_controls_age_residualized",
          88, 70, paths$output_dir)
save_plot(p_positive_age, "extended_data_fig02b_positive_controls_age_included",
          88, 75, paths$output_dir)

# Extended Data Figure 3: longevity and healthspan associations (former main Figure 6).
longevity_domain <- prepare_domain_combined(
  DOMAIN_SETS$Longevity$prs, vae_resid, vae_age, repro_threshold = 50
)
save_plot(
  plot_domain_lollipop(longevity_domain, DOMAIN_SETS$Longevity$center_label, legend_position = "none"),
  "extended_data_fig03_longevity_lollipop", 180, 180, paths$output_dir
)

# Extended Data Figure 4: AoU candidate-set size and directional concordance.
aou <- read_aou_inputs(paths$data_dir)
method_summary <- prepare_aou_method_summary(aou$results)
save_plot(plot_aou_coverage(method_summary$coverage),
          "extended_data_fig04a_aou_coverage_concordance", 88, 80, paths$output_dir)
save_plot(
  plot_aou_domain_counts(
    method_summary$domains, "proteomics",
    "Directionally concordant protein associations by disease domain"
  ),
  "extended_data_fig04b_protein_domains", 88, 88, paths$output_dir
)
save_plot(
  plot_aou_domain_counts(
    method_summary$domains, "chemistries",
    "Directionally concordant clinical associations by disease domain"
  ),
  "extended_data_fig04c_chemistry_domains", 88, 88, paths$output_dir
)
write_tsv(method_summary$coverage,
          file.path(paths$output_dir, "extended_data_fig04a_data.tsv"))
write_tsv(method_summary$domains,
          file.path(paths$output_dir, "extended_data_fig04bc_data.tsv"))

message("Extended Data panels written to: ", paths$output_dir)
