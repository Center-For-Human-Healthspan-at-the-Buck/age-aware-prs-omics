#!/usr/bin/env Rscript

# Reproduce Supplementary Figures 1-9 as individual files.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(file_arg)) else getwd()
repo_guess <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_guess, "R", "common.R"))
paths <- parse_path_args("supplementary")
source(file.path(paths$repo_root, "R", "discovery_plots.R"))

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
D <- load_discovery_data(paths$data_dir)
vae_resid <- D$vae_residualized
vae_age <- D$vae_included

# Supplementary Figures 1-2: cosine similarity of PRS perturbation profiles.
cos_resid <- build_cosine_matrix(vae_resid, stratified = FALSE)
cos_age <- build_cosine_matrix(vae_age, stratified = TRUE)

p_cos_resid <- pheatmap(
  cos_resid,
  main = "Cosine similarity of PRS (age-residualized)",
  cluster_rows = TRUE, cluster_cols = TRUE,
  color = colorRampPalette(c(NATURE_COLORS[["blue"]], "white", EFFECT_COLORS[["Negative"]]))(100),
  fontsize = NATURE$small_pt, fontsize_row = 3.8, fontsize_col = 3.8,
  silent = TRUE
)
p_cos_age <- pheatmap(
  cos_age,
  main = "Cosine similarity of PRS (age-included)",
  cluster_rows = TRUE, cluster_cols = TRUE,
  color = colorRampPalette(c(NATURE_COLORS[["blue"]], "white", EFFECT_COLORS[["Negative"]]))(100),
  fontsize = NATURE$small_pt, fontsize_row = 3.8, fontsize_col = 3.8,
  silent = TRUE
)
save_pheatmap(p_cos_resid, "supp_fig01_prs_cosine_age_residualized", 180, 180, paths$output_dir)
save_pheatmap(p_cos_age, "supp_fig02_prs_cosine_age_included", 180, 180, paths$output_dir)

# Supplementary Figure 3: Mantel comparison of genetic and omic dissimilarity.
p_mantel <- plot_mantel_relationship(D$prs_correlation, cos_resid, cos_age)
save_plot(p_mantel, "supp_fig03_mantel_genetic_omic_dissimilarity", 88, 75, paths$output_dir)

# Supplementary Figures 4-9: complete domain-specific forest plots.
forest_map <- list(
  AD = "supp_fig04_ad_forest",
  PD = "supp_fig05_pd_forest",
  ALS = "supp_fig06_als_forest",
  T2D = "supp_fig07_t2d_forest",
  CVD = "supp_fig08_cvd_forest",
  Longevity = "supp_fig09_longevity_forest"
)

for (domain_key in names(forest_map)) {
  definition <- DOMAIN_SETS[[domain_key]]
  panel_data <- prepare_domain_combined(
    definition$prs, vae_resid, vae_age, repro_threshold = 50
  )
  n_features <- n_distinct(panel_data$Feature)
  axis_pt <- switch(domain_key, T2D = 3.8, Longevity = 2.9, NATURE$small_pt)
  plot <- plot_domain_forest(
    panel_data,
    title = paste0(definition$title, " PRS associations"),
    axis_text_y_pt = axis_pt
  )
  save_plot(
    plot, forest_map[[domain_key]], 180,
    min(NATURE$max_height_mm, max(80, 22 + 3.2 * n_features)),
    paths$output_dir
  )
  write_tsv(
    panel_data %>% mutate(across(where(is.factor), as.character)),
    file.path(paths$output_dir, paste0(forest_map[[domain_key]], "_data.tsv"))
  )
}

message("Supplementary-figure panels written to: ", paths$output_dir)
