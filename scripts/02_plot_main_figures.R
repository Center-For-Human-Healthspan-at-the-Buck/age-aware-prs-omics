#!/usr/bin/env Rscript

# Reproduce the individual panels used in main Figures 1-6.
# Each panel is exported separately for multipanel assembly.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(file_arg)) else getwd()
repo_guess <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_guess, "R", "common.R"))
paths <- parse_path_args("main")
source(file.path(paths$repo_root, "R", "discovery_plots.R"))
source(file.path(paths$repo_root, "R", "aou_plots.R"))

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
D <- load_discovery_data(paths$data_dir)
vae_resid <- D$vae_residualized
vae_age <- D$vae_included

# Figure 1: static workflow schematic.
for (extension in c("pdf", "png")) {
  copy_static_panel(
    file.path(paths$repo_root, "assets", paste0("fig01_workflow_schematic.", extension)),
    file.path(paths$output_dir, paste0("fig01_workflow_schematic.", extension))
  )
}

# Figure 2: age representation and association composition.
fig2_resid_data <- filter_manuscript_prs(vae_resid, age_included = FALSE)
fig2_age_data <- filter_manuscript_prs(vae_age, age_included = TRUE)
chord_resid <- make_chord_data(fig2_resid_data, stratified = FALSE, repro_threshold = 50)
chord_age <- make_chord_data(fig2_age_data, stratified = TRUE, repro_threshold = 50)

save_base_plot(
  function() plot_chord_counts(chord_resid, title = "Age-residualized VAE model"),
  "fig02a_age_residualized_chord", 88, 88, paths$output_dir
)
save_base_plot(
  function() plot_chord_counts(chord_age, title = "Age-included VAE model"),
  "fig02b_age_included_chord", 88, 88, paths$output_dir
)

# The top-feature panels use the full PRS discovery set; the chord panels remain
# restricted to the aging-related domains.
p_fig2c <- plot_top_omic_jitter(
  vae_resid,
  NULL,
  n_features = 10
)
p_fig2d <- plot_top_omic_jitter(
  vae_age,
  NULL,
  n_features = 10
)
save_plot(p_fig2c, "fig02c_age_residualized_top_features", 180, 90, paths$output_dir)
save_plot(p_fig2d, "fig02d_age_included_top_features", 180, 90, paths$output_dir)
write_tsv(chord_resid, file.path(paths$output_dir, "fig02a_chord_counts.tsv"))
write_tsv(chord_age, file.path(paths$output_dir, "fig02b_chord_counts.tsv"))

# Figure 3: shared neurodegeneration fingerprint and complete PD/ALS matrices.
p_neuro <- plot_neuro_fingerprint(vae_resid, vae_age)
save_plot(p_neuro, "fig03a_neurodegeneration_fingerprint", 160, 125, paths$output_dir)

pd_domain <- prepare_domain_combined(
  DOMAIN_SETS$PD$prs, vae_resid, vae_age, repro_threshold = 50
)
als_domain <- prepare_domain_combined(
  DOMAIN_SETS$ALS$prs, vae_resid, vae_age, repro_threshold = 50
)
pd_matrix <- prepare_neuro_all_feature_matrix(pd_domain, "PD")
als_matrix <- prepare_neuro_all_feature_matrix(als_domain, "ALS")
shared_max <- max(c(pd_matrix$data$abs_est, als_matrix$data$abs_est), na.rm = TRUE)
if (!is.finite(shared_max) || shared_max <= 0) shared_max <- 1

p_pd <- plot_neuro_all_feature_matrix(
  pd_matrix$data, "Parkinson's disease", size_limits = c(0, shared_max), show_legend = FALSE,axis_text_y_pt = 6
)
p_als <- plot_neuro_all_feature_matrix(
  als_matrix$data, "Amyotrophic lateral sclerosis", size_limits = c(0, shared_max), show_legend = FALSE,axis_text_y_pt = 6
)
save_plot(p_pd, "fig03b_pd_feature_by_model_matrix", 110,
          min(210, max(125, 28 + 4.15 * pd_matrix$n_features)), paths$output_dir)
save_plot(p_als, "fig03c_als_feature_by_model_matrix", 110,
          min(210, max(125, 28 + 4.15 * als_matrix$n_features)), paths$output_dir)
legend_source <- plot_neuro_all_feature_matrix(
  pd_matrix$data, NULL, size_limits = c(0, shared_max), show_legend = TRUE
)
legend_grob <- extract_plot_legend(legend_source)
if (!is.null(legend_grob)) {
  save_base_plot(function() grid::grid.draw(legend_grob),
                 "fig03bc_shared_legend", 180, 24, paths$output_dir)
}
write_tsv(pd_matrix$data %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "fig03b_panel_data.tsv"))
write_tsv(als_matrix$data %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "fig03c_panel_data.tsv"))

# Figure 4: AD lollipop and age-context panel.
ad_domain <- prepare_domain_combined(
  DOMAIN_SETS$AD$prs, vae_resid, vae_age, repro_threshold = 50
)
p_ad <- plot_domain_lollipop(ad_domain, DOMAIN_SETS$AD$center_label, legend_position = "none")
save_plot(p_ad, "fig04a_ad_lollipop", 180, 180, paths$output_dir)
ad_age <- prepare_ad_age_strata_plot_data(
  vae_age, DOMAIN_SETS$AD$prs, repro_threshold = 50, label_max_chars = 40
)
p_ad_age <- plot_ad_age_strata(ad_age, title = NULL, show_legend = FALSE)
save_plot(p_ad_age, "fig04b_ad_age_context", 76, 126, paths$output_dir)
ad_age_legend_source <- plot_ad_age_strata(
  ad_age, title = NULL, show_legend = TRUE
)
ad_age_legend <- extract_plot_legend(ad_age_legend_source)
if (!is.null(ad_age_legend)) {
  save_base_plot(function() grid::grid.draw(ad_age_legend),
                 "fig04b_legend", 180, 24, paths$output_dir)
}
write_tsv(ad_age %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "fig04b_panel_data.tsv"))

# Figure 5: T2D and CVD domain panels.
t2d_domain <- prepare_domain_combined(
  DOMAIN_SETS$T2D$prs, vae_resid, vae_age, repro_threshold = 50
)
cvd_domain <- prepare_domain_combined(
  DOMAIN_SETS$CVD$prs, vae_resid, vae_age, repro_threshold = 50
)
save_plot(
  plot_domain_lollipop(t2d_domain, DOMAIN_SETS$T2D$center_label, legend_position = "none"),
  "fig05a_t2d_lollipop", 180, 180, paths$output_dir
)
save_plot(
  plot_domain_lollipop(cvd_domain, DOMAIN_SETS$CVD$center_label, legend_position = "none"),
  "fig05b_cvd_lollipop", 180, 180, paths$output_dir
)
cvd_repeated <- make_domain_replicated_table(DOMAIN_SETS$CVD$prs, vae_resid, vae_age)
p_cvd_bubble <- plot_domain_bubble(cvd_repeated, title = "CVD: associations recurring across age models")
save_plot(p_cvd_bubble, "fig05c_cvd_recurring_associations", 88,
          min(125, max(92, 50 + 2.4 * n_distinct(cvd_repeated$omic_clean))), paths$output_dir)
write_tsv(cvd_repeated %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "fig05c_panel_data.tsv"))

# Figure 6: static disease-association workflow, ancestry pies, and VAE-prioritized forests.
for (extension in c("pdf", "png")) {
  copy_static_panel(
    file.path(paths$repo_root, "assets", paste0("fig06a_disease_association_workflow.", extension)),
    file.path(paths$output_dir, paste0("fig06a1_disease_association_workflow.", extension))
  )
}
aou <- read_aou_inputs(paths$data_dir)
aou_pairs <- prepare_aou_vae_pairs(aou$results, aou$qtl)
ancestry <- prepare_ancestry_data(aou$demographics)
save_plot(plot_ancestry_pie(ancestry, "Chemistries", show_legend = FALSE),
          "fig06a2_chemistry_ancestry", 88, 78, paths$output_dir)
save_plot(plot_ancestry_pie(ancestry, "Proteomics", show_legend = FALSE),
          "fig06a3_proteomics_ancestry", 88, 78, paths$output_dir)
draw_ancestry_legend <- function() {
  grid::grid.newpage()

  column_widths <- c(14, rep((75 - 14) / 3, 3))

  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 2,
        ncol = 4,
        widths = grid::unit(column_widths, "mm"),
        heights = grid::unit(c(6.5, 6.5), "mm")
      )
    )
  )

  grid::grid.text(
    "Genetic\nancestry",
    vp = grid::viewport(layout.pos.row = 1:2, layout.pos.col = 1),
    just = "center",
    gp = grid::gpar(
      fontfamily = NATURE$font,
      fontsize = NATURE$small_pt,
      fontface = "bold",
      lineheight = 0.9
    )
  )

  legend_items <- data.frame(
    row = c(1, 2, 1, 2, 1, 2),
    column = c(2, 2, 3, 3, 4, 4),
    label = c(
      "Admixed\nAmerican",
      "African",
      "East\nAsian",
      "European",
      "Middle\nEastern",
      "South\nAsian"
    ),
    colour = c(
      "#CA699D",
      "#E66101",
      "#0571B0",
      "#B2B2B2",
      "#E6AB02",
      "#1B9E77"
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(legend_items))) {
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = legend_items$row[i],
        layout.pos.col = legend_items$column[i]
      )
    )

    grid::grid.rect(
      x = grid::unit(0.8, "mm"),
      y = grid::unit(0.5, "npc"),
      width = grid::unit(3.6, "mm"),
      height = grid::unit(3.6, "mm"),
      just = c("left", "center"),
      gp = grid::gpar(
        fill = legend_items$colour[i],
        col = NA
      )
    )

    grid::grid.text(
      legend_items$label[i],
      x = grid::unit(9.5, "mm"),
      y = grid::unit(0.5, "npc"),
      just = "center",
      gp = grid::gpar(
        fontfamily = NATURE$font,
        fontsize = 4.8,
        lineheight = 0.9
      )
    )

    grid::popViewport()
  }

  grid::popViewport()
}

save_base_plot(
  draw_ancestry_legend,
  "fig06a_ancestry_legend",
  75,
  13,
  paths$output_dir
)

forest_specs <- tribble(
  ~omic, ~relation, ~stub, ~title,
  "chemistries", "same_direction", "fig06b1_clinical_directional", "Clinical features",
  "proteomics", "same_direction", "fig06b2_protein_directional", "Proteins",
  "chemistries", "opposite_direction", "fig06c1_clinical_discordant", "Clinical features",
  "proteomics", "opposite_direction", "fig06c2_protein_discordant", "Proteins"
)
purrr::pwalk(forest_specs, function(omic, relation, stub, title) {
  n_rows <- aou_pairs %>% filter(
    omic_type == omic, modelled, !is.na(replication_q), replication_q < 0.10,
    direction_relation == relation
  ) %>% nrow()
  plot <- plot_aou_forest(aou_pairs, omic, relation, title)
  save_plot(plot, stub, 88, min(210, max(105, 18 + 5.0 * n_rows)), paths$output_dir)
})

legend_plot <- plot_aou_forest(
  aou_pairs, "proteomics", "same_direction", NULL, show_legend = TRUE
)
aou_legend <- extract_plot_legend(legend_plot)
if (!is.null(aou_legend)) {
  save_base_plot(function() grid::grid.draw(aou_legend),
                 "fig06_shared_legend", 180, 28, paths$output_dir)
}
write_tsv(aou_pairs %>% mutate(across(where(is.factor), as.character)),
          file.path(paths$output_dir, "fig06_pair_level_data.tsv"))

message("Main-figure panels written to: ", paths$output_dir)
