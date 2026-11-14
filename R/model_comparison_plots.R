# Extended Data Fig. 1: comparison of significant association sets.

MODEL_META <- tribble(
  ~model_code, ~model_label, ~model_short, ~model_order,
  "vae_age_residualized", "VAE: age-residualized", "VAE\nage-resid.", 1,
  "vae_age_included", "VAE: age-included", "VAE\nage-incl.", 2,
  "ols_age_residualized", "OLS: age-residualized", "OLS\nage-resid.", 3,
  "ols_age_included", "OLS: age-included", "OLS\nage-incl.", 4,
  "linear_age_interaction", "OLS: linear PRS x age", "Linear\nPRS x age", 5,
  "spline_age_interaction", "OLS: PRS x spline(age)", "Spline\nPRS x age", 6
)
MODEL_LEVELS <- MODEL_META$model_code
DOMAIN_LEVELS_COMPARISON <- c("AD", "PD", "ALS", "CVD", "T2D", "Longevity")

read_model_comparison_hits <- function(data_dir) {
  hits <- read_tsv_checked(file.path(data_dir, "manuscript", "model_comparison_hits.tsv"))
  ensure_columns(hits, c(
    "model_code", "domain", "prs_base", "omic_key", "age_group",
    "pair_id_any_age", "pair_id_exact_age"
  ), "model_comparison_hits.tsv")
  hits %>%
    filter(model_code %in% MODEL_LEVELS, domain %in% DOMAIN_LEVELS_COMPARISON) %>%
    mutate(
      model_code = factor(model_code, levels = MODEL_LEVELS),
      domain = factor(domain, levels = DOMAIN_LEVELS_COMPARISON),
      age_group = factor(age_group, levels = c("All", "<45", ">45"))
    )
}

make_model_count_table <- function(hits) {
  hits %>%
    distinct(model_code, domain, pair_id_any_age) %>%
    count(model_code, domain, name = "n_significant_pairs") %>%
    complete(
      model_code = factor(MODEL_LEVELS, levels = MODEL_LEVELS),
      domain = factor(DOMAIN_LEVELS_COMPARISON, levels = DOMAIN_LEVELS_COMPARISON),
      fill = list(n_significant_pairs = 0L)
    ) %>%
    left_join(MODEL_META, by = "model_code") %>%
    mutate(
      model_label = factor(model_label, levels = MODEL_META$model_label),
      model_short = factor(model_short, levels = MODEL_META$model_short)
    ) %>%
    arrange(model_order, domain)
}

make_overlap_table <- function(hits) {
  sets <- lapply(MODEL_LEVELS, function(code) {
    hits %>% filter(as.character(model_code) == code) %>%
      distinct(pair_id_any_age) %>% pull(pair_id_any_age)
  })
  names(sets) <- MODEL_LEVELS

  crossing(model_code_row = MODEL_LEVELS, model_code_col = MODEL_LEVELS) %>%
    rowwise() %>%
    mutate(
      n_row = length(sets[[model_code_row]]),
      n_col = length(sets[[model_code_col]]),
      n_common = length(intersect(sets[[model_code_row]], sets[[model_code_col]])),
      n_union = length(union(sets[[model_code_row]], sets[[model_code_col]])),
      jaccard = if_else(n_union > 0, n_common / n_union, 0)
    ) %>%
    ungroup() %>%
    left_join(
      MODEL_META %>% transmute(
        model_code_row = model_code,
        model_row_short = model_short,
        row_order = model_order
      ), by = "model_code_row"
    ) %>%
    left_join(
      MODEL_META %>% transmute(
        model_code_col = model_code,
        model_col_short = model_short,
        col_order = model_order
      ), by = "model_code_col"
    ) %>%
    mutate(
      model_row_short = factor(model_row_short, levels = rev(MODEL_META$model_short)),
      model_col_short = factor(model_col_short, levels = MODEL_META$model_short)
    )
}

make_age_context_table <- function(hits) {
  hits %>%
    filter(as.character(model_code) == "vae_age_included") %>%
    distinct(domain, age_group, pair_id_exact_age) %>%
    count(domain, age_group, name = "n_significant_pairs") %>%
    complete(
      domain = factor(DOMAIN_LEVELS_COMPARISON, levels = DOMAIN_LEVELS_COMPARISON),
      age_group = factor(c("All", "<45", ">45"), levels = c("All", "<45", ">45")),
      fill = list(n_significant_pairs = 0L)
    )
}

plot_model_counts <- function(counts) {
  ggplot(counts, aes(x = domain, y = fct_rev(model_label), fill = n_significant_pairs)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = n_significant_pairs), size = pt_to_mm(5.5),
              family = NATURE$font, fontface = "bold") +
    scale_fill_gradient(low = "white", high = NATURE_COLORS[["blue"]],
                        name = "Significant\npairs") +
    labs(title = "FDR-significant PRS-omic associations by disease domain",
         subtitle = "Unique base-PRS + omic-feature pairs; age strata collapsed",
         x = NULL, y = NULL) +
    coord_fixed(ratio = 0.72) +
    theme_nature(legend_position = "right", grid = FALSE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold"))
}

plot_model_overlap <- function(overlap) {
  ggplot(overlap, aes(x = model_col_short, y = model_row_short, fill = 100 * jaccard)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = n_common), size = pt_to_mm(5.2),
              family = NATURE$font, fontface = "bold") +
    scale_fill_gradient(low = "white", high = NATURE_COLORS[["green"]],
                        labels = function(x) paste0(signif(x, 2), "%"),
                        name = "Jaccard\nsimilarity") +
    labs(title = "Pairwise overlap of FDR-significant associations",
         subtitle = "Text: common base-PRS-feature pairs; fill: Jaccard similarity",
         x = NULL, y = NULL) +
    coord_fixed() +
    theme_nature(legend_position = "right", grid = FALSE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

plot_age_context_counts <- function(age_counts) {
  ggplot(age_counts, aes(x = domain, y = age_group, fill = n_significant_pairs)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = n_significant_pairs), size = pt_to_mm(5.4),
              family = NATURE$font, fontface = "bold") +
    scale_fill_gradient(low = "white", high = NATURE_COLORS[["vermillion"]],
                        name = "Significant\npairs") +
    labs(title = "Age context of age-included VAE associations",
         subtitle = "The same base pair is counted separately across age contexts",
         x = NULL, y = "Age context") +
    theme_nature(legend_position = "right", grid = FALSE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold"))
}
