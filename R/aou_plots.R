# Main Figure 6 and Extended Data Figure 4 from locked aggregate All of Us outputs.

AOU_DOMAIN_LEVELS <- c(
  "AD_dementia", "Parkinsons", "ALS_MND", "CAD_coronary",
  "Heart_failure", "CVD_composite", "T2D"
)
AOU_DOMAIN_LABELS <- c(
  AD_dementia = "AD/dementia",
  Parkinsons = "PD",
  ALS_MND = "ALS",
  CAD_coronary = "CAD",
  Heart_failure = "Heart failure",
  CVD_composite = "Composite CVD",
  T2D = "T2D"
)
AOU_DOMAIN_COLORS <- c(
  AD_dementia = "#CC79A7", Parkinsons = "#56B4E9", ALS_MND = "#0072B2",
  CAD_coronary = "#D55E00", Heart_failure = "#E69F00",
  CVD_composite = "#8C564B", T2D = "#009E73"
)

first_nonmissing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & str_squish(x) != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}
first_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else x[[1]]
}

base_prs_name <- function(x) strip_age_suffix(x)
canonical_qtl_feature_key <- function(x) canonical_omic_key(x)

collapse_qtl_type <- function(x) {
  x <- unique(na.omit(as.character(x)))
  has_e <- any(x %in% c("eQTL", "both"))
  has_p <- any(x %in% c("pQTL", "both"))
  if (has_e && has_p) return("both")
  if (has_e) return("eQTL")
  if (has_p) return("pQTL")
  if (length(x) > 0 && any(x == "none")) return("none")
  "Not annotated"
}

read_aou_inputs <- function(data_dir) {
  results <- read_tsv_checked(file.path(data_dir, "aou", "replication_summary.tsv"))
  demographics <- readr::read_csv(
    file.path(data_dir, "aou", "demographics_summary.csv"),
    show_col_types = FALSE
  )
  qtl <- readr::read_csv(
    file.path(data_dir, "aou", "qtl_annotations.csv"),
    show_col_types = FALSE
  )
  ensure_columns(results, c(
    "discovery_framework", "discovery_model", "omic_type", "analysis_domain",
    "feature_display", "uniprot", "gene_name", "prs_names", "discovery_direction",
    "replication_test_family", "beta_case", "ci_low", "ci_high",
    "discovery_concordant_p", "discovery_concordant_q", "direction_relation",
    "discovery_concordant_status"
  ), "AoU replication summary")
  list(results = results, demographics = demographics, qtl = qtl)
}

prepare_aou_vae_pairs <- function(results, qtl_table) {
  rows <- results %>%
    filter(
      discovery_framework == "VAE",
      discovery_model %in% c("VAE_age", "VAE_resid"),
      replication_test_family == "disease main effect",
      omic_type %in% c("proteomics", "chemistries"),
      analysis_domain %in% AOU_DOMAIN_LEVELS
    ) %>%
    mutate(
      feature_identity = case_when(
        omic_type == "proteomics" & !is.na(uniprot) & str_squish(uniprot) != "" ~
          paste0("uniprot:", str_squish(uniprot)),
        TRUE ~ paste0("label:", canonical_omic_key(feature_display))
      ),
      beta_case = suppressWarnings(as.numeric(beta_case)),
      ci_low = suppressWarnings(as.numeric(ci_low)),
      ci_high = suppressWarnings(as.numeric(ci_high)),
      replication_p = suppressWarnings(as.numeric(discovery_concordant_p)),
      replication_q = suppressWarnings(as.numeric(discovery_concordant_q))
    )

  pairs <- rows %>%
    group_by(omic_type, analysis_domain, feature_identity) %>%
    summarise(
      feature_display = first_nonmissing(feature_display),
      gene_name = first_nonmissing(gene_name),
      uniprot = first_nonmissing(uniprot),
      prs_names = paste(sort(unique(na.omit(prs_names))), collapse = "; "),
      discovery_models = paste(sort(unique(discovery_model)), collapse = "; "),
      discovery_direction = first_nonmissing(discovery_direction),
      beta_case = first_numeric(beta_case),
      ci_low = first_numeric(ci_low),
      ci_high = first_numeric(ci_high),
      replication_p = first_numeric(replication_p),
      replication_q = first_numeric(replication_q),
      direction_relation = first_nonmissing(direction_relation),
      modelled = any(discovery_concordant_status == "ok", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      domain_label = unname(AOU_DOMAIN_LABELS[analysis_domain]),
      display_feature = case_when(
        omic_type == "proteomics" & !is.na(gene_name) & str_squish(gene_name) != "" ~ gene_name,
        TRUE ~ feature_display
      )
    )

  qtl_lookup <- qtl_table %>%
    transmute(
      prs_key = base_prs_name(feature_a_original_name),
      feature_key = canonical_qtl_feature_key(omic_original_name),
      qtl_type = as.character(qtl_type)
    ) %>%
    filter(!is.na(prs_key), prs_key != "", !is.na(feature_key), feature_key != "")

  annotations <- lapply(seq_len(nrow(pairs)), function(i) {
    row_i <- pairs[i, , drop = FALSE]
    prs_keys <- str_split(as.character(row_i$prs_names), "\\s*;\\s*")[[1]] %>%
      base_prs_name() %>% unique()
    feature_keys <- unique(c(
      canonical_qtl_feature_key(row_i$feature_display),
      canonical_qtl_feature_key(row_i$gene_name)
    ))
    matched <- qtl_lookup %>% filter(prs_key %in% prs_keys, feature_key %in% feature_keys)
    tibble(qtl_type = collapse_qtl_type(matched$qtl_type))
  }) %>% bind_rows()

  bind_cols(pairs, annotations) %>%
    mutate(qtl_type = factor(qtl_type, levels = names(QTL_SHAPES)))
}

parse_percent <- function(x) {
  suppressWarnings(as.numeric(str_extract(as.character(x), "[0-9.]+(?=%)")))
}

prepare_ancestry_data <- function(demographics) {
  ancestry <- demographics %>%
    filter(str_detect(characteristic, "ancestry, n")) %>%
    mutate(ancestry = str_remove(characteristic, ", n \\(.*$")) %>%
    select(ancestry, Chemistries, Proteomics) %>%
    pivot_longer(c(Chemistries, Proteomics), names_to = "modality", values_to = "value") %>%
    mutate(percent = parse_percent(value))

  sample_sizes <- demographics %>%
    filter(characteristic == "Unique participants, N") %>%
    pivot_longer(c(Chemistries, Proteomics), names_to = "modality", values_to = "N") %>%
    mutate(N = suppressWarnings(as.numeric(str_remove_all(N, ","))))

  ancestry %>% left_join(sample_sizes, by = "modality")
}

plot_ancestry_pie <- function(ancestry_df, modality_name, show_legend = TRUE) {
  plot_df <- ancestry_df %>% filter(modality == modality_name)
  n_value <- unique(plot_df$N)
  ggplot(plot_df, aes(x = "", y = percent, fill = ancestry)) +
    geom_col(width = 1, colour = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = c(
        "African ancestry" = "#E66101", "Admixed American ancestry" = "#CA699D",
        "East Asian ancestry" = "#0571B0", "European ancestry" = "#B2B2B2",
        "Middle Eastern ancestry" = "#E6AB02", "South Asian ancestry" = "#1B9E77"
      ),
      labels = c(
        "African ancestry" = "African",
        "Admixed American ancestry" = "Admixed American",
        "East Asian ancestry" = "East Asian",
        "European ancestry" = "European",
        "Middle Eastern ancestry" = "Middle Eastern",
        "South Asian ancestry" = "South Asian"
      ),
      name = "Genetic ancestry"
    ) +
    labs(title = modality_name, subtitle = paste0("N = ", scales::comma(n_value))) +
    theme_void(base_family = NATURE$font) +
    theme(
      plot.title = element_text(size = NATURE$title_pt + 2.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = NATURE$base_pt, hjust = 0.5),
      legend.title = element_text(size = NATURE$small_pt, face = "bold"),
      legend.text = element_text(size = 4.8),
      legend.position = if (isTRUE(show_legend)) "bottom" else "none"
    ) +
    guides(fill = guide_legend(ncol = 3))
}

plot_aou_forest <- function(pair_level, omic_type_value, relation,
                            title, q_cutoff = 0.10, show_legend = FALSE) {
  plot_df <- pair_level %>%
    filter(
      omic_type == omic_type_value,
      modelled,
      !is.na(replication_q), replication_q < q_cutoff,
      direction_relation == relation,
      is.finite(beta_case), is.finite(ci_low), is.finite(ci_high)
    ) %>%
    mutate(
      feature_domain = paste0(display_feature, " – ", domain_label),
      q_label = paste0("q = ", format(replication_q, digits = 2, scientific = TRUE)),
      feature_domain = fct_reorder(feature_domain, beta_case)
    )

  if (nrow(plot_df) == 0) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "No associations") + theme_void())
  }

  x_limits <- pad_continuous_limits(c(plot_df$ci_low, plot_df$ci_high),
                                    lower_mult = 0.10, upper_mult = 0.35,
                                    min_pad = 0.08, include_zero = TRUE)
  q_x <- x_limits[[2]] - 0.02 * diff(x_limits)

  ggplot(plot_df, aes(x = beta_case, y = feature_domain)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.28,
                   colour = "grey45") +
    geom_point(
      aes(shape = qtl_type, fill = beta_case), size = 2.0,
      colour = if (relation == "same_direction") "grey25" else EFFECT_COLORS[["Negative"]],
      stroke = 0.35
    ) +
    geom_text(aes(x = q_x, label = q_label), hjust = 1,
              size = pt_to_mm(4.8), family = NATURE$font) +
    scale_shape_manual(values = QTL_SHAPES, name = "QTL annotation", drop = FALSE) +
    scale_fill_gradient2(low = "#67A9CF", mid = "white", high = "#EF8A62",
                         midpoint = 0, name = "AoU beta") +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(title = title, x = "AoU disease-control beta (95% CI)", y = NULL) +
    theme_nature(legend_position = if (isTRUE(show_legend)) "bottom" else "none", grid = FALSE) +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey92"),
      axis.text.y = element_text(size = 4.8),
      plot.margin = margin(4, 14, 4, 4)
    )
}

prepare_aou_method_summary <- function(results, q_cutoff = 0.10) {
  base <- results %>%
    filter(
      replication_test_family == "disease main effect",
      discovery_framework %in% c("VAE", "Univariate"),
      omic_type %in% c("proteomics", "chemistries"),
      analysis_domain %in% AOU_DOMAIN_LEVELS
    ) %>%
    mutate(
      method = recode(discovery_framework, "Univariate" = "OLS"),
      feature_identity = case_when(
        omic_type == "proteomics" & !is.na(uniprot) & str_squish(uniprot) != "" ~
          paste0("uniprot:", str_squish(uniprot)),
        TRUE ~ paste0("label:", canonical_omic_key(feature_display))
      ),
      replication_q = suppressWarnings(as.numeric(discovery_concordant_q)),
      modelled = discovery_concordant_status == "ok",
      directional = modelled & !is.na(replication_q) & replication_q < q_cutoff &
        direction_relation == "same_direction"
    ) %>%
    group_by(method, omic_type, analysis_domain, feature_identity) %>%
    summarise(
      feature_display = first_nonmissing(feature_display),
      modelled = any(modelled, na.rm = TRUE),
      directional = any(directional, na.rm = TRUE),
      .groups = "drop"
    )

  coverage <- base %>%
    group_by(method, omic_type) %>%
    summarise(
      n_nominated = n(),
      n_modelled = sum(modelled),
      n_directionally_concordant = sum(directional),
      directional_concordance_fraction = if_else(
        n_modelled > 0, n_directionally_concordant / n_modelled, NA_real_
      ),
      .groups = "drop"
    )

  domains <- base %>%
    filter(directional) %>%
    distinct(method, omic_type, analysis_domain, feature_identity) %>%
    count(method, omic_type, analysis_domain, name = "n_directionally_concordant") %>%
    mutate(analysis_domain_label = unname(AOU_DOMAIN_LABELS[analysis_domain]))

  list(pair_level = base, coverage = coverage, domains = domains)
}

plot_aou_coverage <- function(coverage) {
  plot_df <- coverage %>%
    mutate(
      omic_label = recode(omic_type, chemistries = "Chemistries", proteomics = "Proteomics"),
      point_label = paste(method, omic_label, sep = " – ")
    )
  ggplot(plot_df, aes(x = n_modelled, y = directional_concordance_fraction,
                      shape = method, fill = omic_label)) +
    geom_point(size = 3.1, colour = "black", stroke = 0.35) +
    geom_text_repel(
      aes(label = paste0(point_label, "\n", n_directionally_concordant, "/", n_modelled)),
      size = pt_to_mm(5.0), family = NATURE$font, max.overlaps = Inf
    ) +
    scale_shape_manual(values = c(OLS = 22, VAE = 21), name = "Discovery method") +
    scale_fill_manual(values = c(Chemistries = NATURE_COLORS[["vermillion"]],
                                 Proteomics = NATURE_COLORS[["green"]]), name = "Omic layer") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.85)) +
    labs(title = "Candidate-set size and directional concordance",
         x = "AoU-modelled disease-feature pairs",
         y = "Directionally concordant fraction (q < 0.10)") +
    theme_nature(legend_position = "bottom", grid = TRUE)
}

plot_aou_domain_counts <- function(domains, omic_type_value, title) {
  plot_df <- domains %>%
    filter(omic_type == omic_type_value) %>%
    mutate(
      method = factor(method, levels = c("OLS", "VAE")),
      analysis_domain = factor(analysis_domain, levels = AOU_DOMAIN_LEVELS)
    )
  ggplot(plot_df, aes(x = method, y = n_directionally_concordant,
                      fill = analysis_domain)) +
    geom_col(width = 0.65, colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = AOU_DOMAIN_COLORS, labels = AOU_DOMAIN_LABELS,
                      name = "Disease domain", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, x = NULL, y = "Directionally concordant associations (n)") +
    theme_nature(legend_position = "bottom", grid = FALSE)
}
