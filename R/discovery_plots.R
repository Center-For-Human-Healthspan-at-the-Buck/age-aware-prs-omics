# Figure functions for the Arivale discovery analyses.

canonical_pairs <- tribble(
  ~omics_feature,             ~prs_feature,                         ~prs_label,
  "BILIRUBIN, TOTAL",         "Direct_bilirubin_umol_L_PGS000681", "Bilirubin PRS",
  "LDL-CHOL CALCULATION",     "LDL_cholesterol_PGS003873",         "LDL PRS",
  "HDL CHOL DIRECT",          "HDL_cholesterol_PGS002781",         "HDL PRS"
) %>% mutate(omics_key = str_to_lower(omics_feature))

canonical_prs <- c(
  "Direct_bilirubin_umol_L_PGS000681",
  "HDL_cholesterol_PGS002781",
  "LDL_cholesterol_PGS003873",
  "High_cholesterol_PGS002334"
)

PRS_DOMAINS <- list(
  CVD_PRS = c(
    "Heart_failure_PGS001790",
    "Coronary_atherosclerosis_PGS002048",
    "Coronary_artery_disease_PGS004941",
    "Cardiovascular_disease_PGS002316"
  ),
  T2D_PRS = c(
    "Type_2_diabetes_SNPs_insulin_action_PGS000850",
    "Type_2_diabetes_SNPs_insulin_action_secretion_PGS000851",
    "T2D_cases_vs_HbA1c_filtered_controls_PGS000712",
    "Metabolic_syndrome_PGS004928",
    "Type_2_diabetes_T2D_PGS004923"
  ),
  Neurodegeneration_PRS = c(
    "amyotrophic_lateral_sclerosis",
    "Parkinson_disease_PGS000903",
    "Alzheimer_disease_late_onset_PGS000053",
    "Alzheimers_PGS002753"
  ),
  Longevity_PRS = c(
    "Longevity_PGS000906",
    "mvAge",
    "Timmers_healthspan",
    "Timmers_lifespan",
    "Timmers_longevity"
  )
)

MANUSCRIPT_PRS <- unique(unlist(PRS_DOMAINS, use.names = FALSE))
DOMAIN_LOOKUP <- enframe(PRS_DOMAINS, name = "PRS_domain", value = "prs") %>%
  unnest_longer(prs)
SECTOR_ORDER <- c(
  "CVD_PRS", "T2D_PRS", "Neurodegeneration_PRS", "Longevity_PRS",
  "Proteomics", "Metabolomics", "Chemistries"
)
DOMAIN_COLORS <- c(
  CVD_PRS = "#D55E00", T2D_PRS = "#0072B2",
  Neurodegeneration_PRS = "#009E73", Longevity_PRS = "#CC79A7",
  Proteomics = CATEGORY_COLORS[["proteomics"]],
  Metabolomics = CATEGORY_COLORS[["metabolomics"]],
  Chemistries = CATEGORY_COLORS[["chemistries"]]
)

# PRSs displayed in Supplementary Figures 1-3. This is a figure-specific
# selection, not a set of standalone disease-domain plotting modules.
SUPPLEMENTARY_SIMILARITY_PRS <- c(
  "frailty", "gout", "Atrial_fibrillation_PGS000016",
  "Myocardial_infarction_PGS000710", "Hypertrophic_cardiomyopathy_PGS000739",
  "Essential_hypertension_PGS000958", "Coronary_artery_disease_PGS004941",
  "Hypertension_PGS002765", "Alzheimer_disease_late_onset_PGS000053",
  "Alzheimers_PGS002753", "Parkinson_disease_PGS000903",
  "amyotrophic_lateral_sclerosis", "Colon_cancer_PGS000380",
  "All_Cancers_PGS004875", "breast_cancer", "Type_1_diabetes_T1D_PGS004009",
  "Type_2_diabetes_T2D_PGS004923", "Heart_failure_PGS001790",
  "Coronary_atherosclerosis_PGS002048", "Knee_osteoarthritis_PGS002767",
  "Hip_osteoarthritis_PGS002763", "Osteoporosis_PGS002768",
  "Rheumatoid_arthritis_PGS003994", "Age_related_macular_degeneration_PGS004607",
  "Venous_thromboembolism_PGS004854", "Dilated_cardiomyopathy_PGS004861",
  "MultiAncPhenoAge", "MultiAncestryGrimAge", "Timmers_healthspan",
  "Timmers_lifespan", "Timmers_longevity", "mvAge", "systolic_blood_pressure"
)

DOMAIN_SETS <- list(
  AD = list(
    title = "Alzheimer's disease",
    center_label = "Alzheimer's\nDisease",
    prs = c("Alzheimer_disease_late_onset_PGS000053", "Alzheimers_PGS002753")
  ),
  PD = list(
    title = "Parkinson's disease",
    center_label = "Parkinson's\nDisease",
    prs = c("Parkinson_disease_PGS000903")
  ),
  ALS = list(
    title = "Amyotrophic lateral sclerosis",
    center_label = "Amyotrophic\nLateral\nSclerosis",
    prs = c("amyotrophic_lateral_sclerosis")
  ),
  CVD = list(
    title = "Cardiovascular disease",
    center_label = "Cardiovascular\nDisease",
    # The detailed CVD figure/source data include this additional CVD score.
    # Figure 2 uses the four-score PRS_DOMAINS definition retained above to
    # match the manuscript's 493-versus-120 association-count comparison.
    prs = c(PRS_DOMAINS$CVD_PRS, "Coronary_vascular_disease_CVD_PGS004514")
  ),
  T2D = list(
    title = "Type 2 diabetes",
    center_label = "Type 2\nDiabetes\nPRS",
    prs = PRS_DOMAINS$T2D_PRS
  ),
  Longevity = list(
    title = "Longevity",
    center_label = "Longevity",
    prs = PRS_DOMAINS$Longevity_PRS
  )
)

filter_manuscript_prs <- function(df, age_included = FALSE) {
  if (isTRUE(age_included)) {
    df %>% mutate(prs_base = strip_age_suffix(feature_a_original_name)) %>%
      filter(prs_base %in% MANUSCRIPT_PRS)
  } else {
    df %>% filter(feature_a_original_name %in% MANUSCRIPT_PRS)
  }
}
smart_truncate_label <- function(x, max_chars = 34) {
  x <- as.character(x)
  vapply(x, function(label) {
    if (is.na(label) || nchar(label) <= max_chars) return(label)

    # Preserve trailing chemical/lipid formulae when possible.
    trailing_formula <- str_extract(label, "\\([^()]{2,42}\\)\\*?$")
    if (!is.na(trailing_formula) && nchar(trailing_formula) <= max_chars - 8) {
      prefix_n <- max_chars - nchar(trailing_formula) - 4
      return(paste0(str_sub(label, 1, max(6, prefix_n)), "... ", trailing_formula))
    }

    paste0(str_sub(label, 1, max_chars - 3), "...")
  }, character(1), USE.NAMES = FALSE)
}

abbreviate_lollipop_label <- function(x, max_chars = 34) {
  out <- x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_squish()

  # Protein / cytokine names.
  out <- out %>%
    str_replace_all(regex("chemokine\\s*\\(c\\s*-\\s*c\\s*motif\\)\\s*ligand\\s*([0-9a-z]+)", ignore_case = TRUE), "CCL\\1") %>%
    str_replace_all(regex("chemokine\\s*\\(c\\s*-\\s*x\\s*-\\s*c\\s*motif\\)\\s*ligand\\s*([0-9a-z]+)", ignore_case = TRUE), "CXCL\\1") %>%
    str_replace_all(regex("chemokine\\s*\\(c\\s*motif\\)\\s*ligand\\s*([0-9a-z]+)", ignore_case = TRUE), "XCL\\1") %>%
    str_replace_all(regex("interleukin\\s+([0-9]+)([a-z]?)\\s+receptor,\\s*alpha\\s*([0-9]?)", ignore_case = TRUE), "IL-\\1\\2Ra\\3") %>%
    str_replace_all(regex("interleukin\\s+([0-9]+)([a-z]?)\\s+receptor,\\s*beta\\s*([0-9]?)", ignore_case = TRUE), "IL-\\1\\2Rb\\3") %>%
    str_replace_all(regex("interleukin\\s+([0-9]+)([a-z]?)\\s+receptor-like\\s*([0-9]*)", ignore_case = TRUE), "IL-\\1\\2RL\\3") %>%
    str_replace_all(regex("interleukin\\s+([0-9]+)ra([0-9])", ignore_case = TRUE), "IL-\\1Ra\\2") %>%
    str_replace_all(regex("interleukin\\s+([0-9]+)([a-z]?)", ignore_case = TRUE), "IL-\\1\\2") %>%
    str_replace_all(regex("tumor necrosis factor receptor superfamily,\\s*member\\s*([0-9a-z]+)", ignore_case = TRUE), "TNFRSF\\1") %>%
    str_replace_all(regex("tnf receptor superfamily,\\s*member\\s*([0-9a-z]+)", ignore_case = TRUE), "TNFRSF\\1") %>%
    str_replace_all(regex("tumor necrosis factor\\s*\\(ligand\\)\\s*superfamily,\\s*member\\s*([0-9a-z]+)", ignore_case = TRUE), "TNFSF\\1") %>%
    str_replace_all(regex("tumor necrosis factor", ignore_case = TRUE), "TNF") %>%
    str_replace_all(regex("fibroblast growth factor\\s*([0-9a-z]+)", ignore_case = TRUE), "FGF\\1") %>%
    str_replace_all(regex("growth differentiation factor\\s*([0-9a-z]+)", ignore_case = TRUE), "GDF\\1") %>%
    str_replace_all(regex("bone morphogenetic protein\\s*([0-9a-z]+)", ignore_case = TRUE), "BMP\\1") %>%
    str_replace_all(regex("matrix metallopeptidase\\s*([0-9]+).*", ignore_case = TRUE), "MMP\\1") %>%
    str_replace_all(regex("interstitial collagenase", ignore_case = TRUE), "MMP1") %>%
    str_replace_all(regex("glial cell derived neurotrophic factor", ignore_case = TRUE), "GDNF") %>%
    str_replace_all(regex("vascular endothelial growth factor a|vegfa", ignore_case = TRUE), "VEGFA") %>%
    str_replace_all(regex("n-terminal prohorm brain natriuretic pep", ignore_case = TRUE), "NT-proBNP") %>%
    str_replace_all(regex("natriuretic peptide b", ignore_case = TRUE), "BNP") %>%
    str_replace_all(regex("programmed cell death 1 ligand 2", ignore_case = TRUE), "PD-L2") %>%
    str_replace_all(regex("epithelial cell adhesion molecule", ignore_case = TRUE), "EPCAM") %>%
    str_replace_all(regex("plasminogen activator,\\s*urokinase", ignore_case = TRUE), "uPA") %>%
    str_replace_all(regex("plasminogen activator,\\s*tissue", ignore_case = TRUE), "tPA") %>%
    str_replace_all(regex("chitinase 1 \\(chitotriosidase\\)", ignore_case = TRUE), "CHIT1") %>%
    str_replace_all(regex("chitinase 3-like 1", ignore_case = TRUE), "CHI3L1") %>%
    str_replace_all(regex("cub domain containing protein 1", ignore_case = TRUE), "CDCP1") %>%
    str_replace_all(regex("serpin peptidase inhibitor,\\s*clade a,\\s*member 12", ignore_case = TRUE), "SERPINA12") %>%
    str_replace_all(regex("proline/arginine-rich end leucine-rich repeat protein", ignore_case = TRUE), "PRELP") %>%
    str_replace_all(regex("v-set igg domain containing 2", ignore_case = TRUE), "VSIG2") %>%
    str_replace_all(regex("v-src avian sarcoma viral oncogene homolog", ignore_case = TRUE), "v-Src homolog") %>%
    str_replace_all(regex("fms-related tyrosine kinase 3 ligand", ignore_case = TRUE), "FLT3 ligand") %>%
    str_replace_all(regex("tek tyrosine kinase, endothelial", ignore_case = TRUE), "TEK") %>%
    str_replace_all(regex("c-mer proto-oncogene tyrosine kinase", ignore_case = TRUE), "MERTK") %>%
    str_replace_all(regex("axl receptor tyrosine kinase", ignore_case = TRUE), "AXL") %>%
    str_replace_all(regex("leukemia inhibitory factor", ignore_case = TRUE), "LIF")

  # Common chemistry / clinical chemistry abbreviations.
  out <- out %>%
    str_replace_all(regex("mean arterial blood pressure", ignore_case = TRUE), "mean arterial BP") %>%
    str_replace_all(regex("triglyceride hdl ratio", ignore_case = TRUE), "TG/HDL ratio") %>%
    str_replace_all(regex("ldl-chol calculation", ignore_case = TRUE), "LDL-C calc.") %>%
    str_replace_all(regex("hdl chol direct", ignore_case = TRUE), "HDL-C direct") %>%
    str_replace_all(regex("ldl particle number", ignore_case = TRUE), "LDL-P") %>%
    str_replace_all(regex("hdl particle number", ignore_case = TRUE), "HDL-P") %>%
    str_replace_all(regex("asat \\(sgot\\)", ignore_case = TRUE), "AST (SGOT)") %>%
    str_replace_all(regex("alat \\(sgpt\\)", ignore_case = TRUE), "ALT (SGPT)") %>%
    str_replace_all(regex("bun/creat ratio", ignore_case = TRUE), "BUN/creat ratio") %>%
    str_replace_all(regex("crp high sensitivity", ignore_case = TRUE), "hsCRP") %>%
    str_replace_all(regex("carbon dioxide \\(co2\\)", ignore_case = TRUE), "CO2") %>%
    str_replace_all(regex("1,5-anhydroglucitol \\(1,5-ag\\)", ignore_case = TRUE), "1,5-AG") %>%
    str_replace_all(regex("5alpha-preg", ignore_case = TRUE), "5a-preg") %>%
    str_replace_all(regex("pregnanediol-3-glucuronide", ignore_case = TRUE), "pregnanediol-3-gluc.") %>%
    str_replace_all(regex("pregnenetriol disulfate", ignore_case = TRUE), "pregnenetriol disulf.") %>%
    str_replace_all(regex("pregnenediol disulfate", ignore_case = TRUE), "pregnenediol disulf.") %>%
    str_replace_all(regex("glutamine conjugate of c6h10o2", ignore_case = TRUE), "glutamine conj. C6H10O2") %>%
    str_replace_all(regex("sulfate of piperine metabolite c16h19no3", ignore_case = TRUE), "piperine metabolite sulfate") %>%
    str_replace_all(regex("3-bromo-5-chloro-2,6-dihydroxybenzoic acid", ignore_case = TRUE), "bromo-chloro-DHBA") %>%
    str_replace_all(regex("3,5-dichloro-2,6-dihydroxybenzoic acid", ignore_case = TRUE), "dichloro-DHBA")

  # Lipid shorthand: preserve class and acyl composition for very long systematic names.
  out <- vapply(out, function(label) {
    if (is.na(label) || nchar(label) <= max_chars) return(label)

    star <- if_else(str_detect(label, "\\*$"), "*", "")
    formula <- str_extract(label, regex("\\((p-)?[0-9]{1,2}:[0-9][a-z0-9]*/[0-9]{1,2}:[0-9][a-z0-9]*\\)", ignore_case = TRUE))
    lipid_class <- case_when(
      str_detect(label, regex("gpc", ignore_case = TRUE)) ~ "GPC",
      str_detect(label, regex("gpe", ignore_case = TRUE)) ~ "GPE",
      str_detect(label, regex("gpi", ignore_case = TRUE)) ~ "GPI",
      TRUE ~ NA_character_
    )
    if (!is.na(formula) && !is.na(lipid_class)) {
      formula <- str_replace_all(formula, regex("^\\(p-", ignore_case = TRUE), "(P-")
      return(paste0(lipid_class, " ", formula, star))
    }

    sphingo_formula <- str_extract(label, regex("\\((d[0-9]{1,2}:[0-9]/[0-9]{1,2}:[0-9]).*\\)", ignore_case = TRUE))
    sphingo_class <- case_when(
      str_detect(label, regex("glycosyl.*ceramide", ignore_case = TRUE)) ~ "glycosyl ceramide",
      str_detect(label, regex("ceramide", ignore_case = TRUE)) ~ "ceramide",
      str_detect(label, regex("sphingomyelin", ignore_case = TRUE)) ~ "sphingomyelin",
      TRUE ~ NA_character_
    )
    if (!is.na(sphingo_formula) && !is.na(sphingo_class)) {
      sphingo_formula <- str_replace(sphingo_formula, ",.*\\)", ", ...)")
      return(paste0(sphingo_class, " ", sphingo_formula, star))
    }

    label
  }, character(1), USE.NAMES = FALSE)

  out <- out %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()

  smart_truncate_label(out, max_chars = max_chars)
}

resolve_lollipop_label_max_chars <- function(n_labels, requested = NULL) {
  if (!is.null(requested) && !is.na(requested)) return(requested)
  case_when(
    n_labels >= 110 ~ 28,
    n_labels >= 80  ~ 30,
    n_labels >= 55  ~ 34,
    TRUE            ~ 42
  )
}

make_lollipop_label_key <- function(plot_df, label_max_chars = NULL) {
  plot_df <- ensure_omic_name_columns(plot_df)
  full_labels <- plot_df %>%
    transmute(
      omic_original_name_lower,
      omic_category,
      full_label = clean_omic_label(omic_original_name_lower, title_case = FALSE)
    ) %>%
    distinct()

  max_chars <- resolve_lollipop_label_max_chars(nrow(full_labels), label_max_chars)

  full_labels %>%
    mutate(
      plotted_label = abbreviate_lollipop_label(full_label, max_chars = max_chars),
      label_max_chars = max_chars,
      label_was_abbreviated = plotted_label != full_label
    ) %>%
    arrange(desc(label_was_abbreviated), full_label)
}

plot_canonical_benchmark <- function(ols_age_residualized, unique_pairs, canonical_pairs) {
  # The current age-residualized OLS table contains one all-participant estimate
  # per PRS-feature pair. Keep the explicit filter as a defensive check because
  # prepare_ols_results() adds age_group = "All" for shared helper compatibility.
  ols_all <- ols_age_residualized %>% filter(age_group == "All")

  r2_df <- ols_all %>%
    mutate(omics_key = str_to_lower(omics_feature)) %>%
    filter(
      prs_feature %in% canonical_pairs$prs_feature,
      omics_key %in% canonical_pairs$omics_key
    ) %>%
    inner_join(canonical_pairs, by = c("omics_key", "prs_feature")) %>%
    transmute(
      label = paste0(prs_label, "\n", str_wrap(omics_feature.y, width = 20)),
      value = as.numeric(r_sq),
      type = "R² (linear model)",
      ymin = NA_real_,
      ymax = NA_real_
    ) %>%
    distinct(label, type, .keep_all = TRUE)

  if (any(r2_df$value < 0 | r2_df$value > 1, na.rm = TRUE)) {
    stop(
      "Canonical OLS R-squared values fall outside [0, 1]. Check the OLS input table.",
      call. = FALSE
    )
  }

  vae_df <- unique_pairs %>%
    mutate(omics_key = str_to_lower(omic_original_name)) %>%
    filter(
      feature_a_original_name %in% canonical_pairs$prs_feature,
      omics_key %in% canonical_pairs$omics_key,
      percent_reproducibility >= 50
    ) %>%
    inner_join(
      canonical_pairs,
      by = c("feature_a_original_name" = "prs_feature", "omics_key")
    ) %>%
    transmute(
      label = paste0(prs_label, "\n", str_wrap(omics_feature, width = 20)),
      value = as.numeric(norm_bayes_k),
      type = "Normalized Bayes K (VAE)",
      ymin = value - as.numeric(norm_sd_bayes_k),
      ymax = value + as.numeric(norm_sd_bayes_k)
    ) %>%
    distinct(label, type, .keep_all = TRUE)

  combined_df <- bind_rows(r2_df, vae_df) %>%
    mutate(label = fct_inorder(label))

  ggplot(combined_df, aes(x = label)) +
    geom_col(
      data = filter(combined_df, type == "R² (linear model)"),
      aes(y = value, fill = type),
      width = 0.55,
      colour = "black",
      linewidth = 0.15
    ) +
    geom_errorbar(
      data = filter(combined_df, type == "Normalized Bayes K (VAE)"),
      aes(ymin = ymin, ymax = ymax, colour = type),
      width = 0.15,
      linewidth = 0.25
    ) +
    geom_point(
      data = filter(combined_df, type == "Normalized Bayes K (VAE)"),
      aes(y = value, colour = type),
      size = 1.6
    ) +
    scale_fill_manual(
      values = c("R² (linear model)" = NATURE_COLORS[["grey"]]),
      name = NULL
    ) +
    scale_colour_manual(
      values = c("Normalized Bayes K (VAE)" = NATURE_COLORS[["blue"]]),
      name = NULL
    ) +
    labs(x = NULL, y = "Effect size") +
    coord_cartesian(clip = "off") +
    theme_nature(legend_position = "top") +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 1,
        size = 5.0,
        lineheight = 0.9,
        margin = margin(t = 3)
      ),
      plot.margin = margin(4, 5, 12, 5)
    )
}

plot_canonical_benchmark_age <- function(ols_age_included, unique_pairs_age, canonical_pairs) {
  vae_age_df <- unique_pairs_age %>%
    mutate(
      prs_base = str_remove(feature_a_original_name, "_Age_>_45$|_Age_<_45$"),
      age_group = case_when(
        str_detect(feature_a_original_name, "_Age_>_45$") ~ ">45",
        str_detect(feature_a_original_name, "_Age_<_45$") ~ "<45",
        TRUE ~ "All"
      ),
      omics_key = str_to_lower(omic_original_name)
    ) %>%
    filter(
      percent_reproducibility >= 50,
      prs_base %in% canonical_pairs$prs_feature,
      omics_key %in% canonical_pairs$omics_key
    ) %>%
    inner_join(canonical_pairs, by = c("prs_base" = "prs_feature", "omics_key")) %>%
    transmute(
      label = paste0(prs_label, "\n", str_wrap(omics_feature, width = 20)),
      value = as.numeric(norm_bayes_k),
      ymin = value - as.numeric(norm_sd_bayes_k),
      ymax = value + as.numeric(norm_sd_bayes_k),
      type = "Normalized Bayes K (VAE)",
      model = age_group
    ) %>%
    distinct(label, type, model, .keep_all = TRUE)

  # The current age-included OLS table is not age-stratified. It contributes
  # one all-participant R-squared value per canonical pair; the <45 and >45
  # contexts shown in this panel come from the age-included VAE only.
  linear_age_df <- ols_age_included %>%
    filter(age_group == "All") %>%
    mutate(omics_key = str_to_lower(omics_feature)) %>%
    filter(
      prs_feature %in% canonical_pairs$prs_feature,
      omics_key %in% canonical_pairs$omics_key
    ) %>%
    inner_join(canonical_pairs, by = c("omics_key", "prs_feature")) %>%
    transmute(
      label = paste0(prs_label, "\n", str_wrap(omics_feature.y, width = 20)),
      value = as.numeric(r_sq),
      ymin = NA_real_,
      ymax = NA_real_,
      type = "R² (linear model)",
      model = "All"
    ) %>%
    distinct(label, type, model, .keep_all = TRUE)

  if (any(linear_age_df$value < 0 | linear_age_df$value > 1, na.rm = TRUE)) {
    stop(
      "Age-included all-participant OLS R-squared values fall outside [0, 1]. Check the OLS input table.",
      call. = FALSE
    )
  }

  label_levels <- unique(c(linear_age_df$label, vae_age_df$label))
  linear_age_df <- linear_age_df %>%
    mutate(label = factor(label, levels = label_levels))
  vae_age_df <- vae_age_df %>%
    mutate(
      label = factor(label, levels = label_levels),
      model = factor(model, levels = c("All", "<45", ">45"))
    )

  pos <- position_dodge(width = 0.65)

  ggplot() +
    geom_col(
      data = linear_age_df,
      aes(x = label, y = value, fill = type),
      width = 0.55,
      colour = "black",
      linewidth = 0.15
    ) +
    geom_errorbar(
      data = vae_age_df,
      aes(x = label, ymin = ymin, ymax = ymax, group = model, colour = model),
      position = pos,
      width = 0.12,
      linewidth = 0.25
    ) +
    geom_point(
      data = vae_age_df,
      aes(x = label, y = value, group = model, colour = model),
      size = 1.5,
      position = pos
    ) +
    scale_fill_manual(
      values = c("R² (linear model)" = NATURE_COLORS[["grey"]]),
      name = NULL
    ) +
    scale_colour_manual(
      values = MODEL_COLORS[c("All", "<45", ">45")],
      name = "VAE age context",
      drop = FALSE
    ) +
    labs(x = NULL, y = "Effect size") +
    coord_cartesian(clip = "off") +
    theme_nature(legend_position = "top") +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 1,
        size = 5.0,
        lineheight = 0.9,
        margin = margin(t = 3)
      ),
      plot.margin = margin(4, 5, 12, 5)
    )
}

make_chord_data <- function(unique_pairs_df, stratified = FALSE, repro_threshold = 50) {
  df <- ensure_omic_name_columns(unique_pairs_df)
  if (isTRUE(stratified)) {
    df <- df %>% mutate(prs_base = str_remove(feature_a_original_name, "_Age_>_45$|_Age_<_45$"))
  } else {
    df <- df %>% mutate(prs_base = feature_a_original_name)
  }

  if (!is.null(repro_threshold)) {
    df <- df %>% filter(percent_reproducibility >= repro_threshold)
  }

  df %>%
    left_join(DOMAIN_LOOKUP, by = c("prs_base" = "prs")) %>%
    filter(!is.na(PRS_domain)) %>%
    mutate(omic_category = str_to_title(omic_category)) %>%
    count(PRS_domain, omic_category, name = "n") %>%
    rename(from = PRS_domain, to = omic_category)
}

plot_chord_counts <- function(
  chord_data,
  title = NULL,
  sector_order = SECTOR_ORDER,
  grid_col = DOMAIN_COLORS,
  label_cex = 0.62,
  title_cex = 0.78,
  label_offset = 1.08,
  track_height = 0.045,
  directional = TRUE
) {
  sectors_present <- unique(c(chord_data$from, chord_data$to))
  order_use <- sector_order[sector_order %in% sectors_present]
  grid_col_use <- grid_col[names(grid_col) %in% sectors_present]

  gap_after <- rep(2, length(order_use))
  names(gap_after) <- order_use
  if ("Longevity_PRS" %in% names(gap_after)) gap_after[["Longevity_PRS"]] <- 8
  gap_after[length(gap_after)] <- 8

  circos.clear()
  op <- par(mar = c(1, 1, 1, 1), family = NATURE$font)
  on.exit({
    par(op)
    circos.clear()
  }, add = TRUE)

  circos.par(gap.after = gap_after, start.degree = 90, points.overflow.warning = FALSE)

  chordDiagram(
    chord_data,
    grid.col = grid_col_use,
    transparency = 0.38,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = track_height),
    order = order_use,
    directional = as.integer(directional),
    direction.type = if (isTRUE(directional)) "arrows" else "diffHeight",
    link.arr.type = if (isTRUE(directional)) "big.arrow" else NULL,
    link.sort = TRUE,
    link.largest.ontop = TRUE
  )

  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")

    label_col <- ifelse(str_detect(sector_name, "_PRS$"), "#A00000", "black")
    sector_label <- sector_name %>% str_replace_all("_", " ")

    circos.text(
      x = mean(xlim),
      y = label_offset,
      labels = sector_label,
      facing = "bending",
      niceFacing = TRUE,
      cex = label_cex,
      font = 2,
      col = label_col
    )
  }, bg.border = NA)

  if (!is.null(title)) {
    title(title, cex.main = title_cex, font.main = 2, family = NATURE$font)
  }
}

plot_top_omic_jitter <- function(df, title, n_features = 10, repro_threshold = 50, label_n = 2) {
  top_features <- df %>%
    filter(percent_reproducibility >= repro_threshold, !feature_a_original_name %in% canonical_prs) %>%
    count(omic_original_name_lower, omic_category, sort = TRUE) %>%
    slice_head(n = n_features)

  plot_df <- df %>%
    filter(
      percent_reproducibility >= repro_threshold,
      omic_original_name_lower %in% top_features$omic_original_name_lower,
      !feature_a_original_name %in% canonical_prs
    ) %>%
    mutate(
      abs_norm_bayes_k = abs(as.numeric(norm_bayes_k)),
      omic_label = str_wrap(clean_omic_label(omic_original_name_lower, title_case = FALSE), 34),
      prs_clean = clean_prs_label(feature_a_original_name),
      omic_category = factor(omic_category, levels = names(CATEGORY_COLORS))
    )

  top_labeled <- plot_df %>%
    group_by(omic_label) %>%
    arrange(desc(abs_norm_bayes_k)) %>%
    slice_head(n = label_n) %>%
    ungroup()

  ggplot(plot_df, aes(x = fct_reorder(omic_label, abs_norm_bayes_k, .fun = max), y = abs_norm_bayes_k, colour = omic_category)) +
    geom_jitter(width = 0.16, alpha = 0.85, size = 0.8) +
    geom_text_repel(
      data = top_labeled,
      aes(label = prs_clean),
      size = pt_to_mm(5.2),
      family = NATURE$font,
      fontface = "bold",
      box.padding = 0.15,
      max.overlaps = 30,
      segment.size = 0.15,
      min.segment.length = 0
    ) +
    scale_colour_manual(values = CATEGORY_COLORS, name = "Omic category") +
    labs(title = title, x = NULL, y = "|Normalized Bayes K|") +
    theme_nature(legend_position = "top") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
}

clean_matrix <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat
}

build_cosine_matrix <- function(unique_pairs_df, stratified = FALSE, prs_keep = SUPPLEMENTARY_SIMILARITY_PRS, repro_threshold = 50) {
  df <- unique_pairs_df %>%
    filter(percent_reproducibility >= repro_threshold) %>%
    mutate(
      prs_base = if (isTRUE(stratified)) str_remove(feature_a_original_name, "_Age_>_45$|_Age_<_45$") else feature_a_original_name,
      norm_bayes_k = as.numeric(norm_bayes_k),
      omic_original_name_lower = as.character(omic_original_name_lower)
    ) %>%
    filter(prs_base %in% prs_keep) %>%
    distinct(prs_base, omic_original_name_lower, .keep_all = TRUE)

  prs_phenotype_matrix <- df %>%
    select(prs_base, omic_original_name_lower, norm_bayes_k) %>%
    pivot_wider(names_from = omic_original_name_lower, values_from = norm_bayes_k, values_fill = 0)

  row_names <- prs_phenotype_matrix$prs_base
  phenotype_matrix <- prs_phenotype_matrix %>% select(-prs_base) %>% as.matrix()
  rownames(phenotype_matrix) <- row_names

  text2vec::sim2(phenotype_matrix, method = "cosine", norm = "l2")
}

plot_mantel_relationship <- function(prs_cor_matrix, cos_sim_matrix, cos_sim_matrix_age) {
  prs_cor_clean_resid <- clean_matrix(prs_cor_matrix)
  prs_cosine_clean_resid <- clean_matrix(cos_sim_matrix)
  common_resid <- intersect(rownames(prs_cor_clean_resid), rownames(prs_cosine_clean_resid))

  gen_dist_resid <- as.dist(1 - prs_cor_clean_resid[common_resid, common_resid])
  vae_dist_resid <- as.dist(1 - prs_cosine_clean_resid[common_resid, common_resid])
  mantel_resid <- vegan::mantel(gen_dist_resid, vae_dist_resid, method = "spearman", permutations = 999)

  prs_cor_clean_strat <- clean_matrix(prs_cor_matrix)
  prs_cosine_clean_strat <- clean_matrix(cos_sim_matrix_age)
  common_strat <- intersect(rownames(prs_cor_clean_strat), rownames(prs_cosine_clean_strat))

  gen_dist_strat <- as.dist(1 - prs_cor_clean_strat[common_strat, common_strat])
  vae_dist_strat <- as.dist(1 - prs_cosine_clean_strat[common_strat, common_strat])
  mantel_strat <- vegan::mantel(gen_dist_strat, vae_dist_strat, method = "spearman", permutations = 999)

  combined_mantel_df <- bind_rows(
    tibble(Genetic_Dissimilarity = as.vector(gen_dist_resid), Phenotypic_Dissimilarity = as.vector(vae_dist_resid), Dataset = "Residualized"),
    tibble(Genetic_Dissimilarity = as.vector(gen_dist_strat), Phenotypic_Dissimilarity = as.vector(vae_dist_strat), Dataset = "Age-stratified")
  )

  label_resid <- paste0("Residualized: r = ", round(mantel_resid$statistic, 2), ", p = ", signif(mantel_resid$signif, 2))
  label_strat <- paste0("Age-stratified: r = ", round(mantel_strat$statistic, 2), ", p = ", signif(mantel_strat$signif, 2))

  ggplot(combined_mantel_df, aes(x = Genetic_Dissimilarity, y = Phenotypic_Dissimilarity, colour = Dataset)) +
    geom_point(alpha = 0.55, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.3) +
    scale_colour_manual(values = c("Residualized" = NATURE_COLORS[["blue"]], "Age-stratified" = NATURE_COLORS[["vermillion"]]), name = NULL) +
    annotate("text", x = Inf, y = Inf, label = label_resid, hjust = 1.02, vjust = 1.4, size = pt_to_mm(5.2), colour = NATURE_COLORS[["blue"]], family = NATURE$font) +
    annotate("text", x = Inf, y = Inf, label = label_strat, hjust = 1.02, vjust = 2.8, size = pt_to_mm(5.2), colour = NATURE_COLORS[["vermillion"]], family = NATURE$font) +
    labs(x = "Genetic dissimilarity (1 − PRS correlation)", y = "Omic dissimilarity (1 − cosine similarity)") +
    theme_nature(legend_position = "top")
}

prepare_domain_qtl_data <- function(
  df,
  prs_to_plot,
  stratified = FALSE,
  repro_threshold = 50,
  ci_mult = 1.96
) {
  df <- ensure_omic_name_columns(df)
  if (isTRUE(stratified)) {
    df <- df %>%
      mutate(
        prs_base = str_remove(feature_a_original_name, "_Age_>_45$|_Age_<_45$"),
        AgeGroup = case_when(
          str_detect(feature_a_original_name, "_Age_>_45$") ~ ">45",
          str_detect(feature_a_original_name, "_Age_<_45$") ~ "<45",
          TRUE ~ "All"
        )
      )
  } else {
    df <- df %>% mutate(prs_base = feature_a_original_name, AgeGroup = NA_character_)
  }

  df %>%
    filter(prs_base %in% prs_to_plot, percent_reproducibility >= repro_threshold) %>%
    distinct(feature_a_original_name, omic_original_name_lower, .keep_all = TRUE) %>%
    mutate(
      Feature = clean_omic_label(omic_original_name_lower),
      est = as.numeric(norm_bayes_k),
      se = as.numeric(norm_sd_bayes_k),
      ci_lo = est - ci_mult * se,
      ci_hi = est + ci_mult * se,
      has_eQTL = valid_qtl_value(Mapped_gene_e),
      has_pQTL = valid_qtl_value(Mapped_gene_p),
      QTL_type = case_when(
        has_eQTL & has_pQTL ~ "both",
        has_eQTL ~ "eQTL",
        has_pQTL ~ "pQTL",
        TRUE ~ "none"
      ),
      QTL_type = factor(QTL_type, levels = names(QTL_SHAPES)),
      omic_category = factor(str_to_lower(omic_category), levels = names(CATEGORY_COLORS))
    )
}

prepare_domain_combined <- function(
  prs_to_plot,
  unique_pairs,
  unique_pairs_age,
  repro_threshold = 50
) {
  resid <- prepare_domain_qtl_data(unique_pairs, prs_to_plot, stratified = FALSE, repro_threshold = repro_threshold) %>%
    mutate(Model = "Residualized")

  strat <- prepare_domain_qtl_data(unique_pairs_age, prs_to_plot, stratified = TRUE, repro_threshold = repro_threshold) %>%
    mutate(Model = AgeGroup)

  bind_rows(resid, strat) %>%
    mutate(
      Model = factor(Model, levels = names(MODEL_COLORS)),
      ci_width = pmax(ci_hi - ci_lo, 1e-6),
      ball_size = 1 / sqrt(ci_width),
      effect_dir = if_else(est >= 0, "Positive", "Negative")
    )
}

plot_domain_forest <- function(
  plot_df,
  title = NULL,
  axis_text_y_pt = NATURE$small_pt
) {
  feature_order <- plot_df %>%
    group_by(Feature) %>%
    summarise(mean_est = mean(est, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_est) %>%
    pull(Feature) %>%
    as.character()

  model_offsets <- c(
    "Residualized" = -0.27,
    "<45" = -0.09,
    ">45" = 0.09,
    "All" = 0.27
  )

  plot_df <- plot_df %>%
    mutate(
      Feature = as.character(Feature),
      feature_y = match(Feature, feature_order),
      model_offset = unname(model_offsets[as.character(Model)]),
      model_offset = replace_na(model_offset, 0),
      y_plot = feature_y + model_offset
    )

  ggplot(plot_df, aes(x = est, y = y_plot)) +
    # Explicit y/yend coordinates replace position_dodge(). The old use of
    # geom_segment() with position_dodge() could move only one endpoint in some
    # ggplot2 versions, which produced visibly slanted confidence intervals.
    geom_segment(
      aes(x = ci_lo, xend = ci_hi, yend = y_plot),
      linewidth = 0.25,
      colour = "grey35"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey55"
    ) +
    geom_point(
      aes(shape = QTL_type, fill = Model),
      size = 1.6,
      stroke = 0.25,
      colour = "black"
    ) +
    scale_y_continuous(
      breaks = seq_along(feature_order),
      labels = feature_order,
      limits = c(0.5, length(feature_order) + 0.5),
      expand = expansion(add = 0)
    ) +
    scale_shape_manual(values = QTL_SHAPES, name = "QTL type") +
    scale_fill_manual(
      values = MODEL_COLORS,
      name = "Model / age group",
      drop = FALSE
    ) +
    labs(title = title, x = "Normalized Bayes K ± 95% CI", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_nature(legend_position = "top") +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(
        size = axis_text_y_pt,
        colour = "black",
        lineheight = 0.86
      )
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(fill = "white")),
      fill = guide_legend(
        order = 2,
        override.aes = list(shape = 21, colour = "black")
      )
    )
}

plot_domain_lollipop <- function(
  plot_df,
  center_label,
  label_size_pt = 5.6,
  center_size_pt = 7,
  radius_base = 0.15,
  stem_scale = 1.0,
  label_pad = 0.035,
  min_stem = 0.018,
  legend_position = "none",
  label_max_chars = NULL
) {
  if (nrow(plot_df) == 0) return(NULL)

  plot_df <- plot_df %>%
    mutate(
      omic_label_full = clean_omic_label(omic_original_name_lower, title_case = FALSE),
      label_for_group = omic_original_name_lower
    )

  stem_df <- plot_df %>%
    group_by(label_for_group) %>%
    slice_max(abs(est), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(abs(est))) %>%
    mutate(id = row_number())

  if (nrow(stem_df) == 0) return(NULL)

  n_labels <- nrow(stem_df)
  max_chars_use <- resolve_lollipop_label_max_chars(n_labels, label_max_chars)

  stem_df <- stem_df %>%
    mutate(omic_label = abbreviate_lollipop_label(omic_label_full, max_chars = max_chars_use))

  plot_points <- plot_df %>%
    inner_join(stem_df %>% select(label_for_group, id), by = "label_for_group")

  stem_df <- stem_df %>%
    mutate(
      angle_raw = 90 - 360 * ((id - 0.5) / n_labels),
      hjust = if_else(angle_raw < -90, 1, 0),
      angle = if_else(angle_raw < -90, angle_raw + 180, angle_raw),
      stem_y = radius_base + pmax(abs(est) * stem_scale, min_stem),
      # Tuned to keep labels inside the square plotting region while providing
      # slightly more clearance for long labels and large stems.
      label_distance = label_pad + pmin(nchar(omic_label) * 0.0008, 0.035),
      label_y = stem_y + label_distance
    )

  # Shortened labels plus a modest radial buffer keep text inside the PDF page
  # without compressing the circular display further.
  max_y <- max(stem_df$label_y, na.rm = TRUE) * 1.10
  center_x <- as.character(stem_df$id[[1]])

  ggplot() +
    geom_segment(
      data = stem_df,
      aes(
        x = factor(id), xend = factor(id),
        y = radius_base, yend = stem_y,
        colour = effect_dir
      ),
      linewidth = 0.35
    ) +
    scale_colour_manual(name = "Effect direction", values = EFFECT_COLORS) +
    ggnewscale::new_scale_color() +
    geom_text(
      data = stem_df,
      aes(
        x = factor(id), y = label_y,
        label = omic_label,
        angle = angle,
        hjust = hjust,
        colour = omic_category
      ),
      size = pt_to_mm(label_size_pt),
      family = NATURE$font,
      fontface = "bold",
      lineheight = 0.86,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = CATEGORY_COLORS, guide = "none", drop = FALSE) +
    geom_point(
      data = plot_points,
      aes(
        x = factor(id),
        y = radius_base + pmax(abs(est) * stem_scale, min_stem),
        size = ball_size,
        shape = QTL_type,
        fill = Model
      ),
      stroke = 0.25,
      colour = "black",
      alpha = 0.95
    ) +
    annotate(
      "text",
      x = center_x,
      y = 0,
      label = center_label,
      size = pt_to_mm(center_size_pt),
      family = NATURE$font,
      fontface = "bold",
      lineheight = 0.9,
      hjust = 0.5,
      vjust = 0.5
    ) +
    coord_polar(clip = "off") +
    scale_y_continuous(limits = c(0, max_y), expand = expansion(mult = c(0, 0))) +
    scale_size_continuous(range = c(1.2, 3.4), guide = "none") +
    scale_fill_manual(name = "Model / age group", values = MODEL_COLORS, drop = FALSE) +
    scale_shape_manual(name = "QTL type", values = QTL_SHAPES, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = NATURE$base_pt, base_family = NATURE$font) +
    theme(
      plot.margin = margin(14, 14, 14, 14),
      legend.position = legend_position,
      legend.title = element_text(size = NATURE$base_pt, face = "bold"),
      legend.text = element_text(size = NATURE$small_pt)
    )
}

make_domain_replicated_table <- function(prs_to_plot, unique_pairs, unique_pairs_age, repro_threshold = 50) {
  base_df <- prepare_domain_qtl_data(unique_pairs, prs_to_plot, stratified = FALSE, repro_threshold = repro_threshold) %>%
    mutate(
      prs_base = feature_a_original_name,
      omic_clean = clean_omic_label(omic_original_name_lower),
      label = paste(prs_base, omic_original_name_lower, sep = " | "),
      est_base = est,
      se_base = se,
      ci_lo_base = ci_lo,
      ci_hi_base = ci_hi
    )

  age_df <- prepare_domain_qtl_data(unique_pairs_age, prs_to_plot, stratified = TRUE, repro_threshold = repro_threshold) %>%
    mutate(
      omic_clean = clean_omic_label(omic_original_name_lower),
      label = paste(prs_base, omic_original_name_lower, sep = " | "),
      est_age = est,
      se_age = se,
      ci_lo_age = ci_lo,
      ci_hi_age = ci_hi,
      age_group = AgeGroup
    )

  common_labels <- intersect(base_df$label, age_df$label)

  inner_join(
    base_df %>% filter(label %in% common_labels) %>%
      select(label, omic_clean, omic_category, est_base, se_base, ci_lo_base, ci_hi_base, QTL_type, any_of(c("fdr", "FDR_q"))),
    age_df %>% filter(label %in% common_labels) %>%
      select(label, omic_clean, est_age, se_age, ci_lo_age, ci_hi_age, age_group),
    by = c("label", "omic_clean")
  )
}

plot_domain_bubble <- function(
  replicated_df,
  title = NULL,
  label_n = 14,
  label_max_chars = 38
) {
  if (nrow(replicated_df) == 0) return(NULL)

  fdr_values <- if ("fdr" %in% names(replicated_df)) {
    replicated_df$fdr
  } else if ("FDR_q" %in% names(replicated_df)) {
    replicated_df$FDR_q
  } else {
    rep(NA_real_, nrow(replicated_df))
  }

  plot_df <- replicated_df %>%
    mutate(
      fdr_plot = suppressWarnings(as.numeric(fdr_values)),
      fdr_plot = if_else(is.na(fdr_plot) | fdr_plot <= 0, 1e-300, fdr_plot),
      ci_width = pmax(ci_hi_base - ci_lo_base, 1e-6),
      bubble_size = 1 / sqrt(ci_width),
      neg_log10_fdr = -log10(fdr_plot),
      omic_category = factor(str_to_lower(omic_category), levels = names(CATEGORY_COLORS))
    ) %>%
    group_by(omic_clean) %>%
    slice_max(neg_log10_fdr, n = 1, with_ties = FALSE) %>%
    ungroup()

  label_df <- plot_df %>%
    arrange(desc(neg_log10_fdr), desc(abs(est_base))) %>%
    slice_head(n = label_n) %>%
    mutate(bubble_label = abbreviate_lollipop_label(omic_clean, max_chars = label_max_chars))

  x_limits <- pad_continuous_limits(
    plot_df$est_base,
    lower_mult = 0.18,
    upper_mult = 0.18,
    min_pad = 0.035,
    include_zero = TRUE
  )
  y_limits <- pad_continuous_limits(
    plot_df$neg_log10_fdr,
    lower_mult = 0.06,
    upper_mult = 0.18,
    min_pad = 0.18,
    include_zero = FALSE
  )
  y_limits[[1]] <- max(0, y_limits[[1]])

  ggplot(plot_df, aes(x = est_base, y = neg_log10_fdr)) +
    geom_point(
      aes(size = bubble_size, fill = omic_category, shape = QTL_type),
      colour = "black",
      stroke = 0.25,
      alpha = 0.92
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
    geom_text_repel(
      data = label_df,
      aes(label = bubble_label),
      size = pt_to_mm(5.2),
      family = NATURE$font,
      fontface = "bold",
      box.padding = 0.24,
      point.padding = 0.12,
      max.overlaps = Inf,
      segment.size = 0.15,
      min.segment.length = 0,
      force = 1.15,
      force_pull = 0.25,
      xlim = x_limits,
      ylim = y_limits
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    scale_shape_manual(values = QTL_SHAPES, name = "QTL type", drop = FALSE) +
    scale_fill_manual(values = CATEGORY_COLORS, name = "Omic category", drop = FALSE) +
    scale_size_continuous(range = c(1.2, 3.4), guide = "none") +
    labs(title = title, x = "Normalized Bayes K", y = expression(-log[10](FDR))) +
    theme_nature(legend_position = "top") +
    guides(
      fill = guide_legend(
        order = 1,
        nrow = 1,
        title.position = "top",
        override.aes = list(shape = 21, size = 2.1, colour = "black", alpha = 1)
      ),
      shape = guide_legend(
        order = 2,
        nrow = 1,
        title.position = "top",
        override.aes = list(size = 2.1, fill = "white", colour = "black", alpha = 1)
      )
    ) +
    theme(
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = c(0, 0.5),
      legend.direction = "horizontal",
      legend.spacing.y = unit(0.35, "mm"),
      legend.margin = margin(0, 0, 1, 0),
      legend.box.margin = margin(0, 0, 1, 0),
      plot.title = element_text(margin = margin(b = 3)),
      plot.margin = margin(7, 12, 6, 8)
    )
}

prepare_ad_age_strata_plot_data <- function(
  unique_pairs_age,
  prs_to_plot = DOMAIN_SETS$AD$prs,
  repro_threshold = 50,
  label_max_chars = 40
) {
  plot_df <- prepare_domain_qtl_data(
    unique_pairs_age,
    prs_to_plot = prs_to_plot,
    stratified = TRUE,
    repro_threshold = repro_threshold
  ) %>%
    filter(AgeGroup %in% AGE_SET_LEVELS) %>%
    mutate(
      AgeSet = factor(AgeGroup, levels = AGE_SET_LEVELS),
      AgeSetIndex = unname(AGE_SET_X[as.character(AgeSet)]),
      effect_direction = factor(if_else(est < 0, "Negative", "Positive"), levels = c("Negative", "Positive")),
      effect_size = abs(est),
      feature_label_full = clean_omic_label(omic_original_name_lower, title_case = FALSE),
      feature_label = abbreviate_lollipop_label(feature_label_full, max_chars = label_max_chars)
    ) %>%
    group_by(feature_label_full, AgeSet) %>%
    slice_max(effect_size, n = 1, with_ties = FALSE) %>%
    ungroup()

  if (nrow(plot_df) == 0) return(plot_df)

  feature_summary <- plot_df %>%
    group_by(feature_label_full, feature_label) %>%
    summarise(
      age_pattern = paste(AGE_SET_LEVELS[AGE_SET_LEVELS %in% as.character(AgeSet)], collapse = " + "),
      max_effect_size = max(effect_size, na.rm = TRUE),
      mean_est = mean(est, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pattern_priority = case_when(
        age_pattern == ">45" ~ 1L,
        str_detect(age_pattern, "All") &
          (str_detect(age_pattern, fixed("<45")) | str_detect(age_pattern, fixed(">45"))) ~ 2L,
        age_pattern == "All" ~ 3L,
        age_pattern == "<45" ~ 4L,
        TRUE ~ 5L
      )
    ) %>%
    arrange(pattern_priority, desc(max_effect_size), feature_label) %>%
    mutate(
      row_top = row_number(),
      row_plot = n() - row_top + 1
    )

  plot_df %>%
    left_join(feature_summary, by = c("feature_label_full", "feature_label")) %>%
    mutate(
      AgeSet = factor(AgeSet, levels = AGE_SET_LEVELS),
      AgeSetIndex = as.numeric(AgeSetIndex)
    )
}

plot_ad_age_strata <- function(ad_age_df, title = NULL, show_legend = TRUE) {
  if (nrow(ad_age_df) == 0) return(NULL)

  feature_axis <- ad_age_df %>%
    distinct(feature_label_full, feature_label, row_plot) %>%
    arrange(row_plot)

  separator_df <- ad_age_df %>%
    distinct(feature_label_full, pattern_priority, row_plot) %>%
    group_by(pattern_priority) %>%
    summarise(yintercept = min(row_plot) - 0.5, .groups = "drop") %>%
    filter(yintercept > 0.5)

  legend_position_use <- if (isTRUE(show_legend)) "bottom" else "none"

  ggplot(ad_age_df, aes(x = AgeSetIndex, y = row_plot)) +
    geom_hline(
      data = separator_df,
      aes(yintercept = yintercept),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.25,
      colour = "grey70"
    ) +
    geom_point(
      aes(size = effect_size, fill = effect_direction),
      shape = 21,
      colour = "black",
      stroke = 0.22,
      alpha = 0.95
    ) +
    scale_x_continuous(
      breaks = unname(AGE_SET_X),
      labels = names(AGE_SET_X),
      limits = c(0.5, length(AGE_SET_LEVELS) + 0.5),
      expand = expansion(add = 0),
      position = "top"
    ) +
    scale_y_continuous(
      breaks = feature_axis$row_plot,
      labels = feature_axis$feature_label,
      limits = c(0.5, nrow(feature_axis) + 0.5),
      expand = expansion(add = 0)
    ) +
    scale_fill_manual(
      values = EFFECT_COLORS[c("Negative", "Positive")],
      labels = c("Negative", "Positive"),
      name = "Association"
    ) +
    scale_size_continuous(
      range = c(1.1, 3.2),
      name = "|Normalized Bayes K|"
    ) +
    labs(title = title, x = "Age set", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_nature(legend_position = legend_position_use, grid = FALSE) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(size = 2.3, alpha = 1)),
      size = guide_legend(order = 2)
    ) +
    theme(
      panel.background = element_rect(fill = "grey96", colour = NA),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x.top = element_text(size = NATURE$base_pt, face = "bold", margin = margin(b = 2)),
      axis.text.x.top = element_text(size = NATURE$base_pt, face = "bold", colour = "black"),
      axis.text.x.bottom = element_blank(),
      axis.text.y = element_text(size = NATURE$small_pt, face = "bold", colour = "grey25"),
      legend.box = "vertical",
      legend.margin = margin(1, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(7, 5, 5, 5)
    )
}

plot_neuro_fingerprint <- function(unique_pairs, unique_pairs_age) {
  # This panel intentionally uses the full curated neurodegeneration fingerprint
  # from the original analysis script. These are hand-selected associations, not
  # a reproducibility-threshold top-N table, so keep the whitelist complete.

  allowed_age <- tribble(
    ~PRS,                                     ~omic,
    "amyotrophic_lateral_sclerosis",          "FERRITIN",
    "amyotrophic_lateral_sclerosis_Age_>_45", "FERRITIN",
    "amyotrophic_lateral_sclerosis",          "interleukin 5 (colony-stimulating factor, eosinophil)",
    "amyotrophic_lateral_sclerosis",          "interleukin 13",
    "amyotrophic_lateral_sclerosis",          "myeloperoxidase",
    "Alzheimer_disease_late_onset_PGS000053", "OMEGA_6_TOTAL",
    "Alzheimer_disease_late_onset_PGS000053", "chemokine (C-C motif) ligand 13",
    "Alzheimer_disease_late_onset_PGS000053", "chemokine (C-C motif) ligand 24",
    "Alzheimer_disease_late_onset_PGS000053", "heme oxygenase (decycling) 1",
    "Alzheimers_PGS002753_Age_>_45",          "OMEGA_3_TOTAL",
    "Alzheimers_PGS002753_Age_>_45",          "OMEGA_6_TOTAL",
    "Alzheimers_PGS002753",                   "OMEGA_6_TOTAL",
    "Parkinson_disease_PGS000903",            "4-hydroxyphenylpyruvate",
    "Parkinson_disease_PGS000903_Age_<_45",   "DHA",
    "Alzheimers_PGS002753",                   "interleukin 12B (natural killer cell stimulatory factor 2, cytotoxic lymphocyte maturation factor 2, p40)"
  )

  allowed_core <- tribble(
    ~PRS,                                       ~omic,
    "amyotrophic_lateral_sclerosis",            "FERRITIN",
    "amyotrophic_lateral_sclerosis",            "OMEGA_6_TOTAL",
    "amyotrophic_lateral_sclerosis",            "interleukin 5 (colony-stimulating factor, eosinophil)",
    "amyotrophic_lateral_sclerosis",            "interleukin 13",
    "amyotrophic_lateral_sclerosis",            "interleukin 17C",
    "amyotrophic_lateral_sclerosis",            "myeloperoxidase",
    "amyotrophic_lateral_sclerosis",            "interleukin 16",
    "amyotrophic_lateral_sclerosis",            "interleukin 22 receptor, alpha 1",
    "amyotrophic_lateral_sclerosis",            "1-(1-enyl-stearoyl)-2-arachidonoyl-GPE (P-18:0/20:4)*",
    "amyotrophic_lateral_sclerosis",            "1-(1-enyl-stearoyl)-2-docosapentaenoyl-GPE (P-18:0/22:5n3)*",
    "amyotrophic_lateral_sclerosis",            "1-(1-enyl-stearoyl)-2-arachidonoyl-GPC (P-18:0/20:4)",
    "Alzheimer_disease_late_onset_PGS000053",   "chemokine (C-C motif) ligand 2",
    "Alzheimer_disease_late_onset_PGS000053",   "chemokine (C-C motif) ligand 13",
    "Alzheimer_disease_late_onset_PGS000053",   "chemokine (C-C motif) ligand 24",
    "Alzheimer_disease_late_onset_PGS000053",   "interleukin 5 (colony-stimulating factor, eosinophil)",
    "Alzheimer_disease_late_onset_PGS000053",   "2-oleoyl-GPE (18:1)*",
    "Alzheimer_disease_late_onset_PGS000053",   "heme oxygenase (decycling) 1",
    "Alzheimers_PGS002753",                     "OMEGA_6_TOTAL",
    "Alzheimers_PGS002753",                     "interleukin 10",
    "Alzheimers_PGS002753",                     "interleukin 12B",
    "Alzheimers_PGS002753",                     "FERRITIN",
    "Alzheimers_PGS002753",                     "1-stearoyl-2-oleoyl-GPE (18:0/18:1)",
    "Alzheimers_PGS002753",                     "1-oleoyl-2-linoleoyl-GPE (18:1/18:2)*",
    "Alzheimers_PGS002753",                     "1-palmitoyl-2-oleoyl-GPE (16:0/18:1)",
    "Alzheimers_PGS002753",                     "1-(1-enyl-stearoyl)-2-oleoyl-GPE (P-18:0/18:1)",
    "Alzheimers_PGS002753",                     "tumor necrosis factor",
    "Parkinson_disease_PGS000903",              "4-hydroxyphenylpyruvate",
    "Parkinson_disease_PGS000903",              "2-stearoyl-GPE (18:0)*",
    "Parkinson_disease_PGS000903",              "1-palmitoleoyl-2-linolenoyl-GPC (16:1/18:3)*",
    "Parkinson_disease_PGS000903",              "1-(1-enyl-palmitoyl)-2-linoleoyl-GPC (P-16:0/18:2)*",
    "Parkinson_disease_PGS000903",              "chemokine (C-C motif) ligand 16",
    "Parkinson_disease_PGS000903",              "interleukin 12B",
    "Parkinson_disease_PGS000903",              "tumor necrosis factor"
  )

  df_age_filt <- unique_pairs_age %>%
    rename(PRS = feature_a_original_name, omic = omic_original_name) %>%
    inner_join(allowed_age, by = c("PRS", "omic"))

  df_core_filt <- unique_pairs %>%
    rename(PRS = feature_a_original_name, omic = omic_original_name) %>%
    inner_join(allowed_core, by = c("PRS", "omic"))

  priority_feats <- c(
    "OMEGA_6_TOTAL", "OMEGA_3_TOTAL", "DHA",
    "myeloperoxidase", "heme oxygenase (decycling) 1",
    "FERRITIN", "4-hydroxyphenylpyruvate"
  )

  eff_df <- bind_rows(df_age_filt, df_core_filt) %>%
    mutate(
      PRS_label = case_when(
        str_detect(PRS, "Alzheimer") & str_detect(PRS, "Age_>_45") ~ "AD PRS Age >45",
        str_detect(PRS, "Alzheimer") ~ "AD PRS",
        str_detect(PRS, "Parkinson") & str_detect(PRS, "Age_<_45") ~ "PD PRS Age <45",
        str_detect(PRS, "Parkinson") ~ "PD PRS",
        str_detect(PRS, "amyotrophic") & str_detect(PRS, "Age_>_45") ~ "ALS PRS Age >45",
        str_detect(PRS, "amyotrophic") ~ "ALS PRS",
        TRUE ~ PRS
      )
    ) %>%
    distinct(PRS_label, omic, .keep_all = TRUE) %>%
    mutate(
      PRS_label = factor(
        PRS_label,
        levels = c("AD PRS", "AD PRS Age >45", "PD PRS", "PD PRS Age <45", "ALS PRS", "ALS PRS Age >45")
      ),
      omic = fct_rev(fct_relevel(omic, priority_feats))
    )

  ggplot(eff_df, aes(x = PRS_label, y = omic, fill = norm_bayes_k)) +
    geom_tile(width = 0.92, height = 0.90, colour = "white", linewidth = 0.18) +
    scale_fill_gradient2(
      low = "#6A3DFF",
      mid = "white",
      high = "#FF2B14",
      midpoint = 0,
      name = "Normalized\nBayes K",
      guide = guide_colorbar(
        barheight = unit(18, "mm"),
        barwidth = unit(3.2, "mm"),
        ticks.colour = "white",
        frame.colour = NA
      )
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_nature(legend_position = "right", grid = FALSE) +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(linewidth = 0.25, colour = "black"),
      axis.line.y = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text.x = element_text(size = 6.0, face = "bold", angle = 45, hjust = 1, vjust = 1, colour = "black"),
      axis.text.y = element_text(size = 5.8, colour = "black"),
      legend.title = element_text(size = 6.5, face = "bold", lineheight = 0.9),
      legend.text = element_text(size = 5.5, colour = "black"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 1),
      plot.margin = margin(3, 4, 4, 3)
    )
}

prepare_neuro_all_feature_matrix <- function(
  domain_df,
  panel_name,
  label_max_chars = NULL
) {
  if (is.null(domain_df) || nrow(domain_df) == 0) {
    stop("No domain data were supplied for ", panel_name, call. = FALSE)
  }

  model_levels <- c("Residualized", "All", "<45", ">45")

  label_key <- make_lollipop_label_key(
    domain_df,
    label_max_chars = label_max_chars
  ) %>%
    select(
      omic_original_name_lower,
      omic_category,
      full_label,
      plotted_label,
      label_max_chars,
      label_was_abbreviated
    ) %>%
    distinct(omic_original_name_lower, .keep_all = TRUE)

  feature_order <- domain_df %>%
    mutate(
      est_numeric = suppressWarnings(as.numeric(est)),
      abs_est = abs(est_numeric)
    ) %>%
    filter(is.finite(est_numeric)) %>%
    group_by(omic_original_name_lower) %>%
    summarise(
      max_abs_est = max(abs_est, na.rm = TRUE),
      mean_est = mean(est_numeric, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(label_key, by = "omic_original_name_lower") %>%
    arrange(desc(max_abs_est), plotted_label) %>%
    mutate(
      feature_rank = row_number(),
      feature_label = make.unique(as.character(plotted_label), sep = " ")
    )

  matrix_df <- domain_df %>%
    mutate(
      ModelDisplay = factor(as.character(Model), levels = model_levels),
      est = suppressWarnings(as.numeric(est)),
      abs_est = abs(est),
      effect_direction = factor(
        if_else(est < 0, "Negative", "Positive"),
        levels = c("Negative", "Positive")
      ),
      source_feature = as.character(omic_original_name_lower),
      source_prs = as.character(feature_a_original_name)
    ) %>%
    filter(!is.na(ModelDisplay), is.finite(est)) %>%
    group_by(omic_original_name_lower, ModelDisplay) %>%
    slice_max(order_by = abs_est, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(
      feature_order %>%
        select(
          omic_original_name_lower,
          feature_label,
          full_label,
          plotted_label,
          feature_rank,
          max_abs_est
        ),
      by = "omic_original_name_lower"
    ) %>%
    mutate(
      feature_label = factor(
        feature_label,
        levels = rev(feature_order$feature_label)
      ),
      ModelDisplay = factor(ModelDisplay, levels = model_levels),
      QTL_type = factor(as.character(QTL_type), levels = names(QTL_SHAPES)),
      omic_category = factor(
        as.character(omic_category),
        levels = names(CATEGORY_COLORS)
      )
    )

  list(
    data = matrix_df,
    label_key = feature_order,
    n_features = nrow(feature_order)
  )
}

plot_neuro_all_feature_matrix <- function(
  matrix_df,
  title,
  size_limits = NULL,
  show_legend = FALSE,
  axis_text_y_pt = NULL
) {
  if (is.null(matrix_df) || nrow(matrix_df) == 0) return(NULL)

  model_levels <- c("Residualized", "All", "<45", ">45")
  feature_levels <- levels(matrix_df$feature_label)
  n_features <- length(feature_levels)

  if (is.null(axis_text_y_pt)) {
    axis_text_y_pt <- case_when(
      n_features >= 45 ~ 4.0,
      n_features >= 38 ~ 4.3,
      n_features >= 30 ~ 4.6,
      TRUE ~ 5.0
    )
  }

  background_df <- tidyr::expand_grid(
    feature_label = factor(feature_levels, levels = feature_levels),
    ModelDisplay = factor(model_levels, levels = model_levels)
  )

  if (is.null(size_limits)) {
    max_abs <- max(matrix_df$abs_est, na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
    size_limits <- c(0, max_abs)
  }

  legend_position_use <- if (isTRUE(show_legend)) "bottom" else "none"

  ggplot(background_df, aes(x = ModelDisplay, y = feature_label)) +
    geom_tile(
      fill = "grey96",
      colour = "white",
      linewidth = 0.30,
      width = 0.96,
      height = 0.92
    ) +
    geom_point(
      data = matrix_df,
      aes(
        size = abs_est,
        fill = effect_direction,
        shape = QTL_type
      ),
      colour = "black",
      stroke = 0.25,
      alpha = 0.97,
      na.rm = TRUE
    ) +
    scale_x_discrete(
      drop = FALSE,
      position = "top",
      labels = c(
        "Residualized" = "Age-
residualized",
        "All" = "Age-included
all ages",
        "<45" = "Age-included
<45",
        ">45" = "Age-included
>45"
      )
    ) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_manual(
      values = EFFECT_COLORS[c("Negative", "Positive")],
      name = "Association",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = QTL_SHAPES,
      name = "QTL type",
      drop = FALSE
    ) +
    scale_size_continuous(
      limits = size_limits,
      range = c(0.9, 3.4),
      name = "|Normalized Bayes K|"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_nature(legend_position = legend_position_use, grid = FALSE) +
    guides(
      fill = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(
          shape = 21,
          size = 2.3,
          colour = "black",
          alpha = 1
        )
      ),
      shape = guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(
          size = 2.3,
          fill = "white",
          colour = "black",
          alpha = 1
        )
      ),
      size = guide_legend(order = 3, nrow = 1)
    ) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x.top = element_text(
        size = 5.4,
        face = "bold",
        lineheight = 0.88,
        hjust = 0.5,
        vjust = 0,
        margin = margin(b = 3)
      ),
      axis.text.x.bottom = element_blank(),
      axis.text.y = element_text(
        size = axis_text_y_pt,
        face = "bold",
        colour = "grey20",
        lineheight = 0.84,
        margin = margin(r = 2)
      ),
      plot.title = element_text(
        size = NATURE$title_pt,
        face = "bold",
        hjust = 0,
        margin = margin(b = 5)
      ),
      legend.position = legend_position_use,
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = c(0, 0.5),
      legend.spacing.y = unit(0.25, "mm"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

extract_plot_legend <- function(plot_object) {
  plot_grob <- ggplotGrob(plot_object)
  guide_indices <- which(str_detect(plot_grob$layout$name, "^guide-box"))

  for (idx in guide_indices) {
    candidate <- plot_grob$grobs[[idx]]
    if (!inherits(candidate, "zeroGrob")) return(candidate)
  }

  NULL
}
