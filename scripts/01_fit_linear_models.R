#!/usr/bin/env Rscript

# Generic analyte-wise OLS workflow for continuous PRS and continuous omic traits.
# It fits the four model families used in the manuscript:
#   1. age-residualized OLS main effect
#   2. age-included OLS main effect
#   3. linear PRS x age interaction
#   4. joint PRS x natural-spline(age) interaction
#
# Example:
# Rscript scripts/01_fit_linear_models.R \
#   --matrix path/to/analysis_matrix.tsv \
#   --features data/templates/feature_manifest.tsv \
#   --prs data/templates/prs_manifest.tsv \
#   --age-column age \
#   --covariates sex,PC1,PC2 \
#   --output-dir outputs/linear_models

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    if (required) stop("Missing required argument: ", flag, call. = FALSE)
    return(default)
  }
  if (idx == length(args)) stop("Missing value after ", flag, call. = FALSE)
  args[[idx + 1]]
}

matrix_file <- get_arg("--matrix", required = TRUE)
feature_file <- get_arg("--features", required = TRUE)
prs_file <- get_arg("--prs", required = TRUE)
out_dir <- get_arg("--output-dir", "outputs/linear_models")
age_column <- get_arg("--age-column", "age")
id_column <- get_arg("--id-column", "")
covariates_arg <- get_arg("--covariates", "")
covariates <- if (covariates_arg == "") character(0) else str_split(covariates_arg, ",", simplify = TRUE) %>% as.character() %>% str_trim()
min_n <- as.integer(get_arg("--min-n", "100"))
spline_df <- as.integer(get_arg("--spline-df", "4"))
cores <- as.integer(get_arg("--cores", "1"))
fdr_scope <- match.arg(get_arg("--fdr-scope", "global"), c("global", "omics"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  if (str_ends(str_to_lower(path), ".csv")) {
    read_csv(path, show_col_types = FALSE)
  } else {
    read_tsv(path, show_col_types = FALSE)
  }
}

matrix <- read_table(matrix_file)
features <- read_table(feature_file)
prs_manifest <- read_table(prs_file)

required_feature_cols <- c("column", "omics_feature", "omics_dataset")
required_prs_cols <- c("column", "prs_feature")
if (length(setdiff(required_feature_cols, names(features))) > 0) {
  stop("Feature manifest requires columns: ", paste(required_feature_cols, collapse = ", "), call. = FALSE)
}
if (length(setdiff(required_prs_cols, names(prs_manifest))) > 0) {
  stop("PRS manifest requires columns: ", paste(required_prs_cols, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(features$column)) {
  stop("Feature manifest contains duplicate analysis-matrix column names.", call. = FALSE)
}
if (anyDuplicated(prs_manifest$column)) {
  stop("PRS manifest contains duplicate analysis-matrix column names.", call. = FALSE)
}
if (!is.finite(min_n) || min_n < 10) {
  stop("--min-n must be an integer of at least 10.", call. = FALSE)
}
if (!is.finite(spline_df) || spline_df < 2) {
  stop("--spline-df must be an integer of at least 2.", call. = FALSE)
}

id_columns <- if (nzchar(id_column)) id_column else character(0)
required_matrix_cols <- unique(c(id_columns, age_column, covariates, features$column, prs_manifest$column))
missing_matrix_cols <- setdiff(required_matrix_cols, names(matrix))
if (length(missing_matrix_cols) > 0) {
  stop("Analysis matrix is missing columns: ", paste(missing_matrix_cols, collapse = ", "), call. = FALSE)
}

# Convert modeled columns to numeric except explicitly supplied covariates, which
# may be numeric or categorical. Character covariates are converted to factors.
modeled_numeric <- unique(c(age_column, features$column, prs_manifest$column))
matrix <- matrix %>%
  mutate(across(all_of(modeled_numeric), ~suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(covariates), ~if (is.character(.x)) factor(.x) else .x))

bt <- function(x) paste0("`", gsub("`", "", x, fixed = TRUE), "`")
formula_rhs <- function(terms) {
  terms <- terms[!is.na(terms) & terms != ""]
  if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
}
normalize_coefficient_name <- function(x) {
  gsub("`", "", as.character(x), fixed = TRUE)
}
find_coefficient_name <- function(fit, term) {
  coefficient_names <- rownames(summary(fit)$coefficients)
  candidates <- unique(c(term, bt(term), make.names(term)))
  exact <- candidates[candidates %in% coefficient_names]
  if (length(exact) > 0) return(exact[[1]])
  normalized <- normalize_coefficient_name(coefficient_names)
  index <- match(normalize_coefficient_name(term), normalized)
  if (is.na(index)) NA_character_ else coefficient_names[[index]]
}
extract_term <- function(fit, term) {
  coefs <- summary(fit)$coefficients
  coefficient_name <- find_coefficient_name(fit, term)
  if (is.na(coefficient_name)) {
    return(c(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  c(beta = unname(coefs[coefficient_name, "Estimate"]),
    se = unname(coefs[coefficient_name, "Std. Error"]),
    p = unname(coefs[coefficient_name, "Pr(>|t|)"]))
}
fit_metrics <- function(fit) {
  s <- summary(fit)
  c(r_sq = unname(s$r.squared), adj_r_sq = unname(s$adj.r.squared))
}
empty_main <- function(meta, n) {
  tibble(!!!meta, N = n, beta = NA_real_, std_err = NA_real_, p_val = NA_real_,
         r_sq = NA_real_, adj_r_sq = NA_real_, error = NA_character_)
}

fit_one_pair <- function(feature_row, prs_row) {
  feature_col <- feature_row$column
  prs_col <- prs_row$column
  keep <- unique(c(feature_col, prs_col, age_column, covariates))
  dat <- matrix %>% select(all_of(keep)) %>% drop_na()
  n <- nrow(dat)

  meta <- list(
    omics_dataset = as.character(feature_row$omics_dataset),
    omics_feature = as.character(feature_row$omics_feature),
    prs_feature = as.character(prs_row$prs_feature),
    SUB_PATHWAY = if ("SUB_PATHWAY" %in% names(feature_row)) as.character(feature_row$SUB_PATHWAY) else NA_character_,
    SUPER_PATHWAY = if ("SUPER_PATHWAY" %in% names(feature_row)) as.character(feature_row$SUPER_PATHWAY) else NA_character_
  )

  modeled_sd <- c(
    outcome = sd(dat[[feature_col]]),
    prs = sd(dat[[prs_col]]),
    age = sd(dat[[age_column]])
  )
  if (n < min_n || any(!is.finite(modeled_sd) | modeled_sd == 0)) {
    base <- empty_main(meta, n)
    return(list(
      residualized = base,
      age_included = base,
      linear = base %>% mutate(interaction_beta = NA_real_, interaction_std_err = NA_real_, interaction_p_val = NA_real_),
      spline = tibble(
        !!!meta, N = n, prs_predictor = prs_col, spline_df = spline_df, spline_degree = 3,
        linear_beta_prs = NA_real_, linear_p_val_prs = NA_real_,
        linear_interaction_beta = NA_real_, linear_interaction_p_val = NA_real_,
        linear_r_sq = NA_real_, linear_adj_r_sq = NA_real_,
        spline_beta_prs = NA_real_, spline_p_val_prs = NA_real_,
        spline_interaction_F = NA_real_, spline_interaction_p_val = NA_real_,
        spline_interaction_df_num = NA_real_, spline_interaction_df_denom = NA_real_,
        spline_interaction_partial_r2 = NA_real_, spline_r_sq = NA_real_,
        spline_adj_r_sq = NA_real_, error = "insufficient data or zero variance"
      )
    ))
  }

  outcome <- bt(feature_col)
  prs_term <- bt(prs_col)
  age_term <- bt(age_column)
  cov_terms <- bt(covariates)

  result <- tryCatch({
    # 1) Residualize the outcome on age and covariates, then regress residuals on PRS.
    residualization_formula <- as.formula(paste(outcome, "~", formula_rhs(c(age_term, cov_terms))))
    residualization_fit <- lm(residualization_formula, data = dat)
    dat$.age_residualized_outcome <- residuals(residualization_fit)
    residualized_fit <- lm(as.formula(paste(".age_residualized_outcome ~", prs_term)), data = dat)
    resid_effect <- extract_term(residualized_fit, prs_col)
    resid_metrics <- fit_metrics(residualized_fit)

    # 2) Main effect with age retained as a covariate.
    included_formula <- as.formula(paste(outcome, "~", formula_rhs(c(prs_term, age_term, cov_terms))))
    included_fit <- lm(included_formula, data = dat)
    included_effect <- extract_term(included_fit, prs_col)
    included_metrics <- fit_metrics(included_fit)

    # 3) Linear interaction.
    interaction_formula <- as.formula(paste(
      outcome, "~", formula_rhs(c(paste0(prs_term, " * ", age_term), cov_terms))
    ))
    interaction_fit <- lm(interaction_formula, data = dat)
    interaction_name_candidates <- c(
      paste0(prs_col, ":", age_column),
      paste0(age_column, ":", prs_col)
    )
    coefficient_names <- rownames(summary(interaction_fit)$coefficients)
    normalized_names <- normalize_coefficient_name(coefficient_names)
    normalized_candidates <- normalize_coefficient_name(interaction_name_candidates)
    interaction_index <- match(normalized_candidates, normalized_names, nomatch = 0L)
    interaction_index <- interaction_index[interaction_index > 0L]
    interaction_name <- if (length(interaction_index) > 0) {
      coefficient_names[[interaction_index[[1]]]]
    } else {
      interaction_name_candidates[[1]]
    }
    linear_main <- extract_term(interaction_fit, prs_col)
    linear_interaction <- extract_term(interaction_fit, interaction_name)
    linear_metrics <- fit_metrics(interaction_fit)

    # 4) Natural-spline interaction, tested with a nested-model F test.
    ns_term <- paste0("splines::ns(", age_term, ", df = ", spline_df, ")")
    spline_reduced_formula <- as.formula(paste(
      outcome, "~", formula_rhs(c(prs_term, ns_term, cov_terms))
    ))
    spline_full_formula <- as.formula(paste(
      outcome, "~", formula_rhs(c(prs_term, ns_term, paste0(prs_term, ":", ns_term), cov_terms))
    ))
    spline_reduced <- lm(spline_reduced_formula, data = dat)
    spline_full <- lm(spline_full_formula, data = dat)
    spline_test <- anova(spline_reduced, spline_full)
    spline_effect <- extract_term(spline_full, prs_col)
    spline_metrics <- fit_metrics(spline_full)
    f_stat <- suppressWarnings(as.numeric(spline_test$F[[2]]))
    f_p <- suppressWarnings(as.numeric(spline_test$`Pr(>F)`[[2]]))
    df_num <- suppressWarnings(as.numeric(spline_test$Df[[2]]))
    df_denom <- df.residual(spline_full)
    rss_reduced <- deviance(spline_reduced)
    rss_full <- deviance(spline_full)
    partial_r2 <- if (is.finite(rss_reduced) && rss_reduced > 0) (rss_reduced - rss_full) / rss_reduced else NA_real_

    list(
      residualized = tibble(
        !!!meta, N = n, beta = resid_effect[["beta"]], std_err = resid_effect[["se"]],
        p_val = resid_effect[["p"]], r_sq = resid_metrics[["r_sq"]],
        adj_r_sq = resid_metrics[["adj_r_sq"]], error = NA_character_
      ),
      age_included = tibble(
        !!!meta, N = n, beta = included_effect[["beta"]], std_err = included_effect[["se"]],
        p_val = included_effect[["p"]], r_sq = included_metrics[["r_sq"]],
        adj_r_sq = included_metrics[["adj_r_sq"]], error = NA_character_
      ),
      linear = tibble(
        !!!meta, N = n, beta = linear_main[["beta"]], std_err = linear_main[["se"]],
        p_val = linear_main[["p"]], interaction_beta = linear_interaction[["beta"]],
        interaction_std_err = linear_interaction[["se"]],
        interaction_p_val = linear_interaction[["p"]],
        r_sq = linear_metrics[["r_sq"]], adj_r_sq = linear_metrics[["adj_r_sq"]],
        error = NA_character_
      ),
      spline = tibble(
        !!!meta, N = n, prs_predictor = prs_col, spline_df = spline_df, spline_degree = 3,
        linear_beta_prs = linear_main[["beta"]], linear_p_val_prs = linear_main[["p"]],
        linear_interaction_beta = linear_interaction[["beta"]],
        linear_interaction_p_val = linear_interaction[["p"]],
        linear_r_sq = linear_metrics[["r_sq"]], linear_adj_r_sq = linear_metrics[["adj_r_sq"]],
        spline_beta_prs = spline_effect[["beta"]], spline_p_val_prs = spline_effect[["p"]],
        spline_interaction_F = f_stat, spline_interaction_p_val = f_p,
        spline_interaction_df_num = df_num, spline_interaction_df_denom = df_denom,
        spline_interaction_partial_r2 = partial_r2,
        spline_r_sq = spline_metrics[["r_sq"]], spline_adj_r_sq = spline_metrics[["adj_r_sq"]],
        error = NA_character_
      )
    )
  }, error = function(e) {
    base <- empty_main(meta, n) %>% mutate(error = conditionMessage(e))
    list(
      residualized = base,
      age_included = base,
      linear = base %>% mutate(interaction_beta = NA_real_, interaction_std_err = NA_real_, interaction_p_val = NA_real_),
      spline = tibble(
        !!!meta, N = n, prs_predictor = prs_col, spline_df = spline_df, spline_degree = 3,
        linear_beta_prs = NA_real_, linear_p_val_prs = NA_real_,
        linear_interaction_beta = NA_real_, linear_interaction_p_val = NA_real_,
        linear_r_sq = NA_real_, linear_adj_r_sq = NA_real_,
        spline_beta_prs = NA_real_, spline_p_val_prs = NA_real_,
        spline_interaction_F = NA_real_, spline_interaction_p_val = NA_real_,
        spline_interaction_df_num = NA_real_, spline_interaction_df_denom = NA_real_,
        spline_interaction_partial_r2 = NA_real_, spline_r_sq = NA_real_,
        spline_adj_r_sq = NA_real_, error = conditionMessage(e)
      )
    )
  })
  result
}

jobs <- tidyr::crossing(feature_index = seq_len(nrow(features)), prs_index = seq_len(nrow(prs_manifest)))
run_job <- function(i) {
  fit_one_pair(features[jobs$feature_index[[i]], , drop = FALSE],
               prs_manifest[jobs$prs_index[[i]], , drop = FALSE])
}

message("Fitting ", nrow(jobs), " PRS-feature combinations...")
if (cores > 1 && .Platform$OS.type != "windows") {
  fits <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cores)
} else {
  fits <- lapply(seq_len(nrow(jobs)), run_job)
}

residualized <- bind_rows(lapply(fits, function(x) x$residualized))
age_included <- bind_rows(lapply(fits, function(x) x$age_included))
linear <- bind_rows(lapply(fits, function(x) x$linear))
spline <- bind_rows(lapply(fits, function(x) x$spline))

add_fdr <- function(df, p_column, q_column) {
  if (fdr_scope == "omics") {
    df %>% group_by(omics_dataset) %>%
      mutate(!!rlang::sym(q_column) := p.adjust(.data[[p_column]], method = "BH")) %>% ungroup()
  } else {
    df %>% mutate(!!rlang::sym(q_column) := p.adjust(.data[[p_column]], method = "BH"))
  }
}

residualized <- add_fdr(residualized, "p_val", "FDR_q")
age_included <- add_fdr(age_included, "p_val", "FDR_q")
linear <- add_fdr(linear, "p_val", "FDR_q") %>%
  add_fdr("interaction_p_val", "interaction_FDR_q")

spline <- spline %>%
  add_fdr("linear_p_val_prs", "linear_FDR_q_prs") %>%
  add_fdr("linear_interaction_p_val", "linear_interaction_FDR_q") %>%
  add_fdr("spline_p_val_prs", "spline_FDR_q_prs") %>%
  add_fdr("spline_interaction_p_val", "spline_interaction_FDR_q") %>%
  group_by(omics_dataset) %>%
  mutate(
    spline_interaction_FDR_q_within_omics = p.adjust(spline_interaction_p_val, method = "BH"),
    linear_interaction_FDR_q_within_omics = p.adjust(linear_interaction_p_val, method = "BH")
  ) %>% ungroup()

write_tsv(residualized, file.path(out_dir, "ols_main_age_residualized.tsv"), na = "")
write_tsv(age_included, file.path(out_dir, "ols_main_age_included.tsv"), na = "")
write_tsv(linear, file.path(out_dir, "ols_linear_age_interaction.tsv"), na = "")
write_tsv(spline, file.path(out_dir, "ols_spline_age_interaction.tsv"), na = "")

message("Completed. Results written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
