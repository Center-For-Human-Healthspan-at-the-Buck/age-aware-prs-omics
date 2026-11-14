# Shared utilities for Nogal et al. PRS-omics figure reproduction.

required_packages <- c(
  "tidyverse", "readr", "stringr", "forcats", "ggplot2", "ggrepel",
  "circlize", "text2vec", "pheatmap", "vegan", "ggnewscale",
  "patchwork", "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the missing R packages before running the repository scripts: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(circlize)
  library(text2vec)
  library(pheatmap)
  library(vegan)
  library(ggnewscale)
  library(patchwork)
  library(grid)
})

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "README.md")) &&
        dir.exists(file.path(current, "scripts")) &&
        dir.exists(file.path(current, "data"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the repository root from: ", start, call. = FALSE)
}

parse_path_args <- function(default_output_subdir) {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(flag, default) {
    idx <- match(flag, args)
    if (is.na(idx)) return(default)
    if (idx == length(args)) stop("Missing value after ", flag, call. = FALSE)
    args[[idx + 1]]
  }

  root <- get_arg("--repo-root", find_repo_root())
  list(
    repo_root = normalizePath(root, winslash = "/", mustWork = TRUE),
    data_dir = normalizePath(
      get_arg("--data-dir", file.path(root, "data")),
      winslash = "/", mustWork = TRUE
    ),
    output_dir = normalizePath(
      get_arg("--output-dir", file.path(root, "outputs", default_output_subdir)),
      winslash = "/", mustWork = FALSE
    )
  )
}

NATURE <- list(
  font = "Helvetica",
  base_pt = 6.5,
  small_pt = 5.2,
  title_pt = 7.0,
  dpi = 300,
  width_1col_mm = 88,
  width_2col_mm = 180,
  max_width_mm = 180,
  max_height_mm = 230
)

NATURE_COLORS <- c(
  grey = "#4D4D4D", light_grey = "#BDBDBD", blue = "#0072B2",
  sky_blue = "#56B4E9", vermillion = "#D55E00", orange = "#E69F00",
  green = "#009E73", purple = "#CC79A7", black = "#000000"
)

CATEGORY_COLORS <- c(
  chemistries = NATURE_COLORS[["vermillion"]],
  metabolomics = NATURE_COLORS[["green"]],
  proteomics = NATURE_COLORS[["orange"]]
)

MODEL_COLORS <- c(
  Residualized = NATURE_COLORS[["grey"]],
  `<45` = NATURE_COLORS[["blue"]],
  `>45` = NATURE_COLORS[["vermillion"]],
  All = NATURE_COLORS[["orange"]]
)

EFFECT_COLORS <- c(Positive = "#1B7837", Negative = "#D73027")
QTL_SHAPES <- c(none = 21, eQTL = 24, pQTL = 22, both = 23, `Not annotated` = 4)
AGE_SET_LEVELS <- c("<45", "All", ">45")
AGE_SET_X <- setNames(seq_along(AGE_SET_LEVELS), AGE_SET_LEVELS)

pt_to_mm <- function(pt) pt / ggplot2::.pt
mm_to_in <- function(mm) mm / 25.4

theme_nature <- function(
  base_size = NATURE$base_pt,
  base_family = NATURE$font,
  legend_position = "right",
  grid = TRUE
) {
  out <- theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(family = base_family, colour = "black"),
      plot.title = element_text(
        size = NATURE$title_pt, face = "bold", hjust = 0,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(size = base_size, hjust = 0, margin = margin(b = 2)),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = max(NATURE$small_pt, base_size - 0.5), colour = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = max(NATURE$small_pt, base_size - 0.5)),
      strip.text = element_text(size = base_size, face = "bold"),
      line = element_line(linewidth = 0.25),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      legend.key.size = unit(3.2, "mm"),
      legend.box.spacing = unit(1.5, "mm"),
      plot.margin = margin(3, 3, 3, 3),
      legend.position = legend_position
    )
  if (isTRUE(grid)) {
    out <- out + theme(
      panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey90"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  }
  out
}

theme_set(theme_nature())

validate_figure_size <- function(width_mm, height_mm) {
  if (width_mm > NATURE$max_width_mm) {
    stop("Figure width exceeds 180 mm: ", width_mm, call. = FALSE)
  }
  if (height_mm > NATURE$max_height_mm) {
    warning("Figure height exceeds ", NATURE$max_height_mm, " mm.", call. = FALSE)
  }
}

check_written_file <- function(path) {
  info <- if (file.exists(path)) file.info(path) else NULL
  if (is.null(info) || is.na(info$size) || info$size <= 0) {
    stop("Output file was not written or is empty: ", path, call. = FALSE)
  }
  message("Wrote: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(path)
}

save_plot <- function(
  plot,
  file_stub,
  width_mm,
  height_mm,
  out_dir,
  formats = c("pdf", "png"),
  dpi = NATURE$dpi,
  bg = "white"
) {
  validate_figure_size(width_mm, height_mm)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  width_in <- mm_to_in(width_mm)
  height_in <- mm_to_in(height_mm)

  if ("pdf" %in% formats) {
    path <- file.path(out_dir, paste0(file_stub, ".pdf"))
    grDevices::pdf(
      path, width = width_in, height = height_in, family = "Helvetica",
      bg = bg, onefile = FALSE, useDingbats = FALSE
    )
    device_id <- grDevices::dev.cur()
    tryCatch(print(plot), finally = {
      if (grDevices::dev.cur() == device_id) grDevices::dev.off()
    })
    check_written_file(path)
  }

  if ("png" %in% formats) {
    path <- file.path(out_dir, paste0(file_stub, ".png"))
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggsave(path, plot = plot, width = width_in, height = height_in,
             units = "in", dpi = dpi, device = ragg::agg_png, bg = bg)
    } else {
      ggsave(path, plot = plot, width = width_in, height = height_in,
             units = "in", dpi = dpi, bg = bg)
    }
    check_written_file(path)
  }
  invisible(file_stub)
}

save_base_plot <- function(draw_fun, file_stub, width_mm, height_mm, out_dir,
                           formats = c("pdf", "png"), dpi = NATURE$dpi) {
  validate_figure_size(width_mm, height_mm)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if ("pdf" %in% formats) {
    path <- file.path(out_dir, paste0(file_stub, ".pdf"))
    grDevices::pdf(path, width = mm_to_in(width_mm), height = mm_to_in(height_mm),
                   family = "Helvetica", onefile = FALSE, useDingbats = FALSE)
    device_id <- grDevices::dev.cur()
    tryCatch(draw_fun(), finally = {
      if (grDevices::dev.cur() == device_id) grDevices::dev.off()
    })
    check_written_file(path)
  }

  if ("png" %in% formats) {
    path <- file.path(out_dir, paste0(file_stub, ".png"))
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = width_mm, height = height_mm,
                    units = "mm", res = dpi, background = "white")
    } else {
      grDevices::png(path, width = width_mm, height = height_mm,
                     units = "mm", res = dpi, bg = "white")
    }
    device_id <- grDevices::dev.cur()
    tryCatch(draw_fun(), finally = {
      if (grDevices::dev.cur() == device_id) grDevices::dev.off()
    })
    check_written_file(path)
  }
  invisible(file_stub)
}

save_pheatmap <- function(object, file_stub, width_mm, height_mm, out_dir) {
  save_base_plot(
    function() grid::grid.draw(object$gtable),
    file_stub, width_mm, height_mm, out_dir
  )
}

read_tsv_checked <- function(path, comment = NULL) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  if (is.null(comment)) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_tsv(path, comment = comment, show_col_types = FALSE)
  }
}

ensure_columns <- function(df, columns, label) {
  missing <- setdiff(columns, names(df))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(df)
}

ensure_omic_name_columns <- function(df) {
  if (!"omic_original_name" %in% names(df) && "omic_original_name_lower" %in% names(df)) {
    df$omic_original_name <- df$omic_original_name_lower
  }
  if (!"omic_original_name_lower" %in% names(df) && "omic_original_name" %in% names(df)) {
    df$omic_original_name_lower <- str_to_lower(df$omic_original_name)
  }
  df
}

normalise_bayes_k <- function(df) {
  max_value <- max(abs(as.numeric(df$avg_bayes_k)), na.rm = TRUE)
  if (!is.finite(max_value) || max_value <= 0) max_value <- 1
  df %>% mutate(
    avg_bayes_k = as.numeric(avg_bayes_k),
    sd_bayes_k = as.numeric(sd_bayes_k),
    norm_bayes_k = avg_bayes_k / max_value,
    norm_sd_bayes_k = sd_bayes_k / max_value
  )
}

assign_omic_category <- function(df, proteomics_meta, metabolomics_meta) {
  df <- ensure_omic_name_columns(df)
  prot <- str_to_lower(str_trim(as.character(proteomics_meta$gene_description)))
  metab <- str_to_lower(str_trim(as.character(metabolomics_meta$BIOCHEMICAL_NAME)))
  df %>% mutate(
    omic_original_name = str_trim(as.character(omic_original_name)),
    omic_original_name_lower = str_to_lower(str_trim(as.character(omic_original_name_lower))),
    omic_category = case_when(
      omic_original_name_lower %in% prot ~ "proteomics",
      omic_original_name_lower %in% metab ~ "metabolomics",
      TRUE ~ "chemistries"
    )
  )
}

valid_qtl_value <- function(x) {
  value <- str_to_lower(str_trim(as.character(x)))
  !is.na(value) & !value %in% c("", "nan", "na", "none", "null") &
    !str_detect(value, "^(nan|na|none)(;(nan|na|none))*$")
}

clean_prs_label <- function(x) {
  x %>%
    str_remove("_PGS\\d+($|_Age_.*$)") %>%
    str_remove("_Age_[><]_45$") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_title()
}

clean_omic_label <- function(x, title_case = TRUE) {
  out <- x %>% as.character() %>% str_replace_all("_", " ") %>% str_squish()
  if (isTRUE(title_case)) out <- str_to_title(out)
  out
}

pad_continuous_limits <- function(x, lower_mult = 0.08, upper_mult = 0.08,
                                  min_pad = 0.02, include_zero = FALSE) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) == 0) return(c(0, 1))
  rng <- range(x)
  if (isTRUE(include_zero)) rng <- range(c(rng, 0))
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) span <- max(abs(rng), 1) * 0.10
  c(rng[[1]] - max(span * lower_mult, min_pad),
    rng[[2]] + max(span * upper_mult, min_pad))
}

canonical_omic_key <- function(x) {
  x %>% as.character() %>% str_replace_all("[‐‑‒–—―−]", "-") %>%
    str_replace_all("_", " ") %>% str_replace_all("[[:space:]]+", " ") %>%
    str_squish() %>% str_to_lower()
}

strip_age_suffix <- function(x) {
  str_remove(as.character(x), regex("_Age_[<>]_?45$", ignore_case = TRUE))
}

load_discovery_data <- function(data_dir) {
  metabolomics <- read_tsv_checked(file.path(data_dir, "metadata", "metabolomics.tsv"), comment = "#")
  proteomics <- read_tsv_checked(file.path(data_dir, "metadata", "proteomics.tsv"), comment = "#")
  prs_correlation_raw <- read_tsv_checked(
    file.path(data_dir, "metadata", "prs_correlation.tsv")
  )
  ensure_columns(prs_correlation_raw, "prs_feature", "PRS correlation matrix")
  prs_names <- as.character(prs_correlation_raw$prs_feature)
  prs_correlation <- prs_correlation_raw %>%
    select(-prs_feature) %>%
    as.matrix()
  storage.mode(prs_correlation) <- "numeric"
  if (nrow(prs_correlation) != ncol(prs_correlation) ||
      nrow(prs_correlation) != length(prs_names)) {
    stop("PRS correlation input must be a square matrix with one row per prs_feature.", call. = FALSE)
  }
  if (!setequal(colnames(prs_correlation), prs_names)) {
    stop("PRS correlation row and column names do not match.", call. = FALSE)
  }
  prs_correlation <- prs_correlation[, prs_names, drop = FALSE]
  rownames(prs_correlation) <- prs_names

  prepare_vae <- function(path, label) {
    raw <- read_tsv_checked(path)
    ensure_columns(raw, c(
      "fdr", "feature_a_original_name", "omic_original_name", "avg_bayes_k",
      "sd_bayes_k", "percent_reproducibility", "Mapped_gene_e", "Mapped_gene_p"
    ), label)
    raw %>%
      ensure_omic_name_columns() %>%
      mutate(across(any_of(c(
        "fdr", "n_runs", "percent_reproducibility", "avg_bayes_k", "sd_bayes_k"
      )), ~suppressWarnings(as.numeric(.x)))) %>%
      normalise_bayes_k() %>%
      assign_omic_category(proteomics, metabolomics) %>%
      mutate(vae_model = label)
  }

  list(
    vae_residualized = prepare_vae(file.path(data_dir, "vae", "age_residualized.tsv"), "Age-residualized"),
    vae_included = prepare_vae(file.path(data_dir, "vae", "age_included.tsv"), "Age-included"),
    prs_correlation = prs_correlation,
    ols_positive_controls = read_tsv_checked(file.path(data_dir, "manuscript", "ols_positive_controls.tsv")) %>% mutate(age_group = "All"),
    metabolomics_meta = metabolomics,
    proteomics_meta = proteomics
  )
}


copy_static_panel <- function(source, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = TRUE)) {
    stop("Could not copy static panel: ", source, call. = FALSE)
  }
  check_written_file(destination)
}
