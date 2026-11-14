#!/usr/bin/env Rscript
packages <- c(
  "tidyverse", "readr", "stringr", "forcats", "ggplot2", "ggrepel",
  "circlize", "text2vec", "pheatmap", "vegan", "ggnewscale",
  "patchwork", "grid", "ragg"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
message("R dependencies are installed.")
