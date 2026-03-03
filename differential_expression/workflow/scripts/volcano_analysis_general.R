#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

option_list <- list(
  make_option(c("-r", "--results"), type = "character", dest = "results", help = "edgeR results TSV.gz"),
  make_option(c("-m", "--mode"), type = "character", dest = "mode", help = "Mode: plot or stats"),
  make_option(c("-o", "--output"), type = "character", dest = "output", help = "Output path"),
  make_option(c("-l", "--lfc-fold"), type = "double", dest = "lfc_fold", default = 1.5, help = "Fold-change threshold applied as log2(fold)"),
  make_option(c("-f", "--fdr"), type = "double", dest = "fdr_threshold", default = 0.01, help = "FDR significance threshold")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$results) || is.null(opt$mode) || is.null(opt$output)) {
  stop("--results, --mode, and --output are required", call. = FALSE)
}

mode <- match.arg(opt$mode, choices = c("plot", "stats"))

data <- read_tsv(opt$results, progress = FALSE)

required_cols <- c("contrast", "stratum", "condition", "logFC", "padj_global")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in results: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

lfc_fold <- ifelse(is.null(opt$lfc_fold), 1.5, opt$lfc_fold)
if (is.na(lfc_fold) || lfc_fold <= 0) {
  stop("--lfc-fold must be > 0", call. = FALSE)
}
lfc_threshold <- log2(lfc_fold)

fdr_threshold <- ifelse(is.null(opt$fdr_threshold), 0.01, opt$fdr_threshold)
if (is.na(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold >= 1) {
  stop("--fdr must be between 0 and 1", call. = FALSE)
}

data <- data %>%
  mutate(
    direction = case_when(
      padj_global < fdr_threshold & logFC > lfc_threshold ~ "up",
      padj_global < fdr_threshold & logFC < -lfc_threshold ~ "down",
      TRUE ~ "nonsign"
    )
  )

contrast_meta <- data %>% distinct(contrast, stratum, condition)

summary_stats <- data %>%
  group_by(contrast, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(contrast, direction, fill = list(n = 0)) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  left_join(contrast_meta, by = "contrast")

if (mode == "stats") {
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  write_tsv(summary_stats, opt$output)
  quit(save = "no")
}

label_map <- setNames(
  vapply(summary_stats$contrast, function(ct) {
    row <- summary_stats %>% filter(contrast == ct)
    sprintf(
      "%s | %s\nup: %d; down: %d; non-sign: %d",
      row$stratum[[1]],
      row$condition[[1]],
      row$up[[1]],
      row$down[[1]],
      row$nonsign[[1]]
    )
  }, character(1)),
  summary_stats$contrast
)

data$contrast <- factor(data$contrast, levels = summary_stats$contrast)

palette <- c(up = "#d73027", down = "#4575b4", nonsign = "#bdbdbd")

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

volcano <- ggplot(data, aes(x = logFC, y = -log10(padj_global))) +
  geom_point(aes(color = direction), alpha = 0.5, size = 0.6) +
  scale_color_manual(values = palette, name = "Direction") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey40") +
  facet_wrap(~ contrast, ncol = 2, labeller = as_labeller(label_map)) +
  labs(
    title = "Differential expression volcano plots (general)",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10)
  )

ggsave(opt$output, volcano, width = 7, height = 14, bg = "white")
