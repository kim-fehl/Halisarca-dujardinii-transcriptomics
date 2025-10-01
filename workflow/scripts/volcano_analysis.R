#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

option_list <- list(
  make_option("--results", type = "character", help = "edgeR results TSV.gz"),
  make_option("--mode", type = "character", help = "Mode: plot or stats"),
  make_option("--output", type = "character", help = "Output path")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$results) || is.null(opt$mode) || is.null(opt$output)) {
  stop("--results, --mode, and --output are required", call. = FALSE)
}

mode <- match.arg(opt$mode, choices = c("plot", "stats"))

data <- read_tsv(opt$results, progress = FALSE)

lfc_threshold <- log2(2)
padj_threshold <- 0.05

data <- data %>%
  mutate(
    direction = case_when(
      padj_global < padj_threshold & logFC > lfc_threshold ~ "up",
      padj_global < padj_threshold & logFC < -lfc_threshold ~ "down",
      TRUE ~ "nonsign"
    )
  )

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

if (mode == "stats") {
  summary_tbl <- data %>%
    group_by(contrast, direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(direction = factor(direction, levels = c("up", "down", "nonsign"))) %>%
    arrange(contrast, direction) %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  write_tsv(summary_tbl, opt$output)
} else {
  palette <- c(up = "#d73027", down = "#4575b4", nonsign = "#bdbdbd")
  volcano <- ggplot(data, aes(x = logFC, y = -log10(padj_global))) +
    geom_point(aes(color = direction), alpha = 0.4, size = 0.5) +
    scale_color_manual(values = palette, name = "Direction") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey40") +
    facet_wrap(~contrast, scales = "free_y") +
    labs(
      title = "Differential expression volcano plots",
      x = "log2 fold change",
      y = "-log10 adjusted p-value"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(opt$output, volcano, width = 10, height = 12, bg = "white")
}
