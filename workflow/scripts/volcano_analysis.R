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
  make_option(c("-r", "--results"), type = "character", dest = "results", help = "edgeR results TSV.gz"),
  make_option(c("-m", "--mode"), type = "character", dest = "mode", help = "Mode: plot or stats"),
  make_option(c("-o", "--output"), type = "character", dest = "output", help = "Output path")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$results) || is.null(opt$mode) || is.null(opt$output)) {
  stop("--results, --mode, and --output are required", call. = FALSE)
}

mode <- match.arg(opt$mode, choices = c("plot", "stats"))

data <- read_tsv(opt$results, progress = FALSE)

lfc_threshold <- log2(1.5)
fdr_threshold <- 0.001

season_labels <- c(
  Autumn = "Autumn",
  Winter = "Winter",
  Spring = "Spring",
  Summer = "Summer"
)
condition_labels <- c(
  Cells = "Cells",
  Aggregates = "Aggregates"
)

data <- data %>%
  mutate(
    direction = case_when(
      padj_global < fdr_threshold & logFC > lfc_threshold ~ "up",
      padj_global < fdr_threshold & logFC < -lfc_threshold ~ "down",
      TRUE ~ "nonsign"
    )
  )

summary_stats <- data %>%
  group_by(contrast, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(contrast, direction, fill = list(n = 0)) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)

if (mode == "stats") {
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  write_tsv(summary_stats, opt$output)
  quit(save = "no")
}

label_map <- setNames(
  vapply(summary_stats$contrast, function(ct) {
    row <- summary_stats %>% filter(contrast == ct)
    parts <- str_split(ct, "_", simplify = TRUE)
    season <- parts[1]
    condition <- if (ncol(parts) >= 2) parts[2] else ""
    season_name <- if (season %in% names(season_labels)) season_labels[[season]] else season
    condition_name <- if (condition %in% names(condition_labels)) condition_labels[[condition]] else condition
    sprintf(
      "%s %s vs Body\nup: %d; down: %d; non-sign: %d",
      season_name,
      condition_name,
      row$up,
      row$down,
      row$nonsign
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
  facet_wrap(~ contrast, ncol = 4, labeller = as_labeller(label_map)) +
  labs(
    title = "Differential expression volcano plots",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10)
  )

ggsave(opt$output, volcano, width = 14, height = 7, bg = "white")
