#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option(c("-i", "--input-rds"), type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option(c("-c", "--combat-seq"), type = "character", dest = "combat_seq", help = "RDS with ComBat-seq counts"),
  make_option(c("-r", "--combat-ref"), type = "character", dest = "combat_ref", help = "RDS with ComBat-ref counts"),
  make_option(c("-o", "--output"), type = "character", dest = "output", help = "Output PDF path")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input_rds) || is.null(opt$combat_seq) || is.null(opt$combat_ref) || is.null(opt$output)) {
  stop("--input-rds, --combat-seq, --combat-ref, and --output are required", call. = FALSE)
}

de_data <- readRDS(opt$input_rds)
combat_seq_counts <- readRDS(opt$combat_seq)
combat_ref_counts <- readRDS(opt$combat_ref)
counts <- de_data$counts
metadata <- de_data$metadata

stopifnot(all(colnames(counts) == metadata$sample_id))
stopifnot(all(colnames(combat_seq_counts) == metadata$sample_id))
stopifnot(all(colnames(combat_ref_counts) == metadata$sample_id))

logcpm_raw <- cpm(counts, log = TRUE, prior.count = 1)
logcpm_combat_seq <- cpm(combat_seq_counts, log = TRUE, prior.count = 1)
logcpm_combat_ref <- cpm(combat_ref_counts, log = TRUE, prior.count = 1)

calc_pca_df <- function(logcpm, label) {
  pca <- prcomp(t(logcpm), scale. = FALSE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  scores <- as_tibble(pca$x[, 1:2], rownames = "sample_id")
  left_join(scores, metadata, by = "sample_id") %>%
    mutate(
      dataset = label,
      pc1_var = scales::percent(var_expl[1], accuracy = 0.1),
      pc2_var = scales::percent(var_expl[2], accuracy = 0.1)
    )
}

suppressPackageStartupMessages(library(scales))

df_raw <- calc_pca_df(logcpm_raw, "Raw counts")
df_combat_seq <- calc_pca_df(logcpm_combat_seq, "ComBat-seq")
df_combat_ref <- calc_pca_df(logcpm_combat_ref, "ComBat-ref")

season_colors <- c(
  "Winter" = "#1f78b4",
  "Spring" = "#33a02c",
  "Summer" = "#e31a1c",
  "Autumn" = "#ff7f00"
)
condition_shapes <- c("Body" = 16, "Cells" = 17, "Aggregates" = 15)

plot_panel <- function(df) {
  ggplot(df, aes(x = PC1, y = PC2, color = season, shape = condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = season_colors) +
    scale_shape_manual(values = condition_shapes) +
    labs(
      title = unique(df$dataset),
      x = sprintf("PC1 (%s)", unique(df$pc1_var)),
      y = sprintf("PC2 (%s)", unique(df$pc2_var))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

p_raw <- plot_panel(df_raw)
p_combat_seq <- plot_panel(df_combat_seq)
p_combat_ref <- plot_panel(df_combat_ref)

combined <- plot_grid(p_raw + theme(legend.position = "none"),
                      p_combat_seq + theme(legend.position = "none"),
                      p_combat_ref + theme(legend.position = "none"),
                      ncol = 3)
legend <- cowplot::get_legend(p_raw)
final_plot <- plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.12))

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

ggsave(opt$output, final_plot, width = 15, height = 5)
