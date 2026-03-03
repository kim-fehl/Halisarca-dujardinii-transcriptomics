#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input-rds"), type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option(c("-o", "--output"), type = "character", dest = "output", help = "Output PDF path")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds) || is.null(opt$output)) {
  stop("--input-rds and --output are required", call. = FALSE)
}

de_data <- readRDS(opt$input_rds)
counts <- de_data$counts
metadata <- de_data$metadata

if (is.null(counts) || is.null(metadata)) {
  stop("Input RDS must contain 'counts' and 'metadata'", call. = FALSE)
}

if (!all(colnames(counts) %in% metadata$sample_id)) {
  stop("Counts columns and metadata sample_id values are inconsistent", call. = FALSE)
}

metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    season = as.character(season),
    condition = as.character(condition)
  ) %>%
  mutate(order_index = match(sample_id, colnames(counts))) %>%
  filter(!is.na(order_index)) %>%
  arrange(order_index) %>%
  select(-order_index)

counts <- counts[, metadata$sample_id, drop = FALSE]

logcpm <- cpm(counts, log = TRUE, prior.count = 1)
pca <- prcomp(t(logcpm), scale. = FALSE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- as_tibble(pca$x[, 1:2], rownames = "sample_id") %>%
  left_join(metadata, by = "sample_id") %>%
  mutate(
    pc1_var = scales::percent(var_expl[1], accuracy = 0.1),
    pc2_var = scales::percent(var_expl[2], accuracy = 0.1)
  )

season_levels <- unique(pca_df$season)
season_defaults <- c("Winter" = "dodgerblue", "Spring" = "limegreen", "Summer" = "magenta", "Autumn" = "orange")
season_colors <- setNames(rep("grey60", length(season_levels)), season_levels)
known_seasons <- intersect(names(season_defaults), season_levels)
if (length(known_seasons)) {
  season_colors[known_seasons] <- season_defaults[known_seasons]
}
unknown_seasons <- setdiff(season_levels, names(season_defaults))
if (length(unknown_seasons)) {
  season_colors[unknown_seasons] <- hue_pal()(length(unknown_seasons))
}

condition_levels <- unique(pca_df$condition)
shape_defaults <- c("Body" = 2, "Cells" = 0, "Aggregates" = 1)
shape_values <- setNames(rep(16, length(condition_levels)), condition_levels)
known_conditions <- intersect(names(shape_defaults), condition_levels)
if (length(known_conditions)) {
  shape_values[known_conditions] <- shape_defaults[known_conditions]
}
unknown_conditions <- setdiff(condition_levels, names(shape_defaults))
if (length(unknown_conditions)) {
  fallback_shapes <- c(3, 4, 5, 6, 7, 8, 15, 17, 18)
  shape_values[unknown_conditions] <- fallback_shapes[((seq_len(length(unknown_conditions)) - 1) %% length(fallback_shapes)) + 1]
}

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = season, shape = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = season_colors) +
  scale_shape_manual(values = shape_values) +
  labs(
    title = "PCA (raw counts)",
    x = sprintf("PC1 (%s)", unique(pca_df$pc1_var)),
    y = sprintf("PC2 (%s)", unique(pca_df$pc2_var))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
ggsave(opt$output, p, width = 6, height = 5)
