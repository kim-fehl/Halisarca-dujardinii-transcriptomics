#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readxl)
  library(stringr)
  library(scales)
  library(grid)
  library(tibble)
})

option_list <- list(
  make_option("--input-rds", type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option("--geneset", type = "character", dest = "geneset", help = "Path to geneset XLSX"),
  make_option("--baseline-level", type = "character", dest = "baseline_level", default = "Body", help = "Condition used for baseline median CPM annotation"),
  make_option("--output-pdf", type = "character", dest = "output_pdf", help = "Output heatmap PDF")
)

opt <- parse_args(OptionParser(option_list = option_list))
required <- c("input_rds", "geneset", "output_pdf")
missing <- required[!nzchar(trimws(unlist(opt[required])))]
if (length(missing)) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

stopifnot(file.exists(opt$input_rds), file.exists(opt$geneset))

de_data <- readRDS(opt$input_rds)
counts <- de_data$counts
metadata <- de_data$metadata

if (is.null(counts) || is.null(metadata)) {
  stop("Input RDS must contain 'counts' and 'metadata'", call. = FALSE)
}

required_meta <- c("sample_id", "season", "condition")
missing_meta <- setdiff(required_meta, colnames(metadata))
if (length(missing_meta)) {
  stop(sprintf("Metadata missing columns: %s", paste(missing_meta, collapse = ", ")), call. = FALSE)
}

metadata <- metadata %>%
  transmute(
    sample_id = as.character(sample_id),
    stratum = as.character(season),
    condition = as.character(condition),
    replicate = if ("replicate" %in% colnames(de_data$metadata)) as.integer(replicate) else NA_integer_
  ) %>%
  mutate(order_index = match(sample_id, colnames(counts))) %>%
  filter(!is.na(order_index)) %>%
  arrange(order_index) %>%
  select(-order_index)

counts <- counts[, metadata$sample_id, drop = FALSE]

ann_raw <- read_excel(opt$geneset)
if (!("gene_id" %in% colnames(ann_raw)) || !("name" %in% colnames(ann_raw))) {
  stop("Geneset must contain at least 'gene_id' and 'name' columns", call. = FALSE)
}
if (!("group" %in% colnames(ann_raw))) {
  ann_raw$group <- NA_character_
}

ann <- ann_raw %>%
  transmute(
    gene_id = as.character(gene_id),
    name = as.character(name),
    group = as.character(group)
  ) %>%
  filter(!is.na(gene_id) & nzchar(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

go_idx <- ann$gene_id %in% rownames(counts)
if (!any(go_idx)) {
  stop("None of geneset gene_id values are present in counts", call. = FALSE)
}
ann <- ann[go_idx, , drop = FALSE]
ann$name <- make.unique(ann$name, sep = "_dup")

go_ids <- ann$gene_id

logcpm <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
cpm_mat <- edgeR::cpm(counts, log = FALSE)

heat_mat <- logcpm[go_ids, , drop = FALSE]
rownames(heat_mat) <- ann$name

heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z[is.na(heat_mat_z)] <- 0

baseline_samples <- metadata %>%
  filter(condition == opt$baseline_level) %>%
  pull(sample_id)
if (!length(baseline_samples)) {
  stop(sprintf("No samples found for baseline level '%s'", opt$baseline_level), call. = FALSE)
}

median_baseline <- apply(cpm_mat[go_ids, baseline_samples, drop = FALSE], 1, median, na.rm = TRUE)
names(median_baseline) <- ann$name

all_condition_levels <- unique(metadata$condition)
condition_levels <- c(
  opt$baseline_level,
  setdiff(all_condition_levels, opt$baseline_level)
)
stratum_levels <- unique(metadata$stratum)

metadata <- metadata %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    stratum = factor(stratum, levels = stratum_levels)
  ) %>%
  arrange(condition, stratum, replicate, sample_id)

sample_order <- as.character(metadata$sample_id)
heat_mat_z <- heat_mat_z[, sample_order, drop = FALSE]

season_defaults <- c("Winter" = "dodgerblue", "Spring" = "limegreen", "Summer" = "magenta", "Autumn" = "orange")
stratum_colors <- setNames(rep("grey60", length(stratum_levels)), stratum_levels)
known_strata <- intersect(names(season_defaults), stratum_levels)
if (length(known_strata)) {
  stratum_colors[known_strata] <- season_defaults[known_strata]
}
unknown_strata <- setdiff(stratum_levels, names(season_defaults))
if (length(unknown_strata)) {
  stratum_colors[unknown_strata] <- hue_pal()(length(unknown_strata))
}

condition_defaults <- c("Body" = "#7A7A7A", "Cells" = "#5B8FF9", "Aggregates" = "#FFAF00")
condition_colors <- setNames(rep("grey70", length(condition_levels)), condition_levels)
known_conditions <- intersect(names(condition_defaults), condition_levels)
if (length(known_conditions)) {
  condition_colors[known_conditions] <- condition_defaults[known_conditions]
}
unknown_conditions <- setdiff(condition_levels, names(condition_defaults))
if (length(unknown_conditions)) {
  condition_colors[unknown_conditions] <- hue_pal()(length(unknown_conditions))
}

ha_col <- HeatmapAnnotation(
  Condition = metadata$condition,
  Stratum = metadata$stratum,
  show_legend = TRUE,
  col = list(
    Condition = condition_colors,
    Stratum = stratum_colors
  ),
  show_annotation_name = FALSE
)

row_groups <- ann$group
row_groups <- trimws(row_groups)
row_groups[row_groups == ""] <- NA
if (all(is.na(row_groups)) || length(na.omit(unique(row_groups))) <= 1) {
  row_groups <- NULL
} else {
  row_groups[is.na(row_groups)] <- "Ungrouped"
  row_groups <- factor(row_groups, levels = unique(row_groups))
}

row_ha <- rowAnnotation(
  `MedianCPM_baseline` = anno_barplot(
    median_baseline,
    border = FALSE,
    width = unit(2, "cm"),
    gp = gpar(fill = "#208a76", col = NA)
  ),
  annotation_name_rot = 0,
  show_annotation_name = TRUE
)

max_name_len <- max(nchar(rownames(heat_mat_z)), na.rm = TRUE)
plot_width <- 6 + max(0.0, 0.08 * (max_name_len - 10))
plot_height <- 1.56 + 0.174 * nrow(heat_mat_z)

ht <- Heatmap(
  heat_mat_z,
  name = "Z-score",
  top_annotation = ha_col,
  left_annotation = row_ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  cluster_columns = TRUE,
  column_split = metadata$condition,
  cluster_column_slices = FALSE,
  row_names_max_width = unit(0.09 * max_name_len, "in"),
  col = colorRamp2(c(-3, 0, 3), c("#218cf6", "gray95", "#fe1112")),
  column_gap = unit(2, "mm"),
  column_title = "",
  heatmap_legend_param = list(
    title = "Expr. (Z-score logCPM)",
    direction = "horizontal"
  ),
  row_split = row_groups,
  row_title = NA
)

dir.create(dirname(opt$output_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(opt$output_pdf, width = plot_width, height = plot_height)
draw(
  ht,
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()
