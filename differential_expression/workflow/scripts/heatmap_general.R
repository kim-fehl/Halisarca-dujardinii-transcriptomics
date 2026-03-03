#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(circlize)
  library(ComplexHeatmap)
  library(dplyr)
  library(RColorBrewer)
  library(readr)
  library(readxl)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(writexl)
})

option_list <- list(
  make_option("--cpm-rds", type = "character", dest = "cpm_rds"),
  make_option("--sample-stats", type = "character", dest = "sample_stats"),
  make_option("--geneset", type = "character", dest = "geneset"),
  make_option("--baseline-level", type = "character", dest = "baseline_level", default = "Body"),
  make_option("--genome-name", type = "character", dest = "genome_name", default = "default_genome"),
  make_option("--output-pdf", type = "character", dest = "output_pdf"),
  make_option("--output-xlsx", type = "character", dest = "output_xlsx")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_opts <- c("cpm_rds", "sample_stats", "geneset", "output_pdf", "output_xlsx")
missing_required <- required_opts[!nzchar(trimws(unlist(opt[required_opts])))]
if (length(missing_required)) {
  stop(sprintf("Missing required arguments: %s", paste(missing_required, collapse = ", ")), call. = FALSE)
}

stopifnot(file.exists(opt$cpm_rds))
stopifnot(file.exists(opt$sample_stats))
stopifnot(file.exists(opt$geneset))

options(repr.plot.res = 200, warn = -1)

name_column <- "name"
descr_column <- "descr"
gene_id_column <- "gene_id"
geneset_group_source <- "group"
gene_group_column <- "heatmap_group"

heatmap_data <- read_rds(opt$cpm_rds)
if (is.data.frame(heatmap_data)) {
  stop("General heatmap input RDS must be produced by prepare_heatmap_inputs_general.R", call. = FALSE)
}

required_rds_fields <- c("table", "contrast_meta", "baseline_level")
missing_rds_fields <- setdiff(required_rds_fields, names(heatmap_data))
if (length(missing_rds_fields) > 0) {
  stop(sprintf("General heatmap RDS missing fields: %s", paste(missing_rds_fields, collapse = ", ")), call. = FALSE)
}

df_cpm_table <- heatmap_data$table
contrast_meta <- heatmap_data$contrast_meta %>%
  mutate(
    contrast = as.character(contrast),
    stratum = as.character(stratum),
    condition = as.character(condition),
    baseline = as.character(baseline)
  )
baseline_level <- as.character(heatmap_data$baseline_level)
if (nzchar(opt$baseline_level)) {
  baseline_level <- opt$baseline_level
}

df_samples_stats_all <- read_tsv(opt$sample_stats, show_col_types = FALSE)
required_sample_cols <- c("sample", "stratum", "condition")
missing_sample_cols <- setdiff(required_sample_cols, colnames(df_samples_stats_all))
if (length(missing_sample_cols)) {
  stop(sprintf("Sample stats file must contain columns: %s", paste(required_sample_cols, collapse = ", ")), call. = FALSE)
}

if (!"sample_group" %in% colnames(df_samples_stats_all)) {
  df_samples_stats_all <- df_samples_stats_all %>%
    mutate(sample_group = paste(stratum, condition, sep = "."))
}

required_gene_cols <- c(name_column, descr_column, gene_id_column)
df_geneset_raw <- read_excel(opt$geneset)
if (!all(required_gene_cols %in% colnames(df_geneset_raw))) {
  stop(sprintf("Gene set file must contain columns: %s", paste(required_gene_cols, collapse = ", ")), call. = FALSE)
}
if (geneset_group_source %in% colnames(df_geneset_raw)) {
  df_geneset_raw <- df_geneset_raw %>%
    mutate(!!geneset_group_source := as.character(.data[[geneset_group_source]])) %>%
    rename(!!gene_group_column := all_of(geneset_group_source))
} else {
  df_geneset_raw <- df_geneset_raw %>%
    mutate(!!gene_group_column := NA_character_)
}

df_geneset <- df_geneset_raw %>%
  select(all_of(c(required_gene_cols, gene_group_column))) %>%
  mutate(!!name_column := make.unique(as.character(.data[[name_column]]), sep = "_dup")) %>%
  mutate(!!gene_group_column := if_else(is.na(.data[[gene_group_column]]) | !nzchar(.data[[gene_group_column]]),
                                        " ",
                                        as.character(.data[[gene_group_column]])))

if (!"gene_id" %in% colnames(df_cpm_table)) {
  stop("CPM table must contain a 'gene_id' column", call. = FALSE)
}

sample_columns <- intersect(df_samples_stats_all$sample, colnames(df_cpm_table))
if (!length(sample_columns)) {
  stop("No overlapping samples between CPM table and sample stats", call. = FALSE)
}

d <- inner_join(df_geneset, df_cpm_table, by = setNames("gene_id", gene_id_column))
if (!nrow(d)) {
  stop("No overlap between gene set and CPM table", call. = FALSE)
}

sample_long <- d %>%
  select(all_of(c(name_column, gene_group_column, sample_columns))) %>%
  pivot_longer(cols = all_of(sample_columns), names_to = "sample", values_to = "CPM") %>%
  inner_join(df_samples_stats_all, by = "sample")

mean_cpm <- sample_long %>%
  filter(condition == baseline_level) %>%
  group_by_at(c(name_column, gene_group_column, "stratum")) %>%
  summarize(mean_CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = all_of(c(name_column, gene_group_column)), names_from = stratum, names_prefix = "mean_CPM.", values_from = mean_CPM)

if (!nrow(mean_cpm)) {
  stop(sprintf("Unable to derive mean CPM values for baseline '%s'", baseline_level), call. = FALSE)
}

d <- inner_join(d, mean_cpm, by = c(name_column, gene_group_column))

stratum_order <- unique(contrast_meta$stratum)
if (!length(stratum_order)) {
  stop("No strata found in contrast metadata", call. = FALSE)
}

baseline_mean_cols <- paste0("mean_CPM.", stratum_order)
baseline_mean_cols_present <- intersect(baseline_mean_cols, colnames(d))
if (!length(baseline_mean_cols_present)) {
  stop("Mean CPM columns for baseline samples are missing", call. = FALSE)
}

baseline_matrix <- d %>%
  select(all_of(c(name_column, baseline_mean_cols_present))) %>%
  column_to_rownames(name_column)
mx_log10means_base_CPM <- log10(0.01 + as.matrix(baseline_matrix))
mx_log10means_base_CPM <- mx_log10means_base_CPM[, baseline_mean_cols_present, drop = FALSE]

contrast_meta <- contrast_meta %>%
  mutate(
    logFC_col = paste0("logFC_", contrast),
    adjP_col = paste0("adjP_", contrast)
  ) %>%
  filter(logFC_col %in% colnames(d), adjP_col %in% colnames(d))

if (!nrow(contrast_meta)) {
  stop("No logFC/adjP columns detected in heatmap table", call. = FALSE)
}

compare_order <- unique(contrast_meta$condition)
contrast_meta <- contrast_meta %>%
  mutate(
    stratum_factor = factor(stratum, levels = stratum_order),
    condition_factor = factor(condition, levels = compare_order)
  ) %>%
  arrange(stratum_factor, condition_factor)

mx_lfc <- d %>%
  select(all_of(c(name_column, contrast_meta$logFC_col))) %>%
  column_to_rownames(name_column) %>%
  as.matrix()
mx_fdr <- d %>%
  select(all_of(c(name_column, contrast_meta$adjP_col))) %>%
  column_to_rownames(name_column) %>%
  as.matrix()

mx_lfc <- mx_lfc[, contrast_meta$logFC_col, drop = FALSE]
mx_fdr <- mx_fdr[, contrast_meta$adjP_col, drop = FALSE]

colnames(mx_lfc) <- paste0(contrast_meta$condition, " vs ", contrast_meta$baseline)
colnames(mx_fdr) <- colnames(mx_lfc)
col_splits_lfc <- factor(contrast_meta$stratum, levels = stratum_order)

all_strata <- unique(c(stratum_order, str_remove(colnames(mx_log10means_base_CPM), "^mean_CPM\\.")))
season_colors_default <- c("Winter" = "dodgerblue", "Spring" = "limegreen", "Summer" = "magenta", "Autumn" = "orange")
stratum_colors <- setNames(rep("grey60", length(all_strata)), all_strata)
known_strata <- intersect(names(season_colors_default), all_strata)
if (length(known_strata)) {
  stratum_colors[known_strata] <- season_colors_default[known_strata]
}
unknown_strata <- setdiff(all_strata, names(season_colors_default))
if (length(unknown_strata)) {
  stratum_colors[unknown_strata] <- hue_pal()(length(unknown_strata))
}

condition_defaults <- c("Body" = "#7A7A7A", "Cells" = "#5B8FF9", "Aggregates" = "#FFAF00")
condition_levels <- unique(c(as.character(contrast_meta$condition), baseline_level))
condition_colors <- setNames(rep("grey70", length(condition_levels)), condition_levels)
known_conditions <- intersect(names(condition_defaults), condition_levels)
if (length(known_conditions)) {
  condition_colors[known_conditions] <- condition_defaults[known_conditions]
}
unknown_conditions <- setdiff(condition_levels, names(condition_defaults))
if (length(unknown_conditions)) {
  condition_colors[unknown_conditions] <- hue_pal()(length(unknown_conditions))
}

hm_tiss_color_pal <- colorRamp2(c(0, 3, 4), c("#000000", "#17ffc6", "#baffee"))

hm_cell_number_width_mm <- 8
hm_cell_width_mm <- 5
hm_cell_height_mm <- 5
extra_width_mm <- 100
extra_height_mm <- 70

total_width_inch <- min(
  220,
  (ncol(mx_log10means_base_CPM) * hm_cell_width_mm + ncol(mx_lfc) * hm_cell_number_width_mm + extra_width_mm)
) / 25.4
total_height_inch <- min(1500, (nrow(mx_lfc) * hm_cell_height_mm + extra_height_mm)) / 25.4

ann_gene_names <- rowAnnotation(gene = anno_text(rownames(mx_log10means_base_CPM), gp = gpar(fontsize = 8)))

group_factor <- factor(d[[gene_group_column]], levels = str_sort(unique(d[[gene_group_column]]), numeric = TRUE))

mean_strata <- str_remove(colnames(mx_log10means_base_CPM), "^mean_CPM\\.")
ann_means <- HeatmapAnnotation(
  stratum = factor(mean_strata, levels = all_strata),
  col = list(stratum = stratum_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    stratum = list(title = "Stratum", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  )
)

hm_means <- Heatmap(
  mx_log10means_base_CPM,
  name = paste0("meanCPM_", baseline_level),
  top_annotation = ann_means,
  hm_tiss_color_pal,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_split = group_factor,
  row_title_gp = gpar(fontsize = 8),
  column_title = sprintf("Mean expression\nin baseline (%s), CPM", baseline_level),
  column_title_gp = gpar(fontsize = 10, just = "left"),
  row_gap = unit(2, "mm"),
  width = ncol(mx_log10means_base_CPM) * unit(hm_cell_width_mm, "mm"),
  height = nrow(mx_log10means_base_CPM) * unit(hm_cell_height_mm, "mm"),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  column_labels = mean_strata,
  heatmap_legend_param = list(
    title = sprintf("Mean expression\nin baseline (%s), CPM", baseline_level),
    at = c(0, 1, 2, 3, 4),
    labels = c(0, 10, 100, 1000, 10000),
    direction = "horizontal",
    legend_width = unit(35, "mm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 8)
  )
)

ann_lfc <- HeatmapAnnotation(
  stratum = factor(contrast_meta$stratum, levels = stratum_order),
  condition = factor(contrast_meta$condition, levels = compare_order),
  col = list(
    stratum = stratum_colors[stratum_order],
    condition = condition_colors[compare_order]
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    stratum = list(title = "Stratum", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
    condition = list(title = "Condition", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  )
)

hm_lfc <- Heatmap(
  mx_lfc,
  name = "LFC",
  top_annotation = ann_lfc,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  row_split = group_factor,
  row_title_gp = gpar(fontsize = 8),
  row_gap = unit(2, "mm"),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  column_gap = unit(2, "mm"),
  column_split = col_splits_lfc,
  column_title = sprintf("Log2(fold change) of expression\nvs baseline (%s)", baseline_level),
  column_title_gp = gpar(fontsize = 10, just = "left"),
  col = colorRamp2(c(-3, 0, 2), c("#218cf6", "gray95", "#fe1112")),
  width = ncol(mx_lfc) * unit(hm_cell_number_width_mm, "mm"),
  height = nrow(mx_lfc) * unit(hm_cell_height_mm, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(mx_lfc[i, j]) && !is.na(mx_fdr[i, j])) {
      if (abs(mx_lfc[i, j]) >= 0.7 && mx_fdr[i, j] < 0.001) {
        grid.text(sprintf("%.1f", mx_lfc[i, j]), x, y, gp = gpar(fontsize = 7))
      }
    }
  },
  heatmap_legend_param = list(
    title = "\nlog2(fold change)",
    at = c(-3, -1, 0, 1, 3),
    labels = c(-3, -1, 0, 1, 3),
    direction = "horizontal",
    legend_width = unit(35, "mm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 8)
  )
)

for (path in c(opt$output_pdf, opt$output_xlsx)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

pdf(opt$output_pdf, width = total_width_inch, height = total_height_inch)
draw(hm_means + ann_gene_names + hm_lfc,
     main_heatmap = "LFC",
     heatmap_legend_side = "bottom",
     row_sub_title_side = "left",
     annotation_legend_side = "right")
dev.off()

draw(hm_means + ann_gene_names + hm_lfc,
     main_heatmap = "LFC",
     heatmap_legend_side = "bottom",
     row_sub_title_side = "left",
     annotation_legend_side = "right")

write_xlsx(d, opt$output_xlsx)
