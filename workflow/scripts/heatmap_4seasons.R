#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(circlize)
  library(ComplexHeatmap)
  library(dplyr)
  library(ggplot2)
  library(grImport)
  library(RColorBrewer)
  library(readr)
  library(readxl)
  library(reshape2)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(tools)
  library(writexl)
})

option_list <- list(
  make_option("--cpm-rds", type = "character", dest = "cpm_rds"),
  make_option("--sample-stats", type = "character", dest = "sample_stats"),
  make_option("--geneset", type = "character", dest = "geneset"),
  make_option("--set-name", type = "character", dest = "set_name", default = "default_set"),
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

icon_paths <- c(intact = 'resources/icons/sponge_icons_t.pdf', cells = 'resources/icons/sponge_icons_c.pdf', aggregates = 'resources/icons/sponge_icons_a.pdf')
missing_icons <- names(icon_paths)[!nzchar(icon_paths) | !file.exists(icon_paths)]
if (length(missing_icons)) {
  stop(sprintf("Missing icon files for: %s. Check resources/icons/", paste(missing_icons, collapse = ", ")), call. = FALSE)
}
icon_paths <- setNames(
  vapply(icon_paths, function(path) {
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }, character(1), USE.NAMES = FALSE),
  names(icon_paths)
)
gs_path <- Sys.which("gs")
if (!nzchar(gs_path)) {
  message("Ghostscript 'gs' not detected in PATH; icon annotations may fall back to color bars.")
} else {
  current_gs <- Sys.getenv("R_GSCMD", unset = "")
  if (!nzchar(current_gs) || current_gs != gs_path) {
    Sys.setenv(R_GSCMD = gs_path)
  }
}

options(repr.plot.res = 200, warn = -1)
set_name <- opt$set_name

season_colors <- c('Winter' = 'dodgerblue', 'Spring' = 'limegreen', 'Summer' = 'magenta', 'Autumn' = 'orange')
cond_labels   <- c('Body' = "intact body", 'Cells' = "dissociated cells", 'Aggregates' = "aggregates")
cond_shapes   <- c('Body' = 2, 'Cells' = 0, 'Aggregates' = 1)
season_labels <- c('Winter' = 'Winter', 'Spring' = 'Spring', 'Summer' = 'Summer', 'Autumn' = 'Autumn')
season_levels <- names(season_labels)
hm_tiss_color_pal <- colorRamp2(c(0, 3, 4), c("#000000", "#17ffc6", "#baffee"))

name_column <- "name"
descr_column <- "descr"
gene_id_column <- "gene_id"
geneset_group_source <- "group" # source column name in gene set XLSX
gene_group_column <- "heatmap_group"

df_cpm_table <- read_rds(opt$cpm_rds)
df_samples_stats_all <- read_tsv(opt$sample_stats, show_col_types = FALSE)

required_sample_cols <- c("sample", "sample_group")
missing_sample_cols <- setdiff(required_sample_cols, colnames(df_samples_stats_all))
if (length(missing_sample_cols)) {
  stop(sprintf("Sample stats file must contain columns: %s", paste(required_sample_cols, collapse = ", ")), call. = FALSE)
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

sample_long <- d %>%
  select(all_of(c(name_column, gene_group_column, sample_columns))) %>%
  pivot_longer(cols = all_of(sample_columns), names_to = "sample", values_to = "CPM") %>%
  inner_join(df_samples_stats_all, by = "sample") %>%
  mutate(sample_group = if_else(!is.na(sample_group), sample_group, paste(season, condition, sep = ".")))

mean_cpm <- sample_long %>%
  group_by_at(c(name_column, gene_group_column, "sample_group")) %>%
  summarize(mean_CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = all_of(c(name_column, gene_group_column)), names_from = sample_group, names_prefix = "mean_CPM.", values_from = mean_CPM)

if (!nrow(mean_cpm)) {
  stop("Unable to derive mean CPM values", call. = FALSE)
}

d <- inner_join(d, mean_cpm, by = c(name_column, gene_group_column))

body_mean_cols <- paste0("mean_CPM.", season_levels, ".Body")
body_mean_cols_present <- intersect(body_mean_cols, colnames(d))
if (!length(body_mean_cols_present)) {
  stop("Mean CPM columns for intact body samples are missing", call. = FALSE)
}

body_matrix <- d %>%
  select(all_of(c(name_column, body_mean_cols_present))) %>%
  column_to_rownames(name_column)
mx_log10means_tiss_CPM <- log10(0.01 + as.matrix(body_matrix))
mx_log10means_tiss_CPM <- mx_log10means_tiss_CPM[, body_mean_cols_present, drop = FALSE]

lfc_cols <- grep("^logFC_", colnames(d), value = TRUE)
adj_cols <- gsub("^logFC_", "adjP_", lfc_cols)
if (!length(lfc_cols)) {
  stop("No logFC columns detected in CPM table", call. = FALSE)
}
mx_lfc <- d %>% select(all_of(c(name_column, lfc_cols))) %>% column_to_rownames(name_column) %>% as.matrix()
mx_fdr <- d %>% select(all_of(c(name_column, adj_cols))) %>% column_to_rownames(name_column) %>% as.matrix()

lfc_meta <- tibble(column = colnames(mx_lfc)) %>%
  mutate(
    contrast = str_remove(column, "^logFC_"),
    condition = str_match(contrast, "_(Cells|Aggregates)_vs_Body$")[, 2],
    season = str_replace(contrast, "_(Cells|Aggregates)_vs_Body$", "")
  )

if (any(is.na(lfc_meta$condition)) || any(is.na(lfc_meta$season))) {
  stop("Unable to parse season/condition from logFC column names", call. = FALSE)
}

season_order <- season_levels[season_levels %in% lfc_meta$season]
if (!length(season_order)) {
  stop("No recognized seasons found in logFC columns", call. = FALSE)
}
condition_order <- c("Cells", "Aggregates")
lfc_meta <- lfc_meta %>%
  mutate(
    season_factor = factor(season, levels = season_order),
    condition_factor = factor(condition, levels = condition_order)
  ) %>%
  arrange(season_factor, condition_factor)

mx_lfc <- mx_lfc[, lfc_meta$column, drop = FALSE]
mx_fdr <- mx_fdr[, gsub("^logFC_", "adjP_", lfc_meta$column), drop = FALSE]
col_splits_lfc <- factor(as.integer(lfc_meta$season_factor), levels = seq_along(season_order))

hm_cell_number_width_mm  <- 8
hm_cell_width_mm         <- 5
hm_cell_height_mm        <- 5
extra_width_mm           <- 100
extra_height_mm          <- 70

total_width_inch  <- min(220,(ncol(mx_log10means_tiss_CPM)*hm_cell_width_mm + ncol(mx_lfc)*hm_cell_number_width_mm + extra_width_mm)) / 25.4
total_height_inch <- min(1500,(nrow(mx_lfc) * hm_cell_height_mm + extra_height_mm)) / 25.4

ann_gene_names <- rowAnnotation(gene = anno_text(rownames(mx_log10means_tiss_CPM), gp = gpar(fontsize = 8)))

body_color_keys <- str_replace(colnames(mx_log10means_tiss_CPM), "^mean_CPM.([^.]+).*", "\\1")
body_colors <- season_colors[body_color_keys]
names(body_colors) <- colnames(mx_log10means_tiss_CPM)

create_season_bar_annotation <- function(keys) {
  HeatmapAnnotation(
    season = keys,
    col = list(season = season_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Season", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  )
}

ann_img <- tryCatch({
  image_pdf <- rep(icon_paths[["intact"]], length.out = ncol(mx_log10means_tiss_CPM))
  HeatmapAnnotation(
    img = anno_image(image_pdf, border = FALSE, gp = gpar(fill = body_colors, col = NA)),
    show_annotation_name = FALSE
  )
}, error = function(e) {
  message("Icon annotation failed (", conditionMessage(e), "); falling back to color bars.")
  create_season_bar_annotation(body_color_keys)
})

group_factor <- factor(d[[gene_group_column]], levels = str_sort(unique(d[[gene_group_column]]), numeric = TRUE))

hm_means <- Heatmap(
  mx_log10means_tiss_CPM,
  name = "meanCPM_intact_body",
  top_annotation = ann_img,
  hm_tiss_color_pal,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_split = group_factor,
  row_title_gp = gpar(fontsize = 8),
  column_title = "Mean expression\nin intact body, CPM",
  column_title_gp = gpar(fontsize = 10, just = "left"),
  row_gap = unit(2, "mm"),
  width  = ncol(mx_log10means_tiss_CPM) * unit(hm_cell_width_mm, "mm"),
  height = nrow(mx_log10means_tiss_CPM) * unit(hm_cell_height_mm, "mm"),
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Mean expression\nin intact body, CPM",
    at = c(0, 1, 2, 3, 4),
    labels = c(0, 10, 100, 1000, 10000),
    direction = "horizontal",
    legend_width = unit(35, "mm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 8)
  )
)

lfc_image_files <- ifelse(lfc_meta$condition == "Cells", icon_paths[["cells"]], icon_paths[["aggregates"]])
season_colors_lfc <- season_colors[lfc_meta$season]
names(season_colors_lfc) <- lfc_meta$column

ann_img_lfc <- tryCatch({
  HeatmapAnnotation(
    img = anno_image(lfc_image_files, border = FALSE, gp = gpar(fill = season_colors_lfc, col = NA)),
    show_annotation_name = FALSE
  )
}, error = function(e) {
  message("LFC icon annotation failed (", conditionMessage(e), "); using color bars instead.")
  HeatmapAnnotation(
    season = factor(lfc_meta$season, levels = season_levels),
    condition = factor(lfc_meta$condition, levels = c("Cells", "Aggregates")),
    col = list(
      season = season_colors,
      condition = c("Cells" = "#5B8FF9", "Aggregates" = "#FFAF00")
    ),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      season = list(title = "Season", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
      condition = list(title = "Condition", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
    )
  )
})

hm_lfc <- Heatmap(
  mx_lfc,
  name = "LFC",
  top_annotation = ann_img_lfc,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  row_split = group_factor,
  row_title_gp = gpar(fontsize = 8),
  row_gap = unit(2, "mm"),
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  column_split = col_splits_lfc,
  column_title = "Log2(fold change) of expression\nin cell suspensions\nor aggregates vs. intact body",
  column_title_gp = gpar(fontsize = 10, just = "left"),
  col = colorRamp2(c(-3, 0, 2), c("#218cf6", "gray95", "#fe1112")),
  width  = ncol(mx_lfc) * unit(hm_cell_number_width_mm, "mm"),
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
