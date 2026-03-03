#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

option_list <- list(
  make_option("--input-rds", type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option("--edgeR-tsv", type = "character", dest = "edger_tsv", help = "General edgeR results_long TSV.gz"),
  make_option("--baseline-level", type = "character", dest = "baseline_level", default = "Body",
              help = "Reference condition level [default %default]"),
  make_option("--output-rds", type = "character", dest = "output_rds", help = "Output RDS with CPM/logFC/padj and contrast metadata"),
  make_option("--samples-stats", type = "character", dest = "samples_stats", help = "Output TSV with per-sample stats")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("input_rds", "edger_tsv", "output_rds", "samples_stats")
missing_args <- required[!nzchar(trimws(unlist(opt[required])))]
if (length(missing_args) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
}

for (path in c(opt$output_rds, opt$samples_stats)) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

de_data <- readRDS(opt$input_rds)
counts <- de_data$counts
metadata <- de_data$metadata

if (is.null(counts) || is.null(metadata)) {
  stop("Input RDS must contain 'counts' and 'metadata' entries", call. = FALSE)
}

metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    stratum = as.character(season),
    condition = as.character(condition)
  ) %>%
  mutate(order_index = match(sample_id, colnames(counts))) %>%
  filter(!is.na(order_index)) %>%
  arrange(order_index) %>%
  select(-order_index) %>%
  mutate(sample_group = interaction(stratum, condition, sep = ".", drop = TRUE))

counts <- counts[, metadata$sample_id, drop = FALSE]
metadata <- metadata %>%
  mutate(sample_id = factor(sample_id, levels = colnames(counts)))

sample_group <- metadata$sample_group

dgl <- DGEList(counts = counts)
keep <- filterByExpr(dgl, group = sample_group)
dgl <- dgl[keep, , keep.lib.sizes = FALSE]
dgl <- calcNormFactors(dgl)

cpm_table <- as.data.frame(edgeR::cpm(dgl)) %>%
  rownames_to_column("gene_id")

sample_stats <- metadata %>%
  transmute(
    sample = as.character(sample_id),
    stratum = stratum,
    condition = condition,
    sample_group = as.character(sample_group),
    replicate = replicate,
    hpd = hpd
  ) %>%
  left_join(rownames_to_column(as.data.frame(dgl$samples), var = "sample"), by = "sample") %>%
  arrange(sample)

write_tsv(sample_stats, opt$samples_stats)

edgeR_results <- read_tsv(opt$edger_tsv, progress = FALSE, show_col_types = FALSE)
required_edgeR_cols <- c("gene_id", "contrast", "logFC", "padj_global", "stratum", "condition", "baseline")
missing_edgeR_cols <- setdiff(required_edgeR_cols, colnames(edgeR_results))
if (length(missing_edgeR_cols) > 0) {
  stop(
    sprintf("edgeR table is missing columns: %s", paste(missing_edgeR_cols, collapse = ", ")),
    call. = FALSE
  )
}

contrast_meta <- edgeR_results %>%
  select(contrast, stratum, condition, baseline) %>%
  distinct() %>%
  arrange(stratum, condition)

if (nrow(contrast_meta) == 0) {
  stop("No contrasts found in edgeR results", call. = FALSE)
}

baselines <- unique(contrast_meta$baseline)
if (length(baselines) != 1) {
  stop(
    sprintf("Expected a single baseline level in edgeR results, found: %s", paste(baselines, collapse = ", ")),
    call. = FALSE
  )
}
if (nzchar(opt$baseline_level) && baselines[[1]] != opt$baseline_level) {
  warning(
    sprintf(
      "Requested baseline '%s' differs from edgeR baseline '%s'; using edgeR baseline",
      opt$baseline_level,
      baselines[[1]]
    ),
    call. = FALSE
  )
}

logfc_wide <- edgeR_results %>%
  select(gene_id, contrast, logFC) %>%
  distinct() %>%
  pivot_wider(names_from = contrast, values_from = logFC, names_prefix = "logFC_")

padj_wide <- edgeR_results %>%
  select(gene_id, contrast, padj_global) %>%
  distinct() %>%
  pivot_wider(names_from = contrast, values_from = padj_global, names_prefix = "adjP_")

combined_table <- cpm_table %>%
  left_join(logfc_wide, by = "gene_id") %>%
  left_join(padj_wide, by = "gene_id")

out <- list(
  table = combined_table,
  contrast_meta = contrast_meta,
  baseline_level = baselines[[1]]
)

write_rds(out, opt$output_rds, compress = "bz2")
