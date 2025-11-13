#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

option_list <- list(
  make_option("--input-rds", type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option("--edgeR-tsv", type = "character", dest = "edger_tsv", help = "edgeR results_long.tsv.gz"),
  make_option("--set-name", type = "character", dest = "set_name", help = "Label used for output summarizing"),
  make_option("--output-rds", type = "character", dest = "output_rds", help = "Output RDS file with CPM/logFC/padj"),
  make_option("--samples-stats", type = "character", dest = "samples_stats", help = "Output TSV with per-sample stats")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("input_rds", "edger_tsv", "set_name", "output_rds", "samples_stats")
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
    season = factor(as.character(season), levels = c("Autumn", "Winter", "Spring", "Summer")),
    condition = factor(as.character(condition), levels = c("Body", "Cells", "Aggregates"))
  ) %>%
  mutate(order_index = match(sample_id, colnames(counts))) %>%
  filter(!is.na(order_index)) %>%
  arrange(order_index) %>%
  select(-order_index) %>%
  mutate(sample_group = interaction(season, condition, sep = ".", drop = TRUE))

excluded_autumn3 <- metadata %>%
  filter(season == "Autumn", replicate == 3) %>%
  pull(sample_id)

if (length(excluded_autumn3)) {
  message(sprintf("Excluding %d Autumn replicate 3 samples from heatmap inputs: %s",
                  length(excluded_autumn3), paste(excluded_autumn3, collapse = ", ")))
  metadata <- metadata %>%
    filter(!(season == "Autumn" & replicate == 3))
}

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
    season = as.character(season),
    condition = as.character(condition),
    sample_group = as.character(sample_group),
    replicate = replicate,
    hpd = hpd
  ) %>%
  left_join(rownames_to_column(as.data.frame(dgl$samples), var = "sample"), by = "sample") %>%
  arrange(sample)

write_tsv(sample_stats, opt$samples_stats)

edgeR_results <- read_tsv(opt$edger_tsv, progress = FALSE)

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

write_rds(combined_table, opt$output_rds, compress = "bz2")
