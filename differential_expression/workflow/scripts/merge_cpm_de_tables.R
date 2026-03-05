#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
})

option_list <- list(
  make_option("--cpm-tsv", type = "character", dest = "cpm_tsv", help = "All-genes CPM TSV.gz"),
  make_option("--edgeR-tsv", type = "character", dest = "edger_tsv", help = "edgeR long-form results TSV.gz"),
  make_option("--output-tsv", type = "character", dest = "output_tsv", help = "Output merged TSV.gz")
)

opt <- parse_args(OptionParser(option_list = option_list))
required <- c("cpm_tsv", "edger_tsv", "output_tsv")
missing <- required[!nzchar(trimws(unlist(opt[required])))]
if (length(missing)) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

cpm_tbl <- read_tsv(opt$cpm_tsv, progress = FALSE, show_col_types = FALSE)
if (!("gene_id" %in% colnames(cpm_tbl))) {
  stop("CPM table must contain 'gene_id' column", call. = FALSE)
}

edger_tbl <- read_tsv(opt$edger_tsv, progress = FALSE, show_col_types = FALSE)
required_edgeR_cols <- c("gene_id", "contrast", "logFC", "padj_global")
missing_edgeR_cols <- setdiff(required_edgeR_cols, colnames(edger_tbl))
if (length(missing_edgeR_cols) > 0) {
  stop(sprintf("edgeR table missing required columns: %s", paste(missing_edgeR_cols, collapse = ", ")), call. = FALSE)
}

logfc_wide <- edger_tbl %>%
  select(gene_id, contrast, logFC) %>%
  distinct() %>%
  pivot_wider(names_from = contrast, values_from = logFC, names_prefix = "logFC_")

padj_wide <- edger_tbl %>%
  select(gene_id, contrast, padj_global) %>%
  distinct() %>%
  pivot_wider(names_from = contrast, values_from = padj_global, names_prefix = "adjP_")

merged <- cpm_tbl %>%
  left_join(logfc_wide, by = "gene_id") %>%
  left_join(padj_wide, by = "gene_id")

dir.create(dirname(opt$output_tsv), recursive = TRUE, showWarnings = FALSE)
write_tsv(merged, opt$output_tsv)
