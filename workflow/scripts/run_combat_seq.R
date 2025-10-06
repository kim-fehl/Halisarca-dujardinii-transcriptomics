#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(sva)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input-rds"), type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option(c("-b", "--batch-column"), type = "character", dest = "batch_column", default = "batch", help = "Metadata column to use as batch"),
  make_option(c("-o", "--output-rds"), type = "character", dest = "output_rds", help = "Output RDS with ComBat-seq adjusted counts")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input_rds) || is.null(opt$output_rds)) {
  stop("--input-rds and --output-rds are required", call. = FALSE)
}

de_data <- readRDS(opt$input_rds)
counts <- de_data$counts
metadata <- de_data$metadata

if (is.null(counts) || is.null(metadata)) {
  stop("Input RDS must contain 'counts' matrix and 'metadata' tibble", call. = FALSE)
}

batch_col <- opt$batch_column
if (!batch_col %in% colnames(metadata)) {
  stop(sprintf("Batch column '%s' not found in metadata", batch_col), call. = FALSE)
}

batch <- droplevels(metadata[[batch_col]])
if (length(unique(batch)) < 2) {
  stop("Batch correction requires at least two batches", call. = FALSE)
}

counts_matrix <- as.matrix(counts)
if (any(counts_matrix < 0)) {
  stop("Counts must be non-negative for ComBat-seq", call. = FALSE)
}

adjusted_counts <- ComBat_seq(counts = counts_matrix, batch = batch, shrink = FALSE)

dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)

saveRDS(adjusted_counts, file = opt$output_rds)
