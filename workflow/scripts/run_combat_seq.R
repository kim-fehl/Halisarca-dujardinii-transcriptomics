#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(sva)
  library(readr)
})

option_list <- list(
  make_option("--input-rds", type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option("--batch-column", type = "character", dest = "batch_column", default = "season", help = "Metadata column to use as batch"),
  make_option("--output-rds", type = "character", dest = "output_rds", help = "Output RDS with ComBat-seq adjusted counts")
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

if (!opt$batch_column %in% colnames(metadata)) {
  stop(sprintf("Batch column '%s' not found in metadata", opt$batch_column), call. = FALSE)
}

batch <- metadata[[opt$batch_column]]
if (any(is.na(batch))) {
  stop("Batch column contains NA values", call. = FALSE)
}

if (length(unique(batch)) < 2) {
  stop("Batch correction requires at least two batches", call. = FALSE)
}

counts_matrix <- as.matrix(counts)
if (any(counts_matrix < 0)) {
  stop("Counts must be non-negative for ComBat-seq", call. = FALSE)
}

adjusted_counts <- ComBat_seq(counts = counts_matrix, batch = batch)

for (path in c(opt$output_rds)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

saveRDS(adjusted_counts, file = opt$output_rds)
