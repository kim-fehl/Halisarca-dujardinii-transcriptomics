#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input-rds"), type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option(c("-o", "--output-rds"), type = "character", dest = "output_rds", help = "Output RDS with ComBat-ref adjusted counts"),
  make_option(c("-d", "--repo-dir"), type = "character", dest = "repo_dir", help = "Directory containing Combat-ref scripts")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input_rds) || is.null(opt$output_rds)) {
  stop("--input-rds and --output-rds are required", call. = FALSE)
}

required_pkgs <- c(
  "sva", "dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra",
  "scales", "RUVSeq", "ggpubr", "BatchQC", "edgeR"
)
pkg_status <- sapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
})
if (any(!pkg_status)) {
  missing <- paste(required_pkgs[!pkg_status], collapse = ", ")
  stop(sprintf("Missing required R packages for ComBat-ref: %s", missing), call. = FALSE)
}

de_data <- readRDS(opt$input_rds)
counts <- de_data$counts
metadata <- de_data$metadata

if (is.null(counts) || is.null(metadata)) {
  stop("Input RDS must contain 'counts' matrix and 'metadata' tibble", call. = FALSE)
}

if (!"batch" %in% colnames(metadata)) {
  stop("Metadata must contain 'batch' column for ComBat-ref", call. = FALSE)
}

batch <- droplevels(metadata$batch)
if (length(unique(batch)) < 2) {
  stop("ComBat-ref requires at least two batches", call. = FALSE)
}

group <- droplevels(interaction(metadata$season, metadata$condition, sep = "."))
if (length(unique(group)) < 2) {
  group <- metadata$condition
}

repo_dir <- opt$repo_dir
if (is.null(repo_dir) || repo_dir == "") {
  repo_dir <- file.path("workflow", "scripts", "Combat-ref")
}

combat_fn <- NULL
if (dir.exists(repo_dir)) {
  script_dir <- normalizePath(repo_dir, winslash = "/", mustWork = FALSE)
  safe_source <- function(path) {
    if (!file.exists(path)) {
      warning(sprintf("Combat-ref helper missing: %s", path))
      return(FALSE)
    }
    tryCatch({
      source(path)
      TRUE
    }, error = function(e) {
      warning(sprintf("Failed to source %s: %s", path, conditionMessage(e)))
      FALSE
    })
  }
  helper_path <- file.path(repo_dir, "helper_seq.R")
  combat_ref_path <- file.path(repo_dir, "Combat_ref.R")
  combat_seq_path <- file.path(repo_dir, "ComBat_seq.R")
  safe_source(helper_path)
  safe_source(combat_seq_path)
  safe_source(combat_ref_path)
  if (exists("ComBat_ref", mode = "function")) {
    combat_fn <- get("ComBat_ref")
  }
}

if (is.null(combat_fn)) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("ComBat_ref function not found and package 'sva' is unavailable", call. = FALSE)
  }
  combat_fn <- sva::ComBat_ref
}

counts_matrix <- as.matrix(counts)
if (any(counts_matrix < 0)) {
  stop("Counts must be non-negative for ComBat-ref", call. = FALSE)
}

adjusted_counts <- combat_fn(counts = counts_matrix, batch = batch, group = group, genewise.disp = FALSE)

dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)

saveRDS(adjusted_counts, file = opt$output_rds)
