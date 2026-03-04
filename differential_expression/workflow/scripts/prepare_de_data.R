#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

option_list <- list(
  make_option("--counts", type = "character", dest = "counts", help = "featureCounts TSV (gzipped)"),
  make_option("--metadata", type = "character", dest = "metadata", help = "Sample metadata TSV (UTF-8/ASCII)"),
  make_option("--stratum-column", type = "character", dest = "stratum_column", default = "season",
              help = "Metadata column used as stratum/batch (normalized to 'season') [default %default]"),
  make_option("--condition-column", type = "character", dest = "condition_column", default = "aggregation_stage",
              help = "Metadata column used as condition (normalized to 'condition') [default %default]"),
  make_option("--batch-column", type = "character", dest = "batch_column", default = NULL,
              help = "Metadata batch column (normalized to 'batch'); defaults to stratum-column"),
  make_option("--output-rds", type = "character", dest = "output_rds", help = "Output RDS with counts and metadata"),
  make_option("--output-metadata", type = "character", dest = "output_metadata", help = "Output TSV with filtered metadata"),
  make_option("--output-cpm-all", type = "character", dest = "output_cpm_all", default = "",
              help = "Optional output TSV.gz with normalized CPM for all retained genes/samples")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$counts) || is.null(opt$metadata) ||
    is.null(opt$stratum_column) || is.null(opt$condition_column) ||
    is.null(opt$output_rds) || is.null(opt$output_metadata)) {
  stop("All arguments are required", call. = FALSE)
}

for (path in c(opt$output_rds, opt$output_metadata)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
if (!is.null(opt$output_cpm_all) && nzchar(trimws(opt$output_cpm_all))) {
  dir.create(dirname(opt$output_cpm_all), recursive = TRUE, showWarnings = FALSE)
}

read_metadata_table <- function(path) {
  result <- tryCatch(
    read_tsv(path, progress = FALSE, show_col_types = FALSE),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    stop(
      sprintf("Unable to read metadata table as UTF-8/ASCII: %s\n%s", path, conditionMessage(result)),
      call. = FALSE
    )
  }

  result
}

metadata_raw <- read_metadata_table(opt$metadata)

sample_id_candidates <- c(
  "sample_name",
  "sample_id"
)
sample_id_col <- sample_id_candidates[sample_id_candidates %in% names(metadata_raw)]
if (length(sample_id_col) == 0) {
  stop(
    paste0(
      "Metadata must contain one of sample ID columns: ",
      paste(sample_id_candidates, collapse = ", "),
      ". Found columns: ",
      paste(names(metadata_raw), collapse = ", ")
    ),
    call. = FALSE
  )
}
sample_id_col <- sample_id_col[[1]]

required_meta_cols <- c("run_accession", opt$stratum_column, opt$condition_column, "replicate")
missing_meta_cols <- setdiff(unique(required_meta_cols), names(metadata_raw))
if (length(missing_meta_cols) > 0) {
  stop(
    sprintf("Metadata missing required column(s): %s", paste(missing_meta_cols, collapse = ", ")),
    call. = FALSE
  )
}

if (!("hpd" %in% names(metadata_raw))) {
  metadata_raw$hpd <- NA_real_
}

batch_col <- if (!is.null(opt$batch_column) && nzchar(opt$batch_column)) {
  opt$batch_column
} else {
  opt$stratum_column
}

if (!(batch_col %in% names(metadata_raw))) {
  stop(sprintf("Metadata missing batch column: %s", batch_col), call. = FALSE)
}

as_stratum_factor <- function(x, column_name) {
  x_chr <- as.character(x)
  if (identical(column_name, "season")) {
    factor(x_chr, levels = c("Autumn", "Winter", "Spring", "Summer"))
  } else {
    factor(x_chr)
  }
}

as_condition_factor <- function(x, column_name) {
  x_chr <- as.character(x)
  if (identical(column_name, "aggregation_stage")) {
    factor(x_chr, levels = c("Body", "Cells", "Aggregates"))
  } else {
    factor(x_chr)
  }
}

metadata <- metadata_raw %>%
  transmute(
    sample_id = .data[[sample_id_col]],
    run = run_accession,
    season = as_stratum_factor(.data[[opt$stratum_column]], opt$stratum_column),
    condition = as_condition_factor(.data[[opt$condition_column]], opt$condition_column),
    hpd = hpd,
    replicate = as.integer(replicate),
    batch = as.character(.data[[batch_col]])
  ) %>%
  filter(!is.na(sample_id) & !is.na(run)) %>%
  filter(!is.na(season) & !is.na(condition))

if (identical(opt$stratum_column, "season") && identical(opt$condition_column, "aggregation_stage")) {
  metadata <- metadata %>%
    filter(!(season == "Summer" & condition == "Aggregates" & !is.na(hpd) & hpd != 24))
}

metadata <- metadata %>%
  mutate(
    season = droplevels(season),
    condition = droplevels(condition),
    season_condition = factor(paste(season, condition, sep = ".")),
    batch = factor(batch)
  ) %>%
  arrange(season, condition, replicate, sample_id)

counts_df <- read_tsv(
  opt$counts,
  comment = "#",
  progress = FALSE,
  col_types = cols(
    Geneid = col_character(),
    Chr = col_character(),
    Start = col_double(),
    End = col_double(),
    Strand = col_character(),
    Length = col_double(),
    .default = col_double()
  )
)

counts_mat <- as.matrix(counts_df[, -(1:6), drop = FALSE])
rownames(counts_mat) <- counts_df$Geneid
colnames(counts_mat) <- basename(colnames(counts_mat))

metadata <- metadata %>% filter(run %in% colnames(counts_mat))

missing_cols <- setdiff(metadata$run, colnames(counts_mat))
if (length(missing_cols) > 0) {
  stop(sprintf("Counts file missing runs: %s", paste(missing_cols, collapse = ", ")))
}

counts_mat <- counts_mat[, metadata$run, drop = FALSE]
colnames(counts_mat) <- metadata$sample_id

if (any(duplicated(metadata$sample_id))) {
  dupes <- metadata$sample_id[duplicated(metadata$sample_id)]
  stop(sprintf("Duplicate sample_ids found: %s", paste(unique(dupes), collapse = ", ")))
}

if (!is.null(opt$output_cpm_all) && nzchar(trimws(opt$output_cpm_all))) {
  dgl <- DGEList(counts = counts_mat)
  dgl <- calcNormFactors(dgl)
  cpm_all <- as.data.frame(edgeR::cpm(dgl)) %>%
    rownames_to_column("gene_id")
  write_tsv(cpm_all, opt$output_cpm_all)
}

saveRDS(list(counts = counts_mat, metadata = metadata), file = opt$output_rds)
write_tsv(metadata, opt$output_metadata)
