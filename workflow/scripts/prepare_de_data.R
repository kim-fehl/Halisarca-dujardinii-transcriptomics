#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

option_list <- list(
  make_option("--counts", type = "character", help = "featureCounts TSV (gzipped)"),
  make_option("--metadata", type = "character", help = "Sample metadata TSV (UTF-16)"),
  make_option("--output-rds", type = "character", help = "Output RDS with counts and metadata"),
  make_option("--output-metadata", type = "character", help = "Output TSV with filtered metadata")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts) || is.null(opt$metadata) ||
    is.null(opt$output_rds) || is.null(opt$output_metadata)) {
  stop("All arguments are required", call. = FALSE)
}

for (path in c(opt$output_rds, opt$output_metadata)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

metadata <- read_tsv(
  opt$metadata,
  locale = locale(encoding = "UTF-16"),
  progress = FALSE
) %>%
  transmute(
    sample_id = sample_name,
    run = run_accession,
    season = factor(season, levels = c("Autumn", "Winter", "Spring", "Summer")),
    condition = factor(aggregation_stage, levels = c("Body", "Cells", "Aggregates")),
    hpd = hpd,
    replicate = as.integer(replicate)
  ) %>%
  filter(!is.na(sample_id) & !is.na(run)) %>%
  filter(!is.na(season) & !is.na(condition))

metadata <- metadata %>%
  filter(!(season == "Summer" & condition == "Aggregates" & !is.na(hpd) & hpd != 24))

metadata <- metadata %>%
  mutate(
    season = droplevels(season),
    condition = droplevels(condition),
    season_condition = factor(paste(season, condition, sep = ".")),
    batch = season
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

saveRDS(list(counts = counts_mat, metadata = metadata), file = opt$output_rds)
write_tsv(metadata, opt$output_metadata)
