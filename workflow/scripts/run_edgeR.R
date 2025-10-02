#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

option_list <- list(
  make_option(c("-i", "--input-rds"), type = "character", dest = "input_rds", help = "RDS produced by prepare_de_data.R"),
  make_option(c("-o", "--output-tsv"), type = "character", dest = "output_tsv", help = "Output TSV.gz with edgeR results")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input_rds) || is.null(opt$output_tsv)) {
  stop("--input-rds and --output-tsv are required", call. = FALSE)
}

de_data <- readRDS(opt$input_rds)
counts <- de_data$counts
metadata <- de_data$metadata

if (is.null(counts) || is.null(metadata)) {
  stop("Input RDS must contain 'counts' matrix and 'metadata' tibble", call. = FALSE)
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
  select(-order_index)

counts <- counts[, metadata$sample_id, drop = FALSE]
metadata <- metadata %>% mutate(sample_id = factor(sample_id, levels = colnames(counts)))

group <- interaction(metadata$season, metadata$condition, sep = ".", drop = TRUE)

dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, group = group)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

metadata <- metadata %>% mutate(sample_id = factor(as.character(sample_id), levels = colnames(dge)))
group <- droplevels(interaction(metadata$season, metadata$condition, sep = "."))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast_exprs <- c()

seasons <- levels(metadata$season)
for (season in seasons) {
  if (is.na(season)) next
  base <- paste0(season, ".Body")
  if (!(base %in% colnames(design))) next
  cells <- paste0(season, ".Cells")
  aggregates <- paste0(season, ".Aggregates")
  if (cells %in% colnames(design)) {
    contrast_exprs[paste0(season, "_Cells_vs_Body")] <- sprintf("%s - %s", cells, base)
  }
  if (aggregates %in% colnames(design)) {
    contrast_exprs[paste0(season, "_Aggregates_vs_Body")] <- sprintf("%s - %s", aggregates, base)
  }
}

if (length(contrast_exprs) == 0) {
  stop("No valid contrasts could be constructed", call. = FALSE)
}

contrasts <- makeContrasts(contrasts = contrast_exprs, levels = design)

fit <- glmQLFit(dge, design, robust = TRUE)

result_tables <- list()
pvalue_list <- list()

for (i in seq_len(ncol(contrasts))) {
  contrast_name <- colnames(contrasts)[i]
  qlf <- glmQLFTest(fit, contrast = contrasts[, i])
  tbl <- as.data.frame(qlf$table)
  tbl$gene_id <- rownames(tbl)
  tbl <- as_tibble(tbl) %>%
    relocate(gene_id)
  tbl$contrast <- contrast_name
  result_tables[[contrast_name]] <- tbl
  pvalue_list[[contrast_name]] <- tbl$PValue
}

adj_all <- p.adjust(unlist(pvalue_list), method = "BH")
start <- 1
for (name in names(result_tables)) {
  n <- nrow(result_tables[[name]])
  result_tables[[name]]$padj_global <- adj_all[start:(start + n - 1)]
  start <- start + n
}

results_long <- bind_rows(result_tables) %>%
  mutate(
    season = str_replace(contrast, "(_.*)$", ""),
    comparison = str_replace(contrast, "^[^_]+_", ""),
    condition = str_replace(comparison, "_vs_Body", ""),
    logCPM = logCPM,
    logFC = logFC,
    PValue = PValue,
    FDR = FDR
  ) %>%
  select(gene_id, season, condition, contrast, logFC, logCPM, PValue, FDR, padj_global)

for (path in c(opt$output_tsv)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

gz_con <- gzfile(opt$output_tsv, open = "wt")
write_tsv(results_long, gz_con)
close(gz_con)
