# _Halisarca dujardinii_ Differential Expression Pipeline

This repository contains a Snakemake workflow for processing _Halisarca dujardinii_ RNA-seq data end-to-end: download, preprocessing, alignment/quantification, differential expression (edgeR), QA/QC, and heatmap reporting of curated gene sets.

## Dataset description

The dataset (NCBI BioProject PRJNA594150) contains _H. dujardinii_ bulk RNA-seq samples spanning four seasons (Winter, Spring, Summer, Autumn) and three reaggregation stages (intact body, 30-min cell suspension post dissociation, 24-h aggregates). Most stages include three replicates; Autumn samples and the Summer intact-body condition retain two after QC. Because each season forms a separate sequencing batch, season is modeled as a biological covariate rather than fully batch-corrected. ComBat-based corrections are used only for PCA and other exploratory comparisons. Differential expression is evaluated within seasons using a quasi-likelihood GLM for the contrasts “cell suspension vs intact body” and “cell aggregates vs intact body”.

## Pipeline steps
- Data acquisition from NCBI: `sra-tools`
- QC and trimming: `fastp`
- Mapping: `STAR`
- Quantification: `featureCounts`
- Batch-adjusted counts (`ComBat-seq`, `ComBat-ref`) for PCA/correlation checks
- Replicate filtering based on PCA/correlation (removal of Autumn #3)
- Differential expression: `edgeR` (QL GLM)
- Visualization: PCA, volcano plots, heatmaps (`ComplexHeatmap`)

## Prerequisites

- Linux or macOS with [Mamba/Conda](https://mamba.readthedocs.io/) and [Snakemake](https://snakemake.readthedocs.io/) installed
- Sample metadata (`config/samples_metadata.tsv`)
- Genome FASTA/GFF3 resources under `resources/genome/`, as referenced in `config/config.yaml`
- Raw reads accessible to the workflow (downloaded automatically from NCBI)
- Gene-set `*.xlsx` tables whose paths are listed under `heatmap.genesets`; each sheet must provide `group` (may be blank), `gene_id`, `name`, and `descr`

## Quick start

1. Edit `config/config.yaml` and update:
   - `genome.fasta` and `genome.annotation`
   - Each entry under `heatmap.genesets` with the desired `name` (used in output filenames) and Excel `path`
2. Launch the workflow:

   ```bash
   snakemake --cores "$(nproc)" --use-conda
   ```

> **Important:** Always run the pipeline with `--use-conda` so rule-specific environments (especially the R dependencies for DE and heatmaps) are resolved automatically.

This command materializes every declared target under `results`, including QC, DE, and heatmaps. Re-run it whenever inputs or configuration change; Snakemake rebuilds only stale products.

To regenerate heatmaps for all configured gene sets only:

```bash
snakemake --cores 8 --use-conda heatmap_plot
```

## Key outputs

| Stage | Outputs |
| --- | --- |
| Counts & metadata | `results/counts/counts_exons.tsv.gz`, `results/de/data/de_data.rds` |
| QC | `results/qc/multiqc_fastp/*.html`, `results/qc/rseqc/featurecounts_strand.txt` |
| Batch assessment | `results/de/plots/pca_batch_correction.pdf` |
| Volcano plots | `results/de/plots/volcano_plot.png`, `results/de/stats/volcano_counts.tsv` |
| Heatmaps & gene-set tables | `results/de/plots/heatmap_<genome>_<geneset>.pdf`, `results/de/edgeR/heatmap_<genome>_<geneset>.xlsx` |
