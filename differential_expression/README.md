# _Halisarca dujardinii_ Differential Expression Pipeline

This repository contains a Snakemake workflow for processing _Halisarca dujardinii_ RNA-seq data from FASTQ inputs through preprocessing, alignment/quantification, differential expression (edgeR), QA/QC, and heatmap reporting of curated gene sets.

## Dataset description

The dataset (NCBI BioProject PRJNA594150) contains _H. dujardinii_ bulk RNA-seq samples spanning four seasons (Winter, Spring, Summer, Autumn) and three reaggregation stages (intact body, 30-min cell suspension post dissociation, 24-h aggregates). Most stages include three replicates; Autumn samples and the Summer intact-body condition retain two after QC. Because each season forms a separate sequencing batch, season is modeled as a biological covariate rather than fully batch-corrected. ComBat-based corrections are used only for PCA and other exploratory comparisons. Differential expression is evaluated within seasons using a quasi-likelihood GLM for the contrasts “cell suspension vs intact body” and “cell aggregates vs intact body”.

## Pipeline steps
- External FASTQ acquisition/preparation (for example `nf-core/fetchngs`, not embedded in this workflow)
- QC and trimming: `fastp`
- Mapping: `STAR`
- Quantification: `featureCounts`
- Batch-adjusted counts (`ComBat-seq`, `ComBat-ref`) for PCA/correlation checks
- Replicate filtering based on PCA/correlation (removal of Autumn #3)
- Differential expression: `edgeR` (QL GLM)
- Visualization: PCA, volcano plots, heatmaps (`ComplexHeatmap`)

## Prerequisites

- Linux or macOS with [Mamba/Conda](https://mamba.readthedocs.io/) and [Snakemake](https://snakemake.readthedocs.io/) installed
- Sample metadata table (default: `config/samples_metadata.tsv`)
- Samplesheet with FASTQ paths (configured via `project.samplesheet_path`)
- Genome FASTA/GFF3 resources under `resources/genome/`, as referenced in `config/config.yaml`
- Raw reads already available on disk (paths come from the samplesheet; relative paths are resolved against `project.fastq_path`)
- Gene-set `*.xlsx` tables whose paths are listed under `genesets`; each sheet must provide `group` (may be blank), `gene_id`, `name`, and `descr`

## Quick start

1. Edit `config/config.yaml` and update:
   - `project.metadata`
   - `project.samplesheet_path`
   - `project.fastq_path` (base directory for relative FASTQ paths in the samplesheet)
   - `genome.fasta` and `genome.annotation`
   - `de.stratum_column`, `de.condition_column`, `de.baseline_level` (defaults match the current 4-season dataset)
   - `batch_correction.enabled` (`auto`/`true`/`false`; `auto` runs ComBat/PCA only when the stratum has >1 level)
   - `heatmap.set_name` (current 4-season heatmap subset label; used in intermediate filenames)
   - Each entry under `genesets` with the desired `name` (used in output filenames) and Excel `path`
2. Build all pre-DE outputs (FASTQ staging, trimming, mapping, counts, QC, and optional batch-correction PCA when enabled):

   ```bash
   snakemake --cores "$(nproc)" --use-conda pre_de
   ```

3. Build the current 4-season DE/heatmap module (or run `all` for everything):

   ```bash
   snakemake --cores "$(nproc)" --use-conda de_four_seasons
   ```

> **Important:** Always run the pipeline with `--use-conda` so rule-specific environments (especially the R dependencies for DE and heatmaps) are resolved automatically.

Use `snakemake --cores "$(nproc)" --use-conda all` to materialize everything in one run. Re-run any target whenever inputs or configuration change; Snakemake rebuilds only stale products.

To regenerate heatmaps for all configured gene sets only:

```bash
snakemake --cores 8 --use-conda heatmap_plot
```

## Input preparation

### Option A: `nf-core/fetchngs` (recommended for SRA/ENA accessions)

This workflow no longer downloads from NCBI directly. Use `fetchngs` first, then point Snakemake at the generated samplesheet and FASTQs.

Example outline (adjust to your `fetchngs` version and profile):

```bash
nextflow run nf-core/fetchngs \
  --input accessions.csv \
  --outdir fetchngs_out \
  -profile conda
```

Then set in `config/config.yaml`:

```yaml
project:
  metadata: config/samples_metadata.tsv
  samplesheet_path: fetchngs_out/samplesheet/samplesheet.csv
  fastq_path: fetchngs_out
```

Notes:
- `fetchngs` FASTQ columns are often relative paths; this workflow prepends `project.fastq_path` to relative entries.
- One row per `run_accession` is expected.
- If a run has multiple FASTQ chunks/lanes, place them in one cell separated by `;` (the workflow will concatenate them).
- Paired-end runs are supported via `fastq_1` + `fastq_2` (with matching semicolon-separated chunk counts).
- Keep one sequencing layout per run set (all selected runs SE or all selected runs PE); mixed SE/PE input is not supported in one execution.

### Option B: Local FASTQs

You can provide your own samplesheet and FASTQ files without `fetchngs`.

Recommended `samplesheet_path` columns (minimum):

| Column | Required | Notes | 
| --- | --- | --- |
| `run_accession` | Yes | Primary key used to join with metadata |
| `fastq_1` (or `fastq`) | Yes | Relative to `project.fastq_path` or absolute path |
| `fastq_2` | For PE | Required for paired-end runs; leave empty for single-end |
| `single_end` | No | Optional; `true/1` for SE, `false/0` for PE |

For multiple files per run (e.g. lanes), put semicolon-separated lists in `fastq_1` (and `fastq_2` for PE). Keep all selected runs the same layout within one workflow run.

Recommended metadata columns for the current 4-season DE/PCA module (with defaults `de.stratum_column=season`, `de.condition_column=aggregation_stage`):

| Column | Required | Notes |
| --- | --- | --- |
| `run_accession` | Yes | Must match the samplesheet |
| `sample_name` or `sample_id` | Yes | Human-readable sample label |
| `season` | Yes | Expected values: `Autumn`, `Winter`, `Spring`, `Summer` |
| `aggregation_stage` | Yes | Expected values: `Body`, `Cells`, `Aggregates` |
| `replicate` | Yes | Integer replicate index |
| `hpd` | Recommended | Used for the current aggregate filtering logic (24 h aggregates) |

The metadata reader accepts UTF-16 (legacy table) and UTF-8 TSV.

## Targets

- `pre_de`: FASTQ staging, trimming, mapping, counts, QC reports, and PCA batch-correction plot (no edgeR/volcano/heatmaps)
- `de_four_seasons`: current season-specific DE outputs, volcano plot/stat summary, and heatmaps
- `all`: `pre_de` + `de_four_seasons`

Notes:
- `de_four_seasons` is still dataset-specific (expects `de.baseline_level: Body` and the current 4-season downstream scripts).
- `pre_de` includes the ComBat/PCA plot only when batch correction is enabled (or `auto` detects >1 batch level).

## Key outputs

| Stage | Outputs |
| --- | --- |
| Counts & metadata | `results/counts/counts_exons.tsv.gz`, `results/de/data/de_data.rds` |
| QC | `results/qc/multiqc_fastp/*.html`, `results/qc/rseqc/featurecounts_strand.txt` |
| Batch assessment (when batch correction enabled) | `results/de/plots/pca_batch_correction.pdf` |
| Volcano plots | `results/de/plots/volcano_plot.png`, `results/de/stats/volcano_counts.tsv` |
| Heatmaps & gene-set tables | `results/de/plots/heatmap_<genome>_<geneset>.pdf`, `results/de/edgeR/heatmap_<genome>_<geneset>.xlsx` |
