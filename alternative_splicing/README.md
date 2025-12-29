# Alternative Splicing Pipeline (Nextflow)

Nextflow DSL2 scaffold for alternative splicing analyses on genome-aligned paired-end RNA-seq BAMs. The workflow is organized for junction discovery, event-level calls, isoform assembly/quantification, and reporting.

## Layout
- `main.nf` – pipeline entry with stubbed processes for junction extraction (regtools), event callers (rMATS/LeafCutter), assembly/quantification (StringTie/gffcompare/Salmon), and plotting hooks.
- `nextflow.config` – defaults for params, executors, and log/output locations; extend with per-environment profiles under `conf/`.
- `metadata/` – sample sheet lives here (see `metadata/samples.tsv` template).
- `resources/` – reference genome/annotation (FASTA/GTF) and any tool-specific indexes.
- `modules/local/`, `bin/` – place custom Nextflow modules and helper scripts.
- `results/`, `logs/` – pipeline outputs; `work/` is created by Nextflow at runtime.

## Metadata (uses your provided table)
`metadata/samples.tsv` is the real dataset you supplied (columns: `Run`, `BiologicalState`, `Stranded`, `ReadLength`, ...). The pipeline auto-derives:
- `sample_id` ← `Run`
- `condition` ← `BiologicalState`
- `strandedness` ← `Stranded` (expects `RF` or `FR`; falls back to `params.strandedness`)
- `readlen` ← `ReadLength` (used to pick/build STAR indices)
- `bam` (optional) ← `${params.bam_dir}/${Run}.bam` unless an explicit `bam` column is present (default `results/bam`)
- `precomputed GTF` (optional) ← `${params.stringtie_dir}/${Run}.gtf` unless an explicit `gtf` column is present

If you want to skip work for some samples, drop ready BAMs into `results/bam/<Run>.bam` (with `.bai`) and optional per-run StringTie GTFs into `resources/stringtie/<Run>.gtf`. The pipeline will reuse them and only download/align the remaining runs. You can change `--bam_dir` if you prefer another location.

## Running
Update `gtf` and `fasta` (reference) in `nextflow.config` or via `-params-file` and launch:

```bash
cd alternative_splicing
nextflow run main.nf -profile local \
  --samples metadata/samples.tsv \
  --gtf resources/reference.gtf \
  --fasta resources/reference.fasta \
  --outdir results
```

Outputs land under `results/` and `logs/`. Swap `-profile` once you add HPC/cloud settings under `conf/`.

## What the pipeline does now (StringTie.sh parity)
- Downloads SRA runs (from the `Run` column) with `fasterq-dump` → gzip
- QC/trim with `fastp`
- Builds STAR indices per read length automatically under `resources/star_index/`
- Aligns with STAR → sorted BAM + index
- Assembles with StringTie (reuses per-run GTFs if provided)
- Downstream placeholders: regtools junctions, gffcompare merge/compare, LeafCutter/rMATS prep, transcriptome build, Salmon quant

## Requirements
- SRA Toolkit (`prefetch`/`fasterq-dump`), `fastp`, `STAR`, `samtools`, `stringtie`
- Reference FASTA/GTF set via `--fasta` / `--gtf`
