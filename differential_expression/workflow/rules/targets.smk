ruleorder: sort_bam > featurecounts


if IS_PAIRED_END:
    RAW_FASTQ_TARGETS = [
        *expand("resources/raw/{run}_1.fastq.gz", run=RUN_IDS),
        *expand("resources/raw/{run}_2.fastq.gz", run=RUN_IDS),
    ]
else:
    RAW_FASTQ_TARGETS = expand("resources/raw/{run}.fastq.gz", run=RUN_IDS)

SORTED_BAM_TARGETS = expand("results/bam/{run}.sorted.bam", run=RUN_IDS)

PRE_DE_TARGETS = [
    *RAW_FASTQ_TARGETS,
    *SORTED_BAM_TARGETS,
    "results/counts/counts_exons.tsv.gz",
    "results/qc/multiqc_fastp/multiqc_report.html",
    "results/qc/multiqc_fastp_star/multiqc_report.html",
    "results/qc/rseqc/featurecounts_strand.txt",
]

if BATCH_CORR_ENABLED:
    PRE_DE_TARGETS.append(rules.pca_batch_correction.output.pdf)

FOUR_SEASONS_DE_TARGETS = [
    rules.edgeR_results.output.tsv,
    rules.volcano_plot.output.png,
    rules.volcano_stats.output.tsv,
    rules.heatmap_data.output.cpm_lfc,
    rules.heatmap_data.output.sample_stats,
    *HEATMAP_PDF_TARGETS,
    *HEATMAP_XLSX_TARGETS,
]
