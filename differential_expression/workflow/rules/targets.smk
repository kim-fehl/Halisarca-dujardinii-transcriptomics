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
    rules.pca_general.output.pdf,
]

FOUR_SEASONS_DE_TARGETS = [
    rules.edgeR_results_four_seasons.output.tsv,
    rules.volcano_plot_four_seasons.output.png,
    rules.volcano_stats_four_seasons.output.tsv,
    rules.heatmap_data_four_seasons.output.cpm_lfc,
    rules.heatmap_data_four_seasons.output.sample_stats,
    *HEATMAP_PDF_FOUR_SEASONS_TARGETS,
    *HEATMAP_XLSX_FOUR_SEASONS_TARGETS,
    *HEATMAP_ZSCORE_PDF_FOUR_SEASONS_TARGETS,
]
if DE_SAVE_ALL_GENES_CPM_TABLE:
    FOUR_SEASONS_DE_TARGETS.append(rules.all_genes_cpm_with_de_four_seasons.output.tsv)

GENERAL_DE_TARGETS = [
    rules.edgeR_results_general.output.tsv,
    rules.volcano_plot_general.output.png,
    rules.volcano_stats_general.output.tsv,
    rules.heatmap_data_general.output.cpm_lfc,
    rules.heatmap_data_general.output.sample_stats,
    *HEATMAP_PDF_GENERAL_TARGETS,
    *HEATMAP_XLSX_GENERAL_TARGETS,
    *HEATMAP_ZSCORE_PDF_GENERAL_TARGETS,
]
if DE_SAVE_ALL_GENES_CPM_TABLE:
    GENERAL_DE_TARGETS.append(rules.all_genes_cpm_with_de_general.output.tsv)
