DE_CFG = config.get("de", {}) or {}
DE_BASELINE_LEVEL = str(DE_CFG.get("baseline_level", "Body"))
_raw_compare_levels = DE_CFG.get("compare_levels", None)
if _raw_compare_levels is None:
    DE_COMPARE_LEVELS = None
elif isinstance(_raw_compare_levels, (list, tuple)):
    DE_COMPARE_LEVELS = [str(x) for x in _raw_compare_levels]
else:
    DE_COMPARE_LEVELS = [
        str(x).strip() for x in str(_raw_compare_levels).split(",") if str(x).strip()
    ]

HEATMAP_SET_NAME = (config.get("heatmap", {}) or {}).get("set_name", "24h_Autumn12")
HEATMAP_DATA_RDS = f"results/de/data/{HEATMAP_SET_NAME}_cpm_lfc_padj.rds"
HEATMAP_SAMPLE_STATS = f"results/de/data/{HEATMAP_SET_NAME}_samples_stats.edgeR.tsv"

if DE_BASELINE_LEVEL != "Body":
    raise ValueError(
        "The current 4-season DE/heatmap module assumes de.baseline_level='Body'. "
        "Use the default or generalize the downstream R scripts first."
    )
if DE_COMPARE_LEVELS is not None and sorted(DE_COMPARE_LEVELS) != ["Aggregates", "Cells"]:
    raise ValueError(
        "The current 4-season DE/heatmap module assumes de.compare_levels contains "
        "Cells and Aggregates. Use the default or generalize the downstream R scripts first."
    )


rule edgeR_results:
    input:
        rds=rules.prepare_de_data.output.rds
    output:
        tsv="results/de/edgeR/results_long.tsv.gz"
    params:
        baseline_level=lambda wildcards: DE_BASELINE_LEVEL,
        compare_levels=lambda wildcards: ",".join(DE_COMPARE_LEVELS) if DE_COMPARE_LEVELS else ""
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/run_edgeR.R \
            --input-rds {input.rds} \
            --output-tsv {output.tsv} \
            --baseline-level '{params.baseline_level}' \
            --compare-levels '{params.compare_levels}'
        """


rule volcano_plot:
    input:
        results=rules.edgeR_results.output.tsv
    output:
        png="results/de/plots/volcano_plot.png"
    params:
        lfc_fold=lambda wildcards: float(config["processing"].get("volcano_fold_threshold", 1.5)),
        fdr=lambda wildcards: float(config["processing"].get("volcano_fdr_threshold", 0.01))
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/volcano_analysis.R \
            --results {input.results} \
            --mode plot \
            --output {output.png} \
            --lfc-fold {params.lfc_fold} \
            --fdr {params.fdr}
        """


rule volcano_stats:
    input:
        results=rules.edgeR_results.output.tsv
    output:
        tsv="results/de/stats/volcano_counts.tsv"
    params:
        lfc_fold=lambda wildcards: float(config["processing"].get("volcano_fold_threshold", 1.5)),
        fdr=lambda wildcards: float(config["processing"].get("volcano_fdr_threshold", 0.01))
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/volcano_analysis.R \
            --results {input.results} \
            --mode stats \
            --output {output.tsv} \
            --lfc-fold {params.lfc_fold} \
            --fdr {params.fdr}
        """


rule heatmap_data:
    input:
        rds=rules.prepare_de_data.output.rds,
        results=rules.edgeR_results.output.tsv
    output:
        cpm_lfc=HEATMAP_DATA_RDS,
        sample_stats=HEATMAP_SAMPLE_STATS
    params:
        set_name=HEATMAP_SET_NAME
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/prepare_heatmap_inputs.R \
            --input-rds '{input.rds}' \
            --edgeR-tsv '{input.results}' \
            --set-name '{params.set_name}' \
            --output-rds '{output.cpm_lfc}' \
            --samples-stats '{output.sample_stats}'
        """


rule heatmap_plot:
    input:
        cpm=HEATMAP_DATA_RDS,
        sample_stats=HEATMAP_SAMPLE_STATS,
        geneset=lambda wildcards: HEATMAP_GENESET_MAP[wildcards.geneset]["path"]
    output:
        pdf=f"results/de/plots/heatmap_{HEATMAP_GENOME_NAME}_{{geneset}}.pdf",
        xlsx=f"results/de/edgeR/heatmap_{HEATMAP_GENOME_NAME}_{{geneset}}.xlsx"
    params:
        set_name=HEATMAP_SET_NAME,
        genome=HEATMAP_GENOME_NAME
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/heatmap_4seasons.R \
            --cpm-rds '{input.cpm}' \
            --sample-stats '{input.sample_stats}' \
            --geneset '{input.geneset}' \
            --set-name '{params.set_name}' \
            --genome-name '{params.genome}' \
            --output-pdf '{output.pdf}' \
            --output-xlsx '{output.xlsx}'
        """
