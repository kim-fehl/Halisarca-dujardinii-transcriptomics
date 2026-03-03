DE_GENERAL_CFG = config.get("de", {}) or {}
DE_GENERAL_BASELINE_LEVEL = str(DE_GENERAL_CFG.get("baseline_level", "Body"))
_raw_general_compare_levels = DE_GENERAL_CFG.get("compare_levels", None)
if _raw_general_compare_levels is None:
    DE_GENERAL_COMPARE_LEVELS = None
elif isinstance(_raw_general_compare_levels, (list, tuple)):
    DE_GENERAL_COMPARE_LEVELS = [str(x) for x in _raw_general_compare_levels]
else:
    DE_GENERAL_COMPARE_LEVELS = [
        str(x).strip() for x in str(_raw_general_compare_levels).split(",") if str(x).strip()
    ]


rule edgeR_results_general:
    input:
        rds=rules.prepare_de_data.output.rds
    output:
        tsv="results/de/edgeR/results_long_general.tsv.gz"
    params:
        baseline_level=lambda wildcards: DE_GENERAL_BASELINE_LEVEL,
        compare_levels=lambda wildcards: ",".join(DE_GENERAL_COMPARE_LEVELS) if DE_GENERAL_COMPARE_LEVELS else ""
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/run_edgeR_general.R \
            --input-rds {input.rds} \
            --output-tsv {output.tsv} \
            --baseline-level '{params.baseline_level}' \
            --compare-levels '{params.compare_levels}'
        """


rule volcano_plot_general:
    input:
        results=rules.edgeR_results_general.output.tsv
    output:
        png="results/de/plots/volcano_plot_general.png"
    params:
        lfc_fold=lambda wildcards: float(config["processing"].get("volcano_fold_threshold", 1.5)),
        fdr=lambda wildcards: float(config["processing"].get("volcano_fdr_threshold", 0.01))
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/volcano_analysis_general.R \
            --results {input.results} \
            --mode plot \
            --output {output.png} \
            --lfc-fold {params.lfc_fold} \
            --fdr {params.fdr}
        """


rule volcano_stats_general:
    input:
        results=rules.edgeR_results_general.output.tsv
    output:
        tsv="results/de/stats/volcano_counts_general.tsv"
    params:
        lfc_fold=lambda wildcards: float(config["processing"].get("volcano_fold_threshold", 1.5)),
        fdr=lambda wildcards: float(config["processing"].get("volcano_fdr_threshold", 0.01))
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/volcano_analysis_general.R \
            --results {input.results} \
            --mode stats \
            --output {output.tsv} \
            --lfc-fold {params.lfc_fold} \
            --fdr {params.fdr}
        """
