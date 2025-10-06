rule prepare_de_data:
    input:
        counts="results/counts/counts_exons.tsv.gz",
        metadata="config/samples_metadata.tsv"
    output:
        rds="results/de/data/de_data.rds",
        metadata="results/de/data/sample_metadata.tsv"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/prepare_de_data.R \
            --counts {input.counts} \
            --metadata {input.metadata} \
            --output-rds {output.rds} \
            --output-metadata {output.metadata}
        """


rule combat_ref_repo:
    output:
        head="workflow/scripts/Combat-ref/.git/HEAD"
    shell:
        """
        set -euo pipefail
        mkdir -p workflow/scripts
        if [ -d workflow/scripts/Combat-ref ]; then
            rm -rf workflow/scripts/Combat-ref
        fi
        git clone --depth 1 --recurse-submodules https://github.com/xiaoyu12/Combat-ref workflow/scripts/Combat-ref
        git -C workflow/scripts/Combat-ref submodule update --init --recursive
        touch {output.head}
        """


rule combat_seq_counts:
    input:
        rds=rules.prepare_de_data.output.rds
    output:
        rds="results/de/data/combat_seq_counts.rds"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/run_combat_seq.R \
            --input-rds {input.rds} \
            --output-rds {output.rds}
        """


rule combat_ref_counts:
    input:
        rds=rules.prepare_de_data.output.rds,
        repo=rules.combat_ref_repo.output.head
    output:
        rds="results/de/data/combat_ref_counts.rds"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/run_combat_ref.R \
            --input-rds {input.rds} \
            --output-rds {output.rds} \
            --repo-dir workflow/scripts/Combat-ref
        """


rule edgeR_results:
    input:
        rds=rules.prepare_de_data.output.rds
    output:
        tsv="results/de/edgeR/results_long.tsv.gz"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/run_edgeR.R \
            --input-rds {input.rds} \
            --output-tsv {output.tsv}
        """


rule pca_batch_correction:
    input:
        data=rules.prepare_de_data.output.rds,
        combat_seq=rules.combat_seq_counts.output.rds,
        combat_ref=rules.combat_ref_counts.output.rds
    output:
        pdf="results/de/plots/pca_batch_correction.pdf"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/plot_pca_batch.R \
            --input-rds {input.data} \
            --combat-seq {input.combat_seq} \
            --combat-ref {input.combat_ref} \
            --output {output.pdf}
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
