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
        combat=rules.combat_seq_counts.output.rds
    output:
        pdf="results/de/plots/pca_batch_correction.pdf"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/plot_pca_batch.R \
            --input-rds {input.data} \
            --combat-rds {input.combat} \
            --output {output.pdf}
        """


rule volcano_plot:
    input:
        results=rules.edgeR_results.output.tsv
    output:
        png="results/de/plots/volcano_plot.png"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/volcano_analysis.R \
            --results {input.results} \
            --mode plot \
            --output {output.png}
        """


rule volcano_stats:
    input:
        results=rules.edgeR_results.output.tsv
    output:
        tsv="results/de/stats/volcano_counts.tsv"
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/volcano_analysis.R \
            --results {input.results} \
            --mode stats \
            --output {output.tsv}
        """
