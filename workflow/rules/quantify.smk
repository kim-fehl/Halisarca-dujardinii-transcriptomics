from pathlib import Path


rule featurecounts:
    input:
        bam=expand("results/bam/{run}.sorted.bam", run=RUN_IDS),
        annotation=lambda wildcards: str(Path(config["genome"]["annotation"]).expanduser()),
        strand="results/qc/rseqc/featurecounts_strand.txt"
    output:
        counts=f"results/counts/counts_exons_{config['project']['name']}.tsv",
        log=f"results/counts/featureCounts_{config['project']['name']}.log",
        summary=f"results/counts/counts_exons_{config['project']['name']}.tsv.summary"
    params:
        extra=lambda wildcards: config["processing"].get("featurecounts_extra", "")
    threads:
        lambda wildcards: int(config["processing"]["max_threads"])
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p results/counts
        strand=$(cat {input.strand})
        featureCounts \
            -T {threads} \
            -a {input.annotation} \
            -o {output.counts} \
            -s ${strand} \
            {params.extra} \
            {input.bam} \
            2> {output.log}
        mv {output.counts}.summary {output.summary}
        """


rule compress_counts:
    input:
        f"results/counts/counts_exons_{config['project']['name']}.tsv"
    output:
        gz=f"results/counts/counts_exons_{config['project']['name']}.tsv.gz"
    shell:
        """
        sed 's/.sorted.bam//g' {input} | gzip -c > {output.gz}
        """
