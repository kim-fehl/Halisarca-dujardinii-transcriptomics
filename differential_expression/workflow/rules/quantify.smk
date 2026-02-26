from pathlib import Path


rule featurecounts:
    input:
        bam=expand("results/bam/{run}.sorted.bam", run=RUN_IDS),
        annotation=config["genome"]["annotation"],
        strand="results/qc/rseqc/featurecounts_strand.txt"
    output:
        counts=f"results/counts/counts_exons.tsv",
        summary=f"results/counts/counts_exons.summary.tsv"
    log:
        f"logs/featurecounts/featurecounts.log"
    params:
        paired_extra=lambda wildcards: "-p --countReadPairs" if IS_PAIRED_END else "",
        extra=lambda wildcards: config["processing"].get("featurecounts_extra", "")
    threads: MAX_THREADS
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
            {params.paired_extra} \
            -s $strand \
            {params.extra} \
            {input.bam} \
            > {log} 2>&1
        mv {output.counts}.summary {output.summary}
        """


rule compress_counts:
    input:
        f"results/counts/counts_exons.tsv"
    output:
        gz=f"results/counts/counts_exons.tsv.gz"
    shell:
        """
        sed 's/.sorted.bam//g' {input} | gzip -c > {output.gz}
        """
