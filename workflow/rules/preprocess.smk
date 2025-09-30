rule fastp_single_end:
    input:
        fastq=rules.fasterq_dump.output.fastq
    output:
        fastq="results/fastp/{run}.fastp.fastq.gz",
        html="results/fastp/{run}.fastp.html",
        json="results/fastp/{run}.fastp.json"
    params:
        min_length=lambda wildcards: int(config["processing"].get("fastp_min_length", 25)),
        extra=lambda wildcards: config["processing"].get("fastp_extra", ""),
        sample=lambda wildcards: RUN_TO_SAMPLE[wildcards.run]
    threads: AUX_THREADS
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p results/fastp
        fastp \
            -w {threads} \
            -i {input.fastq} \
            -o {output.fastq} \
            -z 9 \
            -l {params.min_length} \
            -h {output.html} \
            -j {output.json} \
            -R {params.sample} \
            {params.extra}
        """
