if IS_PAIRED_END:
    rule fastp_paired_end:
        input:
            r1=rules.stage_fastq_paired_end.output.r1,
            r2=rules.stage_fastq_paired_end.output.r2
        output:
            r1="results/fastp/{run}_1.fastp.fastq.gz",
            r2="results/fastp/{run}_2.fastp.fastq.gz",
            html="results/fastp/{run}.fastp.html",
            json="results/fastp/{run}.fastp.json"
        log:
            "logs/fastp/{run}.log"
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
                -i {input.r1} \
                -I {input.r2} \
                -o {output.r1} \
                -O {output.r2} \
                --detect_adapter_for_pe \
                -z 9 \
                -l {params.min_length} \
                -h {output.html} \
                -j {output.json} \
                -R {params.sample} \
                {params.extra} > {log} 2>&1
            """
else:
    rule fastp_single_end:
        input:
            fastq=rules.stage_fastq_single_end.output.fastq
        output:
            fastq="results/fastp/{run}.fastp.fastq.gz",
            html="results/fastp/{run}.fastp.html",
            json="results/fastp/{run}.fastp.json"
        log:
            "logs/fastp/{run}.log"
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
                {params.extra} > {log} 2>&1
            """
