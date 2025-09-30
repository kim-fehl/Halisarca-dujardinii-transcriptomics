rule prefetch_sra:
    output:
        sra="resources/sra/{run}.sra"
    threads: 1
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p resources/sra
        prefetch -p -v -O resources/sra {wildcards.run}
        """


rule fasterq_dump:
    input:
        sra=rules.prefetch_sra.output.sra
    output:
        fastq="resources/raw/{run}.fastq.gz"
    threads: MAX_THREADS
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p resources/raw
        fasterq-dump --threads {threads} --outdir resources/raw {wildcards.run}
        pigz -p {threads} -f resources/raw/{wildcards.run}.fastq
        """
