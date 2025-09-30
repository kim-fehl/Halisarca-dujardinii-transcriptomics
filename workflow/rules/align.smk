from pathlib import Path

from snakemake.io import directory


rule star_genome_index:
    input:
        fasta=lambda wildcards: str(Path(config["genome"]["fasta"]).expanduser()),
        gff=lambda wildcards: str(Path(config["genome"]["annotation"]).expanduser()),
        read_length="results/qc/read_length.txt"
    output:
        directory("resources/genome/STAR_index")
    log:
        "logs/star/genome_index.log"
    threads: MAX_THREADS
    params:
        star_sa_index_n_bases=lambda wildcards: config["processing"].get("star_sa_index_n_bases", 12)
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p {output}
        ulimit -n 65536
        read_length=$(cat {input.read_length})
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gff} \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang $((read_length - 1)) \
            --genomeSAindexNbases {params.star_sa_index_n_bases} > {log} 2>&1
        """


rule star_align_single_end:
    input:
        fastq="results/fastp/{run}.fastp.fastq.gz",
        index="resources/genome/STAR_index",
        strand="results/qc/rseqc/{run}.strand.txt"
    output:
        bam="results/bam/{run}_Aligned.out.bam"
    params:
        out_prefix=lambda wildcards: f"results/bam/{wildcards.run}_"
    threads: MAX_THREADS
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p results/bam
        ulimit -n 65536
        strand=$(cat {input.strand})
        strand_args=()
        if [[ "$strand" == "UNSTRANDED" ]]; then
            strand_args+=(--outSAMstrandField intronMotif)
        fi

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.fastq} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMattributes NH HI AS nM XS \
            "${{strand_args[@]}}"
        """


rule sort_bam:
    input:
        bam="results/bam/{run}_Aligned.out.bam"
    output:
        sorted="results/bam/{run}.sorted.bam"
    threads: AUX_THREADS
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        sambamba sort -t {threads} -o {output.sorted} {input.bam}
        """
