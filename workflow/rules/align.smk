from pathlib import Path

from snakemake.io import directory


rule star_genome_index:
    input:
        fasta=lambda wildcards: str(Path(config["genome"]["fasta"]).expanduser()),
        gff=lambda wildcards: str(Path(config["genome"]["annotation"]).expanduser()),
        read_length="results/qc/read_length.txt"
    output:
        directory(str(Path(config["genome"]["star_index_dir"]).expanduser()))
    threads: MAX_THREADS
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p {output}
        read_length=$(cat {input.read_length})
        overhang=$((read_length - 1))
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gff} \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang ${overhang} \
            --genomeSAindexNbases {config["processing"].get("star_sa_index_n_bases", 12)}
        """


rule star_align_single_end:
    input:
        fastq=rules.fastp_single_end.output.fastq,
        index=rules.star_genome_index.output,
        strand=rules.infer_strandness_report.output.strand
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
            "${strand_args[@]}"
        """


rule sort_bam:
    input:
        bam=rules.star_align_single_end.output.bam
    output:
        sorted="results/bam/{run}.sorted.bam"
    threads: AUX_THREADS
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        sambamba sort -t {threads} -o {output.sorted} {input.bam}
        """
