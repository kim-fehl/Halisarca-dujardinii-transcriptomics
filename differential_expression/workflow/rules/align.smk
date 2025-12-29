from pathlib import Path

from snakemake.io import directory


def _read_length_minus_one(path):
    length = int(Path(path).read_text(encoding="ascii").strip())
    return max(1, length - 1)


def _star_index_extra(read_length_path, annotation_path):
    extras = [
        f"--sjdbGTFfile {annotation_path}",
        "--sjdbGTFtagExonParentTranscript Parent",
        f"--sjdbOverhang {_read_length_minus_one(read_length_path)}",
        f"--genomeSAindexNbases {config['processing'].get('star_sa_index_n_bases', 12)}",
    ]
    index_extra = config["processing"].get("star_index_extra", "").strip()
    if index_extra:
        extras.append(index_extra)
    return " ".join(extras)


def _star_align_extra(strand_path):
    strand = Path(strand_path).read_text(encoding="ascii").strip()
    extras = [
        "--outSAMtype BAM Unsorted",
        "--outSAMattributes NH HI AS nM XS",
    ]
    if strand == "UNSTRANDED":
        extras.append("--outSAMstrandField intronMotif")
    align_extra = config["processing"].get("star_align_extra", "").strip()
    if align_extra:
        extras.append(align_extra)
    return " ".join(extras)


rule star_genome_index:
    input:
        fasta=lambda wildcards: str(Path(config["genome"]["fasta"]).expanduser()),
        gff=lambda wildcards: str(Path(config["genome"]["annotation"]).expanduser()),
        read_length="results/qc/read_length.txt"
    output:
        directory("resources/genome/STAR_index")
    params:
        extra=lambda wildcards, input: _star_index_extra(input.read_length, input.gff)
    threads: MAX_THREADS
    log:
        "logs/star/genome_index.log"
    wrapper:
        "v7.2.0/bio/star/index"


rule star_align_single_end:
    input:
        fq1="results/fastp/{run}.fastp.fastq.gz",
        idx="resources/genome/STAR_index",
        strand="results/qc/rseqc/{run}.strand.txt"
    output:
        aln="results/bam/{run}_Aligned.out.bam",
        log="results/bam/{run}_Log.out",
        log_final="results/bam/{run}_Log.final.out",
        unmapped="results/bam/{run}_Unmapped.out.mate1"
    params:
        extra=lambda wildcards, input: _star_align_extra(input.strand)
    threads: MAX_THREADS
    log:
        "logs/star/align/{run}.log"
    wrapper:
        "v7.2.0/bio/star/align"


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
