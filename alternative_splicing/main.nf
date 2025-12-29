nextflow.enable.dsl = 2

log.info "Launching alternative splicing backbone (Nextflow DSL2)"

params.outdir       = params.outdir ?: 'results'
params.samples      = params.samples ?: 'metadata/samples.tsv'
params.gtf          = params.gtf ?: null
params.fasta        = params.fasta ?: null
params.strandedness = params.strandedness ?: 'fr-unstranded'
params.bam_dir      = params.bam_dir ?: 'results/bam'
params.stringtie_dir= params.stringtie_dir ?: 'resources/stringtie'
params.raw_dir      = params.raw_dir ?: "${params.outdir}/raw"
params.fastp_dir    = params.fastp_dir ?: "${params.outdir}/fastp"
params.star_index_base = params.star_index_base ?: 'resources/star_index'
params.threads      = params.threads ?: 8

// Channel to load and validate the sample sheet
def load_samples() {
    if (!params.samples) {
        error "Set --samples to a tab-delimited file; this project ships with metadata/samples.tsv (columns include Run, BiologicalState, Stranded, ReadLength)"
    }

    Channel.fromPath(params.samples)
        .ifEmpty { error "Sample sheet not found: ${params.samples}" }
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def run        = (row.Run ?: row.run ?: row.sample_id)?.toString()?.trim()
            def condition  = (row.BiologicalState ?: row.condition ?: 'NA')?.toString()?.trim()
            def stranded   = (row.Stranded ?: row.strandedness ?: params.strandedness).toString()
            def readlen    = row.ReadLength ? row.ReadLength.toString().replaceAll('\\r','').trim() : null

            if (!run) {
                error "Each row must define Run/sample_id: ${row}"
            }

            def bam_guess = params.bam_dir ? file("${params.bam_dir}/${run}.bam") : null
            def bam_path  = row.bam ? file(row.bam) : bam_guess

            def pre_gtf_guess = params.stringtie_dir ? file("${params.stringtie_dir}/${run}.gtf") : null
            def pre_gtf_file  = row.gtf ? file(row.gtf) : pre_gtf_guess
            def pre_gtf       = (pre_gtf_file && pre_gtf_file.exists()) ? pre_gtf_file : null

            tuple(run as String,
                  condition as String,
                  bam_path,
                  stranded,
                  pre_gtf,
                  readlen)
        }
}

workflow {
    samples_ch = load_samples()

    // Split samples into those with precomputed BAMs vs. those needing alignment
    def branched = samples_ch.branch { sid, cond, bam, strand, gtf, readlen ->
        bam && bam.exists() ? 'prealigned' : 'needs_alignment'
    }

    prealigned_ch = branched.prealigned.map { sid, cond, bam, strand, gtf, readlen ->
        tuple(sid, cond, strand, readlen, bam, gtf)
    }

    align_input_ch = branched.needs_alignment.map { sid, cond, bam, strand, gtf, readlen ->
        if (!readlen) {
            error "ReadLength is required for ${sid} when no precomputed BAM is available"
        }
        tuple(sid, cond, bam, strand, gtf, readlen)
    }

    // Build STAR indices per read length for samples that need alignment
    readlen_ch = align_input_ch.map { sid, cond, bam, strand, gtf, readlen -> readlen }
                                .filter { it }
                                .unique()
    star_indices = STAR_INDEX(readlen_ch)

    // Download + fastp for samples needing alignment
    raw_fastq    = DOWNLOAD_SRA(align_input_ch)
    def fastp_step = FASTP(raw_fastq)
    trimmed      = fastp_step.reads

    // Match trimmed reads to the correct STAR index by read length
    trimmed_keyed = trimmed.map { sid, cond, strand, readlen, gtf, r1, r2 ->
        tuple(readlen, sid, cond, strand, gtf, r1, r2)
    }
    aligned      = STAR_ALIGN(trimmed_keyed.join(star_indices)).aligned_bam

    alignment_ch = aligned.mix(prealigned_ch)

    junctions     = REGTOOLS_JUNCTIONS(alignment_ch.map { sid, cond, strand, readlen, bam, gtf -> tuple(sid, cond, bam, strand) })
    assemblies    = STRINGTIE_ASSEMBLE(alignment_ch)
    merged_gtf    = STRINGTIE_MERGE(assemblies.out.assembled_gtf)
    comparison    = GFFCOMPARE(merged_gtf.out.merged_gtf)
    junction_sets = LEAFCUTTER_PREP(alignment_ch.map{ sid, cond, strand, readlen, bam, gtf -> tuple(sid, cond, bam, strand) })
    rmats_events  = RMATS_PREP(alignment_ch.map{ sid, cond, strand, readlen, bam, gtf -> tuple(sid, cond, bam, strand) })

    transcriptome = BUILD_TRANSCRIPTOME(merged_gtf.out.merged_gtf, params.fasta)
    salmon_inputs = alignment_ch.map{ sid, cond, strand, readlen, bam, gtf -> tuple(sid, cond, bam, strand) }
                                .combine(transcriptome.out.transcript_fasta)
    quant         = SALMON_QUANT(salmon_inputs)
}

process DOWNLOAD_SRA {
    tag { sample_id }
    publishDir "${params.raw_dir}", mode: 'copy'
    cpus params.threads

    input:
        tuple val(sample_id), val(condition), val(bam), val(strandedness), path(pre_gtf) optional true, val(readlen)

    output:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path(pre_gtf) optional true, path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    script:
    """
    set -euo pipefail
    r1="${params.raw_dir}/${sample_id}_1.fastq.gz"
    r2="${params.raw_dir}/${sample_id}_2.fastq.gz"

    if [[ -s "\${r1}" && -s "\${r2}" ]]; then
        echo "[INFO] Reusing existing FASTQ for ${sample_id}"
        cp "\${r1}" "${sample_id}_1.fastq.gz"
        cp "\${r2}" "${sample_id}_2.fastq.gz"
        exit 0
    fi

    echo "[INFO] Downloading ${sample_id} to ${params.raw_dir}"
    mkdir -p "${params.raw_dir}"
    fasterq-dump --threads ${task.cpus} --split-files --outdir "${params.raw_dir}" "${sample_id}"
    pigz -p ${task.cpus} -f "${params.raw_dir}/${sample_id}"_1.fastq "${params.raw_dir}/${sample_id}"_2.fastq
    cp "\${r1}" "${sample_id}_1.fastq.gz"
    cp "\${r2}" "${sample_id}_2.fastq.gz"
    """
}

process FASTP {
    tag { sample_id }
    publishDir "${params.fastp_dir}", mode: 'copy'
    cpus params.threads

    input:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path(pre_gtf) optional true, path(r1), path(r2)

    output:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path(pre_gtf) optional true, path("${sample_id}_clean_1.fastq.gz"), path("${sample_id}_clean_2.fastq.gz"), emit: reads
        path "${sample_id}.fastp.html", emit: report_html
        path "${sample_id}.fastp.json", emit: report_json

    script:
    """
    fastp \\
      -i ${r1} -I ${r2} \\
      -o ${sample_id}_clean_1.fastq.gz \\
      -O ${sample_id}_clean_2.fastq.gz \\
      -w ${task.cpus} \\
      --detect_adapter_for_pe \\
      --cut_tail --cut_front --trim_poly_x \\
      --cut_mean_quality 20 \\
      --length_required 50 \\
      --html ${sample_id}.fastp.html \\
      --json ${sample_id}.fastp.json
    """
}

process STAR_INDEX {
    tag { readlen }
    cpus params.threads

    input:
        val readlen

    output:
        tuple val(readlen), path("star_index")

    when:
        params.fasta && params.gtf

    script:
    """
    set -euo pipefail
    overhang=$(( ${readlen} - 1 ))
    idx="${params.star_index_base}/${readlen}"
    mkdir -p "$idx"

    if [[ -s "$idx/SAindex" ]]; then
      echo "[INFO] Reusing STAR index for read length ${readlen} at $idx"
      echo "$idx" > index.path
      ln -s "$idx" .
      exit 0
    fi

    ulimit -n 65536
    echo "[INFO] Building STAR index for read length ${readlen} at $idx"
    STAR \\
      --runThreadN ${task.cpus} \\
      --runMode genomeGenerate \\
      --genomeDir "$idx" \\
      --genomeFastaFiles "${params.fasta}" \\
      --sjdbGTFfile "${params.gtf}" \\
      --sjdbOverhang "$overhang" \\
      --genomeSAindexNbases 12

    ln -s "$idx" star_index
    """
}

process STAR_ALIGN {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'
    cpus params.threads

    input:
        tuple val(readlen), val(sample_id), val(condition), val(strandedness), path(pre_gtf) optional true, path(r1), path(r2), path(index_dir)

    output:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path("${sample_id}.bam"), path(pre_gtf) optional true, emit: aligned_bam
        path "${sample_id}.bam.bai", emit: bam_index

    script:
    """
    set -euo pipefail
    prefix="${sample_id}."
    ulimit -n 65536
    STAR \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${index_dir} \\
      --readFilesIn ${r1} ${r2} \\
      --readFilesCommand zcat \\
      --outFileNamePrefix "$prefix" \\
      --outSAMtype BAM SortedByCoordinate \\
      --outSAMstrandField intronMotif \\
      --alignIntronMax 50000 \\
      --alignMatesGapMax 50000 \\
      --outFilterMismatchNoverLmax 0.04 \\
      --twopassMode Basic
    mv "${prefix}Aligned.sortedByCoord.out.bam" "${sample_id}.bam"
    samtools index -@ ${task.cpus} "${sample_id}.bam"
    """
}

process REGTOOLS_JUNCTIONS {
    tag { sample_id }
    publishDir "${params.outdir}/junctions", mode: 'copy'

    input:
        tuple val(sample_id), val(condition), path(bam), val(strandedness)

    output:
        path "${sample_id}.junctions.bed"

    script:
    """
    # TODO: regtools junctions extract --bam ${bam} > ${sample_id}.junctions.bed
    echo -e "# placeholder junctions for ${sample_id}\\t${condition}\\t${strandedness}" > ${sample_id}.junctions.bed
    """
}

process STRINGTIE_ASSEMBLE {
    tag { sample_id }
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    cpus params.threads

    input:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path(bam), path(pre_gtf) optional true

    output:
        path "${sample_id}.gtf", emit: assembled_gtf

    script:
    def strandFlag = strandedness.toUpperCase().contains('RF') ? '--rf'
                  : strandedness.toUpperCase().contains('FR') ? '--fr'
                  : ''
    """
    if [[ -s "${pre_gtf}" ]]; then
        echo "[INFO] Re-using precomputed GTF for ${sample_id}"
        cp "${pre_gtf}" "${sample_id}.gtf"
    else
        echo "[INFO] Running StringTie for ${sample_id}"
        stringtie "${bam}" -p ${task.cpus} -o "${sample_id}.gtf" -l "${sample_id}" ${strandFlag} ${params.gtf ? "-G ${params.gtf}" : ""}
    fi
    """
}

process STRINGTIE_MERGE {
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
        path(gtfs)

    output:
        path "merged_transcripts.gtf", emit: merged_gtf

    script:
    """
    stringtie --merge -G ${params.gtf} -o merged_transcripts.gtf ${gtfs}
    """
}

process GFFCOMPARE {
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
        path merged_gtf

    output:
        path "gffcompare.stats"

    script:
    """
    # TODO: gffcompare -r ${params.gtf} -o gffcompare ${merged_gtf}
    echo "placeholder stats vs reference ${params.gtf ?: 'unspecified'}" > gffcompare.stats
    """
}

process LEAFCUTTER_PREP {
    tag { sample_id }
    publishDir "${params.outdir}/leafcutter", mode: 'copy'

    input:
        tuple val(sample_id), val(condition), path(bam), val(strandedness)

    output:
        path "${sample_id}.leafcutter.junc"

    script:
    """
    # TODO: leafcutter clustering prep from junctions (see scripts/leafcutter_juncs.R)
    echo -e "${sample_id}\\t${condition}\\t${bam}\\t${strandedness}" > ${sample_id}.leafcutter.junc
    """
}

process RMATS_PREP {
    tag { sample_id }
    publishDir "${params.outdir}/rmats", mode: 'copy'

    input:
        tuple val(sample_id), val(condition), path(bam), val(strandedness)

    output:
        path "${sample_id}.bam.list"

    script:
    """
    # TODO: rmats input lists per condition; collate BAMs and run rmats.py --gtf ${params.gtf}
    echo "${bam}" > ${sample_id}.bam.list
    """
}

process BUILD_TRANSCRIPTOME {
    publishDir "${params.outdir}/transcriptome", mode: 'copy'

    input:
        path merged_gtf
        val genome_fasta

    output:
        path "transcriptome.fa", emit: transcript_fasta

    when:
        genome_fasta

    script:
    """
    # TODO: gffread -w transcriptome.fa -g ${genome_fasta} ${merged_gtf}
    echo ">placeholder_transcript" > transcriptome.fa
    echo "NNNN" >> transcriptome.fa
    """
}

process SALMON_QUANT {
    tag { sample_id }
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
        tuple val(sample_id), val(condition), path(bam), val(strandedness), path(transcript_fasta)

    output:
        path "${sample_id}.quant.sf"

    script:
    """
    # TODO: salmon quant --libType based on ${strandedness} using ${transcript_fasta}
    echo -e "Name\\tLength\\tTPM\\tNumReads" > ${sample_id}.quant.sf
    """
}
