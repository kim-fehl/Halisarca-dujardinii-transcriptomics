nextflow.enable.dsl = 2

log.info "Launching alternative splicing backbone (Nextflow DSL2)"

params.outdir       = params.outdir ?: 'results'
params.samples      = params.samples ?: 'metadata/samples.tsv'
params.samplesheet  = params.samplesheet ?: 'metadata/samplesheet.csv'
params.gtf          = params.gtf ?: null
params.fasta        = params.fasta ?: null
params.strandedness = params.strandedness ?: 'fr-unstranded'
params.raw_dir      = params.raw_dir ?: "${params.outdir}/raw"
params.fastp_dir    = params.fastp_dir ?: "${params.outdir}/fastp"
params.star_index_base = params.star_index_base ?: 'resources/star_index'
params.threads      = params.threads ?: 8

// Channel to load and validate the sample metadata sheet
def load_samples() {
    if (!params.samples) {
        error "Set --samples to a tab-delimited file; this project ships with metadata/samples.tsv (columns include Run, BiologicalState, Stranded, ReadLength)"
    }

    Channel.fromPath(params.samples)
        .ifEmpty { error "Sample sheet not found: ${params.samples}" }
        .splitCsv(header: true, sep: '\t')
        .filter { row ->
            def flag = row.Use?.toString()?.trim()?.toUpperCase()
            flag in ['TRUE', 'YES', '1']
        }
        .map { row ->
            def run        = (row.Run ?: row.run ?: row.sample_id)?.toString()?.trim()
            def experiment = (row.Experiment ?: row.experiment)?.toString()?.trim()
            def condition  = (row.BiologicalState ?: row.condition ?: 'NA')?.toString()?.trim()
            def stranded   = (row.Stranded ?: row.strandedness ?: params.strandedness).toString()
            def readlen    = row.ReadLength ? row.ReadLength.toString().replaceAll('\\r','').trim() : null

            if (!run) {
                error "Each row must define Run/sample_id: ${row}"
            }
            if (!experiment) {
                error "Each row must define Experiment: ${row}"
            }

            tuple(run as String,
                  condition as String,
                  stranded,
                  readlen,
                  experiment as String)
        }
}

// Channel to load nf-core/fetchngs samplesheet with fastq paths
def load_samplesheet() {
    if (!params.samplesheet) {
        error "Set --samplesheet to a fetchngs samplesheet CSV (sample,fastq_1,fastq_2,experiment_accession,run_accession,...)"
    }

    Channel.fromPath(params.samplesheet)
        .ifEmpty { error "Samplesheet not found: ${params.samplesheet}" }
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def exp = row.experiment_accession?.toString()?.trim()
            def r1  = row.fastq_1?.toString()?.trim()
            def r2  = row.fastq_2?.toString()?.trim()
            if (!exp || !r1 || !r2) {
                error "Samplesheet rows must include experiment_accession, fastq_1, fastq_2: ${row}"
            }
            tuple(exp, file(r1), file(r2))
        }
}

workflow {
    samples_ch = load_samples()
    sheet_ch   = load_samplesheet()

    // Join samples metadata to fastq paths using Experiment accession
    samples_with_fastq = samples_ch
        .map { run, cond, strand, readlen, exp -> tuple(exp, run, cond, strand, readlen) }
        .join(sheet_ch)
        .map { exp, run, cond, strand, readlen, r1, r2 -> tuple(run, cond, strand, readlen, r1, r2) }

    align_input_ch = samples_with_fastq.map { sid, cond, strand, readlen, r1, r2 ->
        if (!readlen) {
            error "ReadLength is required for ${sid}"
        }
        tuple(sid, cond, strand, readlen, r1, r2)
    }

    // Build STAR indices per read length for samples that need alignment
    readlen_ch = align_input_ch.map { sid, cond, strand, readlen, r1, r2 -> readlen }
                                .filter { it }
                                .unique()
    star_indices = STAR_INDEX(readlen_ch)

    // fastp for samples with provided FASTQs
    def fastp_step = FASTP(align_input_ch)
    trimmed      = fastp_step.reads

    // Match trimmed reads to the correct STAR index by read length
    trimmed_keyed = trimmed.map { sid, cond, strand, readlen, r1, r2 ->
        tuple(readlen, sid, cond, strand, r1, r2)
    }
    alignment_ch = STAR_ALIGN(trimmed_keyed.join(star_indices)).aligned_bam

    junctions     = REGTOOLS_JUNCTIONS(alignment_ch.map { sid, cond, strand, readlen, bam -> tuple(sid, cond, bam, strand) })
    assemblies    = STRINGTIE_ASSEMBLE(alignment_ch)
    merged_gtf    = STRINGTIE_MERGE(assemblies.out.assembled_gtf)
    comparison    = GFFCOMPARE(merged_gtf.out.merged_gtf)
    junction_sets = LEAFCUTTER_PREP(alignment_ch.map{ sid, cond, strand, readlen, bam -> tuple(sid, cond, bam, strand) })
    rmats_events  = RMATS_PREP(alignment_ch.map{ sid, cond, strand, readlen, bam -> tuple(sid, cond, bam, strand) })

    transcriptome = BUILD_TRANSCRIPTOME(merged_gtf.out.merged_gtf, params.fasta)
    salmon_inputs = alignment_ch.map{ sid, cond, strand, readlen, bam -> tuple(sid, cond, bam, strand) }
                                .combine(transcriptome.out.transcript_fasta)
    quant         = SALMON_QUANT(salmon_inputs)
}

process FASTP {
    tag { sample_id }
    publishDir "${params.fastp_dir}", mode: 'copy'
    cpus params.threads

    input:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path(r1), path(r2)

    output:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path("${sample_id}_clean_1.fastq.gz"), path("${sample_id}_clean_2.fastq.gz"), emit: reads
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
    storeDir "${params.star_index_base}/readlen_${readlen}"

    input:
        val readlen

    output:
        tuple val(readlen), path("star_index")

    when:
        params.fasta && params.gtf

    script:
    def readlen_val = readlen as Integer
    def overhang = readlen_val - 1
    """
    set -euo pipefail
    ulimit -n 65536
    echo "[INFO] Building STAR index for read length ${readlen_val}"
    STAR \\
      --runThreadN ${task.cpus} \\
      --runMode genomeGenerate \\
      --genomeDir "star_index" \\
      --genomeFastaFiles "${params.fasta}" \\
      --sjdbGTFfile "${params.gtf}" \\
      --sjdbOverhang ${overhang} \\
      --genomeSAindexNbases 12
    """
}

process STAR_ALIGN {
    tag { sample_id }
    publishDir "${params.outdir}/bam", mode: 'copy'
    cpus params.threads

    input:
        tuple val(readlen), val(sample_id), val(condition), val(strandedness), path(r1), path(r2), path(index_dir)

    output:
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path("${sample_id}.bam"), emit: aligned_bam
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
        tuple val(sample_id), val(condition), val(strandedness), val(readlen), path(bam)

    output:
        path "${sample_id}.gtf", emit: assembled_gtf

    script:
    def strandFlag = strandedness.toUpperCase().contains('RF') ? '--rf'
                  : strandedness.toUpperCase().contains('FR') ? '--fr'
                  : ''
    """
    echo "[INFO] Running StringTie for ${sample_id}"
    stringtie "${bam}" -p ${task.cpus} -o "${sample_id}.gtf" -l "${sample_id}" ${strandFlag} ${params.gtf ? "-G ${params.gtf}" : ""}
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
