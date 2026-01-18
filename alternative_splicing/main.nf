nextflow.enable.dsl = 2

log.info "Launching alternative splicing backbone (Nextflow DSL2)"

params.outdir       = params.outdir ?: 'results'
params.samples      = params.samples ?: 'metadata/samples.tsv'
params.samplesheet  = params.samplesheet ?: 'metadata/samplesheet.csv'
params.bam_dir      = params.bam_dir ?: null
params.goi          = params.goi ?: null
params.gtf          = params.gtf ?: null
params.fasta        = params.fasta ?: null
params.palette      = params.palette ?: 'conf/palette.tsv'
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
            def condition  = (row.Sample ?: row.condition ?: 'NA')?.toString()?.trim()
            def dataset    = (row.SourceAbbr ?: row.Source ?: 'NA')?.toString()?.trim()
            def stranded   = (row.Stranded ?: row.strandedness ?: params.strandedness).toString()
            def readlen    = row.ReadLength ? row.ReadLength.toString().replaceAll('\\r','').trim() : null

            if (!run) {
                error "Each row must define Run/sample_id: ${row}"
            }
            if (!experiment) {
                error "Each row must define Experiment: ${row}"
            }
            if (!dataset) {
                dataset = 'NA'
            }

            tuple(run as String,
                  condition as String,
                  dataset as String,
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
    if (!params.gtf) {
        error "Parameter --gtf is required for this pipeline."
    }

    samples_ch = load_samples()
    goi_ch     = Channel.empty()

    if (params.goi) {
        goi_ch = RESOLVE_GOI(file(params.goi), file(params.gtf)).goi_resolved
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.gene_id.toString(), row.region.toString(), (row.label ?: row.description ?: row.gene_id).toString()) }
    }

    if (params.bam_dir) {
        alignment_ch = samples_ch.map { sid, cond, dataset, strand, readlen, exp ->
            def bam_path = file("${params.bam_dir}/${sid}.bam")
            if (!bam_path.exists()) {
                error "BAM not found for ${sid} in --bam_dir: ${bam_path}"
            }
            tuple(sid, cond, dataset, strand, readlen, bam_path)
        }
    } else {
        sheet_ch   = load_samplesheet()

        // Join samples metadata to fastq paths using Experiment accession
        samples_with_fastq = samples_ch
            .map { run, cond, dataset, strand, readlen, exp -> tuple(exp, run, cond, dataset, strand, readlen) }
            .join(sheet_ch)
            .map { exp, run, cond, dataset, strand, readlen, r1, r2 -> tuple(run, cond, dataset, strand, readlen, r1, r2) }

        align_input_ch = samples_with_fastq.map { sid, cond, dataset, strand, readlen, r1, r2 ->
            if (!readlen) {
                error "ReadLength is required for ${sid}"
            }
            tuple(sid, cond, dataset, strand, readlen, r1, r2)
        }

        // Build STAR indices per read length for samples that need alignment
        readlen_ch = align_input_ch.map { sid, cond, dataset, strand, readlen, r1, r2 -> readlen }
                                    .filter { it }
                                    .unique()
        star_indices = STAR_INDEX(readlen_ch)

        // fastp for samples with provided FASTQs
        def fastp_step = FASTP(align_input_ch)
        trimmed      = fastp_step.reads

        // Match trimmed reads to the correct STAR index by read length
        trimmed_grouped = trimmed
            .map { sid, cond, dataset, strand, readlen, r1, r2 -> tuple(readlen, tuple(sid, cond, dataset, strand, r1, r2)) }
            .groupTuple()

        aligned_input = trimmed_grouped.join(star_indices)
            .flatMap { readlen, samples, index_dir ->
                samples.collect { sample ->
                    def (sid, cond, dataset, strand, r1, r2) = sample
                    tuple(readlen, sid, cond, dataset, strand, r1, r2, index_dir)
                }
            }

        alignment_ch = STAR_ALIGN(aligned_input).aligned_bam
    }

    // Junction discovery and summaries (regtools)
    junctions = REGTOOLS_JUNCTIONS(alignment_ch).junctions
    junction_manifest = junctions
        .map { sid, cond, dataset, strand, readlen, bed ->
            def published = file("${params.outdir}/junctions/${sid}.junctions.bed").toAbsolutePath()
            "${sid}\t${cond}\t${dataset}\t${strand}\t${readlen}\t${published}"
        }
        .collectFile(name: "junctions_manifest.tsv", newLine: true,
                     sort: { line -> def toks = line.tokenize('\t'); [toks[2], toks[1], toks[0]] },
                     storeDir: "${params.outdir}/junctions")
    junction_summaries = REGTOOLS_SUMMARIZE(junction_manifest)

    // Temporarily disable downstream event/isoform steps while focusing on junction stats and sashimi plots
    // assemblies    = STRINGTIE_ASSEMBLE(alignment_ch)
    // merged_gtf    = STRINGTIE_MERGE(assemblies)
    // comparison    = GFFCOMPARE(merged_gtf)
    // junction_sets = LEAFCUTTER_PREP(alignment_ch)
    // rmats_events  = RMATS_PREP(alignment_ch)
    // transcriptome = BUILD_TRANSCRIPTOME(merged_gtf, params.fasta)
    // salmon_inputs = alignment_ch.combine(transcriptome)
    // quant         = SALMON_QUANT(salmon_inputs)

    if (params.goi && params.gtf) {
        bam_manifest = alignment_ch
            .map { sid, cond, dataset, strand, readlen, bam -> "${sid}\t${bam.toString()}\t${cond}" }
            .collectFile(name: "bam_manifest.tsv", newLine: true, sort: { it -> it.tokenize('\t')[2] },
                         storeDir: "${params.outdir}/sashimi")

        sashimi_inputs = goi_ch.combine(bam_manifest)
                            .map { g, r, l, manifest -> tuple(g, r, l, manifest, file(params.gtf), file(params.palette)) }
        sashimi_plots  = GGSASHIMI_PLOT(sashimi_inputs)
    }
}

process FASTP {
    tag { sample_id }
    publishDir "${params.fastp_dir}", mode: 'copy'
    cpus params.threads

    input:
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path(r1), path(r2)

    output:
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path("${sample_id}_clean_1.fastq.gz"), path("${sample_id}_clean_2.fastq.gz"), emit: reads
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
        tuple val(readlen), val(sample_id), val(condition), val(dataset), val(strandedness), path(r1), path(r2), path(index_dir)

    output:
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path("${sample_id}.bam"), emit: aligned_bam
        path "${sample_id}.bam.bai", emit: bam_index

    script:
    """
    set -euo pipefail
    ulimit -n 65536
    STAR \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${index_dir} \\
      --readFilesIn ${r1} ${r2} \\
      --readFilesCommand zcat \\
      --outFileNamePrefix "${sample_id}." \\
      --outSAMtype BAM SortedByCoordinate \\
      --outSAMstrandField intronMotif \\
      --alignIntronMax 50000 \\
      --alignMatesGapMax 50000 \\
      --outFilterMismatchNoverLmax 0.04 \\
      --twopassMode Basic
    mv "${sample_id}.Aligned.sortedByCoord.out.bam" "${sample_id}.bam"
    samtools index -@ ${task.cpus} "${sample_id}.bam"
    """
}

process REGTOOLS_JUNCTIONS {
    tag { sample_id }
    publishDir "${params.outdir}/junctions", mode: 'copy'
    conda "envs/regtools.yml"

    input:
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path(bam)

    output:
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path("${sample_id}.junctions.bed"), emit: junctions

    script:
        def strandFlag = strandedness.toUpperCase().contains('RF') ? 1
                       : strandedness.toUpperCase().contains('FR') ? 2
                       : 0
        """
        set -euo pipefail
        regtools junctions extract \\
          -a 8 -m 50 -M 500000 \\
          -s ${strandFlag} \\
          -o ${sample_id}.junctions.bed \\
          ${bam}
        """
}

process REGTOOLS_SUMMARIZE {
    tag { "regtools_summary" }
    publishDir "${params.outdir}/junctions", mode: 'copy'
    cpus 2
    conda "envs/regtools.yml"

    input:
        path manifest

    output:
        path "junction_counts_sample.tsv",  emit: counts_long
        path "junction_counts_matrix.tsv",  emit: counts_matrix
        path "junction_counts_dataset.tsv", emit: counts_dataset
        path "junctions_union.bed",         emit: junctions_union

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    import pandas as pd

    manifest_path = Path("${manifest}").resolve()
    if not manifest_path.exists():
        sys.exit(f"Manifest not found: {manifest_path}")

    names = ["sample_id", "condition", "dataset", "strandedness", "readlen", "bed_path"]
    manifest_df = pd.read_csv(manifest_path, sep="\\t", header=None, names=names, comment="#")
    manifest_df = manifest_df.dropna(how="all")
    if manifest_df.empty:
        sys.exit(f"No entries found in manifest: {manifest_path}")

    # Expand junction beds
    long_frames = []
    for _, row in manifest_df.iterrows():
        bed_file = Path(str(row.bed_path)).resolve()
        if not bed_file.exists():
            sys.exit(f"Missing junction bed: {bed_file}")
        bed_df = pd.read_csv(
            bed_file,
            sep="\\t",
            header=None,
            usecols=[0, 1, 2, 4, 5],
            names=["chrom", "start", "end", "score", "strand"],
            comment="#",
        )
        if bed_df.empty:
            continue
        bed_df = bed_df.assign(
            start=bed_df["start"].astype(int),
            end=bed_df["end"].astype(int),
            count=pd.to_numeric(bed_df["score"], errors="coerce").fillna(0).round().astype(int),
        )
        bed_df["strand"] = bed_df["strand"].where(bed_df["strand"].isin(["+", "-", "."]), ".")
        bed_df["junction_id"] = bed_df["chrom"] + ":" + bed_df["start"].astype(str) + "-" + bed_df["end"].astype(str) + ":" + bed_df["strand"]
        bed_df["sample_id"] = row.sample_id
        bed_df["condition"] = row.condition
        bed_df["dataset"] = row.dataset
        bed_df["readlen"] = row.readlen
        long_frames.append(bed_df[["junction_id", "chrom", "start", "end", "strand", "sample_id", "condition", "dataset", "readlen", "count"]])

    if not long_frames:
        sys.exit("No junctions found across beds.")

    long_df = pd.concat(long_frames, ignore_index=True)
    long_df = long_df.sort_values(["dataset", "condition", "sample_id", "junction_id"])
    long_df.to_csv("junction_counts_sample.tsv", sep="\\t", index=False)

    # Wide matrix (junction x sample)
    matrix_df = long_df.pivot_table(index="junction_id", columns="sample_id", values="count", aggfunc="sum", fill_value=0)
    matrix_df.reset_index(inplace=True)
    matrix_df.to_csv("junction_counts_matrix.tsv", sep="\\t", index=False)

    # Dataset-level sums
    dataset_df = (
        long_df.groupby(["junction_id", "dataset"], as_index=False)["count"]
        .sum()
        .sort_values(["dataset", "junction_id"])
    )
    dataset_df.to_csv("junction_counts_dataset.tsv", sep="\\t", index=False)

    # Union BED with total counts across samples
    union_df = (
        long_df.groupby(["junction_id", "chrom", "start", "end", "strand"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "total_count"})
        .sort_values(["chrom", "start", "end", "strand", "junction_id"])
    )
    union_df[["chrom", "start", "end", "junction_id", "total_count", "strand"]].to_csv(
        "junctions_union.bed", sep="\\t", index=False, header=False
    )
    """
}

process STRINGTIE_ASSEMBLE {
    tag { sample_id }
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    cpus params.threads

    input:
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path(bam)

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
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path(bam)

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
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path(bam)

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
        tuple val(sample_id), val(condition), val(dataset), val(strandedness), val(readlen), path(bam), path(transcript_fasta)

    output:
        path "${sample_id}.quant.sf"

    script:
    """
    # TODO: salmon quant --libType based on ${strandedness} using ${transcript_fasta}
    echo -e "Name\\tLength\\tTPM\\tNumReads" > ${sample_id}.quant.sf
    """
}


process BAM_MANIFEST {
    tag { "bam_manifest" }
    publishDir "${params.outdir}/sashimi", mode: 'copy'

    input:
        val alignments

    output:
        path "bam_manifest.tsv", emit: manifest

    script:
    """
    #!/usr/bin/env python
    from pathlib import Path
    alignments = ${alignments.inspect()}
    rows = []
    for entry in alignments:
        sid, cond, dataset, strand, readlen, bam = entry
        bam_path = Path(bam)
        if not bam_path.exists():
            raise SystemExit(f"ERROR: BAM not found in manifest build: {bam_path}")
        rows.append((dataset, cond, sid, bam_path.name))

    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    with open("bam_manifest.tsv", "w") as out:
        for dataset, cond, sid, bam_name in rows:
            out.write(f"{sid}\\t{bam_name}\\t{cond}\\t{dataset}\\n")
    """
}

process GGSASHIMI_PLOT {
    tag { gene_id }
    publishDir "${params.outdir}/sashimi", mode: 'copy'
    container "docker.io/guigolab/ggsashimi:latest"
    containerOptions '--entrypoint ""'

    input:
        tuple val(gene_id), val(region), val(label), path(bam_manifest), path(gtf), path(palette)

    output:
        path "${gene_id}.sashimi.pdf"
        path "${gene_id}.sashimi.grouped.pdf"

    script:
    """
    set -euo pipefail
    python /ggsashimi.py \\
      --bam "${bam_manifest}" \\
      --coordinates "${region}" \\
      --gtf "${gtf}" \\
      --palette "${palette}" \\
      --labels 3 \\
      --color-factor 3 \\
      --base-size 10 \\
      --width 13 \\
      -o "${gene_id}.sashimi"

    python /ggsashimi.py \\
      --bam "${bam_manifest}" \\
      --coordinates "${region}" \\
      --gtf "${gtf}" \\
      --palette "${palette}" \\
      --labels 3 \\
      --overlay 3 \\
      --color-factor 3 \\
      --base-size 10 \\
      --alpha 0.3 \\
      --aggr median_j \\
      --width 13 \\
      --min-coverage 5 \\
      -o "${gene_id}.sashimi.grouped"
    """
}

process RESOLVE_GOI {
    publishDir "${params.outdir}/sashimi", mode: "copy"
    conda "envs/pyranges.yml"

    input:
        path goi_tsv
        path gtf

    output:
        path "goi.resolved.tsv", emit: goi_resolved

    script:
    """
    #!/usr/bin/env python
    import sys
    import pandas as pd
    import pyranges as pr

    goi_path = "${goi_tsv}"
    gtf_path = "${gtf}"

    goi = pd.read_csv(goi_path, sep="\t")
    if "gene_id" not in goi.columns:
        sys.exit("ERROR:missing_gene_id_column")
    gene_ids = goi["gene_id"].dropna().astype(str).str.strip()
    if gene_ids.empty:
        sys.exit("ERROR:no_gene_ids")

    tx = pr.read_gtf(gtf_path)
    tx = tx[tx.Feature == "transcript"]
    tx = tx.query("gene_id in @gene_ids")
    tx_df = tx[["gene_id", "Chromosome", "Strand", "Start", "End"]]
    if tx_df.empty:
        sys.exit("ERROR:no_transcripts_for_goi")

    agg = (
        tx_df.groupby("gene_id")
             .agg({"Chromosome": "first", "Start": "min", "End": "max"})
             .reset_index()
    )
    agg["region"] = agg.apply(lambda r: f"{r.Chromosome}:{int(r.Start)}-{int(r.End)}", axis=1)

    merged = goi.merge(agg[["gene_id", "region"]], on="gene_id", how="left")
    if merged["region"].isnull().any():
        missing = merged[merged["region"].isnull()]["gene_id"].tolist()
        sys.exit(f"ERROR:region_not_found:{','.join(missing)}")

    merged["label"] = merged["description"].fillna("").str.strip()
    merged.loc[merged["label"] == "", "label"] = merged["gene_id"]

    merged[["gene_id", "region", "label"]].to_csv("goi.resolved.tsv", sep="\t", index=False)
    """
}
