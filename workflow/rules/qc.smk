from pathlib import Path


def _infer_star_extra(read_limit):
    extras = [
        "--readFilesCommand zcat",
        "--outSAMtype BAM Unsorted",
        "--outSAMattributes NH HI AS nM XS",
        f"--readMapNumber {read_limit}",
    ]
    custom = config["processing"].get("infer_experiment_star_extra", "").strip()
    if custom:
        extras.append(custom)
    return " ".join(extras)


rule detect_read_length:
    input:
        json=f"results/fastp/{PRIMARY_RUN}.fastp.json"
    output:
        "results/qc/read_length.txt"
    conda:
        "../envs/pipeline.yaml"
    run:
        import json

        out_path = Path(output[0])
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(input.json, "r", encoding="utf-8") as handle:
            data = json.load(handle)

        summary = data.get("summary", {})
        after = summary.get("after_filtering", {})
        read_length = after.get("read1_mean_length")
        if read_length is None:
            total_reads = after.get("total_reads")
            total_bases = after.get("total_bases")
            if total_reads:
                read_length = total_bases / total_reads

        if read_length is None:
            raise ValueError("Unable to determine read length from FASTQ")

        out_path.write_text(f"{int(round(read_length))}\n", encoding="ascii")


rule annotation_bed:
    input:
        lambda wildcards: str(Path(config["genome"]["annotation"]).expanduser())
    output:
        "resources/genome/annotation.bed"
    conda:
        "../envs/pipeline.yaml"
    run:
        dest = Path(output[0])
        dest.parent.mkdir(parents=True, exist_ok=True)

        def parse_attributes(field):
            data = {}
            for item in field.split(";"):
                item = item.strip()
                if not item:
                    continue
                if "=" in item:
                    key, value = item.split("=", 1)
                elif " " in item:
                    key, value = item.split(" ", 1)
                else:
                    continue
                data[key.strip()] = value.strip().strip('"')
            return data

        with open(input[0], "r", encoding="utf-8", errors="ignore") as src, dest.open("w", encoding="ascii") as dst:
            for line in src:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                feature = parts[2].lower()
                if feature != "exon":
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[3]) - 1
                    end = int(parts[4])
                except ValueError:
                    continue
                strand = parts[6]
                attrs = parse_attributes(parts[8])
                name = attrs.get("Parent") or attrs.get("ID") or f"{chrom}:{start}-{end}"
                dst.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


rule infer_strandness_map:
    input:
        fq1="results/fastp/{run}.fastp.fastq.gz",
        idx="resources/genome/STAR_index"
    output:
        aln="results/qc/rseqc/{run}.strand_infer.bam",
        log="results/qc/rseqc/{run}.strand_infer.Log.out",
        log_final="results/qc/rseqc/{run}.strand_infer.Log.final.out",
        unmapped="results/qc/rseqc/{run}.strand_infer.Unmapped.out.mate1",
        sj="results/qc/rseqc/{run}.strand_infer.SJ.out.tab"
    params:
        extra=lambda wildcards: _infer_star_extra(int(config["processing"].get("infer_experiment_read_limit", 200000)))
    threads: AUX_THREADS
    log:
        "logs/star/{run}.infer_strandness_map.log"
    wrapper:
        "v7.2.0/bio/star/align"


rule infer_strandness_report:
    input:
        bam=rules.infer_strandness_map.output.aln,
        annotation=rules.annotation_bed.output
    output:
        report="results/qc/rseqc/{run}.infer_experiment.txt"
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        infer_experiment.py -r {input.annotation} -i {input.bam} > {output.report}
        """


rule infer_strandness_parse:
    input:
        report=rules.infer_strandness_report.output.report
    output:
        strand="results/qc/rseqc/{run}.strand.txt"
    script:
        "../scripts/parse_infer_experiment.py"


rule summarize_strandness:
    input:
        strands=expand("results/qc/rseqc/{run}.strand.txt", run=RUN_IDS)
    output:
        "results/qc/rseqc/featurecounts_strand.txt",
        summary="results/qc/rseqc/strandness.tsv"
    conda:
        "../envs/pipeline.yaml"
    run:
        from collections import Counter

        strands = []
        for path in input.strands:
            path_obj = Path(path)
            data = path_obj.read_text(encoding="utf-8").strip()
            strands.append((path_obj.stem.split('.')[0], data))

        counts = Counter(label for _, label in strands)
        mode, _ = counts.most_common(1)[0]

        strand_flag = "0"
        if mode == "FR_FIRSTSTRAND":
            strand_flag = "2"
        elif mode == "FR_SECONDSTRAND":
            strand_flag = "1"

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text(strand_flag + "\n", encoding="ascii")

        summary_path = Path(output.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", encoding="ascii") as handle:
            handle.write("run\tstrand\n")
            for run, label in strands:
                handle.write(f"{run}\t{label}\n")
            handle.write(f"mode\t{mode}\n")
            handle.write(f"featureCounts_flag\t{strand_flag}\n")


rule multiqc_fastp:
    input:
        json=expand("results/fastp/{run}.fastp.json", run=RUN_IDS)
    output:
        "results/qc/multiqc/multiqc_report.html"
    conda:
        "../envs/pipeline.yaml"
    shell:
        """
        mkdir -p results/qc/multiqc
        multiqc results/fastp --outdir results/qc/multiqc --force
        """
