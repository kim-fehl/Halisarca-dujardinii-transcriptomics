from pathlib import Path


def _stage_fastq_bundle(srcs, out_path, run):
    import gzip
    import os
    import shutil

    srcs = [Path(p) for p in srcs]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for src in srcs:
        if not src.exists():
            raise FileNotFoundError(f"FASTQ input not found for run {run}: {src}")

    # Fast path: preserve original file without copying when there is exactly one gzipped FASTQ.
    if len(srcs) == 1 and str(srcs[0]).endswith(".gz"):
        if out_path.exists() or out_path.is_symlink():
            out_path.unlink()
        os.symlink(srcs[0].resolve(), out_path)
        return

    all_gz = all(str(src).endswith(".gz") for src in srcs)
    if all_gz:
        with out_path.open("wb") as dst:
            for src in srcs:
                with src.open("rb") as handle:
                    shutil.copyfileobj(handle, dst)
        return

    with gzip.open(out_path, "wb") as dst:
        for src in srcs:
            if str(src).endswith(".gz"):
                with gzip.open(src, "rb") as handle:
                    shutil.copyfileobj(handle, dst)
            else:
                with src.open("rb") as handle:
                    shutil.copyfileobj(handle, dst)


if IS_PAIRED_END:
    rule stage_fastq_paired_end:
        input:
            r1=lambda wildcards: RUN_TO_FASTQ_R1[wildcards.run],
            r2=lambda wildcards: RUN_TO_FASTQ_R2[wildcards.run]
        output:
            r1="resources/raw/{run}_1.fastq.gz",
            r2="resources/raw/{run}_2.fastq.gz"
        run:
            _stage_fastq_bundle(input.r1, output.r1, wildcards.run)
            _stage_fastq_bundle(input.r2, output.r2, wildcards.run)
else:
    rule stage_fastq_single_end:
        input:
            srcs=lambda wildcards: RUN_TO_FASTQ_R1[wildcards.run]
        output:
            fastq="resources/raw/{run}.fastq.gz"
        run:
            _stage_fastq_bundle(input.srcs, output.fastq, wildcards.run)
