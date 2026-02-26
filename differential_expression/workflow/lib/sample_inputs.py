import os
import re
from pathlib import Path

import pandas as pd


def _expand_path(path):
    return os.path.expanduser(path) if isinstance(path, str) else path


def _read_delimited_table(path, *, sep=None, label="table"):
    path = Path(path).expanduser()
    attempts = []
    for encoding in ("utf-16", "utf-8-sig", "utf-8"):
        try:
            kwargs = {"encoding": encoding}
            if sep is None:
                kwargs.update({"sep": None, "engine": "python"})
            else:
                kwargs["sep"] = sep
            df = pd.read_csv(path, **kwargs)
            df.columns = [str(col).strip() for col in df.columns]
            return df
        except Exception as exc:  # pandas parser/encoding errors vary
            attempts.append(f"{encoding}: {exc}")
    raise ValueError(
        f"Unable to read {label} at {path}. Tried encodings utf-16/utf-8-sig/utf-8.\n"
        + "\n".join(attempts)
    )


def _resolve_fastq_paths(value, base_path):
    if pd.isna(value):
        return []
    parts = [part.strip() for part in str(value).split(";")]
    parts = [part for part in parts if part and part.lower() != "nan"]
    resolved = []
    for part in parts:
        p = Path(os.path.expanduser(part))
        if not p.is_absolute():
            p = base_path / p
        resolved.append(str(p))
    return resolved


def _as_bool_or_none(value):
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    return None


def load_sample_context(config, cores):
    project_cfg = config.get("project", {}) or {}
    metadata_path = Path(project_cfg["metadata"]).expanduser()
    samplesheet_path = Path(project_cfg["samplesheet_path"]).expanduser()
    fastq_base_path = Path(_expand_path(project_cfg.get("fastq_path", "."))).expanduser()

    metadata_df = _read_delimited_table(metadata_path, sep="\t", label="metadata")
    if "run_accession" not in metadata_df.columns:
        raise ValueError(f"Metadata file '{metadata_path}' must contain a 'run_accession' column")
    metadata_df = metadata_df.dropna(subset=["run_accession"])
    metadata_df["run_accession"] = metadata_df["run_accession"].astype(str).str.strip()
    metadata_df = metadata_df[metadata_df["run_accession"] != ""]

    name_columns = [
        col
        for col in ("library_name", "sample_accession", "biosample")
        if col in metadata_df.columns
    ]

    if name_columns:
        def _derive_sample_name(row):
            for col in name_columns:
                value = row.get(col)
                if isinstance(value, str) and value.strip():
                    return value.strip()
            return row["run_accession"]
    else:
        def _derive_sample_name(row):
            return row["run_accession"]

    metadata_df["sample_id"] = (
        metadata_df.apply(_derive_sample_name, axis=1)
        .map(lambda name: re.sub(r"[^A-Za-z0-9._-]+", "_", name))
    )

    samplesheet_df = _read_delimited_table(samplesheet_path, label="samplesheet")
    run_col = "run_accession" if "run_accession" in samplesheet_df.columns else None
    if run_col is None and "run" in samplesheet_df.columns:
        run_col = "run"
    if run_col is None:
        raise ValueError(
            f"Samplesheet '{samplesheet_path}' must contain 'run_accession' (preferred) or 'run'"
        )

    fastq1_col = next(
        (col for col in ("fastq_1", "fastq", "fq1") if col in samplesheet_df.columns),
        None,
    )
    if fastq1_col is None:
        raise ValueError(
            f"Samplesheet '{samplesheet_path}' must contain one of: fastq_1, fastq, fq1"
        )
    fastq2_col = next((col for col in ("fastq_2", "fq2") if col in samplesheet_df.columns), None)
    single_end_col = "single_end" if "single_end" in samplesheet_df.columns else None

    run_to_fastq_r1 = {}
    run_to_fastq_r2 = {}
    run_to_layout = {}
    for row in samplesheet_df.to_dict("records"):
        run = str(row.get(run_col, "")).strip()
        if not run or run.lower() == "nan":
            continue
        fq1_paths = _resolve_fastq_paths(row.get(fastq1_col), fastq_base_path)
        fq2_paths = _resolve_fastq_paths(row.get(fastq2_col), fastq_base_path) if fastq2_col else []
        single_end = _as_bool_or_none(row.get(single_end_col)) if single_end_col else None
        if not fq1_paths:
            continue
        is_paired = bool(fq2_paths) or single_end is False
        if is_paired and not fq2_paths:
            raise ValueError(
                f"Run '{run}' is marked paired-end in '{samplesheet_path}' but fastq_2/fq2 is empty"
            )
        if single_end is True and fq2_paths:
            raise ValueError(
                f"Run '{run}' in '{samplesheet_path}' has fastq_2/fq2 but single_end=true"
            )
        if fq2_paths and len(fq1_paths) != len(fq2_paths):
            raise ValueError(
                f"Run '{run}' in '{samplesheet_path}' has {len(fq1_paths)} R1 file(s) and "
                f"{len(fq2_paths)} R2 file(s); provide matching semicolon-separated lists"
            )
        if run in run_to_fastq_r1:
            raise ValueError(
                f"Duplicate samplesheet entries for run '{run}'. "
                "Use one row per run and separate multiple FASTQs with ';' in fastq_1/fastq."
            )
        run_to_fastq_r1[run] = fq1_paths
        run_to_fastq_r2[run] = fq2_paths
        run_to_layout[run] = "PE" if is_paired else "SE"

    run_to_sample = dict(zip(metadata_df["run_accession"], metadata_df["sample_id"]))
    run_ids = sorted(run_to_sample.keys())
    if not run_ids:
        raise ValueError("No run_accession entries found in metadata")

    missing_runs = [run for run in run_ids if run not in run_to_fastq_r1]
    if missing_runs:
        preview = ", ".join(missing_runs[:10])
        suffix = "..." if len(missing_runs) > 10 else ""
        raise ValueError(
            f"{len(missing_runs)} metadata run(s) were not found in samplesheet '{samplesheet_path}': "
            f"{preview}{suffix}"
        )

    run_layouts = sorted({run_to_layout[run] for run in run_ids})
    if len(run_layouts) != 1:
        raise ValueError(
            "Mixed sequencing layouts are not supported in one run yet. "
            f"Found: {', '.join(run_layouts)} across metadata-selected runs."
        )

    max_threads = max(1, cores or 1)

    return {
        "metadata_df": metadata_df,
        "RUN_TO_FASTQ_R1": run_to_fastq_r1,
        "RUN_TO_FASTQ_R2": run_to_fastq_r2,
        "RUN_TO_SAMPLE": run_to_sample,
        "RUN_IDS": run_ids,
        "PRIMARY_RUN": run_ids[0],
        "IS_PAIRED_END": run_layouts[0] == "PE",
        "MAX_THREADS": max_threads,
        "AUX_THREADS": max(1, max_threads // 4),
    }
