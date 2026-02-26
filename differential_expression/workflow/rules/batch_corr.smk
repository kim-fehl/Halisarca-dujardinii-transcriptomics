from pathlib import Path


DE_CFG = config.get("de", {}) or {}
DE_STRATUM_COLUMN = DE_CFG.get("stratum_column", "season")
DE_CONDITION_COLUMN = DE_CFG.get("condition_column", "aggregation_stage")

BATCH_CFG = config.get("batch_correction", {}) or {}
BATCH_CORR_MODE = str(BATCH_CFG.get("enabled", "auto")).strip().lower()
if BATCH_CORR_MODE not in {"auto", "true", "false"}:
    raise ValueError("batch_correction.enabled must be one of: auto, true, false")

BATCH_CORR_BATCH_COLUMN = BATCH_CFG.get("batch_column") or DE_STRATUM_COLUMN

_batch_levels = []
if BATCH_CORR_BATCH_COLUMN in metadata_df.columns:
    _batch_series = (
        metadata_df.loc[metadata_df["run_accession"].isin(RUN_IDS), BATCH_CORR_BATCH_COLUMN]
        .dropna()
        .astype(str)
        .str.strip()
    )
    _batch_levels = sorted([x for x in _batch_series.unique().tolist() if x])

if BATCH_CORR_MODE == "true":
    if BATCH_CORR_BATCH_COLUMN not in metadata_df.columns:
        raise ValueError(
            f"batch_correction.batch_column '{BATCH_CORR_BATCH_COLUMN}' not found in metadata"
        )
    BATCH_CORR_ENABLED = True
elif BATCH_CORR_MODE == "false":
    BATCH_CORR_ENABLED = False
else:
    BATCH_CORR_ENABLED = len(_batch_levels) > 1


rule prepare_de_data:
    input:
        counts="results/counts/counts_exons.tsv.gz",
        metadata=lambda wildcards: str(Path(config["project"]["metadata"]).expanduser())
    output:
        rds="results/de/data/de_data.rds",
        metadata="results/de/data/sample_metadata.tsv"
    params:
        stratum_column=lambda wildcards: DE_STRATUM_COLUMN,
        condition_column=lambda wildcards: DE_CONDITION_COLUMN,
        batch_column=lambda wildcards: BATCH_CORR_BATCH_COLUMN
    conda:
        "../envs/r_de.yaml"
    shell:
        """
        Rscript workflow/scripts/prepare_de_data.R \
            --counts {input.counts} \
            --metadata {input.metadata} \
            --stratum-column '{params.stratum_column}' \
            --condition-column '{params.condition_column}' \
            --batch-column '{params.batch_column}' \
            --output-rds {output.rds} \
            --output-metadata {output.metadata}
        """


if BATCH_CORR_ENABLED:
    rule combat_ref_repo:
        output:
            head="workflow/scripts/Combat-ref/.git/HEAD"
        shell:
            """
            set -euo pipefail
            mkdir -p workflow/scripts
            if [ -d workflow/scripts/Combat-ref ]; then
                rm -rf workflow/scripts/Combat-ref
            fi
            git clone --depth 1 --recurse-submodules https://github.com/xiaoyu12/Combat-ref workflow/scripts/Combat-ref
            git -C workflow/scripts/Combat-ref submodule update --init --recursive
            touch {output.head}
            """


    rule combat_seq_counts:
        input:
            rds=rules.prepare_de_data.output.rds
        output:
            rds="results/de/data/combat_seq_counts.rds"
        conda:
            "../envs/r_de.yaml"
        shell:
            """
            Rscript workflow/scripts/run_combat_seq.R \
                --input-rds {input.rds} \
                --output-rds {output.rds}
            """


    rule combat_ref_counts:
        input:
            rds=rules.prepare_de_data.output.rds,
            repo=rules.combat_ref_repo.output.head
        output:
            rds="results/de/data/combat_ref_counts.rds"
        conda:
            "../envs/r_de.yaml"
        shell:
            """
            Rscript workflow/scripts/run_combat_ref.R \
                --input-rds {input.rds} \
                --output-rds {output.rds} \
                --repo-dir workflow/scripts/Combat-ref
            """


    rule pca_batch_correction:
        input:
            data=rules.prepare_de_data.output.rds,
            combat_seq=rules.combat_seq_counts.output.rds,
            combat_ref=rules.combat_ref_counts.output.rds
        output:
            pdf="results/de/plots/pca_batch_correction.pdf"
        conda:
            "../envs/r_de.yaml"
        shell:
            """
            Rscript workflow/scripts/plot_pca_batch.R \
                --input-rds {input.data} \
                --combat-seq {input.combat_seq} \
                --combat-ref {input.combat_ref} \
                --output {output.pdf}
            """
