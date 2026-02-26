import sys
from pathlib import Path


sys.path.insert(0, str(Path(workflow.basedir) / "lib"))
from sample_inputs import load_sample_context


_sample_ctx = load_sample_context(
    config=config,
    cores=workflow.cores,
    workflow_basedir=workflow.basedir,
)

metadata_df = _sample_ctx["metadata_df"]
RUN_TO_FASTQ_R1 = _sample_ctx["RUN_TO_FASTQ_R1"]
RUN_TO_FASTQ_R2 = _sample_ctx["RUN_TO_FASTQ_R2"]
RUN_TO_SAMPLE = _sample_ctx["RUN_TO_SAMPLE"]
RUN_IDS = _sample_ctx["RUN_IDS"]
PRIMARY_RUN = _sample_ctx["PRIMARY_RUN"]
IS_PAIRED_END = _sample_ctx["IS_PAIRED_END"]
MAX_THREADS = _sample_ctx["MAX_THREADS"]
AUX_THREADS = _sample_ctx["AUX_THREADS"]
