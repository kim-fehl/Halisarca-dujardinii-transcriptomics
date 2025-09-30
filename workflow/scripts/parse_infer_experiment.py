#!/usr/bin/env python3
from pathlib import Path
import sys

report_path = Path(snakemake.input["report"])
strand_path = Path(snakemake.output["strand"])

text = report_path.read_text(encoding="utf-8", errors="ignore")

summary_line = None
for line in text.splitlines():
    if line.startswith("This is"):
        continue
    if line.startswith(">>"):
        summary_line = line
        break

same = None
opposite = None
for line in text.splitlines():
    if '"++,--"' in line:
        same = float(line.rsplit(':', 1)[1])
    if '"+-,-+"' in line:
        opposite = float(line.rsplit(':', 1)[1])

strand = "UNSTRANDED"
if summary_line:
    low = summary_line.lower()
    if "first" in low or "read1:forward" in low:
        strand = "FR_SECONDSTRAND"
    elif "second" in low or "read1:reverse" in low:
        strand = "FR_FIRSTSTRAND"

if strand == "UNSTRANDED" and same is not None and opposite is not None:
    if same >= 0.8 and same - opposite >= 0.2:
        strand = "FR_SECONDSTRAND"
    elif opposite >= 0.8 and opposite - same >= 0.2:
        strand = "FR_FIRSTSTRAND"

strand_path.parent.mkdir(parents=True, exist_ok=True)
strand_path.write_text(strand + "\n", encoding="ascii")
