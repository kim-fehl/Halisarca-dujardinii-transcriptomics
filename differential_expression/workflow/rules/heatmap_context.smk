import os
import re
from pathlib import Path


def _slugify(value):
    if value is None:
        return None
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return slug or None


def _expand_path(path):
    if not isinstance(path, str):
        return path
    p = Path(os.path.expanduser(path))
    if not p.is_absolute():
        p = Path(workflow.basedir).resolve().parent / p
    return str(p)


HEATMAP_GENOME_NAME = config.get("genome").get("name")

_raw_genesets = config.get("genesets") or []
if not _raw_genesets:
    raise ValueError("config['genesets'] must contain at least one entry")

HEATMAP_GENESETS = []
for entry in _raw_genesets:
    path = _expand_path(entry.get("path"))
    if not path:
        raise ValueError("Each heatmap geneset requires a 'path'")
    name = entry.get("name") or Path(path).stem
    slug = _slugify(name)
    if not slug:
        raise ValueError(f"Invalid heatmap geneset name derived from '{name}'")
    HEATMAP_GENESETS.append({
        "name": slug,
        "path": path,
        "pdf_four_seasons": f"results/de/plots/heatmap_{HEATMAP_GENOME_NAME}_{slug}_four_seasons.pdf",
        "xlsx_four_seasons": f"results/de/edgeR/heatmap_{HEATMAP_GENOME_NAME}_{slug}_four_seasons.xlsx",
        "pdf_general": f"results/de/plots/heatmap_{HEATMAP_GENOME_NAME}_{slug}_general.pdf",
        "xlsx_general": f"results/de/edgeR/heatmap_{HEATMAP_GENOME_NAME}_{slug}_general.xlsx",
    })

_names = [g["name"] for g in HEATMAP_GENESETS]
if len(set(_names)) != len(_names):
    raise ValueError("Heatmap geneset names must be unique")

HEATMAP_GENESET_MAP = {g["name"]: g for g in HEATMAP_GENESETS}
HEATMAP_GENESET_NAMES = sorted(HEATMAP_GENESET_MAP.keys())
HEATMAP_PDF_FOUR_SEASONS_TARGETS = [HEATMAP_GENESET_MAP[name]["pdf_four_seasons"] for name in HEATMAP_GENESET_NAMES]
HEATMAP_XLSX_FOUR_SEASONS_TARGETS = [HEATMAP_GENESET_MAP[name]["xlsx_four_seasons"] for name in HEATMAP_GENESET_NAMES]
HEATMAP_PDF_GENERAL_TARGETS = [HEATMAP_GENESET_MAP[name]["pdf_general"] for name in HEATMAP_GENESET_NAMES]
HEATMAP_XLSX_GENERAL_TARGETS = [HEATMAP_GENESET_MAP[name]["xlsx_general"] for name in HEATMAP_GENESET_NAMES]
