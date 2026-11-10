"""Shared helpers for the tick microbiome analysis scripts.

Every script under code/analysis/ imports from here so that config loading,
sample-id canonicalisation and taxonomy access behave identically everywhere.
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import yaml

# Repository root = two levels up from this file (code/analysis/common.py)
ROOT = Path(__file__).resolve().parents[2]

CONFIG_PATH = ROOT / "config" / "analysis_params.yaml"
ALIAS_PATH = ROOT / "metadata" / "sample_aliases.tsv"
CONTAMINANT_PATH = ROOT / "metadata" / "reagent_contaminants.tsv"
EXCLUDED_PATH = ROOT / "metadata" / "excluded_samples.tsv"
REASSIGN_PATH = ROOT / "metadata" / "contig_reassignments.tsv"


def load_config(path: Path | str = CONFIG_PATH) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["_root"] = str(ROOT)
    return cfg


def load_aliases(path: Path | str = ALIAS_PATH) -> dict[str, str]:
    """file_stem -> canonical sample_id. Comment lines start with '#'."""
    aliases: dict[str, str] = {}
    if not Path(path).exists():
        return aliases
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "file_stem":
                continue
            aliases[parts[0]] = parts[1]
    return aliases


def canonical(stem: str, aliases: dict[str, str]) -> str:
    """Map a filename stem to its canonical sample_id."""
    return aliases.get(stem, stem)


def strip_suffixes(name: str, suffixes: list[str]) -> str | None:
    """Return `name` with the first matching suffix removed, else None."""
    for suf in suffixes:
        if name.endswith(suf):
            return name[: -len(suf)]
    return None


def reverse_aliases(aliases: dict[str, str]) -> dict[str, list[str]]:
    """canonical sample_id -> the on-disk stems that map to it."""
    rev: dict[str, list[str]] = {}
    for stem, sid in aliases.items():
        rev.setdefault(sid, []).append(stem)
    return rev


def resolve_artefact(sample: str, directory: Path | str, suffix: str,
                     aliases: dict[str, str] | None = None) -> Path | None:
    """Locate <directory>/<stem><suffix> for a canonical sample_id.

    Artefact directories disagree about sample naming (see
    metadata/sample_aliases.tsv), so a file may be stored under an alias stem
    rather than the canonical id. Try the canonical name first, then any alias
    that resolves to it.
    """
    directory = Path(directory)
    p = directory / f"{sample}{suffix}"
    if p.exists():
        return p
    if aliases is None:
        aliases = load_aliases()
    for stem in reverse_aliases(aliases).get(sample, []):
        p = directory / f"{stem}{suffix}"
        if p.exists():
            return p
    return None


def load_excluded(path: Path | str = EXCLUDED_PATH) -> dict[str, str]:
    """sample_id -> reason, for samples removed from every analysis by decision
    (metadata/excluded_samples.tsv). Step 00 turns these into status
    `excluded_manual`; steps that enumerate contigs/ directly must skip them."""
    out: dict[str, str] = {}
    if not Path(path).exists():
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "sample_id":
                continue
            out[parts[0].strip()] = parts[1].strip() if len(parts) > 1 else ""
    return out


def load_contig_reassignments(path: Path | str = REASSIGN_PATH) -> dict[str, dict[str, int]]:
    """sample_id -> {contig_id: new_taxid}, from metadata/contig_reassignments.tsv.

    Manual curation of contigs whose nt LCA is wrong (e.g. tick repeats carried
    in viral reference records; see extracting_pathogens/*/are_these_real.html).
    Applied on top of lca/<s>.blast.lca by steps 03 and 04; the LCA files
    themselves are never edited."""
    out: dict[str, dict[str, int]] = {}
    if not Path(path).exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(r["sample_id"], {})[r["contig_id"]] = int(r["new_taxid"])
    return out


def load_contaminants(path: Path | str = CONTAMINANT_PATH) -> dict[str, str]:
    """genus -> confidence. Annotation only; never used to filter."""
    out: dict[str, str] = {}
    if not Path(path).exists():
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "genus":
                continue
            out[parts[0]] = parts[1]
    return out


# ---------------------------------------------------------------------------
# Read-file discovery
#
# Reads were repaired in a cascade (see code/fastqrepair.sh): QC -> repaired_QC
# -> fastqrepair -> reformat. The most-repaired copy that exists is the correct
# input for bowtie2 paired-end mapping; the earlier copies have pairing defects
# that make bowtie2 fail or silently mis-pair.
# ---------------------------------------------------------------------------
READ_DIR_PRIORITY = [
    ("repaired_QC/fastqrepair/reformat", "{s}.R1.fq.gz", "{s}.R2.fq.gz"),
    ("repaired_QC/fastqrepair", "{s}.R1.fq.gz", "{s}.R2.fq.gz"),
    ("repaired_QC", "{s}.R1.fq.gz", "{s}.R2.fq.gz"),
    ("QC", "{s}.clean.R1.fq.gz", "{s}.clean.R2.fq.gz"),
    ("QC", "{s}.R1.fq.gz", "{s}.R2.fq.gz"),
    ("reads_with_missing_contigs", "{s}_R1.qc.fq.gz", "{s}_R2.qc.fq.gz"),
    ("reads_with_missing_contigs", "{s}_1.clean.fq.gz", "{s}_2.clean.fq.gz"),
]


def find_reads(sample: str, root: Path = ROOT) -> tuple[str | None, str | None, str | None]:
    """Return (r1, r2, source_dir) for the most-repaired read pair available."""
    for subdir, r1pat, r2pat in READ_DIR_PRIORITY:
        r1 = root / subdir / r1pat.format(s=sample)
        r2 = root / subdir / r2pat.format(s=sample)
        if r1.exists() and r2.exists():
            return str(r1), str(r2), subdir
    return None, None, None


def write_tsv(path: Path | str, rows: list[dict], columns: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        # lineterminator is explicit: csv defaults to CRLF, which breaks every
        # unix consumer of these files (cut/awk pick up a trailing \r).
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.relative_to(ROOT)}  ({len(rows)} rows)", file=sys.stderr)
