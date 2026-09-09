"""Shared helpers for the marker-gene workflow (code/analysis/11..16).

Deliberately dependency-free beyond the stdlib and `common.py`: the marker
scripts must run in the base environment as well as in `env/tick_markers.yml`,
and biopython is in neither.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402

# The field separator inside every sequence name. Chosen because it cannot
# appear in a sample id, a tick species name or a SPAdes contig id, so a name
# can always be split back into its parts without ambiguity.
SEP = "|"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                            "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


# --------------------------------------------------------------------------
# sequence basics
# --------------------------------------------------------------------------
def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(seq: str, frame: int = 0) -> str:
    """Translate `seq` from `frame` (0-based). Unknown codons become X."""
    s = seq[frame:]
    s = s[: len(s) - len(s) % 3].upper()
    return "".join(CODON_TABLE.get(s[i:i + 3], "X") for i in range(0, len(s), 3))


def internal_stops(aa: str) -> int:
    """Stop codons that are not the final residue -- the frameshift signal."""
    return aa.rstrip("*").count("*")


def best_frame(seq: str) -> tuple[int, str]:
    """Pick the forward frame with the fewest internal stops.

    Used instead of trusting tblastn's reported `sframe`, because a locus is
    assembled from merged HSPs whose span may start off-frame, and because the
    frame of a minus-strand HSP is reported relative to the whole contig, not
    to the extracted and reverse-complemented span.

    Ties break on the longest uninterrupted stop-free run, which is what
    separates a real reading frame from a shifted one on a short locus.
    """
    best = (10**9, -1, 0, "")
    for f in (0, 1, 2):
        aa = translate(seq, f)
        if not aa:
            continue
        run = max((len(x) for x in aa.split("*")), default=0)
        key = (internal_stops(aa), -run)
        if key < best[:2]:
            best = (key[0], key[1], f, aa)
    return (best[2], best[3]) if best[2] >= 0 else (0, translate(seq, 0))


# --------------------------------------------------------------------------
# FASTA i/o
# --------------------------------------------------------------------------
def read_fasta(path: Path | str) -> dict[str, str]:
    """Name -> sequence, preserving insertion order. Name is the full header."""
    out: dict[str, str] = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def write_fasta(path: Path | str, seqs: dict[str, str], width: int = 70) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


# --------------------------------------------------------------------------
# naming
# --------------------------------------------------------------------------
def slug(text: str, keep_dot: bool = False) -> str:
    """Filesystem- and phylogenetics-software-safe token."""
    pat = r"[^A-Za-z0-9._-]+" if keep_dot else r"[^A-Za-z0-9_-]+"
    return re.sub(pat, "_", str(text)).strip("_") or "NA"


def contig_name(sample_id: str, tick_species: str, target: str, contig_id: str) -> str:
    """Header for a target contig: sample, tick species, target, contig."""
    return SEP.join([sample_id, slug(tick_species), slug(target), contig_id])


def marker_name(sample_id: str, tick_species: str, target: str, marker: str,
                contig_id: str, length: int, copy: int = 1) -> str:
    """Header for an extracted marker sequence.

    Every element the downstream tree needs to be readable is in the name:
    sample id, tick species, target taxon, marker gene, source contig, length,
    and a copy index for the case where one sample yields the locus on more
    than one contig.
    """
    parts = [sample_id, slug(tick_species), slug(target), slug(marker), contig_id,
             f"{length}bp"]
    if copy > 1:
        parts.append(f"copy{copy}")
    return SEP.join(parts)


# --------------------------------------------------------------------------
# misc
# --------------------------------------------------------------------------
def cfg_markers() -> dict:
    return load_config()["markers"]


def out_root() -> Path:
    return ROOT / cfg_markers()["out_dir"]


def require(tool: str) -> str:
    """Absolute path to `tool`, or exit with an actionable message."""
    from shutil import which
    p = which(tool)
    if p is None:
        sys.exit(f"ERROR: `{tool}` not on PATH.\n"
                 f"  conda activate tick_markers   "
                 f"(mamba env create -f env/tick_markers.yml)")
    return p


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"ERROR: command failed ({r.returncode})\n  "
                 + " ".join(cmd) + "\n" + (r.stderr or "")[-3000:])
    return r


def load_panel(path: Path | str | None = None) -> list[dict]:
    """Parse metadata/target_markers.tsv into a list of row dicts."""
    path = Path(path or ROOT / cfg_markers()["panel"])
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            r = dict(zip(header, parts))
            for k in ("n_refs", "min_aln_len"):
                r[k] = int(r[k])
            r["min_pident"] = float(r["min_pident"])
            r["queries"] = [q.strip() for q in r["query"].split("|||") if q.strip()]
            r["key"] = f"{slug(r['target'])}__{slug(r['marker'])}"
            rows.append(r)
    return rows


def tick_species_map() -> dict[str, tuple[str, str]]:
    """sample_id -> (tick_species, source).

    Prefers the curated metadata/sample_metadata.tsv when it exists. Falls back
    to the id-string guess in the template, which is explicitly a GUESS -- the
    source string travels into the report so a name in a tree can always be
    traced back to how the species was assigned.
    """
    import csv as _csv
    for fname, col, default_src in (
        ("sample_metadata.tsv", "tick_species", "curated_metadata"),
        ("sample_metadata_TEMPLATE.tsv", "tick_species_guess", "id_string_guess"),
    ):
        p = ROOT / "metadata" / fname
        if not p.exists():
            continue
        out: dict[str, tuple[str, str]] = {}
        with open(p) as fh:
            for row in _csv.DictReader(fh, delimiter="\t"):
                sid = (row.get("sample_id") or "").strip()
                if not sid:
                    continue
                val = (row.get(col) or "").strip()
                src = (row.get("tick_species_source") or "").strip() or default_src
                if not val:  # curated cell still blank -> fall back to the guess
                    val = (row.get("tick_species_guess") or "").strip()
                    src = "id_string_guess"
                out[sid] = (val or "unknown_species", src or default_src)
        if out:
            return out
    return {}
