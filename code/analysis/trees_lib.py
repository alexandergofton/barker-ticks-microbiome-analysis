"""Shared helpers for the marker-gene phylogeny workflow (code/analysis/17..19).

Sits on top of `markers_lib`, which owns the sequence, FASTA, naming and Entrez
primitives shared with steps 11-16. Everything here is specific to turning a
recovered locus into a tree: the two reference panels, reference naming, and
the HSP-to-locus assembly used to cut a marker out of a downloaded record.

Same dependency rule as markers_lib -- stdlib and `common.py` only.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402
from markers_lib import SEP, load_panel, out_root, slug  # noqa: E402

# Tip-name tags. A reference tip is
#     REF|<accession>|<organism>|<target>|<marker>|<len>bp
# and a sample tip (written by step 14) is
#     <sample_id>|<tick_species>|<target>|<marker>|<contig_id>|<len>bp
# so both put the target at field 2 and the marker at field 3, and field 0
# alone says which kind of tip it is.
REF_TAG = "REF"
OUTGROUP_TAG = "OUTGROUP"

BLAST_FMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
             "evalue bitscore qlen slen")
BLAST_COLS = BLAST_FMT.split()[1:]


def cfg_trees() -> dict:
    return load_config()["trees"]


def trees_root() -> Path:
    return ROOT / cfg_trees()["out_dir"]


def is_reference(name: str) -> bool:
    return name.split(SEP, 1)[0] in (REF_TAG, OUTGROUP_TAG)


def is_outgroup(name: str) -> bool:
    return name.split(SEP, 1)[0] == OUTGROUP_TAG


def ref_name(tag: str, accession: str, organism: str, target: str, marker: str,
             length: int) -> str:
    return SEP.join([tag, slug(accession, keep_dot=True), slug(organism),
                     slug(target), slug(marker), f"{length}bp"])


# ---------------------------------------------------------------------------
# the two reference panels
# ---------------------------------------------------------------------------
def _read_table(path: Path) -> list[dict]:
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def _queries(text: str) -> list[str]:
    return [q.strip() for q in (text or "").split("|||") if q.strip()]


def load_ref_panel(path: Path | str | None = None) -> dict[tuple[str, str], dict]:
    """(target, marker) -> gene-level reference query and its outgroup."""
    path = Path(path or ROOT / cfg_trees()["refs_panel"])
    out = {}
    for r in _read_table(path):
        out[(r["target"], r["marker"])] = {
            "n_refs": int(r["n_refs"]),
            "queries": _queries(r["query"]),
            "outgroup_n": int(r.get("outgroup_n") or 0),
            "outgroup_queries": _queries(r.get("outgroup_query", "")),
        }
    return out


def load_genome_panel(path: Path | str | None = None) -> dict[str, dict]:
    """target -> whole-genome reference query."""
    path = Path(path or ROOT / cfg_trees()["genomes_panel"])
    return {r["target"]: {"n_genomes": int(r["n_genomes"]),
                          "queries": _queries(r["query"])}
            for r in _read_table(path)}


def auto_query(target: str, marker: str, marker_type: str) -> list[str]:
    """Fallback query for a locus with no row in metadata/tree_references.tsv.

    Deliberately crude -- it exists so a newly added marker still gets a tree
    instead of silently getting none. Every locus built this way is marked
    `auto` in the manifest; write a real row when you care about the result.
    """
    if marker_type == "rrna":
        return [f'{target}[Organism] AND "{marker}"[Title] AND 300:12000[SLEN]',
                f"{target}[Organism] AND {marker}[All Fields] AND 300:12000[SLEN]"]
    return [f'{target}[Organism] AND ({marker}[Gene] OR "{marker}"[Title]) AND 300:9000[SLEN]',
            f"{target}[Organism] AND {marker}[All Fields] AND 300:9000[SLEN]"]


# ---------------------------------------------------------------------------
# which loci have something to build a tree from
# ---------------------------------------------------------------------------
def sample_fasta(target: str, marker: str, ext: str = "fna", flagged: bool = False) -> Path:
    stem = f"{slug(target)}__{slug(marker)}{'.flagged' if flagged else ''}"
    return out_root() / "sequences" / slug(target) / f"{stem}.{ext}"


def loci_with_sequences(targets: list[str] | None = None,
                        markers: list[str] | None = None) -> list[dict]:
    """Panel rows whose step-14 nucleotide FASTA exists and is non-empty.

    A locus with no recovered sample sequence has nothing to place in a tree,
    so it is not carried through steps 17-19 at all -- `marker_gaps.tsv` from
    step 16 already explains why it is absent, and inventing a reference-only
    tree for it would say nothing about these ticks.
    """
    rows = []
    for row in load_panel():
        if targets and row["target"] not in set(targets):
            continue
        if markers and row["marker"] not in set(markers):
            continue
        fna = sample_fasta(row["target"], row["marker"])
        if fna.exists() and fna.stat().st_size > 0:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# organism names out of GenBank deflines / esummary titles
# ---------------------------------------------------------------------------
_EPITHET = re.compile(r"^[a-z][a-z-]*\.?$")


def organism_from_title(title: str) -> str:
    """Best-effort binomial from a GenBank title.

    Titles are free text, so this is a heuristic, but it only has to be good
    enough to (a) label a tip and (b) group records of the same organism for
    the per-organism cap. It handles the four shapes that actually occur:

        Rickettsia felis URRWXCal2, complete genome      -> Rickettsia felis
        Candidatus Midichloria mitochondrii IricVA ...   -> Candidatus Midichloria mitochondrii
        Coxiella endosymbiont of Amblyomma americanum .. -> Coxiella endosymbiont of Amblyomma americanum
        Rickettsia sp. strain AUS-1 gltA gene, partial   -> Rickettsia sp. AUS-1
    """
    text = (title or "").split(",")[0].strip()
    # A defline reaching here may still carry its accession; drop a leading
    # token that looks like one rather than reading it as a genus.
    words = text.split()
    if words and re.match(r"^[A-Z]{1,3}_?[A-Z0-9]+\.\d+$", words[0]):
        words = words[1:]
    if not words:
        return "unknown_organism"

    # "Uncultured Candidatus Midichloria sp. clone .." is a name for the record,
    # not for a genus; drop the qualifier so those records group with their
    # genus rather than all becoming one organism called "Uncultured".
    while words and words[0].lower() in ("uncultured", "unidentified", "unclassified"):
        words = words[1:]
    if not words:
        return "unknown_organism"

    # Viruses are not binomials. "Jingmen tick virus isolate GR_2024 segment 4"
    # must not become the genus "Jingmen tick": take the name up to and
    # including "virus", then the isolate, which is what actually distinguishes
    # two records of the same virus -- and the per-organism cap needs them
    # distinguished, or one virus contributes two sequences to a tree of
    # sixty available isolates.
    low = [w.lower().strip(",;()") for w in words[:6]]
    if "virus" in low:
        k = low.index("virus")
        name = words[:k + 1]
        j = k + 1
        while j < len(words) and words[j].lower() in ("isolate", "strain", "clone", "str."):
            j += 1
        if j < len(words) and j > k + 1 and re.match(r"^[A-Za-z0-9][A-Za-z0-9._/-]*$", words[j]):
            name.append(words[j])
        return " ".join(name)

    out: list[str] = []
    i = 0
    if words[0].lower().rstrip(".") in ("candidatus", "ca"):
        out.append("Candidatus")
        i = 1
    if i >= len(words):
        return "unknown_organism"
    out.append(words[i])          # genus
    i += 1

    if i < len(words) and words[i] == "endosymbiont":
        # "<Genus> endosymbiont of <Host genus> <host species>" is a real and
        # important organism name here -- the Coxiella-like and
        # Francisella-like tick endosymbionts are named this way and collapsing
        # them all to "Coxiella endosymbiont" would merge distinct lineages.
        out.extend(words[i:i + 4])
        return " ".join(out)

    if i < len(words) and _EPITHET.match(words[i]):
        out.append(words[i])
        i += 1
        # "sp." / "cf." on their own do not identify anything; keep the strain
        # token after them so two different isolates stay apart.
        if out[-1] in ("sp.", "cf.", "aff.", "sp"):
            # "sp. strain X" / "sp. isolate X" / "sp. clone X" all name the
            # isolate in the token AFTER the collection word.
            while i < len(words) and words[i].lower() in ("strain", "isolate",
                                                          "clone", "str."):
                i += 1
            if i < len(words) and re.match(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", words[i]):
                out.append(words[i])
    return " ".join(out)


# ---------------------------------------------------------------------------
# HSPs -> one locus per reference record
# ---------------------------------------------------------------------------
def parse_blast(path: Path, min_len: int, min_pid: float) -> list[dict]:
    """Tabular BLAST rows clearing the panel floors for this marker."""
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(BLAST_COLS):
                continue
            h = dict(zip(BLAST_COLS, p))
            h["pident"] = float(h["pident"]); h["length"] = int(h["length"])
            h["bitscore"] = float(h["bitscore"])
            h["sstart"] = int(h["sstart"]); h["send"] = int(h["send"])
            h["qlen"] = int(h["qlen"]); h["slen"] = int(h["slen"])
            if h["length"] < min_len or h["pident"] < min_pid:
                continue
            out.append(h)
    return out


def best_locus_per_subject(hits: list[dict], merge_gap: int, unit: int) -> dict[str, dict]:
    """Merge HSPs into candidate loci and keep the best one per record.

    Same merge rule as step 13 -- HSPs of one strand closer than `merge_gap`
    are one locus -- but no cross-marker paralogue resolution, because only one
    marker's queries are ever blasted against a given reference set here.

    Only the best-scoring locus per record survives: a genome carries one copy
    of groEL and twenty of IS1111, and taking the best is the answer in both
    cases -- an rRNA operon repeat or a transposase copy is not extra
    phylogenetic information, it is the same locus again.
    """
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for h in hits:
        groups[(h["sseqid"], "+" if h["sstart"] <= h["send"] else "-")].append(h)

    best: dict[str, dict] = {}
    for (subject, strand), hs in groups.items():
        spans = sorted(((min(h["sstart"], h["send"]), max(h["sstart"], h["send"]), h)
                        for h in hs), key=lambda s: (s[0], s[1]))
        cur, merged = [spans[0]], []
        for s in spans[1:]:
            if s[0] - cur[-1][1] <= merge_gap:
                cur.append(s)
            else:
                merged.append(cur); cur = [s]
        merged.append(cur)

        for grp in merged:
            lo, hi = min(g[0] for g in grp), max(g[1] for g in grp)
            hsps = [g[2] for g in grp]
            top = max(hsps, key=lambda h: h["bitscore"])
            locus = {
                "subject": subject, "strand": strand, "start": lo, "end": hi,
                "len": hi - lo + 1, "n_hsps": len(hsps),
                "sum_bitscore": round(sum(h["bitscore"] for h in hsps), 1),
                "max_bitscore": round(top["bitscore"], 1),
                "best_pident": round(max(h["pident"] for h in hsps), 1),
                "best_query": top["qseqid"], "subject_len": top["slen"],
                "ref_cov_pct": round(min(100.0 * (hi - lo + 1) /
                                         max(1, top["qlen"] * unit), 999.9), 1),
            }
            prev = best.get(subject)
            if prev is None or locus["sum_bitscore"] > prev["sum_bitscore"]:
                best[subject] = locus
    return best
