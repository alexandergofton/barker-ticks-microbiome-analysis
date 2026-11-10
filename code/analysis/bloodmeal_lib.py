"""Shared helpers for the blood-meal workflow (steps 30-35, run_bloodmeal.sh)."""

from __future__ import annotations

import csv
import re
import sys
from array import array
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config  # noqa: E402

RANKS = ["class", "order", "family", "genus", "species"]
VERT_TAXID = 7742  # Vertebrata

_RANK_CODES = {r: i + 1 for i, r in enumerate(
    ["superkingdom", "kingdom", "phylum", "subphylum", "class", "order", "family",
     "genus", "species", "subspecies"])}


def cfg_bloodmeal() -> dict:
    return load_config().get("bloodmeal", {}) or {}


def bloodmeal_root() -> Path:
    return ROOT / cfg_bloodmeal().get("out_dir", "results/bloodmeal")


def read_tsv(path: Path | str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


# ---------------------------------------------------------------------------
# Taxonomy. nodes.dmp has ~2.7M taxa and this node has ~1 GB of free RAM, so
# parents/ranks live in flat arrays indexed by taxid, and scientific names are
# kept only for Vertebrata (every name this workflow ever prints).
# ---------------------------------------------------------------------------
class Taxonomy:
    def __init__(self, taxdump: Path = ROOT / "taxdump"):
        # Two passes over nodes.dmp (size, then fill) so the 2.7M rows are
        # never held as Python objects.
        with open(taxdump / "nodes.dmp") as fh:
            max_id = max(int(line.split("\t|\t", 1)[0]) for line in fh)
        self.parent = array("i", bytes(4 * (max_id + 1)))
        self.rank = bytearray(max_id + 1)
        with open(taxdump / "nodes.dmp") as fh:
            for line in fh:
                p = line.split("\t|\t", 3)
                t = int(p[0])
                self.parent[t], self.rank[t] = int(p[1]), _RANK_CODES.get(p[2], 0)
        # 0 = not yet known, 1 = vertebrate, 2 = not; filled lazily by is_vertebrate
        self._vert = bytearray(max_id + 1)
        self.merged: dict[int, int] = {}
        with open(taxdump / "merged.dmp") as fh:
            for line in fh:
                p = line.split("\t|\t")
                self.merged[int(p[0])] = int(p[1].split("\t")[0])
        # names for Vertebrata plus its own ancestors (Craniata, Chordata, ...),
        # which an LCA across vertebrate + other hits can land on
        above = set(self.path(VERT_TAXID))
        self.names: dict[int, str] = {}
        with open(taxdump / "names.dmp") as fh:
            for line in fh:
                if "scientific name" not in line:
                    continue
                p = line.split("\t|\t")
                t = int(p[0])
                if t in above or self.is_vertebrate(t):
                    self.names[t] = p[1]

    def resolve(self, taxid) -> int:
        try:
            t = int(taxid)
        except (TypeError, ValueError):
            return 0
        t = self.merged.get(t, t)
        return t if 0 < t < len(self.parent) and self.parent[t] else 0

    def path(self, taxid) -> list[int]:
        """Root-ward list of taxids, starting at `taxid`."""
        t = self.resolve(taxid)
        out = []
        while t and t not in out:
            out.append(t)
            if t == 1:
                break
            t = self.parent[t]
        return out

    def is_vertebrate(self, taxid) -> bool:
        t = self.resolve(taxid)
        walked = []
        while t and not self._vert[t]:
            if t == VERT_TAXID:
                self._vert[t] = 1
                break
            if t == 1 or len(walked) > 100:
                self._vert[t] = 2
                break
            walked.append(t)
            t = self.parent[t]
        state = self._vert[t] if t else 2
        for w in walked:
            self._vert[w] = state
        return state == 1

    def lineage(self, taxid) -> dict[str, str]:
        """rank -> name for RANKS, plus taxid/name/rank of the node itself."""
        p = self.path(taxid)
        out = {r: "" for r in RANKS}
        for t in p:
            for r in RANKS:
                if self.rank[t] == _RANK_CODES[r] and not out[r]:
                    out[r] = self.names.get(t, str(t))
        # NCBI gives turtles and crocodiles no class rank. Without one, a
        # Testudinidae contig cannot count against a python at class level, so
        # the order name stands in as the class.
        if not out["class"] and out["order"]:
            out["class"] = out["order"]
        t0 = p[0] if p else 0
        out["taxid"] = str(t0)
        out["name"] = self.names.get(t0, str(t0))
        rc = self.rank[t0] if t0 else 0
        out["rank"] = next((r for r, c in _RANK_CODES.items() if c == rc), "no rank")
        return out

    def lca(self, taxids) -> int:
        paths = [self.path(t) for t in taxids]
        paths = [p for p in paths if p]
        if not paths:
            return 0
        common = set(paths[0])
        for p in paths[1:]:
            common &= set(p)
        return next((t for t in paths[0] if t in common), 1)

    def ancestor_at(self, taxid, rank: str) -> int:
        """The taxid at `rank` on `taxid`'s lineage, or 0 if it has none."""
        code = _RANK_CODES[rank]
        return next((t for t in self.path(taxid) if self.rank[t] == code), 0)

    def cap(self, taxid, max_rank: str) -> int:
        """Deepest node on `taxid`'s lineage no deeper than `max_rank`.

        A species-level LCA whose best identity only supports family is
        returned as that family; a node already above `max_rank` is returned
        unchanged."""
        limit = RANKS.index(max_rank) if max_rank in RANKS else len(RANKS) - 1
        for t in self.path(taxid):
            code = self.rank[t]
            r = next((r for r, c in _RANK_CODES.items() if c == code), None)
            if r in RANKS and RANKS.index(r) <= limit:
                return t
            if r in ("phylum", "subphylum", "kingdom", "superkingdom"):
                return t
        return self.resolve(taxid)


# ---------------------------------------------------------------------------
# What kind of locus an nt subject is. Step 30's subject_type() is a broad
# superset (any title mentioning "mitochondri"), which lets through nuclear
# genes for mitochondrial proteins -- "translocase of inner mitochondrial
# membrane 17A", "lon peptidase 1, mitochondrial" -- in every sample. Those are
# nuclear and must never be read as mtDNA.
# ---------------------------------------------------------------------------
_NUCLEAR_TX = re.compile(r"^(XM|NM|XR|NR)_")
_MT_GENOME = re.compile(r"mitochondrion|mitochondrial (dna|genome)|; mitochondrial$|"
                        r"mitochondrial$|mitochondrial, (complete|partial)", re.I)
_MT_GENE = {
    "COX1": r"cytochrome c oxidase subunit (1|i)\b|cytochrome oxidase subunit (1|i)\b|\b(coi|cox1|co1)\b",
    "COX2": r"cytochrome c oxidase subunit (2|ii)\b|\b(coii|cox2)\b",
    "COX3": r"cytochrome c oxidase subunit (3|iii)\b|\b(coiii|cox3)\b",
    "CYTB": r"cytochrome b(?![\w-])|\b(cytb|cob)\b",
    "12S": r"12s ribosomal|\b12s rrna\b|small subunit ribosomal rna.*mitochondrial",
    "16S": r"16s ribosomal|\b16s rrna\b",
    "ND": r"nadh dehydrogenase subunit|\bnd[1-6]l?\b",
    "ATP": r"atp(ase)? ?(synthase )?(f0 )?(subunit )?[68]\b|\batp[68]\b",
    "CR": r"d-loop|control region",
}
_NUC_RRNA = {
    "18S": r"18s",
    "28S": r"28s",
    "5.8S": r"5\.8s",
    "ITS": r"internal transcribed|\bits[12]?\b",
}


def locus_class(acc: str, title: str) -> tuple[str, str]:
    """(kind, gene) for an nt subject: kind in mtDNA | nuc_rRNA | other.

    gene is a coarse label (COX1, CYTB, 12S, ...; 'genome' for complete or
    partial mitogenomes; 18S/28S/... for nuclear rRNA)."""
    t = title.lower()
    if _NUCLEAR_TX.match(acc) or " mrna" in t or "-like" in t or "predicted" in t:
        return "other", ""
    mt_word = bool(_MT_GENOME.search(title))
    gene = next((g for g, p in _MT_GENE.items() if re.search(p, t)), "")
    if mt_word or gene:
        if re.search(r"mitochondrion,? (complete|partial) genome|complete mitochondrial genome|"
                     r"mitochondrial (dna|genome),? complete", t):
            return "mtDNA", "genome"
        if mt_word or gene in ("COX1", "COX2", "COX3", "CYTB", "12S", "16S", "CR", "ND", "ATP"):
            return "mtDNA", gene or "mito"
    if re.search(r"ribosomal|rrna|rdna|internal transcribed", t):
        g = next((g for g, p in _NUC_RRNA.items() if re.search(p, t)), "")
        if g:
            return "nuc_rRNA", g
    return "other", ""


def identity_rank(pident: float, caps: dict[str, float]) -> str:
    """Deepest rank supported by `pident` (caps: rank -> minimum identity)."""
    for r in reversed(RANKS):
        if r in caps and pident >= float(caps[r]):
            return r
    return "class"


def contig_cov(contig: str) -> float:
    m = re.search(r"_cov_([0-9.]+)$", contig or "")
    return float(m.group(1)) if m else 0.0


def contig_len(contig: str) -> int:
    m = re.search(r"_length_(\d+)_", contig or "")
    return int(m.group(1)) if m else 0
