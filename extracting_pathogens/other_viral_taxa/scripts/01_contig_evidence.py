#!/usr/bin/env python3
"""Per-contig nucleotide evidence for every viral-LCA contig (other_viral_taxa).

Reads the extraction (lca_assignments.tsv, blast/*.viral.blast.tsv) and the
local taxdump, and writes into are_these_real/:

  viral_taxa.tsv      one row per viral LCA taxon: lineage names, virus group
                      (bacteriophage / eukaryotic virus / unresolved), counts
  allhits.tsv         every BLAST hit, sample-prefixed, with the subject's
                      organism group (virus, tick, arthropod, bacteria, ...)
  contig_evidence.tsv one row per contig: length, k-mer coverage, bp covered by
                      viral vs non-viral hits, best scores, co-hit groups
  all_viral.fa        every contig, id = "<sample>|<contig>"

Remember: nt BLAST kept only 10 subjects per contig, so a MISSING non-viral
hit is weak evidence, a PRESENT one is informative.
"""

from __future__ import annotations

import collections
import csv
import glob
import os
import re
import sys
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
ROOT = OV.parents[1]
W = OV / "are_these_real"
sys.path.insert(0, str(ROOT / "code" / "analysis"))
from bloodmeal_lib import Taxonomy  # noqa: E402
from common import canonical, load_aliases  # noqa: E402

# Organism groups for BLAST subjects, tested in this order on the lineage.
GROUPS = [(10239, "virus"), (6935, "tick"), (6656, "other arthropod"),
          (7742, "vertebrate"), (33208, "other animal"), (2, "bacteria"),
          (2157, "archaea"), (2759, "other eukaryote")]
# A viral lineage passing through any of these is a bacteriophage / archaeal virus.
PHAGE_CLADES = {"Caudoviricetes", "Faserviricetes", "Malgrandaviricetes", "Leviviricetes",
                "Tectiliviricetes", "Vinavirales", "Plasmaviridae", "Cystoviridae",
                "Microviridae", "Inoviridae", "Tectiviridae", "Corticoviridae",
                "unclassified bacterial viruses", "Crassvirales"}
LIN_RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus"]


def merge(iv):
    """Merged intervals. Coverage is kept as interval lists, not position
    sets: 42,789 contigs x 2.0M hits as Python int sets would need ~2 GB."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def span(iv):
    return sum(e - s + 1 for s, e in merge(iv))


def minus(a, b):
    """bp of intervals `a` not covered by intervals `b`."""
    a, b = merge(a), merge(b)
    return span(a) - sum(max(0, min(e, be) - max(s, bs) + 1) for s, e in a for bs, be in b)


def main() -> None:
    tx = Taxonomy(ROOT / "taxdump")
    group_cache: dict[int, str] = {}

    def group(taxid: str) -> str:
        t = tx.resolve(taxid.split(";")[0])
        if t not in group_cache:
            p = set(tx.path(t))
            group_cache[t] = next((g for m, g in GROUPS if m in p), "unclassified")
        return group_cache[t]

    # ---- viral taxa: lineage names -----------------------------------------
    vt = {r["taxid"]: r for r in csv.DictReader(open(OV / "viral_taxids.tsv"), delimiter="\t",
                                                 fieldnames=["taxid", "name", "rank", "family", "genus"])}
    need = {}
    for t in vt:
        need[t] = tx.path(t)
    want = {x for p in need.values() for x in p}
    names, ranks = {}, {}
    with open(ROOT / "taxdump" / "names.dmp") as fh:
        for line in fh:
            if "scientific name" in line:
                p = line.split("\t|\t")
                if int(p[0]) in want:
                    names[int(p[0])] = p[1]
    with open(ROOT / "taxdump" / "nodes.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t", 3)
            if int(p[0]) in want:
                ranks[int(p[0])] = p[2]

    lca = list(csv.DictReader(open(OV / "lca_assignments.tsv"), delimiter="\t"))
    aliases = load_aliases()
    md = {r["sample_id"]: r for r in csv.DictReader(open(ROOT / "metadata" / "sample_metadata.tsv"),
                                                    delimiter="\t")}
    cnt, smp = collections.Counter(), collections.defaultdict(set)
    for r in lca:
        cnt[r["lca_taxid"]] += 1
        smp[r["lca_taxid"]].add(r["sample"])
    vrow = {}
    for t in vt:
        lin = {k: "" for k in LIN_RANKS}
        for x in need[t]:
            rk = ranks.get(x, "")
            if rk in lin and not lin[rk]:
                lin[rk] = names.get(x, "")
        clades = {names.get(x, "") for x in need[t]}
        nm = names.get(tx.resolve(t), vt[t]["name"])
        if clades & PHAGE_CLADES or re.search(r"phage|crass", nm, re.I):
            vg = "bacteriophage"
        elif not lin["family"] and not lin["genus"] and vt[t]["rank"] not in ("species", "no rank", "serotype", "clade") \
                or nm in ("Viruses", "Riboviria", "unclassified viruses"):
            vg = "unresolved"
        else:
            vg = "eukaryotic virus"
        vrow[t] = {"taxid": t, "name": nm, "rank": vt[t]["rank"], **lin, "virus_group": vg,
                   "n_contigs": cnt[t], "n_samples": len(smp[t])}
    with open(W / "viral_taxa.tsv", "w") as o:
        cols = ["taxid", "name", "rank", *LIN_RANKS, "virus_group", "n_contigs", "n_samples"]
        o.write("\t".join(cols) + "\n")
        for r in sorted(vrow.values(), key=lambda r: -r["n_contigs"]):
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- hits ---------------------------------------------------------------
    q = collections.defaultdict(lambda: {"v": [], "o": [], "g": collections.Counter(),
                                         "gbp": collections.defaultdict(list),
                                         "bv": (0.0, 0.0, "", 0, 0), "bo": (0.0, "")})
    with open(W / "allhits.tsv", "w") as out:
        out.write("sample\tcontig\tsaccver\tssciname\tstaxid\tsgroup\tpident\tlength\t"
                  "qstart\tqend\tsstart\tsend\tevalue\tbitscore\n")
        for f in sorted(glob.glob(str(OV / "blast" / "*.viral.blast.tsv"))):
            s = os.path.basename(f)[: -len(".viral.blast.tsv")]
            with open(f) as fh:
                next(fh)
                for ln in fh:
                    a = ln.rstrip("\n").split("\t")
                    if len(a) < 17:
                        continue
                    g = group(a[6])
                    qs, qe = sorted((int(a[11]), int(a[12])))
                    bits = float(a[16])
                    d = q[(s, a[0])]
                    if g == "virus":
                        d["v"].append((qs, qe))
                        if bits > d["bv"][0]:
                            d["bv"] = (bits, float(a[7]), a[1], int(a[13]), int(a[14]))
                    else:
                        d["o"].append((qs, qe))
                        d["gbp"][g].append((qs, qe))
                        d["g"][g] += 1
                        if bits > d["bo"][0]:
                            d["bo"] = (bits, g)
                    out.write("\t".join([s, a[0], a[1], a[4], a[6], g, a[7], a[8], str(qs), str(qe),
                                         a[13], a[14], a[15], a[16]]) + "\n")

    # ---- contigs ------------------------------------------------------------
    cols = ["sample", "sample_id", "status", "tick_genus", "tick_species", "contig", "len", "kmer_cov",
            "lca_taxid", "lca_name", "virus_group", "viral_family", "n_hits", "viral_bp", "viral_frac",
            "best_viral_bits", "best_viral_pident", "top_viral_subject", "top_viral_sstart",
            "top_viral_send", "nonviral_bp", "nonviral_only_bp", "best_nonviral_bits",
            "best_nonviral_group", "tick_bp", "arthropod_bp", "bacteria_bp", "vertebrate_bp",
            "nonviral_groups"]
    with open(W / "contig_evidence.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in lca:
            s, c = r["sample"], r["query"]
            d = q.get((s, c))
            L = int(re.search(r"length_(\d+)", c).group(1))
            cov = float(re.search(r"cov_([\d.]+)", c).group(1))
            sid = canonical(s, aliases)
            m = md.get(sid, {})
            if d is None:
                d = {"v": [], "o": [], "g": {}, "gbp": {}, "bv": (0, 0, "", 0, 0), "bo": (0, "")}
            vb = span(d["v"])
            gbp = d["gbp"]
            row = [s, sid, m.get("status", ""), m.get("tick_genus", "") or "(unrecorded)", m.get("tick_species", ""), c, L, cov,
                   r["lca_taxid"], vrow[r["lca_taxid"]]["name"], vrow[r["lca_taxid"]]["virus_group"],
                   vrow[r["lca_taxid"]]["family"], r["n_hits"], vb, round(vb / L, 4),
                   d["bv"][0], d["bv"][1], d["bv"][2], d["bv"][3], d["bv"][4],
                   span(d["o"]), minus(d["o"], d["v"]), d["bo"][0], d["bo"][1],
                   span(gbp.get("tick", [])), span(gbp.get("other arthropod", [])),
                   span(gbp.get("bacteria", [])), span(gbp.get("vertebrate", [])),
                   ";".join(f"{k}:{v}" for k, v in sorted(d["g"].items(), key=lambda x: -x[1]))]
            o.write("\t".join(map(str, row)) + "\n")

    # ---- sequences, sample-prefixed -------------------------------------------
    with open(W / "all_viral.fa", "w") as o:
        for f in sorted(glob.glob(str(OV / "fasta" / "*.viral.fasta"))):
            s = os.path.basename(f)[: -len(".viral.fasta")]
            for ln in open(f):
                o.write(f">{s}|{ln[1:].split()[0]}\n" if ln.startswith(">") else ln)
    print(f"{len(lca):,} contigs, {len(vrow)} taxa", file=sys.stderr)


if __name__ == "__main__":
    main()
