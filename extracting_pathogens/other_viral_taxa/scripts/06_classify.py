#!/usr/bin/env python3
"""Combine nucleotide, protein, adapter and hotspot evidence into one verdict
per contig (contig_verdicts.tsv) and per viral taxon (taxon_verdicts.tsv).

Contig classes, first match wins:

  adapter artefact          no viral protein, and every viral nt hit overlaps
                            Illumina/Nextera adapter or a poly-G/C run
  viral (protein-supported) viral protein beats every cellular one; no host flank
  phage (protein-supported) same, for bacteriophage taxa, no bacterial flank
  prophage in bacterial contig  phage protein + a bacterial gene/hit beside it
  EVE candidate (host flank)    viral protein + a eukaryotic (tick, arthropod,
                            vertebrate, ...) protein or nt hit OUTSIDE it
  viral ORF + long flank    viral protein, >= FLANK bp of the contig with no
                            viral evidence and no host evidence either
  cellular gene             a cellular protein beats every viral one
  hotspot repeat            no protein; contig is in a reference hotspot
  tick sequence             no protein; tick/arthropod nt co-hit
  nt-only (CDS)             no protein; viral hit inside a reference CDS
  nt-only (non-coding)      no protein; viral hit outside any annotated CDS

Taxon overrides, applied to protein-supported contigs (the viral and EVE
classes above), because their protein can be viral without the virus being
present:

  endogenous retrovirus     Retroviridae: vertebrate ERVs carried in blood-meal DNA
  endogenous adintovirus    Eupolintoviridae: polinton-like elements integrated
                            in host genomes
  human papillomavirus contaminant  Papillomaviridae protein >= 95 % aa identity
                            (human HPV types, i.e. reagent/handling contamination)

The classes then roll up into the five categories the audit asks about:
genuine, EVE, misclassified genomic, technical artefact, unverified.
"""

from __future__ import annotations

import collections
import csv
import re
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
FLANK = 500
HOST = {"tick", "other arthropod", "vertebrate", "other animal", "other eukaryote"}
# Nextera mosaic end, Nextera read 1/2 primer tails, TruSeq adapter, P5, P7 --
# both strands. Short enough to hit partial copies at contig ends.
ADAPTERS = ["CTGTCTCTTATACACATCT", "TCGTCGGCAGCGTCAGATG", "GTCTCGTGGGCTCGG",
            "AGATCGGAAGAGC", "ATCTCGTATGCCGTCTTCTGCTTG", "AATGATACGGCGACCACCGA"]
COMP = str.maketrans("ACGT", "TGCA")
PATS = [re.compile(p) for a in ADAPTERS for p in (a, a.translate(COMP)[::-1])] + \
       [re.compile(r"G{20,}"), re.compile(r"C{20,}")]

ROLLUP = {
    "viral (protein-supported)": "genuine", "phage (protein-supported)": "genuine",
    "prophage in bacterial contig": "genuine (prophage)",
    "EVE candidate (host flank)": "EVE candidate", "viral ORF + long flank": "EVE candidate",
    "cellular gene": "misclassified genomic", "hotspot repeat": "misclassified genomic",
    "tick sequence": "misclassified genomic",
    "endogenous retrovirus": "EVE candidate", "endogenous adintovirus": "EVE candidate",
    "adapter artefact": "technical artefact", "human papillomavirus contaminant": "technical artefact",
    "nt-only (CDS)": "unverified", "nt-only (non-coding)": "unverified",
}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def ranges(txt):
    return [tuple(map(int, x.split("-"))) for x in txt.split(";") if x]


def main() -> None:
    seqs, name = {}, None
    for ln in open(W / "all_viral.fa"):
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(ln.strip().upper())
    adapt = {}
    for k, v in seqs.items():
        s = "".join(v)
        adapt[k] = merge([(m.start() + 1, m.end()) for p in PATS for m in p.finditer(s)])

    ce = {r["sample"] + "|" + r["contig"]: r for r in csv.DictReader(open(W / "contig_evidence.tsv"), delimiter="\t")}
    pe = {r["id"]: r for r in csv.DictReader(open(W / "protein_evidence.tsv"), delimiter="\t")}
    hot = collections.defaultdict(set)
    for r in csv.DictReader(open(W / "hotspot_members.tsv"), delimiter="\t"):
        hot[r["id"]].add(f"{r['saccver']}#{r['peak']}")
    vhits = collections.defaultdict(list)
    for r in csv.DictReader(open(W / "hits_ref.tsv"), delimiter="\t"):
        vhits[r["id"]].append(r)
    qiv = collections.defaultdict(list)          # contig coordinates of viral nt hits
    for r in csv.DictReader(open(W / "allhits.tsv"), delimiter="\t"):
        if r["sgroup"] == "virus":
            qiv[r["sample"] + "|" + r["contig"]].append((int(r["qstart"]), int(r["qend"])))

    out_cols = ["id", "sample_id", "status", "tick_genus", "lca_taxid", "lca_name", "virus_group",
                "viral_family", "len", "kmer_cov", "class", "category", "prot_call", "viral_title",
                "viral_pident", "viral_prot_bp", "cell_group", "cell_title", "cell_only_bp",
                "viral_nt_bp", "best_viral_pident", "flank_bp", "adapter_bp", "hit_in_cds",
                "hotspots", "nonviral_groups"]
    rows = []
    for k, c in ce.items():
        p = pe.get(k, {"prot_call": "none", "viral_ranges": "", "cell_ranges": "", "cell_group": "",
                       "cell_only_bp": "0", "viral_prot_bp": "0", "viral_title": "", "cell_title": "",
                       "viral_pident": "0"})
        L = int(c["len"])
        ad = adapt.get(k, [])
        nt = qiv.get(k, [])
        # adapter-explained: every viral nt hit touches an adapter/homopolymer run (+-20 bp)
        def touches(iv):
            return any(min(iv[1], b + 20) >= max(iv[0], a - 20) for a, b in ad)
        adapter_only = bool(nt) and all(touches(iv) for iv in nt)
        vr = ranges(p["viral_ranges"])
        viral_region = merge(vr + nt)
        flank = L - sum(e - s + 1 for s, e in viral_region)
        cell_only = int(p["cell_only_bp"])
        cgrp = p["cell_group"]
        host_flank = (cell_only >= 150 and cgrp in HOST) or \
                     (int(c["nonviral_only_bp"]) >= 50 and c["best_nonviral_group"] in HOST)
        bact_flank = (cell_only >= 150 and cgrp == "bacteria") or int(c["bacteria_bp"]) > 0
        incds = [int(h["in_cds"]) for h in vhits.get(k, []) if h["ref_has_cds"] == "1"]
        call = p["prot_call"]
        phage = c["virus_group"] == "bacteriophage"

        if call != "viral" and adapter_only:
            cls = "adapter artefact"
        elif call == "viral" and int(p["viral_prot_bp"]) >= 90:
            if phage:
                cls = "prophage in bacterial contig" if bact_flank else "phage (protein-supported)"
            elif host_flank:
                cls = "EVE candidate (host flank)"
            elif flank >= FLANK:
                cls = "viral ORF + long flank"
            else:
                cls = "viral (protein-supported)"
        elif call == "cellular":
            cls = "cellular gene"
        elif k in hot:
            cls = "hotspot repeat"
        elif c["best_nonviral_group"] in ("tick", "other arthropod"):
            cls = "tick sequence"
        elif incds and max(incds):
            cls = "nt-only (CDS)"
        else:
            cls = "nt-only (non-coding)"
        if ROLLUP[cls] in ("genuine", "EVE candidate") and not phage:
            fam = c["viral_family"]
            if fam == "Retroviridae":
                cls = "endogenous retrovirus"
            elif fam == "Eupolintoviridae":
                cls = "endogenous adintovirus"
            elif fam == "Papillomaviridae" and float(p["viral_pident"] or 0) >= 95:
                cls = "human papillomavirus contaminant"
        rows.append({
            "id": k, "sample_id": c["sample_id"], "status": c["status"], "tick_genus": c["tick_genus"],
            "lca_taxid": c["lca_taxid"], "lca_name": c["lca_name"], "virus_group": c["virus_group"],
            "viral_family": c["viral_family"], "len": L, "kmer_cov": c["kmer_cov"], "class": cls,
            "category": ROLLUP[cls], "prot_call": call, "viral_title": p["viral_title"],
            "viral_pident": p["viral_pident"], "viral_prot_bp": p["viral_prot_bp"],
            "cell_group": cgrp, "cell_title": p["cell_title"], "cell_only_bp": cell_only,
            "viral_nt_bp": c["viral_bp"], "best_viral_pident": c["best_viral_pident"], "flank_bp": flank,
            "adapter_bp": sum(b - a + 1 for a, b in ad),
            "hit_in_cds": "" if not incds else round(sum(incds) / len(incds), 2),
            "hotspots": ";".join(sorted(hot.get(k, ()))), "nonviral_groups": c["nonviral_groups"]})
    with open(W / "contig_verdicts.tsv", "w") as o:
        o.write("\t".join(out_cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in out_cols) + "\n")

    # ---- per taxon ----------------------------------------------------------
    cats = ["genuine", "genuine (prophage)", "EVE candidate", "misclassified genomic",
            "technical artefact", "unverified"]
    by = collections.defaultdict(list)
    for r in rows:
        by[r["lca_taxid"]].append(r)
    vt = {r["taxid"]: r for r in csv.DictReader(open(W / "viral_taxa.tsv"), delimiter="\t")}
    with open(W / "taxon_verdicts.tsv", "w") as o:
        cols = ["taxid", "name", "rank", "virus_group", "family", "n_contigs", "n_samples", "total_bp",
                *cats, "genuine_samples", "eve_samples", "main_class", "verdict"]
        o.write("\t".join(cols) + "\n")
        for t, rs in sorted(by.items(), key=lambda x: -len(x[1])):
            n = len(rs)
            cc = collections.Counter(r["category"] for r in rs)
            cl = collections.Counter(r["class"] for r in rs)
            gen = cc["genuine"] + cc["genuine (prophage)"]
            # Verdict: what the MAJORITY of contigs are, but genuine and EVE
            # contigs are never hidden by a majority of artefacts -- a taxon
            # with 9 real contigs among 10,000 artefacts is "mostly artefact,
            # some genuine", exactly the North shore virus result.
            frac = {k: v / n for k, v in cc.items()}
            top = max(cats, key=lambda k: cc[k])
            if frac.get("genuine", 0) + frac.get("genuine (prophage)", 0) >= 0.5:
                v = "genuine" if cc["genuine"] >= cc["genuine (prophage)"] else "genuine (prophage)"
            elif frac.get("EVE candidate", 0) >= 0.5:
                v = "EVE"
            else:
                v = {"misclassified genomic": "misclassified genomic", "technical artefact": "technical artefact",
                     "unverified": "unverified", "EVE candidate": "EVE",
                     "genuine": "genuine", "genuine (prophage)": "genuine (prophage)"}[top]
                extra = [f"{cc[x]} {x}" for x in ("genuine", "genuine (prophage)", "EVE candidate")
                         if cc[x] and x != top]
                if extra:
                    v += " + " + ", ".join(extra)
            m = vt.get(t, {})
            o.write("\t".join(map(str, [
                t, m.get("name", rs[0]["lca_name"]), m.get("rank", ""), rs[0]["virus_group"], m.get("family", ""),
                n, len({r["sample_id"] for r in rs}), sum(r["len"] for r in rs), *[cc[k] for k in cats],
                len({r["sample_id"] for r in rs if r["category"].startswith("genuine")}),
                len({r["sample_id"] for r in rs if r["category"] == "EVE candidate"}),
                cl.most_common(1)[0][0], v])) + "\n")


if __name__ == "__main__":
    main()
