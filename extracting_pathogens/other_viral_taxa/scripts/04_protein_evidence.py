#!/usr/bin/env python3
"""DIAMOND blastx of every viral-LCA contig against protdb (step 03), summarised
per contig into are_these_real/protein_evidence.tsv.

Per contig:
  best_group / best_bits / best_title     best protein hit overall
  viral_bits, viral_pident, viral_aa, viral_title   best VIRAL protein hit
  cell_bits, cell_group, cell_title       best CELLULAR protein hit
  prot_call   viral     best viral hit >= best cellular hit
              cellular  a cellular protein beats every viral one
              none      no protein hit at e <= 1e-5
  viral_prot_bp   contig bp covered by viral HSPs that are not beaten by a
                  cellular HSP overlapping them by >= 50%
  cell_prot_bp    contig bp covered by cellular HSPs (any group)
  cell_only_bp    cellular-covered bp outside every viral HSP (a host gene
                  BESIDE the viral one -- the EVE flank signal)
  viral_ranges / cell_ranges   merged contig intervals ("s-e;s-e")

DIAMOND is run once (cached in diamond.tsv); delete that file to re-run.
Block size 0.2 / 1 index chunk keeps DIAMOND near 1 GB; 0.4 was OOM-killed
on this 3 GB node once the hit-record proteins were added.
"""

from __future__ import annotations

import collections
import csv
import subprocess
import sys
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen".split()


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def span(iv):
    return sum(e - s + 1 for s, e in merge(iv))


def main() -> None:
    dm = W / "diamond.tsv"
    if not dm.exists() or dm.stat().st_size == 0:
        tmp = W / "diamond_tmp"
        tmp.mkdir(exist_ok=True)
        subprocess.run(["diamond", "blastx", "-q", str(W / "all_viral.fa"), "-d", str(W / "protdb"),
                        "-o", str(dm) + ".part", "--sensitive", "-e", "1e-5", "-k", "25",
                        "-b", "0.2", "-c", "1", "-p", "2", "--quiet", "-t", str(tmp),
                        "-f", "6", *FIELDS], check=True)
        Path(str(dm) + ".part").rename(dm)

    grp = {}
    with open(W / "protdb_groups.tsv") as fh:
        next(fh)
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            grp[a[0]] = (a[1], a[3])

    hits = collections.defaultdict(list)
    with open(dm) as fh:
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            qs, qe = sorted((int(a[4]), int(a[5])))
            g, t = grp.get(a[1], ("unclassified", ""))
            hits[a[0]].append((float(a[9]), g, t, float(a[2]), int(a[3]), qs, qe, a[1]))

    cols = ["id", "prot_call", "best_group", "best_bits", "best_title", "viral_bits", "viral_pident",
            "viral_aa", "viral_sseqid", "viral_title", "cell_bits", "cell_group", "cell_title",
            "viral_prot_bp", "cell_prot_bp", "cell_only_bp", "viral_ranges", "cell_ranges"]
    with open(W / "protein_evidence.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for ln in open(W / "all_viral.fa"):
            if not ln.startswith(">"):
                continue
            q = ln[1:].split()[0]
            h = sorted(hits.get(q, []), reverse=True)
            vir = [x for x in h if x[1] == "virus"]
            cel = [x for x in h if x[1] != "virus"]
            b = h[0] if h else None
            v = vir[0] if vir else None
            c = cel[0] if cel else None
            call = "none" if not h else ("viral" if v and (not c or v[0] >= c[0]) else "cellular")
            # a viral HSP only counts if no overlapping cellular HSP outscores it
            vr = []
            for x in vir:
                beaten = any(y[0] > x[0] and min(x[6], y[6]) - max(x[5], y[5]) + 1 >= 0.5 * (x[6] - x[5] + 1)
                             for y in cel)
                if not beaten:
                    vr.append((x[5], x[6]))
            cr = [(y[5], y[6]) for y in cel]
            cell_only = span(cr) - span([(max(s, vs), min(e, ve)) for s, e in merge(cr)
                                         for vs, ve in merge(vr) if min(e, ve) >= max(s, vs)])
            row = [q, call, b[1] if b else "", b[0] if b else 0, b[2] if b else "",
                   v[0] if v else 0, v[3] if v else 0, v[4] if v else 0, v[7] if v else "",
                   v[2] if v else "", c[0] if c else 0, c[1] if c else "", c[2] if c else "",
                   span(vr), span(cr), cell_only,
                   ";".join(f"{s}-{e}" for s, e in merge(vr)),
                   ";".join(f"{s}-{e}" for s, e in merge(cr))]
            o.write("\t".join(map(str, row)) + "\n")
    n = collections.Counter(r["prot_call"] for r in csv.DictReader(open(W / "protein_evidence.tsv"),
                                                                      delimiter="\t"))
    print(dict(n), file=sys.stderr)


if __name__ == "__main__":
    main()
