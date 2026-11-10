#!/usr/bin/env python3
"""Where do the viral nt hits land on the reference records?

  hits_ref.tsv   every viral hit with the reference's CDS overlap (from the
                 feature tables, ncbi/hitrec.ft): in_cds = hit midpoint inside
                 an annotated CDS
  hotspots.tsv   windows of a reference that many contigs stack on
  hotspot_members.tsv  saccver, peak, contig id

Hotspot search, per reference record: slide a WINDOW-bp window, take the window
covered by the most DISTINCT contigs, record it, drop those contigs and repeat
while a window still holds >= MIN_CONTIGS contigs. This is the signature both
earlier audits found -- Cheeloo Jingmen-like virus and North shore virus
contigs piled on ~70-300 bp at the end of the records -- and a genuine viral
genome sampled by many contigs does not produce it (its contigs tile the whole
record instead). `share` = the hotspot's contigs / all contigs hitting the record.
"""

from __future__ import annotations

import collections
import csv
import re
import sys
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
WINDOW, MIN_CONTIGS, MAX_PEAKS = 300, 10, 5


def read_ft(path: Path) -> dict[str, list[tuple[int, int]]]:
    """accession -> CDS intervals (1-based, inclusive)."""
    cds = collections.defaultdict(list)
    acc, key = None, None
    if not path.exists():
        return cds
    for ln in open(path):
        if ln.startswith(">Feature"):
            m = re.search(r"\|([^|]+)\|", ln)
            acc = m.group(1) if m else ln.split()[1]
            key = None
            continue
        a = ln.rstrip("\n").split("\t")
        if len(a) >= 2 and a[0] and a[1]:
            s, e = (int(re.sub(r"[<>^]", "", x) or 0) for x in a[:2])
            if len(a) >= 3 and a[2]:
                key = a[2]
            if key == "CDS":
                cds[acc].append(tuple(sorted((s, e))))
        elif len(a) >= 3 and a[2]:
            key = a[2]
    return cds


def main() -> None:
    cds = read_ft(W / "ncbi" / "hitrec.ft")
    ce = {r["sample"] + "|" + r["contig"]: r for r in csv.DictReader(open(W / "contig_evidence.tsv"), delimiter="\t")}
    pe = {r["id"]: r["prot_call"] for r in csv.DictReader(open(W / "protein_evidence.tsv"), delimiter="\t")} \
        if (W / "protein_evidence.tsv").exists() else {}

    by_ref = collections.defaultdict(list)
    with open(W / "allhits.tsv") as fh, open(W / "hits_ref.tsv", "w") as o:
        o.write("id\tlca_name\tvirus_group\tsaccver\tssciname\tpident\tlength\tsstart\tsend\t"
                "ref_has_cds\tin_cds\n")
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sgroup"] != "virus":
                continue
            k = r["sample"] + "|" + r["contig"]
            s, e = sorted((int(r["sstart"]), int(r["send"])))
            mid = (s + e) / 2
            iv = cds.get(r["saccver"], [])
            inc = int(any(a <= mid <= b for a, b in iv))
            c = ce.get(k, {})
            o.write("\t".join([k, c.get("lca_name", ""), c.get("virus_group", ""), r["saccver"],
                               r["ssciname"], r["pident"], r["length"], str(s), str(e),
                               str(int(bool(iv))), str(inc)]) + "\n")
            by_ref[r["saccver"]].append((s, e, k, inc, r["ssciname"]))

    cols = ["saccver", "ssciname", "peak", "start", "end", "n_contigs", "n_samples", "share",
            "frac_in_cds", "ref_has_cds", "frac_prot_viral", "median_hit_len", "tick_genera",
            "lca_taxa", "contigs_hitting_ref"]
    rows, members = [], []
    for acc, hits in by_ref.items():
        alln = len({h[2] for h in hits})
        if alln < MIN_CONTIGS:
            continue
        left = list(hits)
        for peak in range(1, MAX_PEAKS + 1):
            starts = sorted({h[0] for h in left})
            best, bw = set(), None
            for w0 in starts:   # windows anchored at hit starts are enough
                w1 = w0 + WINDOW
                ks = {h[2] for h in left if h[0] >= w0 and h[1] <= w1 + 50}
                if len(ks) > len(best):
                    best, bw = ks, (w0, w1)
            if len(best) < MIN_CONTIGS:
                break
            inw = [h for h in left if h[2] in best and h[0] >= bw[0] and h[1] <= bw[1] + 50]
            st = sorted(h[0] for h in inw)
            en = sorted(h[1] for h in inw)
            lens = sorted(h[1] - h[0] + 1 for h in inw)
            gen = collections.Counter(ce[k]["tick_genus"] for k in best if k in ce)
            tax = collections.Counter(ce[k]["lca_name"] for k in best if k in ce)
            rows.append({
                "saccver": acc, "ssciname": hits[0][4], "peak": peak,
                "start": st[len(st) // 10], "end": en[len(en) * 9 // 10],   # 10-90% bounds
                "n_contigs": len(best), "n_samples": len({k.split("|")[0] for k in best}),
                "share": round(len(best) / alln, 3),
                "frac_in_cds": round(sum(h[3] for h in inw) / len(inw), 3),
                "ref_has_cds": int(bool(cds.get(acc))),
                "frac_prot_viral": round(sum(pe.get(k) == "viral" for k in best) / len(best), 3) if pe else "",
                "median_hit_len": lens[len(lens) // 2],
                "tick_genera": ";".join(f"{g}:{n}" for g, n in gen.most_common()),
                "lca_taxa": ";".join(f"{g}:{n}" for g, n in tax.most_common(4)),
                "contigs_hitting_ref": alln})
            members += [(acc, peak, k) for k in sorted(best)]
            left = [h for h in left if h[2] not in best]
    rows.sort(key=lambda r: -r["n_contigs"])
    with open(W / "hotspots.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(W / "hotspot_members.tsv", "w") as o:
        o.write("saccver\tpeak\tid\n")
        for m in members:
            o.write("\t".join(map(str, m)) + "\n")
    print(f"{len(by_ref):,} reference records hit; {len(rows)} hotspots >= {MIN_CONTIGS} contigs",
          file=sys.stderr)


if __name__ == "__main__":
    main()
