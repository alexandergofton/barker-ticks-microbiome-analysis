#!/usr/bin/env python3
"""Characterise the reference hotspots (step 05) as sequence elements.

For every hotspot window:
  * cut the window out of the reference record (ncbi/hitrec.fa)
  * flag Illumina/Nextera adapter or poly-G/C content (technical sequence)
  * flag a tRNA B-box (GGTTCGANTCC, either strand) -- the hallmark of
    tRNA-derived SINEs, as found at the 5' end of the North shore virus records
Then cluster the non-adapter windows into ELEMENTS (all-vs-all blastn,
single linkage) and, per element (longest window as probe):
  * copies across every viral-LCA contig (blastn vs a db of all_viral.fa)
  * flank diversity: the 150 bp beside each copy, clustered at 90% (cd-hit-est)
    -- an interspersed genomic repeat sits in many different contexts; a
    genuine viral sequence sits in one (the rest of the virus)
  * similarity to the Cheeloo Jingmen-like virus and North shore virus elements

Writes hotspot_windows.fa, elements.tsv, element_copies.tsv, element_flanks.tsv.
"""

from __future__ import annotations

import collections
import csv
import re
import subprocess
import sys
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
EP = OV.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
ADAPT = __import__("06_classify").PATS
BBOX = re.compile(r"GGTTCGA.TCC|GGA.TCGAACC")
OTHER = {"Cheeloo element": EP / "Cheeloo_Jingmen-like_virus/are_these_real/elem.fa",
         "NSV 5' element": EP / "north_shore_virus/are_these_real/elem5.fa"}


def read_fa(path):
    seqs, n = {}, None
    for ln in open(path):
        if ln.startswith(">"):
            n = ln[1:].split()[0]
            seqs[n] = []
        else:
            seqs[n].append(ln.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True, **kw).stdout


def main() -> None:
    ref = read_fa(W / "ncbi" / "hitrec.fa")
    hs = list(csv.DictReader(open(W / "hotspots.tsv"), delimiter="\t"))
    win = {}
    with open(W / "hotspot_windows.fa", "w") as o:
        for h in hs:
            s = ref.get(h["saccver"])
            if not s:
                continue
            a, b = int(h["start"]), int(h["end"])
            # widen to the full stacked span: 10-90% bounds trim the tails
            a, b = max(1, a - 20), min(len(s), b + 20)
            k = f"{h['saccver']}#{h['peak']}"
            win[k] = s[a - 1:b]
            o.write(f">{k}\n{win[k]}\n")

    flags = {}
    for k, s in win.items():
        ad = sum(len(m.group(0)) for p in ADAPT for m in p.finditer(s))
        flags[k] = {"adapter_bp": ad, "bbox": int(bool(BBOX.search(s))), "len": len(s),
                    "gc": round((s.count("G") + s.count("C")) / max(len(s), 1), 3)}

    # ---- cluster the biological (non-adapter) windows into elements ----------
    bio = [k for k in win if flags[k]["adapter_bp"] < 15]
    with open(W / "bio_windows.fa", "w") as o:
        for k in bio:
            o.write(f">{k}\n{win[k]}\n")
    parent = {k: k for k in bio}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    if bio:
        for ln in sh(f"blastn -task blastn -query {W}/bio_windows.fa -subject {W}/bio_windows.fa "
                     f"-evalue 1e-10 -outfmt '6 qseqid sseqid pident length'").splitlines():
            q, s, pid, L = ln.split("\t")
            if q != s and int(L) >= 40:
                parent[find(q)] = find(s)
    groups = collections.defaultdict(list)
    for k in bio:
        groups[find(k)].append(k)
    hsd = {f"{h['saccver']}#{h['peak']}": h for h in hs}
    elem_of, elems = {}, []
    for i, (root, ks) in enumerate(sorted(groups.items(), key=lambda g: -sum(int(hsd[k]["n_contigs"]) for k in g[1])), 1):
        eid = f"E{i}"
        for k in ks:
            elem_of[k] = eid
        elems.append((eid, ks, max(ks, key=lambda k: len(win[k]))))
    for k in win:
        elem_of.setdefault(k, "adapter")

    # ---- copies and flanks per element -----------------------------------------
    db = W / "allviraldb"
    if not Path(str(db) + ".nsq").exists() and not Path(str(db) + ".00.nsq").exists():
        sh(f"makeblastdb -in {W}/all_viral.fa -dbtype nucl -out {db}")
    allseq = read_fa(W / "all_viral.fa")
    comp = str.maketrans("ACGTN", "TGCAN")
    gen = {r["id"]: r["tick_genus"] for r in csv.DictReader(open(W / "contig_verdicts.tsv"), delimiter="\t")} \
        if (W / "contig_verdicts.tsv").exists() else {}
    erows, crow, frow = [], [], []
    for eid, ks, probe in elems:
        with open(W / f"{eid}.fa", "w") as o:
            o.write(f">{eid}|{probe}\n{win[probe]}\n")
        hits = sh(f"blastn -task blastn -query {W}/{eid}.fa -db {db} -evalue 1e-5 -max_target_seqs 100000 "
                  f"-outfmt '6 sseqid slen sstart send qstart qend pident length'").splitlines()
        seen, flanks = set(), {"up": [], "down": []}
        g = collections.Counter()
        for ln in hits:
            c, slen, ss, se, qs, qe, pid, L = ln.split("\t")
            slen, ss, se = int(slen), int(ss), int(se)
            crow.append((eid, c, ss, se, pid, L))
            if c in seen:
                continue
            seen.add(c)
            g[gen.get(c, "?")] += 1
            s = allseq[c]
            plus = ss < se
            lo, hi = min(ss, se), max(ss, se)
            up = s[max(0, lo - 151):lo - 1] if plus else s[hi:hi + 150].translate(comp)[::-1]
            dn = s[hi:hi + 150] if plus else s[max(0, lo - 151):lo - 1].translate(comp)[::-1]
            if len(up) == 150:
                flanks["up"].append((c, up))
            if len(dn) == 150:
                flanks["down"].append((c, dn))
        fl = {}
        for side, fs in flanks.items():
            if len(fs) < 2:
                fl[side] = (len(fs), len(fs), 0)
                continue
            fa = W / f"{eid}_{side}.fa"
            fa.write_text("".join(f">{c}\n{x}\n" for c, x in fs))
            sh(f"cd-hit-est -i {fa} -o {W}/{eid}_{side}_c90 -c 0.9 -n 8 -T 2 -M 1500 -d 0")
            sizes = collections.Counter()
            cl = None
            for ln in open(f"{W}/{eid}_{side}_c90.clstr"):
                if ln.startswith(">Cluster"):
                    cl = ln
                else:
                    sizes[cl] += 1
            fl[side] = (len(fs), len(sizes), max(sizes.values()))
        sim = []
        for nm, path in OTHER.items():
            if path.exists():
                o = sh(f"blastn -task blastn -query {W}/{eid}.fa -subject {path} -evalue 1e-3 "
                       f"-outfmt '6 pident length'").splitlines()
                if o:
                    sim.append(f"{nm} {o[0].split()[0]}%/{o[0].split()[1]}bp")
        taxa = collections.Counter()
        for k in ks:
            for t in hsd[k]["lca_taxa"].split(";"):
                if ":" in t:
                    n_, v = t.rsplit(":", 1)
                    taxa[n_] += int(v)
        erows.append({"element": eid, "probe": probe, "probe_len": len(win[probe]),
                      "n_windows": len(ks), "records": ";".join(sorted({k.split("#")[0] for k in ks})),
                      "ref_species": ";".join(sorted({hsd[k]["ssciname"] for k in ks})),
                      "hotspot_contigs": sum(int(hsd[k]["n_contigs"]) for k in ks),
                      "copies": len(hits), "contigs_with_copy": len(seen),
                      "genera": ";".join(f"{a}:{b}" for a, b in g.most_common()),
                      "up_flanks": fl["up"][0], "up_clusters": fl["up"][1], "up_largest": fl["up"][2],
                      "down_flanks": fl["down"][0], "down_clusters": fl["down"][1], "down_largest": fl["down"][2],
                      "bbox": max(flags[k]["bbox"] for k in ks), "gc": flags[probe]["gc"],
                      "frac_in_cds": round(sum(float(hsd[k]["frac_in_cds"]) * int(hsd[k]["n_contigs"]) for k in ks)
                                           / max(1, sum(int(hsd[k]["n_contigs"]) for k in ks)), 3),
                      "lca_taxa": ";".join(f"{a}:{b}" for a, b in taxa.most_common(5)),
                      "similar_to": "; ".join(sim), "sequence": win[probe]})
    cols = list(erows[0].keys()) if erows else ["element"]
    with open(W / "elements.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in erows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(W / "element_copies.tsv", "w") as o:
        o.write("element\tid\tsstart\tsend\tpident\tlength\n")
        for r in crow:
            o.write("\t".join(map(str, r)) + "\n")
    with open(W / "hotspot_windows.tsv", "w") as o:
        o.write("hotspot\telement\tlen\tadapter_bp\tbbox\tgc\tsequence\n")
        for k, s in win.items():
            f = flags[k]
            o.write(f"{k}\t{elem_of[k]}\t{f['len']}\t{f['adapter_bp']}\t{f['bbox']}\t{f['gc']}\t{s}\n")
    print(f"{len(win)} hotspot windows: {len(win) - len(bio)} adapter, {len(bio)} biological "
          f"in {len(elems)} elements", file=sys.stderr)


if __name__ == "__main__":
    main()
