#!/usr/bin/env python3
"""Give each credible tick virus from the audit its own directory beside
../north_shore_virus, so it can be carried forward (trees, prevalence).

Per virus, in extracting_pathogens/<name>/:

  contigs.tsv      every contig with that LCA, all categories (from contig_verdicts.tsv)
  genuine.fasta    genuine contigs
  eve.fasta        EVE candidates (only if there are any)
  prevalence.tsv   per sample: tick genus, genuine/EVE contigs, genuine bp,
                   longest genuine contig, median k-mer coverage

The list is by LCA name, i.e. contigs assigned only to a higher rank
(Nairoviridae sp., Riboviria, ...) are not included.
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
EP = OV.parent
VIRUSES = ["Dali Nairo tick virus 1", "Densovirinae sp.", "Hulunbuir Totiv tick virus 2",
           "Butlers creek virus", "Cheeloo toti-like virus 2", "Kaisodi virus", "Haseki tick virus",
           "Bole tick virus 4", "Alongshan virus"]


def main() -> None:
    seqs, name = {}, None
    for ln in open(W / "all_viral.fa"):
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(ln.strip())
    rows = list(csv.DictReader(open(W / "contig_verdicts.tsv"), delimiter="\t"))
    cols = list(rows[0])
    print("virus\tcontigs\tgenuine\tEVE\tgenuine_samples")
    for v in VIRUSES:
        rs = [r for r in rows if r["lca_name"] == v]
        d = EP / v.replace(" ", "_").rstrip(".")
        d.mkdir(exist_ok=True)
        with open(d / "contigs.tsv", "w") as o:
            o.write("\t".join(cols) + "\n")
            for r in rs:
                o.write("\t".join(r[c] for c in cols) + "\n")
        for cat, fn in (("genuine", "genuine.fasta"), ("EVE candidate", "eve.fasta")):
            sel = [r for r in rs if r["category"] == cat]
            if sel:
                with open(d / fn, "w") as o:
                    for r in sel:
                        o.write(f">{r['id']}\n" + "\n".join(seqs[r["id"]]) + "\n")
        by = {}
        for r in rs:
            by.setdefault(r["sample_id"], []).append(r)
        with open(d / "prevalence.tsv", "w") as o:
            o.write("sample_id\ttick_genus\tgenuine_contigs\teve_contigs\tgenuine_bp\tmax_len\tmedian_cov\n")
            for s, sr in sorted(by.items()):
                g = [r for r in sr if r["category"] == "genuine"]
                e = [r for r in sr if r["category"] == "EVE candidate"]
                if not g and not e:
                    continue
                o.write("\t".join(map(str, [
                    s, sr[0]["tick_genus"], len(g), len(e), sum(int(r["len"]) for r in g),
                    max((int(r["len"]) for r in g), default=0),
                    round(statistics.median(float(r["kmer_cov"]) for r in g), 2) if g else ""])) + "\n")
        ng = sum(r["category"] == "genuine" for r in rs)
        ne = sum(r["category"] == "EVE candidate" for r in rs)
        print(f"{v}\t{len(rs)}\t{ng}\t{ne}\t{len({r['sample_id'] for r in rs if r['category'] == 'genuine'})}")


if __name__ == "__main__":
    main()
