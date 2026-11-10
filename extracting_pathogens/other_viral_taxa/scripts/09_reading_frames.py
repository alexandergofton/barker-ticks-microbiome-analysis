#!/usr/bin/env python3
"""Reading-frame integrity of the protein-supported eukaryotic-virus contigs.

For every contig of a eukaryotic or unresolved virus whose class is
viral (protein-supported), EVE candidate (host flank), viral ORF + long
flank, endogenous retrovirus or endogenous adintovirus, its best viral protein (step 04) is aligned back to it with tblastn
(-seg no). Per contig:

  stops     stop codons inside the aligned translation (best HSP)
  frames    distinct reading frames among HSPs to that protein covering
            different parts of it -- > 1 means a frameshift
  aa_cov    fraction of the protein covered

Genuine RNA virus carried over into a DNA extract should be intact; an old
endogenous copy usually carries stops and frameshifts. A young EVE can be
intact, so intact frames do not exclude integration -- flanks decide that.

Writes reading_frames.tsv.
"""

from __future__ import annotations

import collections
import csv
import subprocess
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
CLASSES = {"viral (protein-supported)", "EVE candidate (host flank)", "viral ORF + long flank",
           "endogenous retrovirus", "endogenous adintovirus"}


def main() -> None:
    cv = [r for r in csv.DictReader(open(W / "contig_verdicts.tsv"), delimiter="\t")
          if r["class"] in CLASSES and r["virus_group"] != "bacteriophage"]
    pe = {r["id"]: r["viral_sseqid"] for r in csv.DictReader(open(W / "protein_evidence.tsv"), delimiter="\t")}
    pairs = {r["id"]: pe[r["id"]] for r in cv if pe.get(r["id"])}
    (W / "rf_prot.ids").write_text("\n".join(sorted(set(pairs.values()))) + "\n")
    (W / "rf_contig.ids").write_text("\n".join(sorted(pairs)) + "\n")
    subprocess.run(f"seqkit grep -f {W}/rf_prot.ids {W}/protdb.fa > {W}/rf_prot.fa 2>/dev/null; "
                   f"seqkit grep -f {W}/rf_contig.ids {W}/all_viral.fa > {W}/rf_contigs.fa 2>/dev/null",
                   shell=True, check=True)
    out = subprocess.run(
        ["tblastn", "-query", str(W / "rf_prot.fa"), "-subject", str(W / "rf_contigs.fa"), "-evalue", "1e-3",
         "-seg", "no", "-db_gencode", "1", "-max_target_seqs", "100000",
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send sframe bitscore qlen sseq"],
        check=True, text=True, capture_output=True).stdout
    h = collections.defaultdict(list)
    for ln in out.splitlines():
        a = ln.split("\t")
        if pairs.get(a[1]) == a[0]:
            h[a[1]].append(a)
    with open(W / "reading_frames.tsv", "w") as o:
        o.write("id\tprotein\tpident\taln_aa\taa_cov\tstops\tframes\n")
        for c, p in sorted(pairs.items()):
            hs = sorted(h.get(c, []), key=lambda a: -float(a[9]))
            if not hs:
                o.write(f"{c}\t{p}\t\t0\t0\t\t\n")
                continue
            best = hs[0]
            # HSPs that cover DIFFERENT parts of the protein (<50% overlap with a
            # better one) and sit in another frame = frameshift
            kept = []
            for a in hs:
                s, e = int(a[4]), int(a[5])
                if all(min(e, int(b[5])) - max(s, int(b[4])) < 0.5 * (e - s + 1) for b in kept):
                    kept.append(a)
            cov = set()
            for a in kept:
                cov.update(range(int(a[4]), int(a[5]) + 1))
            o.write("\t".join(map(str, [c, p, best[2], best[3], round(len(cov) / int(best[10]), 3),
                                        best[11].count("*"), len({a[8] for a in kept})])) + "\n")


if __name__ == "__main__":
    main()
