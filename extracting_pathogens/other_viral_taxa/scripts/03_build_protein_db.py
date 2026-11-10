#!/usr/bin/env python3
"""Build the competitive protein database for the blastx-style test.

A contig is only "viral at the protein level" if its best protein hit is viral
AND beats every cellular protein it also hits. A herpesvirus dUTPase or a
phage integrase has cellular homologues, so a viral-only database would
over-call viruses. Components (are_these_real/protdb_download/, fetched once):

  RefSeq viral proteins   viral.1.protein.faa.gz           group "virus"
  UniProt Swiss-Prot      uniprot_sprot.fasta.gz           group from OX= taxid
  I. scapularis proteome  GCF_016920785.2 protein.faa.gz   group "tick"
  hit-record proteins     ncbi/hitrec_prot.fa (step 02)     group "virus"
      CDS translations of the viral nt records our contigs hit, many of which
      (tick RNA viruses) are not in RefSeq.

Identical sequences are written once (first source wins, in the order above):
most hit-record phage proteins repeat RefSeq ones, and without this the
database (1.8M proteins) gets DIAMOND OOM-killed on this 3 GB node.

Writes protdb.fa (ids prefixed rv|, sp|, is|, hr|), protdb_groups.tsv
(id -> group, taxid, title) and the DIAMOND database protdb.dmnd.
"""

from __future__ import annotations

import gzip
import re
import subprocess
import sys
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
ROOT = OV.parents[1]
W = OV / "are_these_real"
D = W / "protdb_download"
sys.path.insert(0, str(ROOT / "code" / "analysis"))
from bloodmeal_lib import Taxonomy  # noqa: E402

GROUPS = [(10239, "virus"), (6935, "tick"), (6656, "other arthropod"),
          (7742, "vertebrate"), (33208, "other animal"), (2, "bacteria"),
          (2157, "archaea"), (2759, "other eukaryote")]


def fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    h, s = None, []
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if h:
                    yield h, "".join(s)
                h, s = ln[1:].rstrip(), []
            else:
                s.append(ln.strip())
    if h:
        yield h, "".join(s)


def main() -> None:
    tx = Taxonomy(ROOT / "taxdump")

    def group(t):
        p = set(tx.path(t))
        return next((g for m, g in GROUPS if m in p), "unclassified")

    n, seen_seq = {}, set()
    with open(W / "protdb.fa", "w") as fa, open(W / "protdb_groups.tsv", "w") as gt:
        gt.write("id\tgroup\ttaxid\ttitle\n")

        def put(pid, g, t, title, seq, src):
            h = hash(seq)
            if h in seen_seq:
                n["(duplicate sequences skipped)"] = n.get("(duplicate sequences skipped)", 0) + 1
                return
            seen_seq.add(h)
            fa.write(f">{pid}\n{seq}\n")
            gt.write(f"{pid}\t{g}\t{t}\t{title}\n")
            n[src] = n.get(src, 0) + 1

        for h, s in fasta(D / "viral.1.protein.faa.gz"):
            acc, _, title = h.partition(" ")
            put(f"rv|{acc}", "virus", "", title, s, "RefSeq viral")
        for h, s in fasta(D / "uniprot_sprot.fasta.gz"):
            m = re.search(r"OX=(\d+)", h)
            t = int(m.group(1)) if m else 0
            acc = h.split("|")[1]
            title = h.split(" ", 1)[1].split(" OX=")[0]
            put(f"sp|{acc}", group(t), t, title, s, "Swiss-Prot")
        for h, s in fasta(D / "iscap.protein.faa.gz"):
            acc, _, title = h.partition(" ")
            put(f"is|{acc}", "tick", 6945, title, s, "I. scapularis")
        hr = W / "ncbi" / "hitrec_prot.fa"
        if hr.exists():
            seen = set()
            for h, s in fasta(hr):
                # >lcl|OL452128.1_prot_UMO75525.1_1 [protein=RdRp] ...
                m = re.match(r"lcl\|(\S+?)_prot_(\S+)", h)
                pid = f"hr|{m.group(1)}|{m.group(2)}" if m else f"hr|{h.split()[0]}"
                if pid in seen or not s:
                    continue
                seen.add(pid)
                pm = re.search(r"\[protein=([^\]]+)\]", h)
                put(pid, "virus", "", pm.group(1) if pm else "", s, "hit records")
    for k, v in n.items():
        print(f"  {k:<14} {v:>9,}", file=sys.stderr)
    subprocess.run(["diamond", "makedb", "--in", str(W / "protdb.fa"), "-d", str(W / "protdb"),
                    "--threads", "2", "--quiet"], check=True)


if __name__ == "__main__":
    main()
