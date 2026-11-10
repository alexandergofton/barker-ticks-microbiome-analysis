#!/usr/bin/env python3
"""Write metadata/contig_reassignments.tsv from the "are these real?" audits.

The nt LCA assigned tens of thousands of contigs to two tick RNA viruses whose
GenBank records carry a tick genomic repeat (see are_these_real.qmd in each
directory). Those contigs are tick genomic DNA, so they are re-assigned to
Ixodida and steps 03/04 then treat them like any other tick contig (scope
`host`, out of the microbial matrices). The lca/ files are never edited.

  Cheeloo Jingmen-like virus (3064049)  every contig       -> Ixodida
  North shore virus (2485875)           all but 9 contigs  -> Ixodida
      the 9 kept are are_these_real/cds.ids (RdRp protein support)

Sample names in the extraction tables are LCA file stems; they are written
as canonical sample_ids (metadata/sample_aliases.tsv).

    python3 extracting_pathogens/curate_tick_repeat_contigs.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "code" / "analysis"))
from common import canonical, load_aliases, write_tsv  # noqa: E402

NEW_TAXID, NEW_NAME = 6935, "Ixodida"
SOURCES = [
    # (directory, original taxid, file of genuine "sample|contig" ids or None, reason)
    ("Cheeloo_Jingmen-like_virus", 3064049, None,
     "tick genomic repeat in the 5' end of the Cheeloo Jingmen-like virus references; "
     "no jingmenvirus protein (Cheeloo_Jingmen-like_virus/are_these_real.html)"),
    ("north_shore_virus", 2485875, "are_these_real/cds.ids",
     "tick SINE-like repeat at the ends of the North shore virus references; "
     "no RdRp protein support (north_shore_virus/are_these_real.html)"),
]


def main() -> None:
    aliases = load_aliases()
    rows = []
    for d, taxid, keep_file, reason in SOURCES:
        keep = set()
        if keep_file:
            keep = {l.strip() for l in open(HERE / d / keep_file) if l.strip()}
        n_in = n_kept = 0
        with open(HERE / d / "lca_assignments.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if int(r["lca_taxid"]) != taxid:
                    continue
                n_in += 1
                if f"{r['sample']}|{r['query']}" in keep:
                    n_kept += 1
                    continue
                rows.append({"sample_id": canonical(r["sample"], aliases),
                             "contig_id": r["query"],
                             "original_taxid": taxid, "original_name": r["lca_name"],
                             "new_taxid": NEW_TAXID, "new_name": NEW_NAME,
                             "source": d, "reason": reason})
        if n_kept != len(keep):
            sys.exit(f"{d}: {len(keep)} genuine ids listed but {n_kept} found")
        print(f"  {d}: {n_in:,} contigs, {n_kept} kept as viral, "
              f"{n_in - n_kept:,} -> {NEW_NAME}", file=sys.stderr)
    write_tsv(ROOT / "metadata" / "contig_reassignments.tsv", rows,
              ["sample_id", "contig_id", "original_taxid", "original_name",
               "new_taxid", "new_name", "source", "reason"])


if __name__ == "__main__":
    main()
