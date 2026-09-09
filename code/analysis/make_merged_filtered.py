#!/usr/bin/env python3
"""Build taxdump/merged.filtered.dmp.

The project's nodes.dmp/names.dmp are the April 2026 NCBI dump and the existing
LCA calls were made against them, so they must not be replaced -- swapping them
would make backfilled samples incomparable with the 120 already processed.

merged.dmp is taken from the *current* dump purely to rescue staxids that NCBI
has merged since April 2026 (previously these fell through to "unclassified").
But ~2% of its targets are taxa created after April 2026 and so are absent from
the local nodes.dmp; taxopy raises KeyError on those during load. This script
drops exactly those rows and keeps the rest.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TAX = ROOT / "taxdump"

nodes = set()
with open(TAX / "nodes.dmp") as fh:
    for line in fh:
        nodes.add(int(line.split("\t", 1)[0]))

kept = dropped = 0
with open(TAX / "merged.dmp") as fh, open(TAX / "merged.filtered.dmp", "w") as out:
    for line in fh:
        f = [x.strip() for x in line.split("|")]
        try:
            old, new = int(f[0]), int(f[1])
        except (ValueError, IndexError):
            continue
        if new in nodes:
            out.write(f"{old}\t|\t{new}\t|\n")
            kept += 1
        else:
            dropped += 1

print(f"nodes.dmp taxids       : {len(nodes):,}", file=sys.stderr)
print(f"merged rows kept       : {kept:,}", file=sys.stderr)
print(f"merged rows dropped    : {dropped:,}  (target postdates local nodes.dmp)",
      file=sys.stderr)
print(f"wrote {TAX / 'merged.filtered.dmp'}", file=sys.stderr)
