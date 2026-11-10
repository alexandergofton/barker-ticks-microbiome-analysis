#!/usr/bin/env python3
"""Fetch, from NCBI E-utilities, the viral nucleotide records our contigs hit.

For every viral subject in are_these_real/allhits.tsv (step 01):
  ncbi/hitrec.fa        nucleotide sequence
  ncbi/hitrec_prot.fa   CDS translations (fasta_cds_aa)
  ncbi/hitrec.ft        feature tables (CDS coordinates; hit-in-CDS test)

The nt snapshot used for the original BLAST is not mounted here, so the public
records are fetched instead. Cached: accessions already present are skipped.
Records longer than MAX_LEN (whole bacterial chromosomes mislabelled as viral,
giant viruses) are fetched as feature tables only.
"""

from __future__ import annotations

import csv
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
N = W / "ncbi"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BATCH = 100
MAX_LEN = 400_000


def efetch(ids: list[str], rettype: str, db: str = "nuccore") -> str:
    data = {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": "text",
            "tool": "tick_viral_audit", "email": "alexander.gofton@gmail.com"}
    if os.environ.get("NCBI_API_KEY"):
        data["api_key"] = os.environ["NCBI_API_KEY"]
    for attempt in range(5):
        try:
            req = urllib.request.Request(E + "efetch.fcgi", urllib.parse.urlencode(data).encode())
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read().decode()
        except Exception as e:  # noqa: BLE001 -- network; retry then give up loudly
            print(f"    retry {attempt + 1} ({rettype}): {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"efetch failed for {ids[:3]}...")


def have(path: Path, kind: str) -> set[str]:
    if not path.exists():
        return set()
    out = set()
    for ln in open(path):
        if kind == "fa" and ln.startswith(">"):
            out.add(ln[1:].split()[0])
        elif kind == "prot" and ln.startswith(">lcl|"):
            out.add(ln[5:].split("_prot_")[0].split("_cds_")[0])
        elif kind == "ft" and ln.startswith(">Feature"):
            out.add(ln.split("|")[1] if "|" in ln else ln.split()[1])
    return out


def main() -> None:
    N.mkdir(exist_ok=True)
    slen: dict[str, int] = {}
    with open(W / "allhits.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sgroup"] == "virus":
                slen[r["saccver"]] = max(slen.get(r["saccver"], 0), int(r["sstart"]), int(r["send"]))
    accs = sorted(slen)
    print(f"{len(accs):,} viral subject records", file=sys.stderr)

    # sequence length from esummary would cost a call per batch; the furthest
    # hit coordinate is a lower bound, good enough to keep megabase records out
    small = [a for a in accs if slen[a] <= MAX_LEN]
    jobs = [("hitrec.ft", "ft", "ft", accs),
            ("hitrec.fa", "fasta", "fa", small),
            ("hitrec_prot.fa", "fasta_cds_aa", "prot", small)]
    for fname, rettype, kind, want in jobs:
        done = have(N / fname, kind)
        todo = [a for a in want if a not in done]
        print(f"  {fname}: {len(done):,} cached, {len(todo):,} to fetch", file=sys.stderr)
        with open(N / fname, "a") as o:
            for i in range(0, len(todo), BATCH):
                o.write(efetch(todo[i:i + BATCH], rettype))
                o.flush()
                time.sleep(0.4)
                if (i // BATCH) % 10 == 0:
                    print(f"    {i + BATCH:,}/{len(todo):,}", file=sys.stderr)
    # the CDS-protein header names the record it came from; records without a
    # CDS annotation simply contribute nothing
    print("done", file=sys.stderr)


if __name__ == "__main__":
    main()
