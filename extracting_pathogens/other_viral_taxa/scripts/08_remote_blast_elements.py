#!/usr/bin/env python3
"""Does each element (step 07) occur in a CELLULAR genome?

The element probes are subsequences of public GenBank virus records (not our
contigs). They are sent once to NCBI BLAST (megablast, core_nt) with viruses
excluded (NOT txid10239), so any hit is a non-viral record carrying the same
sequence -- a tick genome assembly hit is direct evidence of a tick repeat.
The local nt snapshot is not mounted on this node, hence the web service.

Writes remote_blast.json (cached raw result) and remote_blast.tsv
(element, hit rank, accession, title, organism, taxid, identity, aln length,
query coverage, evalue).
"""

from __future__ import annotations

import csv
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

OV = Path(__file__).resolve().parents[1]
W = OV / "are_these_real"
URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
HITS = 10


def call(params: dict, post: bool = False) -> bytes:
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=data) if post else urllib.request.Request(URL + "?" + data.decode())
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read()


def main() -> None:
    cache = W / "remote_blast.json"
    elems = list(csv.DictReader(open(W / "elements.tsv"), delimiter="\t"))
    if not cache.exists():
        query = "".join(f">{e['element']}\n{e['sequence']}\n" for e in elems)
        put = call({"CMD": "Put", "PROGRAM": "blastn", "MEGABLAST": "on", "DATABASE": "core_nt",
                    "QUERY": query, "ENTREZ_QUERY": "NOT txid10239[Organism:exp]",
                    "HITLIST_SIZE": HITS, "EXPECT": "1e-5", "TOOL": "tick_viral_audit",
                    "EMAIL": "alexander.gofton@gmail.com"}, post=True).decode()
        rid = re.search(r"RID = (\S+)", put).group(1)
        print(f"RID {rid}", file=sys.stderr)
        time.sleep(30)
        while True:
            st = call({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}).decode()
            s = re.search(r"Status=(\w+)", st).group(1)
            if s == "READY":
                break
            if s in ("FAILED", "UNKNOWN"):
                sys.exit(f"remote BLAST {s} (RID {rid})")
            time.sleep(60)
        raw = call({"CMD": "Get", "FORMAT_TYPE": "JSON2", "RID": rid})
        out = []
        with zipfile.ZipFile(io.BytesIO(raw)) as z:   # JSON2 = zip, one file per query
            for n in z.namelist():
                if n.endswith(".json") and n != f"{rid}.json":
                    out.append(json.loads(z.read(n)))
        cache.write_text(json.dumps(out))
    res = json.loads(cache.read_text())
    with open(W / "remote_blast.tsv", "w") as o:
        o.write("element\trank\taccession\ttitle\torganism\ttaxid\tpident\taln_len\tquery_cov\tevalue\n")
        for r in res:
            s = r["BlastOutput2"]["report"]["results"]["search"]
            q, qlen = s["query_title"], s["query_len"]
            hits = s.get("hits", [])
            if not hits:
                o.write(f"{q}\t0\t\tno non-viral hit\t\t\t\t\t\t\n")
            for i, h in enumerate(hits, 1):
                d, hsp = h["description"][0], h["hsps"][0]
                o.write("\t".join(map(str, [
                    q, i, d.get("accession", ""), d.get("title", ""), d.get("sciname", ""), d.get("taxid", ""),
                    round(100 * hsp["identity"] / hsp["align_len"], 1), hsp["align_len"],
                    round(100 * (hsp["query_to"] - hsp["query_from"] + 1) / qlen, 1), hsp["evalue"]])) + "\n")


if __name__ == "__main__":
    main()
