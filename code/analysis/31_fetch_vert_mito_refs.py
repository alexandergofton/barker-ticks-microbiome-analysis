#!/usr/bin/env python3
"""Step 31 -- RefSeq vertebrate mitogenomes: FASTA, per-gene feature table and
a local BLAST database. Network step (run with PYTHON_NET); cached.

The blood-meal call itself is made from the samples' nt BLAST (step 30 cache +
contig LCA), because nt is competitive across all of life and RefSeq alone is
not. This curated set is used for two things nt cannot give:

  1. WHICH mitochondrial genes a host contig spans. Step 33 blasts each
     vertebrate-mito contig against these genomes and projects the hit onto the
     subject's annotated features (COX1, CYTB, 12S, ...).
  2. A second, curated best hit per contig (RefSeq mitogenomes are reviewed
     and complete, unlike most nt mito records).

Source: https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/
        mitochondrion.*.genomic.gbff.gz  (all organelles, ~400 MB; streamed and
        reduced to Vertebrata, taxid 7742, using taxdump/)

Outputs (results/bloodmeal/refs/)
  vert_mito_refseq.fa            one record per vertebrate mitogenome (accession)
  vert_mito_refseq.tsv           accession, taxid, organism, length
  vert_mito_features.tsv         accession, gene (canonical), type, start, end, strand
  blastdb/vert_mito_refseq.*     makeblastdb nucleotide database
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv  # noqa: E402
from markers_lib import require, run  # noqa: E402

import importlib
scan = importlib.import_module("30_scan_vertebrate_hits")

FTP = "https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/"

# Canonical vertebrate mitochondrial gene names. Checked in order; first match
# wins. Applied to the gene qualifier, then product, then note.
GENE_PATTERNS = [
    ("ND4L", r"\b(nd4l|nad4l|nadh4l)\b|nadh dehydrogenase subunit 4l"),
    ("COX1", r"\b(cox1|coi|co1|coxi|cox-?i)\b|cytochrome c oxidase subunit (1|i)\b"),
    ("COX2", r"\b(cox2|coii|co2|coxii)\b|cytochrome c oxidase subunit (2|ii)\b"),
    ("COX3", r"\b(cox3|coiii|co3|coxiii)\b|cytochrome c oxidase subunit (3|iii)\b"),
    ("ATP6", r"\b(atp6|atpase6|atpase ?6)\b|atp synthase f0 subunit 6|atp synthase subunit a\b"),
    ("ATP8", r"\b(atp8|atpase8|atpase ?8)\b|atp synthase f0 subunit 8|atp synthase protein 8"),
    ("ND1", r"\b(nd1|nad1|nadh1)\b|nadh dehydrogenase subunit 1\b"),
    ("ND2", r"\b(nd2|nad2|nadh2)\b|nadh dehydrogenase subunit 2\b"),
    ("ND3", r"\b(nd3|nad3|nadh3)\b|nadh dehydrogenase subunit 3\b"),
    ("ND4", r"\b(nd4|nad4|nadh4)\b|nadh dehydrogenase subunit 4\b"),
    ("ND5", r"\b(nd5|nad5|nadh5)\b|nadh dehydrogenase subunit 5\b"),
    ("ND6", r"\b(nd6|nad6|nadh6)\b|nadh dehydrogenase subunit 6\b"),
    ("CYTB", r"\b(cytb|cob|cyt ?b)\b|cytochrome b\b"),
    ("12S", r"\b(rrns|12s|s-rrna|rnr1|mt-rnr1)\b|small subunit ribosomal|12s ribosomal"),
    ("16S", r"\b(rrnl|16s|l-rrna|rnr2|mt-rnr2)\b|large subunit ribosomal|16s ribosomal"),
]
FEATURE_TYPES = {"CDS", "rRNA", "D-loop", "tRNA"}


def canonical_gene(ftype: str, quals: dict[str, str]) -> str:
    if ftype == "D-loop":
        return "CR"
    if ftype == "tRNA":
        return "tRNA"
    for key in ("gene", "product", "note"):
        text = quals.get(key, "").lower()
        if not text:
            continue
        for name, pat in GENE_PATTERNS:
            if re.search(pat, text):
                return name
    return ""


def parse_location(loc: str) -> tuple[int, int, str] | None:
    """Outer span of a GenBank location (join/complement/<>); None if remote."""
    if ":" in loc:  # feature on another record
        return None
    nums = [int(n) for n in re.findall(r"\d+", loc)]
    if not nums:
        return None
    return min(nums), max(nums), "-" if "complement" in loc else "+"


def parse_gbff(fh, keep_taxids: set[str]):
    """Yield (accession, taxid, organism, seq, features) for vertebrate records."""
    acc = org = taxid = ""
    feats: list[tuple[str, int, int, str, str]] = []
    seq_parts: list[str] = []
    section = ""
    cur_type = cur_loc = ""
    cur_quals: dict[str, str] = {}
    last_q = ""

    def flush_feature():
        nonlocal cur_type, cur_loc, cur_quals
        if cur_type in FEATURE_TYPES:
            span = parse_location(cur_loc)
            g = canonical_gene(cur_type, cur_quals)
            if span and g:
                feats.append((g, span[0], span[1], span[2], cur_type))
        cur_type, cur_loc, cur_quals = "", "", {}

    for line in fh:
        if line.startswith("LOCUS"):
            acc = org = taxid = ""
            feats, seq_parts, section = [], [], ""
        elif line.startswith("VERSION"):
            acc = line.split()[1]
        elif line.startswith("  ORGANISM"):
            org = line[12:].strip()
        elif line.startswith("FEATURES"):
            section = "feat"
        elif line.startswith("ORIGIN"):
            flush_feature()
            section = "seq"
        elif line.startswith("//"):
            if taxid in keep_taxids and seq_parts:
                yield acc, taxid, org, "".join(seq_parts).upper(), feats
            section = ""
        elif section == "feat":
            if len(line) > 5 and line[5] != " ":
                flush_feature()
                cur_type, cur_loc = line[5:21].strip(), line[21:].strip()
                last_q = ""
            elif line[21:22] == "/":
                m = re.match(r"/(\w+)=?\"?(.*?)\"?$", line[21:].rstrip())
                if m:
                    last_q = m.group(1)
                    cur_quals[last_q] = m.group(2)
                    if last_q == "db_xref" and m.group(2).startswith("taxon:"):
                        taxid = m.group(2).split(":")[1]
            elif not last_q and cur_type:
                cur_loc += line.strip()  # wrapped location
            elif last_q in ("product", "note", "gene"):
                cur_quals[last_q] += " " + line.strip().rstrip('"')
        elif section == "seq":
            seq_parts.append("".join(line.split()[1:]))


def list_gbff() -> list[str]:
    html = urllib.request.urlopen(FTP, timeout=60).read().decode()
    return sorted(set(re.findall(r'(mitochondrion\.\d+\.genomic\.gbff\.gz)', html)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="re-download and rebuild")
    args = ap.parse_args()

    out = scan.bloodmeal_root() / "refs"
    out.mkdir(parents=True, exist_ok=True)
    fa = out / "vert_mito_refseq.fa"
    db = out / "blastdb" / "vert_mito_refseq"
    if fa.exists() and Path(f"{db}.nsq").exists() and not args.force:
        print(f"cached: {fa}", file=sys.stderr)
        return

    files = list_gbff()
    if not files:
        sys.exit(f"ERROR: no mitochondrion.*.genomic.gbff.gz listed at {FTP}")
    for f in files:
        dest = out / f
        if dest.exists() and not args.force:
            continue
        print(f"downloading {FTP}{f}", file=sys.stderr)
        tmp = dest.with_name(dest.name + ".part")
        with urllib.request.urlopen(FTP + f, timeout=300) as r, open(tmp, "wb") as w:
            shutil.copyfileobj(r, w, 1 << 20)
        tmp.rename(dest)

    vert = scan.subtree(scan.VERT_TAXID)
    meta, feat_rows = [], []
    with open(fa, "w") as fo:
        for f in files:
            with gzip.open(out / f, "rt") as fh:
                for acc, taxid, org, seq, feats in parse_gbff(fh, vert):
                    fo.write(f">{acc}\n")
                    for i in range(0, len(seq), 80):
                        fo.write(seq[i:i + 80] + "\n")
                    meta.append({"accession": acc, "taxid": taxid, "organism": org,
                                 "length": len(seq)})
                    feat_rows += [{"accession": acc, "gene": g, "type": t, "start": s,
                                   "end": e, "strand": st} for g, s, e, st, t in feats]
    write_tsv(out / "vert_mito_refseq.tsv", meta, ["accession", "taxid", "organism", "length"])
    write_tsv(out / "vert_mito_features.tsv", feat_rows,
              ["accession", "gene", "type", "start", "end", "strand"])
    print(f"{len(meta)} vertebrate mitogenomes", file=sys.stderr)

    db.parent.mkdir(exist_ok=True)
    run([require("makeblastdb"), "-in", str(fa), "-dbtype", "nucl", "-out", str(db),
         "-title", "RefSeq vertebrate mitogenomes"])
    print(f"wrote {db}", file=sys.stderr)


if __name__ == "__main__":
    main()
