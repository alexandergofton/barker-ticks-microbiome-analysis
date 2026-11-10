#!/usr/bin/env python3
"""Step 30 -- one cached pass over each sample's raw nt BLAST, keeping only the
HSPs whose subject is a vertebrate (NCBI Vertebrata, taxid 7742).

Blood-meal identification (steps 30-35, run_bloodmeal.sh) needs the vertebrate
hits of every contig, not just the contig's LCA: the LCA cannot say whether the
contig is a mitochondrial segment or a piece of nuclear genome, and it carries
no identity, so a 90% match to the only Varanus mitogenome in nt would read as a
species-level call. The raw BLAST files are up to ~160 GB each (464 GB in all),
so they are streamed once with awk -- the same approach as step 26c -- and
reduced to a small per-sample cache:

  mito / nuc_rRNA subjects   EVERY vertebrate HSP is kept (these are the
                             contigs that will be named; the full hit list is
                             needed for a per-contig LCA and identity cap)
  nuclear subjects           dropped. The nuclear vertebrate background is
                             summarised from the contig LCA files instead
                             (step 31); holding a best hit per contig here
                             runs out of memory on host-rich samples
                             (>600k vertebrate contigs).

Subject type comes from the nt record title, as in step 26c
(`subject_type`): "mitochondri"/cytochrome/COX/COI/12S/16S -> mito,
18S/28S/5.8S/ITS/rRNA -> nuc_rRNA, else nuclear.

Output
  results/bloodmeal/vert_hits/<sample>.tsv.gz   cache, one per sample
    contig  subject_type  staxid  pident  length  qstart  qend  sstart  send
    evalue  bitscore  saccver  stitle

Runtime is dominated by disk read (~50-150 MB/s per file on /data). Already
cached samples are skipped unless --force.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_aliases, load_config, resolve_artefact  # noqa: E402

VERT_TAXID = "7742"  # Vertebrata
HIT_COLS = ["contig", "subject_type", "staxid", "pident", "length", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore", "saccver", "stitle"]


def subject_type(title: str) -> str:
    """Same rules as step 26c, so 'mito' means the same thing in both audits."""
    t = title.lower()
    if "mitochond" in t or re.search(r"cytochrome|\bcox[123]\b|\bcoi\b|12s|16s ribosomal", t):
        return "mito"
    if re.search(r"18s|28s|5\.8s|internal transcribed|\bits[12]?\b|ribosomal rna|rrna", t):
        return "nuc_rRNA"
    return "nuclear"


def subtree(root: str) -> set[str]:
    """Every taxid under `root` in taxdump/nodes.dmp, plus merged (old) taxids
    that now point into it -- the nt BLAST was run against an older taxonomy."""
    kids: dict[str, list[str]] = defaultdict(list)
    with open(ROOT / "taxdump" / "nodes.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t", 2)
            if p[0] != p[1]:
                kids[p[1]].append(p[0])
    out, stack = set(), [root]
    while stack:
        t = stack.pop()
        out.add(t)
        stack.extend(kids.get(t, []))
    with open(ROOT / "taxdump" / "merged.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t")
            if p[1].split("\t")[0] in out:
                out.add(p[0])
    return out


def bloodmeal_root() -> Path:
    cfg = load_config().get("bloodmeal", {}) or {}
    return ROOT / cfg.get("out_dir", "results/bloodmeal")


def analysed_samples() -> list[str]:
    """Every sample step 00 did not exclude (same population as the TSE)."""
    with open(ROOT / "results" / "audit" / "sample_completeness.tsv") as fh:
        return sorted(r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")
                      if not r["status"].startswith("excluded") and r["has_blast"] == "1")


def extract(blast: Path, cache: Path, taxid_file: Path) -> tuple[int, int]:
    # mawk is ~2x faster than gawk on these files and the scan is the whole
    # cost of the step. Field 7 can hold several ';'-separated taxids.
    awk_bin = "/usr/bin/mawk" if Path("/usr/bin/mawk").exists() else "awk"
    prog = ("NR==FNR{v[$1]=1; next} {n=split($7,t,\";\"); for(i=1;i<=n;i++) if(t[i] in v)"
            "{print $1\"\\t\"t[i]\"\\t\"$8\"\\t\"$9\"\\t\"$12\"\\t\"$13\"\\t\"$14\"\\t\"$15"
            "\"\\t\"$16\"\\t\"$17\"\\t\"$2\"\\t\"$3; break}}")
    proc = subprocess.Popen([awk_bin, "-F", "\t", prog, str(taxid_file), str(blast)],
                            stdout=subprocess.PIPE, text=True,
                            env={**os.environ, "LC_ALL": "C"})
    n_kept = n_nuc = 0
    tmp = cache.with_name(cache.name + ".tmp")
    with gzip.open(tmp, "wt") as out:
        out.write("\t".join(HIT_COLS) + "\n")
        for line in proc.stdout:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            st = subject_type(p[11])
            if st == "nuclear":
                n_nuc += 1
                continue
            out.write("\t".join([p[0], st] + p[1:]) + "\n")
            n_kept += 1
    if proc.wait() != 0:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"awk failed on {blast}")
    tmp.rename(cache)
    return n_kept, n_nuc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", action="append", help="limit to this sample (repeatable)")
    ap.add_argument("--jobs", type=int, default=2, help="parallel BLAST-file scans (I/O bound)")
    ap.add_argument("--force", action="store_true", help="re-scan already cached samples")
    args = ap.parse_args()

    out = bloodmeal_root() / "vert_hits"
    out.mkdir(parents=True, exist_ok=True)
    aliases = load_aliases()
    samples = analysed_samples()
    if args.sample:
        samples = [s for s in samples if s in set(args.sample)]

    todo = []
    for s in samples:
        cache = out / f"{s}.tsv.gz"
        blast = resolve_artefact(s, ROOT / "blast", ".blast.txt", aliases)
        if blast is None:
            print(f"  {s}: no BLAST file, skipped", file=sys.stderr)
        elif args.force or not cache.exists():
            todo.append((s, blast, cache))
    # Smallest first: most samples are done early and partial results are usable.
    todo.sort(key=lambda x: x[1].stat().st_size)
    if not todo:
        print("all vertebrate-hit caches present", file=sys.stderr)
        return

    taxid_file = out / "_vertebrata_taxids.txt"
    with open(taxid_file, "w") as fh:
        fh.writelines(f"{t}\n" for t in sorted(subtree(VERT_TAXID)))
    gb = sum(b.stat().st_size for _, b, _ in todo) / 1e9
    print(f"scanning {len(todo)} raw BLAST files, {gb:.0f} GB ({args.jobs} at a time)",
          file=sys.stderr)
    with ThreadPoolExecutor(args.jobs) as ex:
        futs = {ex.submit(extract, b, c, taxid_file): s for s, b, c in todo}
        for i, f in enumerate(futs, 1):
            nk, nn = f.result()
            print(f"  [{i}/{len(todo)}] {futs[f]}: {nk} mito/rRNA HSPs, "
                  f"{nn} nuclear HSPs dropped", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
