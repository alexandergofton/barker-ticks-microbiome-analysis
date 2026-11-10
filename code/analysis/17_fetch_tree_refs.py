#!/usr/bin/env python3
"""Step 17 -- download the reference sequences that go INTO the trees.

Steps 11-16 stop at an alignment of sample sequences only, on the principle
that choosing a reference set is a judgement call. This step makes that call
from two curated panels instead of by hand:

  metadata/tree_references.tsv    per target x marker -- named species and
                                  published tick-derived sequences for the
                                  locus, plus an outgroup where one is safe
  metadata/tree_ref_genomes.tsv   per target -- complete genomes, which are the
                                  only source of several loci (Rickettsia sca4
                                  and 17 kDa, Coxiella icd/dnaK/IS1111,
                                  Francisella sdhA, all of Rickettsiella and
                                  Midichloria). Step 18 cuts the locus out of
                                  them exactly as step 13 cuts it out of a
                                  contig, so one genome set serves every locus
                                  of its target.

Everything is fetched from nuccore as NUCLEOTIDE, because the trees are built
on nucleotides.

CHOOSING BEFORE DOWNLOADING. Candidates are screened on esummary -- accession,
title, length -- and only the chosen ones are fetched. Two things fall out of
that: a query matching 400 near-identical Rickettsia felis records costs one
summary call rather than 400 downloads, and the quota can be spent on BREADTH,
capping how many records any one organism contributes
(`trees.curate.max_per_organism`) and then filling round-robin across
organisms. A tree of forty R. felis sequences has no topology worth reading.

Nothing here decides whether a record really carries the locus. Step 18 blasts
each one and drops it if it does not, so a loose query costs download time, not
correctness.

Cached: a locus whose FASTA already exists is skipped. --refresh re-fetches,
--target/--marker limit the scope.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import Entrez, dedupe_refs, read_fasta, slug, write_fasta  # noqa: E402
from trees_lib import (auto_query, cfg_trees, load_genome_panel, load_ref_panel,  # noqa: E402
                       loci_with_sequences, organism_from_title, trees_root)

MANIFEST_COLS = ["target", "marker", "tier", "n_records", "n_organisms",
                 "median_len", "query_source", "query_used", "path"]


def spread_by_organism(cands: list[dict], max_per_org: int, want: int) -> list[dict]:
    """Fill the quota round-robin across organisms, relevance order within each.

    Breadth over depth: one record from every species first, then a second from
    each, and so on until the quota is full or `max_per_org` is reached. NCBI's
    relevance sort puts the best-studied organism's many near-identical records
    at the top, which is exactly the wrong reference set for a tree.
    """
    by_org: dict[str, list[dict]] = defaultdict(list)
    for c in cands:
        by_org[organism_from_title(c["title"])].append(c)
    order = list(by_org)          # first-seen organism order = relevance order
    picked: list[dict] = []
    for rnd in range(max_per_org):
        for org in order:
            if len(picked) >= want:
                return picked
            if rnd < len(by_org[org]):
                picked.append({**by_org[org][rnd], "organism": org})
    return picked


def select(ez: Entrez, queries: list[str], want: int, max_len: int,
           max_per_org: int, multiplier: int) -> tuple[list[dict], str]:
    """Run the tiered queries until one yields records; return the chosen set."""
    for q in queries:
        try:
            ids = ez.esearch("nuccore", q, want * multiplier)
        except RuntimeError as e:
            print(f"    esearch error -- {e}", file=sys.stderr)
            continue
        if not ids:
            continue
        try:
            sums = [s for s in ez.esummary("nuccore", ids) if 0 < s["slen"] <= max_len]
        except RuntimeError as e:
            print(f"    esummary error -- {e}", file=sys.stderr)
            continue
        picked = spread_by_organism(sums, max_per_org, want)
        if picked:
            return picked, q
    return [], ""


def fetch(ez: Entrez, picked: list[dict], batch: int, max_len: int) -> dict[str, str]:
    """Download the chosen records. A failure raises, and the caller then writes
    NOTHING for that locus -- a half-downloaded reference set cached as complete
    is worse than no reference set, because a re-run would skip it."""
    seqs: dict[str, str] = {}
    ids = [p["id"] for p in picked]
    for i in range(0, len(ids), batch):
        fasta = ez.efetch_fasta("nuccore", ids[i:i + batch])
        cur: dict[str, str] = {}
        name, buf = None, []
        for line in fasta.splitlines():
            if line.startswith(">"):
                if name:
                    cur[name] = "".join(buf)
                name, buf = line[1:], []
            elif name:
                buf.append(line.strip())
        if name:
            cur[name] = "".join(buf)
        seqs.update(cur)
        print(f"      fetched {min(i + batch, len(ids))}/{len(ids)}", file=sys.stderr)
    return dedupe_refs(seqs, max_len)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    ap.add_argument("--marker", action="append", help="limit to this marker (repeatable)")
    ap.add_argument("--refresh", action="store_true", help="re-fetch cached reference sets")
    ap.add_argument("--skip-genomes", action="store_true",
                    help="gene-level references only (the genome tier is the slow part)")
    args = ap.parse_args()

    tc = cfg_trees()
    fc = tc["fetch"]
    cc = tc["curate"]
    root = trees_root() / "refs"
    for sub in ("gene", "outgroup", "genome"):
        (root / sub).mkdir(parents=True, exist_ok=True)

    api_key = os.environ.get(fc.get("api_key_env", "NCBI_API_KEY")) or None
    ez = Entrez(fc["email"], api_key)
    print(f"Entrez as {fc['email']}"
          f"{' with API key' if api_key else ' (no API key: 3 req/s)'}", file=sys.stderr)

    ref_panel = load_ref_panel()
    genome_panel = load_genome_panel()
    loci = loci_with_sequences(args.target, args.marker)
    if not loci:
        sys.exit("ERROR: no target x marker has sample sequences -- run run_markers.sh first.")
    print(f"{len(loci)} loci with recovered sample sequences", file=sys.stderr)

    manifest: list[dict] = []

    def record(target: str, marker: str, tier: str, path: Path, source: str, query: str):
        seqs = read_fasta(path) if path.exists() else {}
        lens = sorted(len(s) for s in seqs.values())
        orgs = {organism_from_title(n) for n in seqs}
        manifest.append({
            "target": target, "marker": marker, "tier": tier,
            "n_records": len(seqs), "n_organisms": len(orgs),
            "median_len": lens[len(lens) // 2] if lens else 0,
            "query_source": source, "query_used": query,
            "path": str(path.relative_to(ROOT)),
        })

    # ---- gene-level and outgroup references, per locus ---------------------
    for i, row in enumerate(loci, 1):
        target, marker = row["target"], row["marker"]
        key = f"{slug(target)}__{slug(marker)}"
        tag = f"[{i:>3}/{len(loci)}] {target} {marker}"
        cfg = ref_panel.get((target, marker))
        source = "panel" if cfg else "auto"
        if cfg is None:
            cfg = {"n_refs": cc["max_refs_per_locus"],
                   "queries": auto_query(target, marker, row["marker_type"]),
                   "outgroup_n": 0, "outgroup_queries": []}

        gene_fa = root / "gene" / f"{key}.fna"
        if gene_fa.exists() and not args.refresh:
            print(f"{tag}: cached, {len(read_fasta(gene_fa))} gene records", file=sys.stderr)
            record(target, marker, "gene", gene_fa, source, "(cached)")
        else:
            print(f"{tag}: searching ({source})", file=sys.stderr)
            picked, used = select(ez, cfg["queries"], cfg["n_refs"],
                                  fc["max_gene_seq_len"], cc["max_per_organism"],
                                  fc["search_multiplier"])
            if picked:
                try:
                    seqs = fetch(ez, picked, fc["batch_size"], fc["max_gene_seq_len"])
                except RuntimeError as e:
                    print(f"{tag}: FETCH FAILED -- {e}\n  nothing cached; re-run "
                          f"this step to retry the locus", file=sys.stderr)
                    seqs = {}
                if seqs:
                    write_fasta(gene_fa, seqs)
                    # Count organisms in what was actually WRITTEN, not in what
                    # was picked: dedupe_refs drops exact duplicate sequences,
                    # so the two numbers differ whenever a locus is deposited
                    # identically under several accessions.
                    print(f"{tag}: {len(seqs)} gene records, "
                          f"{len({organism_from_title(n) for n in seqs})} organisms",
                          file=sys.stderr)
                    record(target, marker, "gene", gene_fa, source, used)
            else:
                print(f"{tag}: no gene-level records from "
                      f"{len(cfg['queries'])} quer{'y' if len(cfg['queries']) == 1 else 'ies'}"
                      f" -- the genome tier may still supply this locus", file=sys.stderr)
                record(target, marker, "gene", gene_fa, source, "(none)")

        out_fa = root / "outgroup" / f"{key}.fna"
        if not (cfg["outgroup_n"] and cfg["outgroup_queries"]):
            continue
        if out_fa.exists() and not args.refresh:
            record(target, marker, "outgroup", out_fa, source, "(cached)")
            continue
        picked, used = select(ez, cfg["outgroup_queries"], cfg["outgroup_n"],
                              fc["max_gene_seq_len"], 1, fc["search_multiplier"])
        if picked:
            try:
                write_fasta(out_fa, fetch(ez, picked, fc["batch_size"],
                                          fc["max_gene_seq_len"]))
                print(f"{tag}: {len(picked)} outgroup records", file=sys.stderr)
                record(target, marker, "outgroup", out_fa, source, used)
            except RuntimeError as e:
                print(f"{tag}: outgroup FETCH FAILED -- {e}", file=sys.stderr)
        else:
            # Not fatal, and not silent: the tree is midpoint-rooted instead,
            # which the figure states on its face.
            print(f"{tag}: outgroup query returned nothing -- tree will be "
                  f"midpoint-rooted", file=sys.stderr)
            record(target, marker, "outgroup", out_fa, source, "(none)")

    # ---- whole genomes, per target ----------------------------------------
    targets = sorted({r["target"] for r in loci})
    if args.skip_genomes:
        print("\n--skip-genomes: genome tier not fetched", file=sys.stderr)
    else:
        print(f"\ngenome tier for {len(targets)} targets", file=sys.stderr)
        for target in targets:
            cfg = genome_panel.get(target)
            fa = root / "genome" / f"{slug(target)}.fna"
            if cfg is None:
                print(f"  {target}: no row in metadata/tree_ref_genomes.tsv "
                      f"-- gene-level references only", file=sys.stderr)
                continue
            if fa.exists() and not args.refresh:
                print(f"  {target}: cached, {len(read_fasta(fa))} genomes", file=sys.stderr)
                record(target, "*", "genome", fa, "panel", "(cached)")
                continue
            picked, used = select(ez, cfg["queries"], cfg["n_genomes"],
                                  fc["max_genome_seq_len"], cc["max_per_organism"],
                                  fc["search_multiplier"])
            if not picked:
                print(f"  {target}: genome query returned nothing", file=sys.stderr)
                record(target, "*", "genome", fa, "panel", "(none)")
                continue
            mb = sum(p["slen"] for p in picked) / 1e6
            print(f"  {target}: fetching {len(picked)} records "
                  f"({mb:.1f} Mb, {len({p['organism'] for p in picked})} organisms)",
                  file=sys.stderr)   # pre-fetch: the picked set is all we know
            try:
                seqs = fetch(ez, picked, fc["genome_batch_size"], fc["max_genome_seq_len"])
            except RuntimeError as e:
                print(f"  {target}: genome FETCH FAILED -- {e}\n"
                      f"    nothing cached; re-run this step to retry", file=sys.stderr)
                continue
            write_fasta(fa, seqs)
            record(target, "*", "genome", fa, "panel", used)

    write_tsv(trees_root() / "tree_reference_manifest.tsv", manifest, MANIFEST_COLS)
    empty = [m for m in manifest if m["tier"] == "gene" and m["n_records"] == 0]
    if empty:
        print(f"\n{len(empty)} loci have no gene-level references: "
              + ", ".join(f"{m['target']}/{m['marker']}" for m in empty[:12])
              + (" ..." if len(empty) > 12 else "")
              + "\n  They fall back to the genome tier. Fix the query in "
                "metadata/tree_references.tsv and re-run with --refresh if that "
                "is not enough.", file=sys.stderr)
    print(f"\nreferences under {(trees_root() / 'refs').relative_to(ROOT)}/", file=sys.stderr)


if __name__ == "__main__":
    main()
