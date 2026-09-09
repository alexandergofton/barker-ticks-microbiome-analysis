#!/usr/bin/env python3
"""Streaming LCA from BLAST outfmt6 -- a drop-in replacement for code/lca.py.

WHY THIS EXISTS
    code/lca.py does `pd.read_csv(blast_file)`, loading the entire BLAST output
    into a DataFrame. Seven of this project's BLAST files are 10-150 GB, so they
    OOM and were never processed. This version streams the file, holding only
    one query's hits at a time, so peak memory is independent of input size.

WHAT IS UNCHANGED
    The classification logic is identical to code/lca.py: filter on e-value and
    percent identity, keep hits within --top-percent of the best remaining
    bitscore, expand semicolon-separated staxids, drop taxids absent from the
    taxonomy, then taxopy.find_lca(). Output columns and ordering are the same.

WHAT IS ADDED
    * merged.dmp support -- BLAST staxids that NCBI has since merged were
      previously silently dropped to "unclassified". Remapped ids are counted
      and reported.
    * Memoisation of Taxon construction and of whole LCA lookups. Taxid sets
      repeat heavily across contigs; this is the difference between hours and
      days on the large files.
    * An input-ordering check. Streaming assumes BLAST emits all hits for a
      query contiguously (it does). If that assumption is ever violated the run
      fails loudly rather than producing a silently wrong LCA.

USAGE
    lca_stream.py -i blast/S.blast.txt -o lca/S.blast.lca \
        --nodes taxdump/nodes.dmp --names taxdump/names.dmp \
        --merged taxdump/merged.dmp
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from itertools import groupby

import taxopy

# outfmt6 column indices as produced by code/blast.slurm.sh
Q, STAXID, PIDENT, EVALUE, BITSCORE = 0, 6, 7, 15, 16
N_COLS = 17


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="BLAST outfmt6 results file")
    p.add_argument("-o", "--output", required=True, help="Output TSV")
    p.add_argument("--nodes", required=True)
    p.add_argument("--names", required=True)
    p.add_argument("--merged", default=None,
                   help="merged.dmp -- rescues staxids NCBI has since merged")
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--pident", type=float, default=0.0)
    p.add_argument("--top-percent", type=float, default=10.0,
                   help="keep hits within this %% of the top bitscore per query")
    p.add_argument("--no-order-check", action="store_true",
                   help="skip the duplicate-query guard (saves memory on huge files)")
    p.add_argument("--progress-every", type=int, default=1_000_000,
                   help="log progress every N input lines")
    return p.parse_args()


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_merged(path: str | None) -> set[int]:
    """Set of retired taxids that merged.dmp remaps. Used only for reporting --
    taxopy applies the remap internally and keeps Taxon.taxid as requested, so
    the rescue is otherwise invisible."""
    if not path:
        return set()
    out = set()
    with open(path) as fh:
        for line in fh:
            f = line.split("|", 1)[0].strip()
            if f.isdigit():
                out.add(int(f))
    return out


class Resolver:
    """Taxid -> Taxon with memoisation, merged-id rescue and LCA caching."""

    def __init__(self, taxdb: taxopy.TaxDb, merged: set[int] | None = None):
        self.db = taxdb
        self.merged = merged or set()
        self._taxon: dict[int, taxopy.Taxon | None] = {}
        self._lca: dict[tuple[int, ...], tuple[int, str, str]] = {}
        self.n_invalid = 0
        self.n_rescued = 0

    def taxon(self, taxid: int):
        """None if the taxid is unknown to the taxonomy."""
        if taxid in self._taxon:
            return self._taxon[taxid]
        t = None
        try:
            t = taxopy.Taxon(taxid, self.db)
            if taxid in self.merged:
                # would have been dropped to "unclassified" without merged.dmp
                self.n_rescued += 1
        except Exception:
            self.n_invalid += 1
        self._taxon[taxid] = t
        return t

    def lca(self, taxids: list[int]) -> tuple[int, str, str]:
        """(taxid, name, rank) for the LCA of valid members of `taxids`."""
        taxa = [t for t in (self.taxon(x) for x in taxids) if t is not None]
        if not taxa:
            return 0, "unclassified", "no rank"
        key = tuple(sorted({t.taxid for t in taxa}))
        hit = self._lca.get(key)
        if hit is not None:
            return hit
        if len(key) == 1:
            t = taxa[0]
            out = (t.taxid, t.name, t.rank)
        else:
            lca = taxopy.find_lca(taxa, self.db)
            out = (lca.taxid, lca.name, lca.rank)
        self._lca[key] = out
        return out


def expand_taxids(field: str) -> list[int]:
    """Handle the semicolon-separated staxids BLAST sometimes returns."""
    if not field or field == "N/A":
        return []
    if ";" not in field:
        return [int(field)] if field.isdigit() else []
    return [int(t) for t in field.split(";") if t.strip().isdigit()]


def main() -> None:
    args = parse_args()
    t0 = time.time()

    print("Loading taxonomy database...", file=sys.stderr)
    taxdb = taxopy.TaxDb(nodes_dmp=args.nodes, names_dmp=args.names,
                         merged_dmp=args.merged)
    resolver = Resolver(taxdb, load_merged(args.merged))

    print(f"Streaming BLAST results from {args.input}...", file=sys.stderr)
    seen: set[int] | None = None if args.no_order_check else set()

    n_lines = n_queries = n_malformed = n_unclassified = 0
    tmp_out = args.output + ".partial"

    def records(fh):
        """Yield (query, [(staxid_str, pident, evalue, bitscore), ...])."""
        nonlocal n_lines, n_malformed

        def rows():
            nonlocal n_lines, n_malformed
            for line in fh:
                n_lines += 1
                if args.progress_every and n_lines % args.progress_every == 0:
                    el = time.time() - t0
                    print(f"  {n_lines:,} lines  {n_queries:,} queries  "
                          f"{el:,.0f}s  ({n_lines/max(el,1e-9):,.0f} lines/s)",
                          file=sys.stderr)
                f = line.rstrip("\n").split("\t")
                if len(f) < N_COLS:
                    n_malformed += 1
                    continue
                yield f

        for query, group in groupby(rows(), key=lambda f: f[Q]):
            yield query, list(group)

    with opener(args.input) as fh, open(tmp_out, "w") as out:
        out.write("query\tlca_taxid\tlca_name\tlca_rank\tn_hits\ttop_bitscore\n")

        for query, hits in records(fh):
            n_queries += 1
            if seen is not None:
                h = hash(query)
                if h in seen:
                    sys.exit(
                        f"ERROR: query {query!r} reappears after its block closed. "
                        "The BLAST file is not grouped by query, so streaming would "
                        "split its hits and give a wrong LCA. Sort by column 1 first, "
                        "or rerun with the non-streaming code/lca.py."
                    )
                seen.add(h)

            # --- same filter order as code/lca.py ---------------------------
            kept = []
            for f in hits:
                try:
                    ev = float(f[EVALUE]); pid = float(f[PIDENT]); bs = float(f[BITSCORE])
                except ValueError:
                    n_malformed += 1
                    continue
                if ev <= args.evalue and pid >= args.pident:
                    kept.append((f[STAXID], bs))
            if not kept:
                continue

            best = max(bs for _, bs in kept)
            window = best * (1 - args.top_percent / 100)
            kept = [(sid, bs) for sid, bs in kept if bs >= window]

            taxids: list[int] = []
            for sid, _ in kept:
                taxids.extend(expand_taxids(sid))

            taxid, name, rank = resolver.lca(taxids)
            if taxid == 0:
                n_unclassified += 1
            out.write(f"{query}\t{taxid}\t{name}\t{rank}\t{len(kept)}\t{best}\n")

    os.replace(tmp_out, args.output)

    el = time.time() - t0
    print("\n".join([
        "",
        f"Done in {el/60:.1f} min -> {args.output}",
        f"  input lines        : {n_lines:,}",
        f"  queries            : {n_queries:,}",
        f"  unclassified       : {n_unclassified:,}",
        f"  malformed lines    : {n_malformed:,}",
        f"  distinct taxids    : {len(resolver._taxon):,}",
        f"  taxids not in db   : {resolver.n_invalid:,}",
        f"  distinct taxids rescued via merged.dmp : {resolver.n_rescued:,}",
        f"  cached LCA lookups : {len(resolver._lca):,}",
    ]), file=sys.stderr)


if __name__ == "__main__":
    main()
