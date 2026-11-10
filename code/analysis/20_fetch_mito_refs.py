#!/usr/bin/env python3
"""Step 20 -- fetch the reference probes for the mitogenome workflow.

For every row of metadata/mito_genes.tsv (13 protein-coding genes + 2 rRNAs)
this runs the Entrez query, takes up to `n_refs` records and caches them under
results/mito/refs/. Every probe is fetched at ORDER level (Ixodida = all ticks)
because the sample set spans Ixodidae and Argasidae -- see the header of
metadata/mito_genes.tsv for why that matters.

These are search probes, exactly as in the marker-gene workflow (step 11):
they locate and orient each gene inside a sample's own contigs. They are not
part of the final alignments or tree. The named-species references that ARE in
the tree are whole published mitogenomes, fetched and annotated by step 20b.

Re-running is cheap: a gene whose FASTA already exists is skipped. Use
--refresh to force a re-fetch, --gene to limit the scope.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (Entrez, cfg_mito, dedupe_refs, load_mito_panel, mito_root,  # noqa: E402
                         read_fasta, write_fasta)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", action="append", help="limit to this gene (repeatable)")
    ap.add_argument("--refresh", action="store_true", help="re-fetch cached genes")
    args = ap.parse_args()

    mc = cfg_mito()
    rc = mc["refs"]
    refs_dir = mito_root() / "refs"
    refs_dir.mkdir(parents=True, exist_ok=True)

    api_key = os.environ.get(rc.get("api_key_env", "NCBI_API_KEY")) or None
    ez = Entrez(rc["email"], api_key, tool="barker_ticks_mito")
    print(f"Entrez as {rc['email']}"
          f"{' with API key' if api_key else ' (no API key: 3 req/s)'}", file=sys.stderr)

    full_panel = load_mito_panel()
    panel = full_panel
    if args.gene:
        panel = [r for r in panel if r["gene"] in set(args.gene)]

    used_query: dict[str, str] = {}
    for i, row in enumerate(panel, 1):
        db = "protein" if row["type"] == "pcg" else "nuccore"
        ext = "faa" if row["type"] == "pcg" else "fna"
        path = refs_dir / f"{row['key']}.{ext}"
        tag = f"[{i:>2}/{len(panel)}] {row['gene']}"

        if path.exists() and not args.refresh:
            print(f"{tag}: cached, {len(read_fasta(path))} refs", file=sys.stderr)
            continue

        got: dict[str, str] = {}
        used = ""
        for q in row["queries"]:
            try:
                ids = ez.esearch(db, q, row["n_refs"] * 3)
            except RuntimeError as e:
                print(f"{tag}: esearch error -- {e}", file=sys.stderr)
                continue
            if not ids:
                continue
            fasta = ez.efetch_fasta(db, ids[: row["n_refs"] * 3])
            tmp = refs_dir / f".{row['key']}.tmp"
            tmp.write_text(fasta)
            got = dedupe_refs(read_fasta(tmp), rc["max_seq_len"])
            tmp.unlink(missing_ok=True)
            got = dict(list(got.items())[: row["n_refs"]])
            if got:
                used = q
                break

        if len(got) < rc["min_refs_ok"]:
            print(f"{tag}: UNRESOLVED -- no usable sequences from "
                  f"{len(row['queries'])} quer{'y' if len(row['queries'])==1 else 'ies'}",
                  file=sys.stderr)
            continue

        write_fasta(path, got)
        print(f"{tag}: {len(got)} refs", file=sys.stderr)
        used_query[row["key"]] = used

    manifest, unresolved = [], []
    for row in full_panel:
        ext = "faa" if row["type"] == "pcg" else "fna"
        path = refs_dir / f"{row['key']}.{ext}"
        if path.exists():
            seqs = read_fasta(path)
            lens = sorted(len(s) for s in seqs.values())
            manifest.append({"gene": row["gene"], "type": row["type"],
                             "n_refs_fetched": len(seqs),
                             "median_ref_len": lens[len(lens) // 2] if lens else 0,
                             "path": str(path.relative_to(ROOT)),
                             "query_used": used_query.get(row["key"], "(cached)")})
        else:
            unresolved.append({"gene": row["gene"], "type": row["type"], "query": row["query"]})

    write_tsv(refs_dir / "reference_manifest.tsv", manifest,
              ["gene", "type", "n_refs_fetched", "median_ref_len", "path", "query_used"])
    write_tsv(refs_dir / "unresolved_queries.tsv", unresolved, ["gene", "type", "query"])
    if unresolved:
        print(f"\n{len(unresolved)} genes returned nothing -- see "
              f"results/mito/refs/unresolved_queries.tsv.", file=sys.stderr)
    print(f"\n{len(manifest)}/{len(full_panel)} mito genes have reference sequences",
          file=sys.stderr)


if __name__ == "__main__":
    main()
