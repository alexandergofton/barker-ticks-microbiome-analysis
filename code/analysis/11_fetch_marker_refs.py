#!/usr/bin/env python3
"""Step 11 -- fetch the reference sequences used to LOCATE each marker locus.

For every row of metadata/target_markers.tsv this runs the Entrez query, takes
up to `n_refs` records and caches them under results/markers/refs/.

WHAT THESE REFERENCES ARE FOR, AND WHAT THEY ARE NOT
  They are search probes. blastn (rRNA) and tblastn (protein) use them to find
  and orient the locus inside the assembled contigs. They are NOT the reference
  set for your phylogeny -- you add that by hand to the alignments produced by
  step 15, choosing type strains and named species yourself. Nothing downstream
  of step 14 uses these sequences.

Consequences of that framing:
  - A loose query is tolerable. Extra or redundant probes cost search time, not
    correctness, because every hit is re-scored against the panel thresholds.
  - Whole-genome records are dropped (`refs.max_seq_len`); a 2 Mb genome as a
    blastn query would hit every contig of the genus.
  - A query returning nothing is reported as UNRESOLVED and the row is skipped
    downstream rather than silently producing no hits.

Re-running is cheap: a locus whose FASTA already exists is skipped. Use
--refresh to force a re-fetch, --target/--marker to limit the scope.

Set NCBI_API_KEY in the environment to lift the rate limit from 3/s to 10/s.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (Entrez, cfg_markers, dedupe_refs, load_panel, out_root,  # noqa: E402
                         read_fasta, write_fasta)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    ap.add_argument("--marker", action="append", help="limit to this marker (repeatable)")
    ap.add_argument("--refresh", action="store_true", help="re-fetch cached loci")
    args = ap.parse_args()

    mc = cfg_markers()
    rc = mc["refs"]
    refs_dir = out_root() / "refs"
    refs_dir.mkdir(parents=True, exist_ok=True)

    api_key = os.environ.get(rc.get("api_key_env", "NCBI_API_KEY")) or None
    ez = Entrez(rc["email"], api_key)
    print(f"Entrez as {rc['email']}"
          f"{' with API key' if api_key else ' (no API key: 3 req/s)'}", file=sys.stderr)

    full_panel = load_panel()
    panel = full_panel
    if args.target:
        panel = [r for r in panel if r["target"] in set(args.target)]
    if args.marker:
        panel = [r for r in panel if r["marker"] in set(args.marker)]

    used_query: dict[str, str] = {}
    for i, row in enumerate(panel, 1):
        ext = "fna" if row["marker_type"] == "rrna" else "faa"
        path = refs_dir / f"{row['key']}.{ext}"
        tag = f"[{i:>3}/{len(panel)}] {row['target']} {row['marker']}"

        if path.exists() and not args.refresh:
            print(f"{tag}: cached, {len(read_fasta(path))} refs", file=sys.stderr)
            continue

        got: dict[str, str] = {}
        used = ""
        for q in row["queries"]:
            try:
                ids = ez.esearch(row["entrez_db"], q, row["n_refs"] * 3)
            except RuntimeError as e:
                print(f"{tag}: esearch error -- {e}", file=sys.stderr)
                continue
            if not ids:
                continue
            fasta = ez.efetch_fasta(row["entrez_db"], ids[: row["n_refs"] * 3])
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

    # The manifest and the unresolved list always describe the WHOLE panel,
    # rebuilt from what is on disk rather than from what this invocation
    # happened to touch. A --target run, or a run interrupted part way, must not
    # leave a manifest that silently claims the other loci do not exist.
    manifest, unresolved = [], []
    for row in full_panel:
        ext = "fna" if row["marker_type"] == "rrna" else "faa"
        path = refs_dir / f"{row['key']}.{ext}"
        if path.exists():
            manifest.append({**_mrow(row, path, read_fasta(path)),
                             "query_used": used_query.get(row["key"], "(cached)")})
        else:
            unresolved.append({"target": row["target"], "marker": row["marker"],
                               "entrez_db": row["entrez_db"], "query": row["query"]})

    write_tsv(refs_dir / "reference_manifest.tsv", manifest,
              ["target", "marker", "marker_type", "entrez_db", "n_refs_fetched",
               "median_ref_len", "path", "query_used"])
    write_tsv(refs_dir / "unresolved_queries.tsv", unresolved,
              ["target", "marker", "entrez_db", "query"])
    if unresolved:
        print(f"\n{len(unresolved)} panel rows returned nothing -- see "
              f"results/markers/refs/unresolved_queries.tsv. Fix the query in "
              f"metadata/target_markers.tsv and re-run; these loci are simply "
              f"absent from the search until then.", file=sys.stderr)
    print(f"\n{len(manifest)}/{len(full_panel)} panel loci have reference sequences",
          file=sys.stderr)


def _mrow(row: dict, path: Path, seqs: dict[str, str]) -> dict:
    lens = sorted(len(s) for s in seqs.values())
    return {"target": row["target"], "marker": row["marker"],
            "marker_type": row["marker_type"], "entrez_db": row["entrez_db"],
            "n_refs_fetched": len(seqs),
            "median_ref_len": lens[len(lens) // 2] if lens else 0,
            "path": str(path.relative_to(ROOT))}


if __name__ == "__main__":
    main()
