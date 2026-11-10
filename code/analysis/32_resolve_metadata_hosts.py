#!/usr/bin/env python3
"""Step 32 -- resolve each metadata host to an NCBI taxon, and record whether
the reference databases can name that host at all.

metadata host_species is free text, "Common name (Scientific name)", with a
handful of misspellings, genus-only entries and one common/scientific
conflict. Corrections live in metadata/host_name_fixes.tsv (config
`bloodmeal.host_fixes`); everything else is resolved from its parenthesised
scientific name against taxdump/names.dmp (scientific names and synonyms, so
Canis familiaris and Macropus rufus resolve to their current NCBI names).

Reference availability is the key to reading a species-level mismatch later:
if nt holds no mitochondrial sequence for the recorded host, the blood meal can
only ever be named as its closest sequenced relative, and a "different species,
same genus" result is the expected outcome, not a conflict.

  refseq_mito_species/genus   RefSeq vertebrate mitogenomes (step 31) in the
                              host's species / genus -- offline
  nt_mito_species/genus       nt records for the host species / genus with a
                              mitochondrial title (Entrez count; network, cached
                              in refs/nt_mito_counts.tsv; skipped with --offline)

Output
  results/bloodmeal/tables/metadata_hosts.tsv   one row per distinct host_species
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from bloodmeal_lib import RANKS, Taxonomy, bloodmeal_root, cfg_bloodmeal, read_tsv  # noqa: E402

COLS = ["host_species", "n_samples", "resolved_name", "fix_note", "meta_taxid", "meta_rank",
        *[f"meta_{r}" for r in RANKS], "status",
        "refseq_mito_species", "refseq_mito_genus", "nt_mito_species", "nt_mito_genus"]

MITO_TERM = ("(mitochondrion[Title] OR mitochondrial[Title] OR \"cytochrome b\"[Title] "
             "OR \"cytochrome oxidase\"[Title] OR \"cytochrome c oxidase\"[Title] "
             "OR 12S[Title] OR 16S[Title])")


def lookup_names(wanted: set[str]) -> dict[str, list[tuple[str, str]]]:
    """lower-case name -> [(taxid, name class)], scientific names first."""
    out: dict[str, list[tuple[str, str]]] = {}
    with open(ROOT / "taxdump" / "names.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t")
            key = p[1].lower()
            if key in wanted:
                cls = p[3].rstrip("\t|\n")
                if cls in ("scientific name", "synonym", "equivalent name", "genbank synonym"):
                    out.setdefault(key, []).append((p[0], cls))
    for v in out.values():
        v.sort(key=lambda x: x[1] != "scientific name")
    return out


def nt_counts(taxa: dict[str, str], cache: Path, offline: bool) -> dict[str, str]:
    """taxid -> nt mito record count (string; '' when unknown)."""
    known = {r["taxid"]: r["count"] for r in read_tsv(cache)} if cache.exists() else {}
    todo = [t for t in taxa if t and t not in known]
    if todo and not offline:
        from markers_lib import Entrez, cfg_mito
        rc = cfg_mito()["refs"]  # same Entrez contact / API key as the mito workflow
        ez = Entrez(rc["email"], os.environ.get(rc.get("api_key_env", "NCBI_API_KEY")) or None,
                    tool="barker_ticks_bloodmeal")
        for i, t in enumerate(todo, 1):
            xml = ez._call("esearch.fcgi", {"db": "nuccore", "rettype": "count",
                                            "term": f"txid{t}[Organism:exp] AND {MITO_TERM}"})
            m = re.search(r"<Count>(\d+)</Count>", xml)
            known[t] = m.group(1) if m else ""
            print(f"  [{i}/{len(todo)}] {taxa[t]}: {known[t]} nt mito records",
                  file=sys.stderr)
        write_tsv(cache, [{"taxid": k, "count": v} for k, v in sorted(known.items())],
                  ["taxid", "count"])
    return known


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--offline", action="store_true",
                    help="skip the Entrez nt counts (cached values are still used)")
    args = ap.parse_args()

    cfg = cfg_bloodmeal()
    root = bloodmeal_root()
    meta = read_tsv(ROOT / "metadata" / "sample_metadata.tsv")
    fixes = {r["host_species"]: r for r in read_tsv(ROOT / cfg.get("host_fixes",
                                                                   "metadata/host_name_fixes.tsv"))}
    hosts: dict[str, int] = {}
    for r in meta:
        if not r["status"].startswith("excluded"):
            hosts[r["host_species"].strip()] = hosts.get(r["host_species"].strip(), 0) + 1

    def sci_name(h: str) -> tuple[str, str]:
        if h in fixes:
            return fixes[h]["resolved_name"].strip(), fixes[h]["note"]
        m = re.search(r"\(([^)]*)\)", h)
        name = m.group(1).strip() if m else h
        if name.lower() == "human":  # "Homo sapiens (Human)": the parenthesis is the common name
            name = h.split("(")[0].strip()
        return name, ""

    resolved = {h: sci_name(h) for h in hosts}
    names = lookup_names({n.lower() for n, _ in resolved.values() if n})
    tax = Taxonomy()

    rows = []
    for h, n in sorted(hosts.items(), key=lambda x: (-x[1], x[0])):
        name, note = resolved[h]
        row = {"host_species": h, "n_samples": n, "resolved_name": name, "fix_note": note}
        if not h:
            row["status"] = "no_host_recorded"
        elif not name:
            row["status"] = "no_host"  # questing
        elif name.lower() not in names:
            row["status"] = "UNRESOLVED_add_to_host_fixes"
        else:
            taxid, cls = names[name.lower()][0]
            lin = tax.lineage(taxid)
            row.update({"meta_taxid": lin["taxid"], "meta_rank": lin["rank"],
                        **{f"meta_{r}": lin[r] for r in RANKS},
                        "status": "resolved" if cls == "scientific name" else f"resolved_via_{cls}"})
        rows.append(row)

    # ---- reference availability -------------------------------------------
    refseq = read_tsv(root / "refs" / "vert_mito_refseq.tsv")
    sp_n, ge_n = {}, {}
    for r in refseq:
        for rank, d in (("species", sp_n), ("genus", ge_n)):
            a = tax.ancestor_at(r["taxid"], rank)
            if a:
                d[a] = d.get(a, 0) + 1
    want: dict[str, str] = {}
    for row in rows:
        if not row.get("meta_taxid"):
            continue
        s, g = tax.ancestor_at(row["meta_taxid"], "species"), tax.ancestor_at(row["meta_taxid"], "genus")
        row["refseq_mito_species"] = sp_n.get(s, 0) if s else ""
        row["refseq_mito_genus"] = ge_n.get(g, 0) if g else ""
        row["_s"], row["_g"] = str(s or ""), str(g or "")
        for t in (s, g):
            if t:
                want[str(t)] = tax.names.get(t, str(t))
    counts = nt_counts(want, root / "refs" / "nt_mito_counts.tsv", args.offline)
    for row in rows:
        row["nt_mito_species"] = counts.get(row.pop("_s", ""), "")
        row["nt_mito_genus"] = counts.get(row.pop("_g", ""), "")

    write_tsv(root / "tables" / "metadata_hosts.tsv", rows, COLS)
    bad = [r["host_species"] for r in rows if r["status"].startswith("UNRESOLVED")]
    if bad:
        print("UNRESOLVED host names (add them to metadata/host_name_fixes.tsv):\n  "
              + "\n  ".join(bad), file=sys.stderr)


if __name__ == "__main__":
    main()
