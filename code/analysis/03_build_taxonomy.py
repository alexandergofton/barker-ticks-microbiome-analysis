#!/usr/bin/env python3
"""Step 3 -- build the taxon lineage / annotation table (rowData for the TSE).

Collects every distinct lca_taxid across lca/*.blast.lca, resolves a full
lineage for each via taxopy, and classifies it into the analysis scope.

CLASSIFICATION RULES (verified against the local taxdump, not assumed)
  * Viruses have an 'acellular root' key and no 'domain' key.
  * Cellular taxa carry domain in {Bacteria, Archaea, Eukaryota}.
  * NCBI's current 'kingdom' for bacteria is a clade name (Pseudomonadati,
    Bacillati) and is useless for host/microbe splitting -- only the eukaryotic
    kingdom is used, where Metazoa / Viridiplantae mark non-microbial.
  * Microbial eukaryotes such as Babesia and Trypanosoma have NO kingdom key at
    all, so "Eukaryota and kingdom not in {Metazoa, Viridiplantae}" captures
    them correctly.
  * An LCA that sits above domain (root, 'cellular organisms') is neither host
    nor microbial -- it is genuinely ambiguous and is binned as such rather
    than silently assigned.

Output: results/taxonomy/taxon_lineage.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import taxopy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_config, load_contaminants, write_tsv  # noqa: E402

RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def collect_taxids(lca_dir: Path) -> dict[int, tuple[str, str]]:
    """taxid -> (name, rank) as recorded in the LCA files themselves."""
    seen: dict[int, tuple[str, str]] = {}
    files = sorted(p for p in lca_dir.glob("*.blast.lca") if not p.name.endswith(".kreport"))
    for path in files:
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                try:
                    t = int(f[1])
                except ValueError:
                    continue
                if t not in seen:
                    seen[t] = (f[2], f[3])
    print(f"  {len(files)} LCA files, {len(seen):,} distinct taxids", file=sys.stderr)
    return seen


def main() -> None:
    cfg = load_config()
    scope = cfg["scope"]
    host_phyla = set(scope["host_phyla"])
    blood_phyla = set(scope["bloodmeal_phyla"])
    excl_kingdoms = set(scope["exclude_euk_kingdoms"])
    contaminants = load_contaminants()

    print("Loading taxonomy database...", file=sys.stderr)
    db = taxopy.TaxDb(
        nodes_dmp=str(ROOT / "taxdump" / "nodes.dmp"),
        names_dmp=str(ROOT / "taxdump" / "names.dmp"),
        merged_dmp=str(ROOT / "taxdump" / "merged.filtered.dmp"),
    )

    print("Collecting taxids from LCA files...", file=sys.stderr)
    seen = collect_taxids(ROOT / "lca")

    rows = []
    n_unresolved = 0
    for taxid, (lca_name, lca_rank) in sorted(seen.items()):
        row = {"taxid": taxid, "lca_name": lca_name, "lca_rank": lca_rank}

        if taxid == 0:
            d = {}
            is_virus = False
        else:
            try:
                tx = taxopy.Taxon(taxid, db)
                d = tx.rank_name_dictionary
                # taxopy reports the name as-recorded; prefer it over the LCA file
                row["lca_name"] = tx.name
                row["lca_rank"] = tx.rank
            except Exception:
                d = {}
                n_unresolved += 1
            is_virus = "acellular root" in d

        for r in RANKS:
            row[r] = d.get(r, "")

        domain = row["domain"]
        if is_virus:
            group = "Viruses"
        elif domain in ("Bacteria", "Archaea", "Eukaryota"):
            group = domain
        else:
            group = "unclassified"
        row["superkingdom_group"] = group

        phylum, kingdom = row["phylum"], row["kingdom"]
        row["is_host"] = int(phylum in host_phyla)
        row["is_bloodmeal"] = int(phylum in blood_phyla)
        row["is_microbial"] = int(
            group in ("Bacteria", "Archaea", "Viruses")
            or (group == "Eukaryota" and kingdom not in excl_kingdoms)
        )
        # Everything that is neither host, bloodmeal nor microbial splits into
        # two genuinely different cases, which must not be conflated:
        #   other_eukaryote -- resolved, but out of scope (plants, algae, other
        #                      Metazoa). Real signal: environmental/dietary carryover.
        #   ambiguous       -- the LCA sits above the level where the host/microbe
        #                      distinction can be made at all (root, 'cellular
        #                      organisms', Pancrustacea). Not a taxon assignment.
        row["is_other_euk"] = int(
            group == "Eukaryota"
            and not row["is_host"] and not row["is_bloodmeal"] and not row["is_microbial"]
        )
        row["is_ambiguous"] = int(
            not row["is_host"] and not row["is_bloodmeal"]
            and not row["is_microbial"] and not row["is_other_euk"]
        )

        # Single mutually-exclusive label -- easier to group on than four booleans.
        row["scope_class"] = (
            "microbial" if row["is_microbial"] else
            "host" if row["is_host"] else
            "bloodmeal" if row["is_bloodmeal"] else
            "other_eukaryote" if row["is_other_euk"] else
            "ambiguous"
        )

        conf = contaminants.get(row["genus"], "")
        row["contaminant_flag"] = int(bool(conf))
        row["contaminant_confidence"] = conf

        rows.append(row)

    cols = (["taxid", "lca_name", "lca_rank"] + RANKS
            + ["superkingdom_group", "scope_class", "is_host", "is_bloodmeal",
               "is_microbial", "is_other_euk", "is_ambiguous",
               "contaminant_flag", "contaminant_confidence"])
    write_tsv(ROOT / "results" / "taxonomy" / "taxon_lineage.tsv", rows, cols)

    # ---- homonym report -----------------------------------------------------
    # NCBI reuses names across kingdoms: Bacillus is both a bacterial genus and
    # a stick-insect genus; Planococcus is both a bacterium and a mealybug.
    # Any code that keys on a rank NAME rather than a taxid can silently pick
    # the wrong lineage, so the collisions are written out and asserted against
    # downstream.
    from collections import Counter, defaultdict
    homonyms = []
    for rank in RANKS:
        by_name = defaultdict(set)
        for r in rows:
            if r[rank]:
                by_name[r[rank]].add(r["superkingdom_group"])
        for name, grps in by_name.items():
            if len(grps) > 1:
                homonyms.append({"rank": rank, "name": name,
                                 "superkingdom_groups": ";".join(sorted(grps))})
    write_tsv(ROOT / "results" / "taxonomy" / "homonyms.tsv",
              homonyms or [{"rank": "", "name": "", "superkingdom_groups": "none"}],
              ["rank", "name", "superkingdom_groups"])
    if homonyms:
        print(f"\ncross-kingdom name collisions: {len(homonyms)}", file=sys.stderr)
        for h in homonyms:
            print(f"  {h['rank']}: {h['name']} -> {h['superkingdom_groups']}", file=sys.stderr)

    # ---- summary ------------------------------------------------------------
    grp = Counter(r["superkingdom_group"] for r in rows)
    print("\ndistinct taxa by group", file=sys.stderr)
    for g, n in grp.most_common():
        print(f"  {g:<16} {n:>7,}", file=sys.stderr)
    print("\ndistinct taxa by scope_class", file=sys.stderr)
    for g, n in Counter(r["scope_class"] for r in rows).most_common():
        print(f"  {g:<18} {n:>7,}", file=sys.stderr)
    print(f"\n  contaminant-flagged {sum(r['contaminant_flag'] for r in rows):>7,}",
          file=sys.stderr)
    if n_unresolved:
        print(f"  taxids not in taxdump {n_unresolved:>7,}", file=sys.stderr)


if __name__ == "__main__":
    main()
