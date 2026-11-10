#!/usr/bin/env python3
"""Step 34 -- one blood-meal call per sample, and how it compares with the
host recorded in the metadata.

Evidence (step 33, status ok only; mtDNA contigs carry step 33c's refined
name where it has one), used hierarchically:
  tier 1  vertebrate mitochondrial genome / segment contigs
  tier 2  ONLY if tier 1 is empty: vertebrate 18S rRNA contigs and short COI
          fragments
A sample with neither is `no_bloodmeal_detected`. When contigs exist but no
taxon clears the read floor, the sub-floor taxa are used and the sample is
marked `evidence = weak` (strong / weak / none).

Host taxa: the tier's contigs are grouped by their identity-capped name
(call_taxid); each taxon must clear `bloodmeal.call.min_contigs` and
`min_reads`. Taxa in `lab_contaminants` (human, chicken, cow, pig) and taxa
flagged as lane cross-talk (`cross_talk_ratio`: a lane-mate carries >= 100x
the reads of the same family) are reported but left out of the primary call
unless they are the only vertebrate signal -- human is also a genuine
recorded host for 4 samples, so it is never discarded.

Primary call: walking species -> genus -> family -> order -> class, the first
rank at which a single taxon holds >= `dominance` of the read weight resolved
to that rank plus the weight of coarser contigs that contradict it. A contig named only to family therefore never blocks a species
call made by other contigs, and two genuinely different species sharing the
weight push the call up to their common genus/family.

Other hosts (`hosts`): each taxon's reads plus those of every taxon nested
under it; those holding >= `secondary_min_share` are listed, keeping only the
deepest of a nested chain. `multiple_families` = yes when two or more non-lab families each pass
that share: a mixed or sequential blood meal, a pooled tube of ticks fed on
different hosts, or cross-sample contamination.

Comparison with metadata (resolved by step 32), rank by rank from class down:
  MATCH_<rank>     agree at every rank both resolve; <rank> is the deepest
                   rank compared (MATCH_family = call only reached family)
  MISMATCH_<rank>  first rank at which they differ
  lab_or_crosstalk_only  the only signal is lab DNA or lane cross-talk and
                   it disagrees with the recorded host (not evidence against it)
  no_host_recorded / questing / no_bloodmeal_detected
`metadata_host_in_hosts` says whether the recorded host is consistent with ANY
reported host taxon, not just the primary.

Outputs (results/bloodmeal/tables/)
  bloodmeal_host_taxa.tsv     sample x host taxon, all tiers, with floors applied
  bloodmeal_calls.tsv         one row per sample
  bloodmeal_concordance.tsv   counts of comparison outcome, overall and by
                              metadata host class / tick genus / life stage
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from bloodmeal_lib import RANKS, bloodmeal_root, cfg_bloodmeal, read_tsv  # noqa: E402

TIER_ORDER = [("tier1", {"tier1_mito"}), ("tier2", {"tier2_COI", "tier2_18S"})]

TAXA_COLS = ["sample_id", "tier", "call_taxid", "call_rank", "call_name",
             *[f"call_{r}" for r in RANKS], "n_contigs", "mapped_reads", "share",
             "best_pident", "max_contig_len", "genes", "lab_contaminant", "cross_talk", "passes_floor"]

CALL_COLS = ["sample_id", "tick_genus", "tick_species", "life_stage", "n_ticks_pooled",
             "host_type", "host_species", "meta_rank", *[f"meta_{r}" for r in RANKS],
             "meta_nt_mito_species",
             "evidence", "tier_used", "n_tier1_contigs", "n_tier2_contigs", "tier_reads",
             "primary_rank", "primary_name", *[f"primary_{r}" for r in RANKS],
             "primary_share", "primary_best_pident", "primary_genes",
             "hosts", "n_host_families", "multiple_families", "lab_taxa_detected",
             "cross_talk_taxa",
             "comparison", "agree_to_rank", "metadata_host_in_hosts", "note"]


def lineage_of(r: dict, prefix: str) -> dict[str, str]:
    return {k: r.get(f"{prefix}_{k}", "") for k in RANKS}


def consistent(a: dict[str, str], b: dict[str, str]) -> bool:
    """No rank at which both are resolved and differ."""
    return all(not a[r] or not b[r] or a[r] == b[r] for r in RANKS)


def compare(call: dict[str, str], meta: dict[str, str]) -> tuple[str, str]:
    deepest = ""
    for r in RANKS:
        if call[r] and meta[r]:
            if call[r] != meta[r]:
                return f"MISMATCH_{r}", deepest
            deepest = r
    return (f"MATCH_{deepest}", deepest) if deepest else ("not_comparable", "")


def primary(taxa: list[dict], dominance: float) -> tuple[dict | None, float]:
    """(lineage of the primary call, its share) -- see module docstring.

    At each rank the leading taxon X is weighed against every contig resolved
    to that rank AND every coarser contig that contradicts X at a rank it does
    resolve (a Peramelidae contig counts against a Varanus species call; a
    Macropodidae contig does not count against Macropus)."""
    for r in reversed(RANKS):
        w: Counter = Counter()
        for t in taxa:
            if t[f"call_{r}"]:
                w[t[f"call_{r}"]] += t["mapped_reads"]
        if not w:
            continue
        name, top = w.most_common(1)[0]
        src = next(t for t in taxa if t[f"call_{r}"] == name)
        coarser = RANKS[:RANKS.index(r)]
        against = sum(t["mapped_reads"] for t in taxa if not t[f"call_{r}"] and any(
            t[f"call_{q}"] and t[f"call_{q}"] != src[f"call_{q}"] for q in coarser))
        tot = sum(w.values()) + against
        if tot and top / tot >= dominance:
            lin = {k: src[f"call_{k}"] for k in RANKS[:RANKS.index(r) + 1]}
            lin.update({k: "" for k in RANKS[RANKS.index(r) + 1:]})
            lin["rank"], lin["name"] = r, name
            return lin, top / tot
    return None, 0.0


def is_ancestor(a: dict, b: dict) -> bool:
    """Taxon a lies on taxon b's lineage, above it. A node above class
    (Craniata, Vertebrata) is an ancestor of every ranked taxon."""
    if a is b:
        return False
    if a["call_rank"] not in RANKS:
        return b["call_rank"] in RANKS
    if b["call_rank"] not in RANKS:
        return False
    return (RANKS.index(b["call_rank"]) > RANKS.index(a["call_rank"])
            and b[f"call_{a['call_rank']}"] == a["call_name"])


def host_groups(taxa: list[dict], min_share: float) -> list[tuple[dict, int]]:
    """Reported hosts: each taxon carries its own reads plus those of every
    taxon nested under it; those >= min_share are kept, and of a nested chain
    only the deepest survives. (Liasis 41k reads + one L. mackloti contig 2k
    -> "Liasis", not a 2% "L. mackloti".)"""
    tot = sum(t["mapped_reads"] for t in taxa) or 1
    cum = {id(t): t["mapped_reads"] + sum(o["mapped_reads"] for o in taxa if is_ancestor(t, o))
           for t in taxa}
    sel = [t for t in taxa if cum[id(t)] / tot >= min_share]
    sel = [t for t in sel if not any(is_ancestor(t, o) for o in sel)]
    # coarser contigs that agree with exactly one reported host count towards it
    # (Peramelidae contigs behind an "Isoodon" call); shared ancestors (Mammalia
    # behind both Isoodon and Macropus) count towards neither
    out = []
    for t in sel:
        own = sum(a["mapped_reads"] for a in taxa if is_ancestor(a, t)
                  and sum(is_ancestor(a, o) for o in sel) == 1)
        out.append((t, cum[id(t)] + own))
    return sorted(out, key=lambda x: -x[1])


def main() -> None:
    cfg = cfg_bloodmeal()
    cc = cfg["call"]
    lab = set(cc.get("lab_contaminants", []))
    root = bloodmeal_root()
    contigs = read_tsv(root / "tables" / "vertebrate_contigs.tsv")
    hosts_meta = {r["host_species"]: r for r in read_tsv(root / "tables" / "metadata_hosts.tsv")}
    meta = {r["sample_id"]: r for r in read_tsv(ROOT / "metadata" / "sample_metadata.tsv")}
    scanned = {p.name[:-len(".tsv.gz")] for p in (root / "vert_hits").glob("*.tsv.gz")}

    # Step 33c's identity-vote name replaces the first-pass nt-LCA name for
    # every mtDNA contig it could re-name (see 33c for why).
    refined_p = root / "tables" / "refined_contig_calls.tsv"
    refined = {(r["sample_id"], r["contig"]): r for r in read_tsv(refined_p)
               if r.get("refined_taxid")} if refined_p.exists() else {}
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for c in contigs:
        if c["status"] == "ok":
            c["mapped_reads"] = int(c["mapped_reads"])
            rf = refined.get((c["sample_id"], c["contig"]))
            c["naming"] = "nt_lca"
            if rf and c["kind"] == "mtDNA":
                c.update({"call_taxid": rf["refined_taxid"], "call_rank": rf["refined_rank"],
                          "call_name": rf["refined_name"],
                          **{f"call_{r}": rf[f"refined_{r}"] for r in RANKS},
                          "best_pident": rf["refined_pident"], "naming": "refined"})
            by_sample[c["sample_id"]].append(c)

    # ---- lane cross-talk: reads per (lane, family) in every sample -----------
    def tkey(c: dict) -> str:
        return c["call_family"] or c["call_name"]
    lane_reads: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for sid, cs in by_sample.items():
        lane = meta.get(sid, {}).get("seq_run", "")
        if lane:
            for c in cs:
                lane_reads[(lane, tkey(c))][sid] += c["mapped_reads"]
    ratio = float(cc.get("cross_talk_ratio", 0) or 0)

    def cross_talk(sid: str, key: str, reads: int) -> str:
        lane = meta.get(sid, {}).get("seq_run", "")
        if not ratio or not lane:
            return ""
        others = [(n, s) for s, n in lane_reads.get((lane, key), {}).items() if s != sid]
        if not others:
            return ""
        n, s = max(others)
        return f"{s} ({n} reads)" if n >= ratio * max(reads, 1) else ""

    taxa_rows, call_rows = [], []
    for sid in sorted(scanned):
        cs = by_sample.get(sid, [])
        m = meta.get(sid, {})
        hm = hosts_meta.get(m.get("host_species", "").strip(), {})
        row = {k: m.get(k, "") for k in ("tick_genus", "tick_species", "life_stage",
                                         "n_ticks_pooled", "host_type", "host_species")}
        row.update({"sample_id": sid, "meta_rank": hm.get("meta_rank", ""),
                    **{f"meta_{r}": hm.get(f"meta_{r}", "") for r in RANKS},
                    "meta_nt_mito_species": hm.get("nt_mito_species", ""),
                    "n_tier1_contigs": sum(c["tier"] == "tier1_mito" for c in cs),
                    "n_tier2_contigs": sum(c["tier"] in ("tier2_COI", "tier2_18S") for c in cs)})

        # ---- host taxa, per tier ------------------------------------------
        tier_taxa: dict[str, list[dict]] = {}
        tier_all: dict[str, list[dict]] = {}
        for tier, tiers in TIER_ORDER:
            g: dict[str, list[dict]] = defaultdict(list)
            for c in cs:
                if c["tier"] in tiers:
                    g[c["call_taxid"]].append(c)
            tot = sum(c["mapped_reads"] for v in g.values() for c in v) or 1
            out = []
            for taxid, v in g.items():
                c0 = v[0]
                genes = {x for c in v for x in (c.get("refseq_genes") or c["genes_nt"]).split(",") if x}
                t = {"sample_id": sid, "tier": tier, "call_taxid": taxid,
                     "call_rank": c0["call_rank"], "call_name": c0["call_name"],
                     **{f"call_{r}": c0[f"call_{r}"] for r in RANKS},
                     "n_contigs": len(v), "mapped_reads": sum(c["mapped_reads"] for c in v),
                     "best_pident": max(float(c["best_pident"]) for c in v),
                     "max_contig_len": max(int(c["contig_len"]) for c in v),
                     "genes": ",".join(sorted(genes)),
                     "lab_contaminant": "yes" if c0["call_species"] in lab else ""}
                t["cross_talk"] = cross_talk(sid, tkey(c0), t["mapped_reads"])
                t["share"] = round(t["mapped_reads"] / tot, 4)
                t["passes_floor"] = "yes" if (t["n_contigs"] >= cc["min_contigs"]
                                              and t["mapped_reads"] >= cc["min_reads"]) else ""
                out.append(t)
            out.sort(key=lambda t: -t["mapped_reads"])
            taxa_rows += out
            tier_taxa[tier] = [t for t in out if t["passes_floor"]]
            tier_all[tier] = out

        used = next((t for t, _ in TIER_ORDER if tier_taxa[t]), "")
        row["evidence"] = "strong" if used else "none"
        taxa = tier_taxa.get(used, [])
        if not used:
            # Nothing clears the read floor: fall back to the sub-floor taxa of
            # the best tier, reported as a WEAK call rather than "not detected"
            # (on this cohort most of them agree with the recorded host).
            used = next((t for t, _ in TIER_ORDER if tier_all[t]), "")
            taxa = tier_all.get(used, [])
            row["evidence"] = "weak" if used else "none"
        row["tier_used"] = used or "none"
        non_lab = [t for t in taxa if not t["lab_contaminant"] and not t["cross_talk"]]
        basis = non_lab or taxa
        row["tier_reads"] = sum(t["mapped_reads"] for t in taxa)
        row["lab_taxa_detected"] = ",".join(t["call_name"] for t in taxa if t["lab_contaminant"])
        row["cross_talk_taxa"] = "; ".join(f"{t['call_name']} <- {t['cross_talk']}"
                                           for t in taxa if t["cross_talk"])

        # ---- primary + other hosts ----------------------------------------
        prim, share = primary(basis, float(cc["dominance"]))
        tot = sum(t["mapped_reads"] for t in taxa) or 1
        groups = host_groups(taxa, float(cc["secondary_min_share"]))
        listed = [t for t, _ in groups]
        row["hosts"] = "; ".join(f"{t['call_name']} ({t['call_rank']}, {n / tot:.0%}"
                                 f", {t['n_contigs']} ctg{', lab?' if t['lab_contaminant'] else ''}"
                                 f"{', cross-talk?' if t['cross_talk'] else ''})"
                                 for t, n in groups)
        fams = {t["call_family"] for t in listed
                if t["call_family"] and not t["lab_contaminant"] and not t["cross_talk"]}
        row["n_host_families"] = len(fams)
        row["multiple_families"] = "yes" if len(fams) > 1 else ""
        if prim:
            sup = [t for t in basis if t[f"call_{prim['rank']}"] == prim["name"]]
            genes = sorted({g for t in sup for g in t["genes"].split(",") if g})
            row.update({"primary_rank": prim["rank"], "primary_name": prim["name"],
                        **{f"primary_{r}": prim[r] for r in RANKS},
                        "primary_share": round(share, 3),
                        "primary_best_pident": max(t["best_pident"] for t in sup),
                        "primary_genes": ",".join(genes)})

        # ---- comparison with metadata -------------------------------------
        mlin = lineage_of(row, "meta")
        notes = []
        if not m.get("host_species", "").strip():
            row["comparison"] = "no_host_recorded"
        elif hm.get("status") == "no_host":
            row["comparison"] = "questing"
        elif not prim:
            row["comparison"] = "no_bloodmeal_detected" if not taxa else "no_dominant_host"
        else:
            row["comparison"], row["agree_to_rank"] = compare(lineage_of(row, "primary"), mlin)
            # Lab DNA / lane bleed is not evidence AGAINST the recorded host.
            if not non_lab and not row["comparison"].startswith("MATCH_"):
                row["comparison"], row["agree_to_rank"] = "lab_or_crosstalk_only", ""
        if any(mlin.values()) and taxa:
            row["metadata_host_in_hosts"] = "yes" if any(
                t["call_class"] and consistent(lineage_of(t, "call"), mlin)
                for t in listed or taxa) else "no"
        if row.get("comparison", "").startswith("MISMATCH_") and row["comparison"] in (
                "MISMATCH_species", "MISMATCH_genus"):
            n = hm.get("nt_mito_species", "")
            if n in ("0", "1", "2", "3", "4", "5"):
                notes.append(f"recorded host has only {n} nt mito record(s): the blood meal may "
                             "be matching its nearest sequenced relative")
        if row["evidence"] == "weak":
            notes.append(f"weak: no host taxon reaches {cc['min_reads']} mapped reads")
        if used == "tier2":
            notes.append("tier 2 only (18S/COI fragments): low taxonomic resolution")
        if taxa and not non_lab:
            notes.append("only lab-contaminant / lane cross-talk taxa detected")
        if hm.get("fix_note", "").startswith("CHECK"):
            notes.append("metadata host name inconsistent: " + hm["fix_note"])
        row["note"] = "; ".join(notes)
        call_rows.append(row)

    out = root / "tables"
    write_tsv(out / "bloodmeal_host_taxa.tsv", taxa_rows, TAXA_COLS)
    write_tsv(out / "bloodmeal_calls.tsv", call_rows, CALL_COLS)

    summ = []
    for by in ("all", "evidence", "meta_class", "tick_genus", "life_stage", "tier_used"):
        cnt = Counter((("all" if by == "all" else r.get(by) or "(blank)"), r["comparison"])
                      for r in call_rows)
        tot = Counter(k for k, _ in cnt.elements())
        summ += [{"grouping": by, "group": g, "comparison": c, "n": n, "n_group": tot[g]}
                 for (g, c), n in sorted(cnt.items())]
    write_tsv(out / "bloodmeal_concordance.tsv", summ,
              ["grouping", "group", "comparison", "n", "n_group"])
    c = Counter(r["comparison"] for r in call_rows)
    print("comparison outcomes:\n  " + "\n  ".join(f"{k}: {v}" for k, v in c.most_common()),
          file=sys.stderr)


if __name__ == "__main__":
    main()
