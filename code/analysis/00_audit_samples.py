#!/usr/bin/env python3
"""Step 0 -- reconcile every sample across every pipeline artefact.

The pipeline was run incrementally and the artefact directories disagree about
which samples exist. This script produces the authoritative inventory that all
later steps consume, so nothing downstream has to guess.

Outputs
  results/audit/sample_completeness.tsv   one row per canonical sample
  results/audit/read_dir_priority.tsv     resolved read paths for samples needing mapping
  results/audit/audit_summary.txt         human-readable reconciliation
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, canonical, find_reads, load_aliases, write_tsv  # noqa: E402


def stems(directory: Path, suffixes: list[str]) -> dict[str, str]:
    """Map canonical sample_id -> the on-disk stem that produced it."""
    out: dict[str, str] = {}
    if not directory.exists():
        return out
    for entry in directory.iterdir():
        for suf in suffixes:
            if entry.name.endswith(suf):
                out[entry.name[: -len(suf)]] = entry.name
                break
    return out


def main() -> None:
    aliases = load_aliases()
    audit_dir = ROOT / "results" / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    # ---- inventory each artefact directory by raw file stem -----------------
    raw = {
        "contigs":       stems(ROOT / "contigs", [".contigs.fa"]),
        "contigs_noread": stems(ROOT / "contigs_with_missing_reads", [".contigs.fa"]),
        "blast":         stems(ROOT / "blast", [".blast.txt"]),
        "lca":           stems(ROOT / "lca", [".blast.lca"]),
        "kreport":       stems(ROOT / "lca", [".blast.lca.kreport"]),
        "coverage":      stems(ROOT / "mapping", [".coverage.tsv"]),
        "bam":           stems(ROOT / "mapping", [".bam"]),
        "bowtie2log":    stems(ROOT / "mapping", [".bowtie2.log"]),
    }
    # `.blast.lca` also matches the `.blast.lca.kreport` prefix strip above; drop overlap
    raw["lca"] = {k: v for k, v in raw["lca"].items() if not k.endswith(".kreport")}

    # ---- manifest ----------------------------------------------------------
    manifest: set[str] = set()
    man_path = ROOT / "final_sam_man.csv"
    if man_path.exists():
        with open(man_path, encoding="utf-8-sig") as fh:
            next(fh, None)
            for line in fh:
                sid = line.split(",")[0].strip()
                if sid:
                    manifest.add(canonical(sid, aliases))

    # ---- collapse to canonical ids -----------------------------------------
    canon: dict[str, dict[str, str]] = {}
    for artefact, mapping in raw.items():
        for stem in mapping:
            sid = canonical(stem, aliases)
            canon.setdefault(sid, {})[artefact] = stem

    all_ids = sorted(set(canon) | manifest)

    rows = []
    for sid in all_ids:
        got = canon.get(sid, {})
        has = {k: k in got for k in raw}
        blast_bytes = 0
        if has["blast"]:
            blast_bytes = (ROOT / "blast" / f"{got['blast']}.blast.txt").stat().st_size

        r1, r2, read_src = find_reads(sid)
        # samples whose contigs live in contigs_with_missing_reads/ are known read-less
        has_contigs = has["contigs"] or has["contigs_noread"]

        # ---- verdict --------------------------------------------------------
        if has["lca"] and has["coverage"]:
            status, action = "complete", "none"
        elif has["lca"] and not has["coverage"]:
            if r1:
                status, action = "needs_mapping", "run 02_backfill_mapping"
            else:
                status, action = "excluded_no_reads", "terminal: reads absent, cannot map"
        elif not has["lca"] and has["blast"]:
            status = "needs_lca"
            action = ("run 01_backfill_lca (large BLAST, streaming required)"
                      if blast_bytes > 10 * 1024**3 else "run 01_backfill_lca")
        elif not has["lca"] and not has["blast"]:
            status, action = "excluded_no_blast", "terminal: no BLAST output"
        else:
            status, action = "unknown", "review manually"

        rows.append({
            "sample_id": sid,
            "in_manifest": int(sid in manifest),
            "has_contigs": int(has_contigs),
            "contigs_have_reads": int(has["contigs"]),
            "has_blast": int(has["blast"]),
            "blast_gb": round(blast_bytes / 1024**3, 3),
            "has_lca": int(has["lca"]),
            "has_kreport": int(has["kreport"]),
            "has_coverage": int(has["coverage"]),
            "has_bam": int(has["bam"]),
            "reads_found": int(r1 is not None),
            "read_source_dir": read_src or "",
            "alias_applied": int(any(s != sid for s in got.values())),
            "status": status,
            "action": action,
        })

    cols = list(rows[0].keys())
    write_tsv(audit_dir / "sample_completeness.tsv", rows, cols)

    # ---- read paths for the mapping backfill -------------------------------
    map_rows = []
    for r in rows:
        if r["status"] != "needs_mapping":
            continue
        sid = r["sample_id"]
        r1, r2, src = find_reads(sid)
        stem = canon[sid].get("contigs") or canon[sid].get("contigs_noread")
        sub = "contigs" if canon[sid].get("contigs") else "contigs_with_missing_reads"
        # Paths are written RELATIVE to the project root. The same filesystem
        # is mounted at different absolute paths on the login node
        # (/data/hb-austick/...) and the cluster (/datasets/work/hb-austicks/...),
        # so absolute paths here would break whichever one did not write them.
        # The backfill scripts cd to the project root before reading this file.
        map_rows.append({
            "sample_id": sid,
            "contigs": f"{sub}/{stem}.contigs.fa",
            "r1": str(Path(r1).relative_to(ROOT)),
            "r2": str(Path(r2).relative_to(ROOT)),
            "read_source_dir": src,
        })
    write_tsv(audit_dir / "read_dir_priority.tsv", map_rows,
              ["sample_id", "contigs", "r1", "r2", "read_source_dir"])

    # ---- samples needing LCA ------------------------------------------------
    # Relative for the same reason as above.
    lca_rows = [{
        "sample_id": r["sample_id"],
        "blast": f"blast/{canon[r['sample_id']]['blast']}.blast.txt",
        "blast_gb": r["blast_gb"],
    } for r in rows if r["status"] == "needs_lca"]
    lca_rows.sort(key=lambda x: -x["blast_gb"])
    write_tsv(audit_dir / "lca_backfill_queue.tsv", lca_rows,
              ["sample_id", "blast", "blast_gb"])

    # ---- summary ------------------------------------------------------------
    counts = Counter(r["status"] for r in rows)
    lines = [
        "Sample reconciliation audit",
        "=" * 60,
        f"canonical samples resolved : {len(rows)}",
        f"present in final_sam_man   : {len(manifest)}",
        "",
        "status breakdown",
        "-" * 60,
    ]
    for status, n in counts.most_common():
        lines.append(f"  {status:<22} {n:>4}")
    lines += [
        "",
        f"projected N after backfill : "
        f"{counts['complete'] + counts['needs_mapping'] + counts['needs_lca']}",
        f"terminal exclusions        : "
        f"{counts['excluded_no_reads'] + counts['excluded_no_blast']}",
        "",
        "aliases applied",
        "-" * 60,
    ]
    for stem, sid in aliases.items():
        lines.append(f"  {stem}  ->  {sid}")
    lines += ["", "terminal exclusions (no abundance data possible)", "-" * 60]
    for r in rows:
        if r["status"].startswith("excluded"):
            lines.append(f"  {r['sample_id']:<34} {r['status']}")
    lines += ["", "LCA backfill queue (largest BLAST first)", "-" * 60]
    for r in lca_rows:
        lines.append(f"  {r['sample_id']:<34} {r['blast_gb']:>8.1f} GB")

    text = "\n".join(lines)
    (audit_dir / "audit_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
