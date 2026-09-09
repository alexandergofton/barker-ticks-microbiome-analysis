#!/usr/bin/env python3
"""Step 13 -- locate each marker locus inside the target contigs.

rRNA loci are found with blastn, protein-coding loci with tblastn, using the
reference probes cached by step 11. Raw HSPs are then assembled into loci:

  1. HSPs failing the panel's min_aln_len / min_pident are dropped.
  2. Remaining HSPs for one marker on one contig and one strand are merged into
     a locus when they lie closer together than `merge_gap_*`. tblastn routinely
     breaks a single divergent gene into several HSPs, so without this a gene
     would be reported as three short fragments.
  3. Loci of DIFFERENT markers overlapping on the same contig are resolved
     competitively -- highest summed bitscore takes the region, the loser is
     recorded as `paralog_suppressed` with the marker that beat it.

Step 3 is the one that matters biologically. Rickettsia ompA, ompB and sca4 are
all Sca-family autotransporters and cross-hit each other's probes; unresolved,
the same stretch of contig would be exported as three different "genes" and
each would produce a tree. The same applies to the Anaplasma msp2/p44 paralogue
family and to the viral polyprotein vs its own NS5 sub-region.

Outputs
  results/markers/hits/<target>__<marker>.blast.tsv   raw HSPs, kept for audit
  results/markers/marker_loci.tsv                     assembled loci, all statuses
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import cfg_markers, load_panel, out_root, require, run, slug  # noqa: E402

FMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
       "evalue bitscore qlen slen")
COLS = FMT.split()[1:]

LOCI_COLS = ["target", "marker", "marker_type", "seq_name", "contig_len",
             "locus_start", "locus_end", "strand", "locus_len", "n_hsps",
             "sum_bitscore", "max_bitscore", "best_pident", "best_ref",
             "ref_len", "ref_cov_pct", "status", "suppressed_by"]


def parse_hits(path: Path, min_len: int, min_pid: float) -> list[dict]:
    """Read a BLAST tabular file, keeping HSPs that clear the panel floors."""
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(COLS):
                continue
            h = dict(zip(COLS, p))
            h["pident"] = float(h["pident"]); h["length"] = int(h["length"])
            h["bitscore"] = float(h["bitscore"]); h["evalue"] = float(h["evalue"])
            h["sstart"] = int(h["sstart"]); h["send"] = int(h["send"])
            h["qlen"] = int(h["qlen"]); h["slen"] = int(h["slen"])
            if h["length"] < min_len or h["pident"] < min_pid:
                continue
            out.append(h)
    return out


def assemble_loci(hits: list[dict], merge_gap: int, max_per_contig: int,
                  row: dict) -> list[dict]:
    """Merge HSPs into loci, per contig and strand."""
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for h in hits:
        strand = "+" if h["sstart"] <= h["send"] else "-"
        groups[(h["sseqid"], strand)].append(h)

    per_contig: dict[str, list[dict]] = defaultdict(list)
    for (contig, strand), hs in groups.items():
        spans = sorted(((min(h["sstart"], h["send"]), max(h["sstart"], h["send"]), h)
                        for h in hs), key=lambda s: (s[0], s[1]))
        cur = [spans[0]]
        merged = []
        for s in spans[1:]:
            if s[0] - cur[-1][1] <= merge_gap:
                cur.append(s)
            else:
                merged.append(cur); cur = [s]
        merged.append(cur)

        for grp in merged:
            lo = min(g[0] for g in grp)
            hi = max(g[1] for g in grp)
            hsps = [g[2] for g in grp]
            top = max(hsps, key=lambda h: h["bitscore"])
            # Reference coverage uses the merged span against the probe's own
            # length, scaled by 3 for tblastn where the subject is nucleotide and
            # the query is protein. It answers "how much of the gene is here?",
            # which is what decides whether a locus is worth a tree.
            unit = 3 if row["marker_type"] == "protein" else 1
            cov = 100.0 * (hi - lo + 1) / max(1, top["qlen"] * unit)
            per_contig[contig].append({
                "target": row["target"], "marker": row["marker"],
                "marker_type": row["marker_type"], "seq_name": contig,
                "contig_len": top["slen"], "locus_start": lo, "locus_end": hi,
                "strand": strand, "locus_len": hi - lo + 1, "n_hsps": len(hsps),
                "sum_bitscore": round(sum(h["bitscore"] for h in hsps), 1),
                "max_bitscore": round(top["bitscore"], 1),
                "best_pident": round(max(h["pident"] for h in hsps), 1),
                "best_ref": top["qseqid"], "ref_len": top["qlen"],
                "ref_cov_pct": round(min(cov, 999.9), 1),
                "status": "kept", "suppressed_by": "",
            })

    out = []
    for contig, loci in per_contig.items():
        loci.sort(key=lambda x: -x["sum_bitscore"])
        out.extend(loci[:max_per_contig])
    return out


def resolve_paralogs(loci: list[dict], min_overlap: float) -> list[dict]:
    """Greedy, highest-bitscore-first claim of each contig region."""
    by_contig: dict[str, list[dict]] = defaultdict(list)
    for l in loci:
        by_contig[l["seq_name"]].append(l)

    for contig, ls in by_contig.items():
        ls.sort(key=lambda x: (-x["max_bitscore"], -x["sum_bitscore"]))
        claimed: list[dict] = []
        for l in ls:
            for c in claimed:
                ov = min(l["locus_end"], c["locus_end"]) - max(l["locus_start"], c["locus_start"]) + 1
                if ov <= 0:
                    continue
                # Fraction of the SHORTER locus: a short probe hit buried inside
                # a long one is the paralogue case we are trying to catch.
                if ov / min(l["locus_len"], c["locus_len"]) >= min_overlap:
                    if l["marker"] != c["marker"]:
                        l["status"] = "paralog_suppressed"
                        l["suppressed_by"] = f"{c['marker']}@bits{c['max_bitscore']}"
                    else:
                        l["status"] = "duplicate_suppressed"
                        l["suppressed_by"] = c["marker"]
                    break
            if l["status"] == "kept":
                claimed.append(l)
    return loci


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    ap.add_argument("--force", action="store_true", help="re-run BLAST for cached hits")
    args = ap.parse_args()

    mc = cfg_markers()
    sc = mc["search"]
    root = out_root()
    contigs_dir, hits_dir, db_dir = root / "contigs", root / "hits", root / "blastdb"
    hits_dir.mkdir(parents=True, exist_ok=True)
    db_dir.mkdir(parents=True, exist_ok=True)

    makeblastdb, blastn, tblastn = require("makeblastdb"), require("blastn"), require("tblastn")

    panel = load_panel()
    if args.target:
        panel = [r for r in panel if r["target"] in set(args.target)]

    by_target: dict[str, list[dict]] = defaultdict(list)
    for r in panel:
        by_target[r["target"]].append(r)

    all_loci: list[dict] = []
    for target, rows in sorted(by_target.items()):
        fa = contigs_dir / f"{slug(target)}.contigs.fa"
        if not fa.exists():
            print(f"{target}: no contigs (step 12 found none) -- skipped", file=sys.stderr)
            continue

        db = db_dir / slug(target)
        if args.force or not (db.with_suffix(".nin").exists() or Path(f"{db}.nin").exists()):
            run([makeblastdb, "-in", str(fa), "-dbtype", "nucl", "-out", str(db)])

        target_loci: list[dict] = []
        for row in rows:
            ext = "fna" if row["marker_type"] == "rrna" else "faa"
            ref = root / "refs" / f"{row['key']}.{ext}"
            if not ref.exists():
                continue
            hit_file = hits_dir / f"{row['key']}.blast.tsv"
            if args.force or not hit_file.exists():
                prog = blastn if row["marker_type"] == "rrna" else tblastn
                cmd = [prog, "-query", str(ref), "-db", str(db),
                       "-out", str(hit_file), "-outfmt", FMT,
                       "-evalue", str(sc["evalue_rrna"] if row["marker_type"] == "rrna"
                                      else sc["evalue_protein"]),
                       "-num_threads", str(sc["threads"]),
                       "-max_target_seqs", "100000"]
                if row["marker_type"] == "rrna":
                    cmd += ["-task", "blastn"]   # not megablast: cross-species rRNA
                run(cmd)

            hits = parse_hits(hit_file, row["min_aln_len"], row["min_pident"])
            gap = sc["merge_gap_rrna"] if row["marker_type"] == "rrna" else sc["merge_gap_protein"]
            loci = assemble_loci(hits, gap, sc["max_loci_per_contig_marker"], row)
            target_loci.extend(loci)
            if loci:
                print(f"  {target:<18} {row['marker']:<14} "
                      f"{len(hits):>7,} HSPs -> {len(loci):>5,} loci", file=sys.stderr)

        target_loci = resolve_paralogs(target_loci, sc["paralog_overlap"])
        kept = sum(1 for l in target_loci if l["status"] == "kept")
        sup = len(target_loci) - kept
        print(f"{target:<18} {kept:>6,} loci kept, {sup:,} suppressed", file=sys.stderr)
        all_loci.extend(target_loci)

    # As in step 12: a --target run must leave the other targets' loci in place,
    # because step 14 reads this table whole.
    loci_path = root / "marker_loci.tsv"
    if args.target and loci_path.exists():
        scope = set(args.target)
        with open(loci_path) as fh:
            others = [r for r in csv.DictReader(fh, delimiter="\t")
                      if r["target"] not in scope]
        for r in others:
            r["sum_bitscore"] = float(r["sum_bitscore"])
        print(f"  merging with {len(others):,} existing loci from other targets",
              file=sys.stderr)
        all_loci = others + all_loci

    all_loci.sort(key=lambda l: (l["target"], l["marker"], l["seq_name"],
                                 -float(l["sum_bitscore"])))
    write_tsv(loci_path, all_loci, LOCI_COLS)
    kept = sum(1 for l in all_loci if l["status"] == "kept")
    print(f"\n{kept:,} loci kept of {len(all_loci):,} candidates", file=sys.stderr)


if __name__ == "__main__":
    main()
