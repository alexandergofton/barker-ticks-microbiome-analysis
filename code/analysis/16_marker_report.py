#!/usr/bin/env python3
"""Step 16 -- what was recovered, and just as importantly, what was not.

  marker_summary.tsv     one row per target x marker: sequences, samples,
                         length and identity distribution, alignment path
  marker_by_sample.tsv   sample x target x marker, long form, with the tick
                         species and the step-10 presence call attached
  marker_gaps.tsv        every target x marker in the panel that produced
                         NOTHING, split by whether the target had contigs at all

marker_gaps.tsv is the point of this step. A marker gene missing from a sample
is the normal case, not a failure: whether gltA is present depends on what the
assembly happened to reconstruct. The table separates the two reasons a locus is
absent -- the target had no contigs in that sample, or it had contigs and the
locus was not among them -- because only the second says anything biological.
"""

from __future__ import annotations

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import load_panel, out_root, read_fasta, slug  # noqa: E402


def read_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def med(xs: list[float]) -> float:
    return round(statistics.median(xs), 1) if xs else 0.0


def main() -> None:
    root = out_root()
    seqs = read_rows(root / "marker_sequences.tsv")
    contig_meta = read_rows(root / "contigs" / "contig_metadata.tsv")
    if not seqs:
        sys.exit("ERROR: results/markers/marker_sequences.tsv is empty -- run step 14.")

    panel = load_panel()
    alndir = root / "alignments"

    # ---- per target x marker ------------------------------------------------
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in seqs:
        grouped[(r["target"], r["marker"])].append(r)

    summary = []
    for row in panel:
        g = grouped.get((row["target"], row["marker"]), [])
        good = [r for r in g if not r["flagged"]]
        lens = [int(r["seq_len_nt"]) for r in good]
        aln = alndir / slug(row["target"]) / f"{slug(row['target'])}__{slug(row['marker'])}.aln.fasta"
        cols = 0
        if aln.exists():
            v = read_fasta(aln)
            cols = len(next(iter(v.values()))) if v else 0
        summary.append({
            "target": row["target"], "marker": row["marker"],
            "marker_type": row["marker_type"],
            "n_sequences": len(good), "n_samples": len({r["sample_id"] for r in good}),
            "n_flagged_internal_stops": len(g) - len(good),
            "median_len_nt": med(lens), "min_len_nt": min(lens) if lens else 0,
            "max_len_nt": max(lens) if lens else 0,
            "median_pident_to_ref": med([float(r["best_pident"]) for r in good]),
            "median_ref_cov_pct": med([float(r["ref_cov_pct"]) for r in good]),
            "alignment_cols": cols,
            "alignment": str(aln.relative_to(root)) if aln.exists() else "",
        })
    summary.sort(key=lambda r: (-r["n_sequences"], r["target"], r["marker"]))
    write_tsv(root / "marker_summary.tsv", summary,
              ["target", "marker", "marker_type", "n_sequences", "n_samples",
               "n_flagged_internal_stops", "median_len_nt", "min_len_nt",
               "max_len_nt", "median_pident_to_ref", "median_ref_cov_pct",
               "alignment_cols", "alignment"])

    # ---- per sample ---------------------------------------------------------
    by_sample = sorted(seqs, key=lambda r: (r["sample_id"], r["target"], r["marker"]))
    write_tsv(root / "marker_by_sample.tsv",
              [{k: r[k] for k in ("sample_id", "tick_species", "tick_species_source",
                                  "target", "marker", "seq_id", "seq_len_nt",
                                  "best_pident", "ref_cov_pct", "n_internal_stops",
                                  "target_present_in_sample", "target_evidence",
                                  "lca_name", "flagged")} for r in by_sample],
              ["sample_id", "tick_species", "tick_species_source", "target", "marker",
               "seq_id", "seq_len_nt", "best_pident", "ref_cov_pct",
               "n_internal_stops", "target_present_in_sample", "target_evidence",
               "lca_name", "flagged"])

    # ---- what is missing, and why -------------------------------------------
    contigs_per_target: dict[str, int] = defaultdict(int)
    samples_per_target: dict[str, set] = defaultdict(set)
    for r in contig_meta:
        contigs_per_target[r["target"]] += 1
        samples_per_target[r["target"]].add(r["sample_id"])

    gaps = []
    for row in panel:
        t, mk = row["target"], row["marker"]
        if grouped.get((t, mk)):
            continue
        n_contigs = contigs_per_target.get(t, 0)
        gaps.append({
            "target": t, "marker": mk, "marker_type": row["marker_type"],
            "n_target_contigs": n_contigs,
            "n_samples_with_target": len(samples_per_target.get(t, ())),
            "reason": ("no contigs assigned to this target at all"
                       if n_contigs == 0 else
                       "target has contigs, but this locus was not among them"),
        })
    gaps.sort(key=lambda r: (-r["n_target_contigs"], r["target"], r["marker"]))
    write_tsv(root / "marker_gaps.tsv", gaps,
              ["target", "marker", "marker_type", "n_target_contigs",
               "n_samples_with_target", "reason"])

    # ---- console ------------------------------------------------------------
    recovered = [s for s in summary if s["n_sequences"]]
    print(f"\n{len(recovered)} of {len(panel)} panel loci recovered at least one "
          f"sequence\n", file=sys.stderr)
    print(f"  {'target':<20} {'marker':<14} {'seqs':>5} {'samples':>8} "
          f"{'med bp':>8} {'med %id':>8} {'aln cols':>9}", file=sys.stderr)
    for s in recovered:
        print(f"  {s['target']:<20} {s['marker']:<14} {s['n_sequences']:>5} "
              f"{s['n_samples']:>8} {s['median_len_nt']:>8.0f} "
              f"{s['median_pident_to_ref']:>8.1f} {s['alignment_cols']:>9}",
              file=sys.stderr)
    n_flag = sum(s["n_flagged_internal_stops"] for s in summary)
    if n_flag:
        print(f"\n{n_flag} sequences carry internal stop codons and are held in "
              f"*.flagged.fna, out of the alignments. Inspect before use: they are "
              f"usually assembly frameshifts or pseudogenes.", file=sys.stderr)
    has_contigs = sum(1 for g in gaps if g["n_target_contigs"])
    print(f"\n{len(gaps)} panel loci produced nothing "
          f"({has_contigs} where the target DID have contigs -- see marker_gaps.tsv)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
