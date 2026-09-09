#!/usr/bin/env python3
"""Step 12 -- pull every target contig out of the per-sample assemblies.

Input   results/targets/target_contig_index.tsv   (written by step 10)
Output  results/markers/contigs/<target>.contigs.fa
        results/markers/contigs/contig_metadata.tsv

Each sequence is renamed on the way out to

    <sample_id>|<tick_species>|<target>|<contig_id>

so the sample and the tick it came from stay attached to the sequence through
every later step. `contig_metadata.tsv` carries the rest -- read evidence, TPM,
the LCA name, the presence call and evidence tier from step 10, and how the
tick species was assigned (curated metadata vs id-string guess).

SCOPE: every contig in the index is taken, including those from samples below
the presence call. A recovered marker gene is itself detection evidence, so
filtering on read counts here would throw away real detections; the call is
carried as annotation instead. See `markers:` in config/analysis_params.yaml.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import contig_name, out_root, slug, tick_species_map, write_fasta  # noqa: E402

INDEX = ROOT / "results" / "targets" / "target_contig_index.tsv"
PREVALENCE = ROOT / "results" / "targets" / "target_prevalence.tsv"


def load_prevalence() -> dict[tuple[str, str], tuple[str, str]]:
    """(target, sample_id) -> (present, evidence tier) from step 10."""
    out: dict[tuple[str, str], tuple[str, str]] = {}
    if not PREVALENCE.exists():
        return out
    with open(PREVALENCE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["target"], r["sample_id"])] = (r.get("present", ""),
                                                  r.get("evidence", ""))
    return out


def stream_needed(fasta: Path, wanted: set[str]) -> dict[str, str]:
    """Read `fasta` once, returning only the contigs in `wanted`.

    Assemblies run to gigabytes, so this streams and never holds a whole file:
    it keeps only the sequences asked for and stops as soon as it has them all.
    """
    out: dict[str, str] = {}
    keep, name, buf = False, None, []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and name is not None:
                    out[name] = "".join(buf)
                    if len(out) == len(wanted):
                        return out
                name = line[1:].split()[0]
                keep, buf = name in wanted, []
            elif keep:
                buf.append(line.strip())
    if keep and name is not None:
        out[name] = "".join(buf)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    args = ap.parse_args()

    if not INDEX.exists():
        sys.exit(f"ERROR: {INDEX.relative_to(ROOT)} not found -- run "
                 f"`python3 code/analysis/10_export_targets.py` first.")

    species = tick_species_map()
    prev = load_prevalence()
    outdir = out_root() / "contigs"
    outdir.mkdir(parents=True, exist_ok=True)

    # sample -> fasta path, and sample -> {contig_id -> [row, ...]}
    rows_by_sample: dict[str, list[dict]] = defaultdict(list)
    fasta_of: dict[str, str] = {}
    keep_targets = set(args.target) if args.target else None
    with open(INDEX) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if keep_targets and r["target"] not in keep_targets:
                continue
            rows_by_sample[r["sample_id"]].append(r)
            if r.get("contigs_fasta"):
                fasta_of[r["sample_id"]] = r["contigs_fasta"]

    n_rows = sum(len(v) for v in rows_by_sample.values())
    print(f"{n_rows:,} target contigs across {len(rows_by_sample)} samples",
          file=sys.stderr)

    by_target: dict[str, dict[str, str]] = defaultdict(dict)
    meta: list[dict] = []
    missing_fasta, missing_contig = [], 0

    for i, (sample, rows) in enumerate(sorted(rows_by_sample.items()), 1):
        path = fasta_of.get(sample, "")
        if not path or not Path(path).exists():
            missing_fasta.append(sample)
            continue
        wanted = {r["contig_id"] for r in rows}
        seqs = stream_needed(Path(path), wanted)
        missing_contig += len(wanted) - len(seqs)
        sp, sp_src = species.get(sample, ("unknown_species", "not_in_metadata"))

        for r in rows:
            seq = seqs.get(r["contig_id"])
            if seq is None:
                continue
            name = contig_name(sample, sp, r["target"], r["contig_id"])
            by_target[r["target"]][name] = seq
            present, evidence = prev.get((r["target"], sample), ("", ""))
            meta.append({
                "seq_name": name, "sample_id": sample, "tick_species": sp,
                "tick_species_source": sp_src, "target": r["target"],
                "contig_id": r["contig_id"], "contig_len": len(seq),
                "mapped_reads": r.get("mapped_reads", ""),
                "tpm_microbial": r.get("tpm_microbial", ""),
                "lca_taxid": r.get("lca_taxid", ""), "lca_name": r.get("lca_name", ""),
                "lca_genus": r.get("genus", ""), "lca_species": r.get("species", ""),
                "top_bitscore": r.get("top_bitscore", ""),
                "target_present_in_sample": present, "target_evidence": evidence,
            })
        if i % 25 == 0 or i == len(rows_by_sample):
            print(f"  [{i}/{len(rows_by_sample)}] {sample}", file=sys.stderr)

    for target, seqs in sorted(by_target.items()):
        write_fasta(outdir / f"{slug(target)}.contigs.fa", seqs)
        print(f"  {target:<20} {len(seqs):>6,} contigs", file=sys.stderr)

    # A --target run must not drop the other targets' rows: contig_metadata.tsv
    # is read whole by steps 14 and 16, so keep what this run did not touch.
    meta_path = outdir / "contig_metadata.tsv"
    if keep_targets and meta_path.exists():
        with open(meta_path) as fh:
            others = [r for r in csv.DictReader(fh, delimiter="\t")
                      if r["target"] not in keep_targets]
        print(f"  merging with {len(others):,} existing rows from "
              f"{len({r['target'] for r in others})} other targets", file=sys.stderr)
        meta = others + meta
        meta.sort(key=lambda r: (r["target"], r["sample_id"], r["contig_id"]))

    write_tsv(meta_path, meta,
              ["seq_name", "sample_id", "tick_species", "tick_species_source",
               "target", "contig_id", "contig_len", "mapped_reads", "tpm_microbial",
               "lca_taxid", "lca_name", "lca_genus", "lca_species", "top_bitscore",
               "target_present_in_sample", "target_evidence"])

    # These two are reported, never swallowed: a silently missing assembly would
    # look identical to a target that simply carries no marker genes.
    if missing_fasta:
        print(f"\nWARNING: {len(missing_fasta)} samples had no readable contig "
              f"FASTA and were skipped: {', '.join(sorted(missing_fasta)[:8])}"
              f"{' ...' if len(missing_fasta) > 8 else ''}", file=sys.stderr)
    if missing_contig:
        print(f"WARNING: {missing_contig} contig ids in the index were not found "
              f"in their sample's FASTA.", file=sys.stderr)
    # Keyed on the VALUE, not the source: a sample can be present in the metadata
    # template and still have no species, because its id carries no abbreviation
    # the guesser recognises. Both cases put `unknown_species` in the tip label,
    # so both have to be reported.
    n_samples = len({m["sample_id"] for m in meta})
    unknown = sorted({m["sample_id"] for m in meta
                      if m["tick_species"] in ("", "unknown_species")})
    if unknown:
        print(f"\nWARNING: {len(unknown)} of {n_samples} samples have no tick "
              f"species and are labelled `unknown_species` in every sequence "
              f"name.\n  Fill in metadata/sample_metadata.tsv and re-run steps "
              f"12-15 to relabel. Affected: "
              f"{', '.join(unknown[:6])}{' ...' if len(unknown) > 6 else ''}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
