#!/usr/bin/env python3
"""Step 14 -- cut the marker sequences out of the contigs and label them.

Reads the loci assembled by step 13 and writes, per target x marker:

  results/markers/sequences/<target>/<target>__<marker>.fna    nucleotide
  results/markers/sequences/<target>/<target>__<marker>.faa    protein loci only

Every sequence is named

  <sample_id>|<tick_species>|<target>|<marker>|<contig_id>|<len>bp[|copyN]

so a tip label in the finished tree states the sample it came from, the tick it
was found in, and which gene it is, with no lookup required. `copyN` appears
only when one sample yields the same locus on more than one contig.

READING FRAME. A protein-coding locus is trimmed to the forward frame with the
fewest internal stop codons, so the exported nucleotide sequence is in frame and
a multiple of three -- which is what makes the codon-aware alignment in step 15
possible. The frame is recovered from the sequence rather than taken from
tblastn's `sframe`, because a locus assembled from merged HSPs need not start on
a codon boundary, and a minus-strand frame is reported relative to the whole
contig rather than the extracted span.

Internal stop codons that survive frame selection are counted, not removed:
they usually mean an assembly frameshift or a pseudogene, and that is something
you must see before putting the sequence in a tree, not something to paper over.
Sequences with any internal stop are written to a parallel `*.flagged.fna`
instead of the main file.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (SEP, best_frame, cfg_markers, internal_stops, marker_name,  # noqa: E402
                         out_root, read_fasta, revcomp, slug, translate, write_fasta)

SEQ_COLS = ["seq_id", "sample_id", "tick_species", "tick_species_source", "target",
            "marker", "marker_type", "contig_id", "contig_len", "locus_start",
            "locus_end", "strand", "seq_len_nt", "aa_len", "n_internal_stops",
            "frame", "best_pident", "max_bitscore", "ref_cov_pct", "best_ref",
            "mapped_reads", "tpm_microbial", "lca_name", "lca_species",
            "target_present_in_sample", "target_evidence", "flagged", "fasta"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    args = ap.parse_args()

    mc = cfg_markers()
    ec = mc["extract"]
    root = out_root()
    loci_path = root / "marker_loci.tsv"
    if not loci_path.exists():
        sys.exit("ERROR: results/markers/marker_loci.tsv not found -- run step 13 first.")

    meta: dict[str, dict] = {}
    with open(root / "contigs" / "contig_metadata.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["seq_name"]] = r

    loci: list[dict] = []
    with open(loci_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] != "kept":
                continue
            if args.target and r["target"] not in set(args.target):
                continue
            loci.append(r)
    print(f"{len(loci):,} kept loci", file=sys.stderr)

    contig_cache: dict[str, dict[str, str]] = {}
    rows: list[dict] = []
    # (sample, target, marker) -> running copy number, assigned longest first so
    # `copy1` is always the best sequence for that sample.
    loci.sort(key=lambda r: (r["target"], r["marker"], r["sample_id"] if "sample_id" in r
                             else r["seq_name"], -int(r["locus_len"])))
    copy_no: dict[tuple[str, str, str], int] = defaultdict(int)
    by_file: dict[tuple[str, str, bool], dict[str, str]] = defaultdict(dict)
    aa_by_file: dict[tuple[str, str, bool], dict[str, str]] = defaultdict(dict)
    n_short = n_flagged = 0

    for l in loci:
        target = l["target"]
        if target not in contig_cache:
            fa = root / "contigs" / f"{slug(target)}.contigs.fa"
            contig_cache[target] = read_fasta(fa) if fa.exists() else {}
        seq_full = contig_cache[target].get(l["seq_name"])
        if seq_full is None:
            continue

        is_prot = l["marker_type"] == "protein"
        flank = 0 if is_prot else int(ec["flank_rrna"])
        start = max(1, int(l["locus_start"]) - flank)
        end = min(len(seq_full), int(l["locus_end"]) + flank)
        sub = seq_full[start - 1:end]
        if l["strand"] == "-":
            sub = revcomp(sub)

        frame, aa = (0, "")
        if is_prot:
            frame, aa = best_frame(sub)
            sub = sub[frame:]
            sub = sub[: len(sub) - len(sub) % 3]   # keep it a whole number of codons
            aa = translate(sub, 0)
        if len(sub) < int(ec["min_locus_len_nt"]):
            n_short += 1
            continue

        m = meta.get(l["seq_name"], {})
        sample = m.get("sample_id", l["seq_name"].split(SEP)[0])
        species = m.get("tick_species", l["seq_name"].split(SEP)[1] if SEP in l["seq_name"] else "unknown")
        contig_id = m.get("contig_id", l["seq_name"].rsplit(SEP, 1)[-1])

        key = (sample, target, l["marker"])
        copy_no[key] += 1
        copy = copy_no[key]
        name = marker_name(sample, species, target, l["marker"], contig_id, len(sub), copy)

        stops = internal_stops(aa) if is_prot else 0
        flagged = stops > 0
        n_flagged += int(flagged)

        fkey = (target, l["marker"], flagged)
        by_file[fkey][name] = sub
        if is_prot:
            aa_by_file[fkey][name] = aa

        suffix = ".flagged" if flagged else ""
        rows.append({
            "seq_id": name, "sample_id": sample, "tick_species": species,
            "tick_species_source": m.get("tick_species_source", ""),
            "target": target, "marker": l["marker"], "marker_type": l["marker_type"],
            "contig_id": contig_id, "contig_len": l["contig_len"],
            "locus_start": start, "locus_end": end, "strand": l["strand"],
            "seq_len_nt": len(sub), "aa_len": len(aa) if is_prot else "",
            "n_internal_stops": stops if is_prot else "", "frame": frame if is_prot else "",
            "best_pident": l["best_pident"], "max_bitscore": l["max_bitscore"],
            "ref_cov_pct": l["ref_cov_pct"], "best_ref": l["best_ref"],
            "mapped_reads": m.get("mapped_reads", ""),
            "tpm_microbial": m.get("tpm_microbial", ""),
            "lca_name": m.get("lca_name", ""), "lca_species": m.get("lca_species", ""),
            "target_present_in_sample": m.get("target_present_in_sample", ""),
            "target_evidence": m.get("target_evidence", ""),
            "flagged": "internal_stop_codons" if flagged else "",
            "fasta": f"sequences/{slug(target)}/{slug(target)}__{slug(l['marker'])}{suffix}.fna",
        })

    for (target, marker, flagged), seqs in sorted(by_file.items()):
        d = root / "sequences" / slug(target)
        stem = f"{slug(target)}__{slug(marker)}{'.flagged' if flagged else ''}"
        write_fasta(d / f"{stem}.fna", seqs)
        if aa_by_file.get((target, marker, flagged)):
            write_fasta(d / f"{stem}.faa", aa_by_file[(target, marker, flagged)])

    # As in steps 12 and 13: keep the other targets' sequences in the table, so
    # a --target re-run does not blank the step-16 report for everything else.
    seq_path = root / "marker_sequences.tsv"
    if args.target and seq_path.exists():
        scope = set(args.target)
        with open(seq_path) as fh:
            others = [r for r in csv.DictReader(fh, delimiter="\t")
                      if r["target"] not in scope]
        print(f"  merging with {len(others):,} existing sequences from other targets",
              file=sys.stderr)
        rows = sorted(others + rows, key=lambda r: (r["target"], r["marker"],
                                                    r["sample_id"]))
    write_tsv(seq_path, rows, SEQ_COLS)

    per = defaultdict(lambda: [0, set()])
    for r in rows:
        if r["flagged"]:
            continue
        k = (r["target"], r["marker"])
        per[k][0] += 1
        per[k][1].add(r["sample_id"])
    print(f"\n{'target':<20} {'marker':<14} {'seqs':>6} {'samples':>8}", file=sys.stderr)
    for (t, mk), (n, s) in sorted(per.items()):
        print(f"{t:<20} {mk:<14} {n:>6,} {len(s):>8}", file=sys.stderr)
    print(f"\n{len(rows):,} sequences written "
          f"({n_flagged:,} flagged for internal stop codons, "
          f"{n_short:,} dropped below {ec['min_locus_len_nt']} bp)", file=sys.stderr)


if __name__ == "__main__":
    main()
