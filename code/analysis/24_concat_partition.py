#!/usr/bin/env python3
"""Step 24 -- concatenate the 13 gene alignments into one partitioned supermatrix.

Only the 13 protein-coding genes go into the tree (the phylogeny requested is
built from PCGs; rrnL/rrnS stay as separate alignments in results/mito/alignments/
for anyone who wants them, but are not part of this matrix). Genes are
concatenated in the fixed order below so partition boundaries are reproducible.

A sample missing a gene gets an all-gap block of that gene's aligned length,
not a dropped row -- IQ-TREE handles missing data as gaps correctly under ML,
and dropping the sample instead would silently shrink every OTHER partition's
taxon set to whatever that sample happened to have. A sample recovering too
few genes overall (`tree.min_genes_for_tree`) is excluded entirely instead --
see the config comment for why.

Published reference mitogenomes (step 20b) are added as extra rows, ids
prefixed `REF_`, when `mito.ref_genomes.include_in_tree` is on; the same
min_genes_for_tree rule applies to them.

Outputs
  results/mito/concat/mito_concat.fasta        one row per included sample
  results/mito/concat/mito_partitions.txt       IQ-TREE -p partition file
  results/mito/concat/mito_partition_summary.tsv
  results/mito/concat/excluded_samples.tsv      samples left out, and why
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import cfg_mito, mito_root, read_fasta, write_fasta  # noqa: E402

PCG_GENES = ["COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L",
            "ND5", "ND6", "ATP6", "ATP8"]


def main() -> None:
    mc = cfg_mito()
    tc = mc["tree"]
    root = mito_root()
    aln_dir, concat_dir = root / "alignments", root / "concat"
    concat_dir.mkdir(parents=True, exist_ok=True)

    with open(root / "mito_selection.tsv") as fh:
        selection = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")
                    if r["status"] == "ok"}
    if mc.get("ref_genomes", {}).get("include_in_tree", False):
        with open(root / "ref_genomes" / "ref_genome_summary.tsv") as fh:
            refs = {r["ref_id"]: {"tick_species": r["species"]}
                    for r in csv.DictReader(fh, delimiter="\t") if r["status"] == "ok"}
        print(f"{len(refs)} reference mitogenomes added as rows", file=sys.stderr)
        selection.update(refs)

    per_gene_seqs: dict[str, dict[str, str]] = {}
    gene_len: dict[str, int] = {}
    skipped_genes = []
    for gene in PCG_GENES:
        f = aln_dir / f"{gene}.aln.fasta"
        if not f.exists():
            skipped_genes.append(gene)
            continue
        seqs = read_fasta(f)
        by_sample: dict[str, str] = {}
        for name, seq in seqs.items():
            sid = name.split("|", 1)[0]
            if sid in by_sample:
                print(f"WARNING: {gene} has more than one sequence for {sid}; "
                      f"keeping the first", file=sys.stderr)
                continue
            by_sample[sid] = seq
        per_gene_seqs[gene] = by_sample
        gene_len[gene] = len(next(iter(seqs.values())))

    used_genes = [g for g in PCG_GENES if g in per_gene_seqs]
    if not used_genes:
        sys.exit("ERROR: no gene alignments found -- run steps 21-23 first.")
    if skipped_genes:
        print(f"genes excluded from the matrix (no alignment; too few sequences): "
              f"{', '.join(skipped_genes)}", file=sys.stderr)

    concat: dict[str, str] = {}
    excluded: list[dict] = []
    for sid, sel in sorted(selection.items()):
        n_present = sum(1 for g in used_genes if sid in per_gene_seqs[g])
        if n_present < tc["min_genes_for_tree"]:
            excluded.append({"sample_id": sid, "tick_species": sel["tick_species"],
                            "reason": f"only {n_present}/{len(used_genes)} genes "
                                     f"recovered (< min_genes_for_tree="
                                     f"{tc['min_genes_for_tree']})"})
            continue
        parts = [per_gene_seqs[g].get(sid, "-" * gene_len[g]) for g in used_genes]
        concat[f"{sid}|{sel['tick_species'].replace(' ', '_')}"] = "".join(parts)

    write_fasta(concat_dir / "mito_concat.fasta", concat)
    write_tsv(concat_dir / "excluded_samples.tsv", excluded, ["sample_id", "tick_species", "reason"])

    part_lines = []
    summary = []
    offset = 1
    for g in used_genes:
        length = gene_len[g]
        start, end = offset, offset + length - 1
        part_lines.append(f"DNA, {g} = {start}-{end}")
        ids = [sid.split("|", 1)[0] for sid in concat]
        n_with = sum(1 for sid in ids if sid in per_gene_seqs[g])
        n_ref = sum(1 for sid in ids if sid in per_gene_seqs[g] and sid.startswith("REF_"))
        summary.append({"gene": g, "start": start, "end": end, "aln_len": length,
                        "n_samples_with_gene": n_with - n_ref, "n_refs_with_gene": n_ref,
                        "pct_missing": round(100 * (1 - n_with / max(1, len(concat))), 1)})
        offset = end + 1
    (concat_dir / "mito_partitions.txt").write_text("\n".join(part_lines) + "\n")
    write_tsv(concat_dir / "mito_partition_summary.tsv", summary,
              ["gene", "start", "end", "aln_len", "n_samples_with_gene", "n_refs_with_gene",
               "pct_missing"])

    total_len = offset - 1
    print(f"\n{len(concat)} rows (samples + references) in the matrix ({len(excluded)} excluded, "
          f"see excluded_samples.tsv), {len(used_genes)} gene partitions, "
          f"{total_len:,} total columns", file=sys.stderr)
    for s in summary:
        print(f"  {s['gene']:<8} {s['start']:>6}-{s['end']:<6} ({s['aln_len']:>5} cols)  "
              f"{s['n_samples_with_gene'] + s['n_refs_with_gene']:>4}/{len(concat)} rows "
              f"({100 - s['pct_missing']:.0f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
