#!/usr/bin/env python3
"""Step 23 -- align each mitochondrial gene separately across all samples.

Protein-coding genes are aligned as amino acids (MAFFT) then back-translated
onto their in-frame nucleotides, giving a codon alignment where every column
triplet is one codon and no gap falls mid-codon -- required for step 24 to
assign a sensible per-codon-position partition later if you want one, and for
IQ-TREE's codon/DNA models to see real codon columns. The two rRNAs are
aligned directly as nucleotide.

If `mito.ref_genomes.include_in_tree` is on, each gene's published reference
sequences (step 20b, results/mito/ref_genomes/genes/) are aligned together
with the samples' -- the merged input is written beside the alignment as
<GENE>.input.fna/.faa so the exact aligned set is on disk.

Genes recovered from fewer than `align.min_seqs_to_align` samples are copied
through unaligned rather than dropped -- COX1 typically has one of the largest
sample counts and the rarer genes (ATP8, ND4L) the fewest; step 24 needs to
know which samples are missing which gene either way.

Outputs
  results/mito/alignments/<GENE>.aln.fasta       nucleotide, codon-aware for PCGs
  results/mito/alignments/<GENE>.aa.aln.fasta     protein, PCGs only
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from markers_lib import align, backtranslate, cfg_mito, mito_root, read_fasta, require, write_fasta  # noqa: E402

PCG_GENES = ["COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L",
            "ND5", "ND6", "ATP6", "ATP8"]
RRNA_GENES = ["rrnL", "rrnS"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", action="append", help="limit to this gene (repeatable)")
    ap.add_argument("--include-flagged", action="store_true",
                    help="also align the *.flagged sets (internal stop codons)")
    args = ap.parse_args()

    mc = cfg_mito()
    ac = mc["align"]
    tool = ac["tool"]
    exe = require(tool)
    root = mito_root()
    gene_dir, aln_dir = root / "genes", root / "alignments"
    aln_dir.mkdir(parents=True, exist_ok=True)

    use_refs = bool(mc.get("ref_genomes", {}).get("include_in_tree", False))
    ref_gene_dir = root / "ref_genomes" / "genes"
    if use_refs and not ref_gene_dir.is_dir():
        sys.exit("ERROR: mito.ref_genomes.include_in_tree is on but "
                 "results/mito/ref_genomes/genes/ is missing -- run step 20b first.")

    genes = PCG_GENES + RRNA_GENES
    if args.gene:
        genes = [g for g in genes if g in set(args.gene)]

    n_aln = n_skip = 0
    for gene in genes:
        is_pcg = gene in PCG_GENES
        candidates = [gene_dir / f"{gene}.fna"]
        if args.include_flagged:
            candidates.append(gene_dir / f"{gene}.flagged.fna")
        fna = candidates[0]
        if not fna.exists():
            print(f"  {gene:<8} no sequences recovered (see mito_gaps in the report)",
                  file=sys.stderr)
            continue

        nt = read_fasta(fna)
        if args.include_flagged and candidates[1].exists():
            nt.update(read_fasta(candidates[1]))
        n_samples = len(nt)
        aa_in = gene_dir / f"{gene}.faa"
        if use_refs and (ref_gene_dir / f"{gene}.fna").exists():
            nt.update(read_fasta(ref_gene_dir / f"{gene}.fna"))
            fna = aln_dir / f"{gene}.input.fna"
            write_fasta(fna, nt)
            if is_pcg:
                aa = read_fasta(gene_dir / f"{gene}.faa")
                aa.update(read_fasta(ref_gene_dir / f"{gene}.faa"))
                aa_in = aln_dir / f"{gene}.input.faa"
                write_fasta(aa_in, aa)

        if len(nt) < int(ac["min_seqs_to_align"]):
            shutil.copyfile(fna, aln_dir / f"{gene}.unaligned.fasta")
            print(f"  {gene:<8} {len(nt):>4} seq -- too few to align, copied through",
                  file=sys.stderr)
            n_skip += 1
            continue

        if is_pcg:
            aa_out = aln_dir / f"{gene}.aa.aln.fasta"
            align(tool, ac["mafft_opts"].split(), aa_in, aa_out, exe)
            codon = backtranslate(read_fasta(aa_out), nt)
            write_fasta(aln_dir / f"{gene}.aln.fasta", codon)
            kind = "codon"
        else:
            align(tool, (ac["mafft_opts"] + " --adjustdirection").split(), fna,
                  aln_dir / f"{gene}.aln.fasta", exe)
            kind = "nt"

        alen = len(next(iter(read_fasta(aln_dir / f"{gene}.aln.fasta").values())))
        print(f"  {gene:<8} {n_samples:>4} samples + {len(nt) - n_samples:>3} refs  "
              f"{alen:>6} cols  ({kind})", file=sys.stderr)
        n_aln += 1

    print(f"\n{n_aln} gene alignments written to results/mito/alignments/, "
          f"{n_skip} genes had too few sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
