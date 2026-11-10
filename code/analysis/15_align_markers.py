#!/usr/bin/env python3
"""Step 15 -- align each target x marker set, ready for reference sequences.

  rRNA loci        MAFFT --adjustdirection on the nucleotides. Contigs assemble
                   in arbitrary orientation and blastn hits on the minus strand
                   are already reverse-complemented by step 14, but
                   --adjustdirection is kept as a safety net; MAFFT marks any
                   sequence it flips with an `_R_` name prefix.

  protein-coding   MAFFT on the translations, then the protein alignment is
                   back-translated onto the in-frame nucleotides. The result is
                   a codon alignment: gaps fall on codon boundaries and column
                   3n+1..3n+3 is always one codon. Aligning these genes as raw
                   nucleotide would put gaps mid-codon and produce a tree built
                   on an alignment no substitution model describes.

Written per target x marker:
  <target>__<marker>.aln.fasta        nucleotide (codon-aware for CDS loci)
  <target>__<marker>.aa.aln.fasta     protein, CDS loci only

NOW ADD YOUR REFERENCES. These alignments contain only sequences recovered from
the samples. Add type strains and named species to the FASTA and realign, or
align your references in as a profile:

  mafft --add my_refs.fasta --keeplength --auto <target>__<marker>.aln.fasta

Sets with fewer than `align.min_seqs_to_align` sequences are copied through
unaligned -- one sequence is not an alignment, and it is better to see the file
than to wonder why the locus vanished.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from markers_lib import (align, backtranslate, cfg_markers, out_root, read_fasta,  # noqa: E402
                         require, write_fasta)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    ap.add_argument("--include-flagged", action="store_true",
                    help="also align the *.flagged sets (internal stop codons)")
    args = ap.parse_args()

    mc = cfg_markers()
    ac = mc["align"]
    tool = ac["tool"]
    exe = require(tool)
    root = out_root()
    seqdir, alndir = root / "sequences", root / "alignments"
    if not seqdir.exists():
        sys.exit("ERROR: results/markers/sequences/ not found -- run step 14 first.")

    files = sorted(seqdir.glob("*/*.fna"))
    if args.target:
        keep = {t.replace(" ", "_") for t in args.target}
        files = [f for f in files if f.parent.name in keep]
    if not args.include_flagged:
        files = [f for f in files if ".flagged" not in f.name]

    n_aln = n_skip = 0
    for fna in files:
        target_dir = alndir / fna.parent.name
        target_dir.mkdir(parents=True, exist_ok=True)
        stem = fna.name[: -len(".fna")]
        faa = fna.with_suffix(".faa")
        is_prot = faa.exists()
        nt = read_fasta(fna)

        if len(nt) < int(ac["min_seqs_to_align"]):
            shutil.copyfile(fna, target_dir / f"{stem}.unaligned.fasta")
            print(f"  {stem:<40} {len(nt):>4} seq -- too few to align, copied through",
                  file=sys.stderr)
            n_skip += 1
            continue

        if is_prot and ac.get("codon_aware", True):
            aa_out = target_dir / f"{stem}.aa.aln.fasta"
            align(tool, ac["mafft_opts"].split(), faa, aa_out, exe)
            codon = backtranslate(read_fasta(aa_out), nt)
            write_fasta(target_dir / f"{stem}.aln.fasta", codon)
            kind = "codon"
        else:
            opts = (ac["mafft_rrna_opts"] if not is_prot else ac["mafft_opts"]).split()
            align(tool, opts, fna, target_dir / f"{stem}.aln.fasta", exe)
            kind = "nt"

        alen = len(next(iter(read_fasta(target_dir / f"{stem}.aln.fasta").values())))
        print(f"  {stem:<40} {len(nt):>4} seq  {alen:>6} cols  ({kind})", file=sys.stderr)
        n_aln += 1

    print(f"\n{n_aln} alignments written to results/markers/alignments/, "
          f"{n_skip} sets had too few sequences", file=sys.stderr)
    print("Add your reference sequences next -- see the header of this script "
          "or README_markers.md", file=sys.stderr)


if __name__ == "__main__":
    main()
