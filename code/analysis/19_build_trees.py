#!/usr/bin/env python3
"""Step 19 -- align samples with their references and build one tree per locus.

For every target x marker that recovered a sample sequence, this combines the
step-14 sample sequences with the step-18 references, aligns them, and runs
FastTree.

  rRNA loci        MAFFT --adjustdirection on the nucleotides.
  protein-coding   MAFFT on the translations, back-translated onto the in-frame
                   nucleotides. Both samples and references were trimmed to a
                   whole number of codons in the forward frame (steps 14 and
                   18), so the result is a true codon alignment -- gaps on codon
                   boundaries, column 3n+1..3n+3 always one codon. Aligning
                   these genes as raw nucleotide would put gaps mid-codon and
                   hand the substitution model an alignment it does not describe.

FastTree, not IQ-TREE, is deliberate: approximately-ML with SH-like local
support, seconds per locus over 47 loci. It is the screening tree. Re-run the
loci that matter under a model-tested partitioned analysis -- `iqtree` is in the
same environment, and step 25 shows the pattern.

Written per locus, all three in one directory so a locus travels as a unit:

  <target>/<target>__<marker>.aln.fna    the alignment the tree was built from
  <target>/<target>__<marker>.nwk        FastTree output, UNROOTED
  <target>/<target>__<marker>.aa.aln.faa protein alignment, CDS loci only

Step 19b roots each tree and draws it. A locus with fewer than
`trees.tree.min_tips` sequences still gets its alignment written -- one tree
fewer is not a reason to lose the alignment.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (align, backtranslate, read_fasta, require, run, slug,  # noqa: E402
                         write_fasta)
from trees_lib import (cfg_trees, is_outgroup, is_reference, loci_with_sequences,  # noqa: E402
                       sample_fasta, trees_root)

SUM_COLS = ["target", "marker", "marker_type", "n_tips", "n_samples", "n_refs",
            "n_outgroups", "n_organisms", "aln_cols", "median_ungapped_nt",
            "pct_gaps", "codon_aware", "trimmed", "tree_built", "note",
            "alignment", "tree"]


def strip_flip_prefix(seqs: dict[str, str]) -> tuple[dict[str, str], int]:
    """Undo MAFFT's `_R_` name prefix on sequences it reverse-complemented.

    --adjustdirection renames what it flips. The sequence content in the
    alignment is already the flipped one, so only the name has to come back --
    and it must, because a tip called `_R_REF|...` no longer parses as a
    reference and would be drawn as a sample.
    """
    out, n = {}, 0
    for name, seq in seqs.items():
        if name.startswith("_R_"):
            name, n = name[3:], n + 1
        out[name] = seq
    return out, n


def trim(spec: str, infile: Path, outfile: Path) -> bool:
    """Apply trimal if `trees.align.trim` asks for it. Returns whether it ran."""
    spec = (spec or "none").strip().lower()
    if spec in ("", "none", "off", "false"):
        return False
    exe = require("trimal")
    cmd = [exe, "-in", str(infile), "-out", str(outfile)]
    if spec == "gappyout":
        cmd.append("-gappyout")
    elif spec.startswith("gt:"):
        cmd += ["-gt", spec.split(":", 1)[1]]
    else:
        sys.exit(f"ERROR: trees.align.trim = {spec!r} -- expected none | gappyout | gt:<0-1>")
    run(cmd)
    return True


def stats(aln: dict[str, str]) -> tuple[int, int, float]:
    cols = len(next(iter(aln.values()))) if aln else 0
    lens = sorted(len(s) - s.count("-") for s in aln.values())
    med = lens[len(lens) // 2] if lens else 0
    total = cols * len(aln)
    gaps = sum(s.count("-") for s in aln.values())
    return cols, med, round(100.0 * gaps / total, 1) if total else 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    ap.add_argument("--marker", action="append", help="limit to this marker (repeatable)")
    ap.add_argument("--include-flagged", action="store_true",
                    help="also include the *.flagged sample sequences (internal stops)")
    ap.add_argument("--force", action="store_true", help="realign and rebuild cached trees")
    args = ap.parse_args()

    tc = cfg_trees()
    ac, trc = tc["align"], tc["tree"]
    root = trees_root()
    sets_dir = root / "refs" / "sets"
    work = root / "work"
    work.mkdir(parents=True, exist_ok=True)

    aligner = require(ac["tool"])
    fasttree = require(trc["tool"])

    loci = loci_with_sequences(args.target, args.marker)
    if not loci:
        sys.exit("ERROR: no target x marker has sample sequences -- run run_markers.sh first.")

    rows: list[dict] = []
    n_trees = 0
    for i, row in enumerate(loci, 1):
        target, marker = row["target"], row["marker"]
        key = f"{slug(target)}__{slug(marker)}"
        is_prot = row["marker_type"] == "protein"
        tag = f"[{i:>3}/{len(loci)}] {key}"
        outdir = root / slug(target)
        outdir.mkdir(parents=True, exist_ok=True)
        aln_path = outdir / f"{key}.aln.fna"
        nwk_path = outdir / f"{key}.nwk"

        nt = dict(read_fasta(sample_fasta(target, marker, "fna")))
        aa = dict(read_fasta(sample_fasta(target, marker, "faa"))) if is_prot else {}
        if args.include_flagged:
            nt.update(read_fasta(sample_fasta(target, marker, "fna", flagged=True))
                      if sample_fasta(target, marker, "fna", flagged=True).exists() else {})
            if is_prot and sample_fasta(target, marker, "faa", flagged=True).exists():
                aa.update(read_fasta(sample_fasta(target, marker, "faa", flagged=True)))
        n_samples = len(nt)

        ref_fna = sets_dir / f"{key}.refs.fna"
        if ref_fna.exists():
            nt.update(read_fasta(ref_fna))
            if is_prot and (sets_dir / f"{key}.refs.faa").exists():
                aa.update(read_fasta(sets_dir / f"{key}.refs.faa"))
        n_refs = sum(1 for n in nt if is_reference(n))
        n_out = sum(1 for n in nt if is_outgroup(n))
        n_org = len({n.split("|")[2] for n in nt if is_reference(n)})

        note = ""
        if n_refs == 0:
            note = "no references: sample sequences only"
        if len(nt) < int(trc["min_tips"]):
            shutil.copyfile(sample_fasta(target, marker, "fna"),
                            outdir / f"{key}.unaligned.fna")
            print(f"{tag}: {len(nt)} seqs -- under min_tips, no tree", file=sys.stderr)
            rows.append({"target": target, "marker": marker,
                         "marker_type": row["marker_type"], "n_tips": len(nt),
                         "n_samples": n_samples, "n_refs": n_refs, "n_outgroups": n_out,
                         "n_organisms": n_org, "aln_cols": 0, "median_ungapped_nt": 0,
                         "pct_gaps": 0.0, "codon_aware": "", "trimmed": "",
                         "tree_built": "no",
                         "note": (note + "; " if note else "") + "fewer tips than trees.tree.min_tips",
                         "alignment": "", "tree": ""})
            continue

        # ---- align -----------------------------------------------------------
        # Rebuild when an input is newer than the alignment. Without this, a
        # step-18 re-run that adds references to a locus would leave step 19
        # happily reusing the alignment built before they existed, and the new
        # references would silently never reach the tree.
        inputs = [sample_fasta(target, marker, "fna"), ref_fna]
        newest = max((p.stat().st_mtime for p in inputs if p.exists()), default=0.0)
        stale = aln_path.exists() and aln_path.stat().st_mtime < newest
        if args.force or not aln_path.exists() or stale:
            if is_prot and ac.get("codon_aware", True):
                missing = [n for n in nt if n not in aa]
                if missing:
                    sys.exit(f"ERROR: {key}: {len(missing)} sequences have no translation "
                             f"(first: {missing[0]}). Re-run step 14 / step 18.")
                aa_in = work / f"{key}.faa"
                write_fasta(aa_in, aa)
                aa_out = outdir / f"{key}.aa.aln.faa"
                align(ac["tool"], ac["mafft_opts"].split(), aa_in, aa_out, aligner)
                aa_aln, _ = strip_flip_prefix(read_fasta(aa_out))
                write_fasta(aa_out, aa_aln)
                write_fasta(aln_path, backtranslate(aa_aln, nt))
                kind = "codon"
            else:
                nt_in = work / f"{key}.fna"
                write_fasta(nt_in, nt)
                align(ac["tool"], ac["mafft_rrna_opts"].split(), nt_in, aln_path, aligner)
                flipped_aln, n_flip = strip_flip_prefix(read_fasta(aln_path))
                write_fasta(aln_path, flipped_aln)
                if n_flip:
                    note = (note + "; " if note else "") + f"{n_flip} sequences reoriented by MAFFT"
                kind = "nt"
            if trim(ac.get("trim", "none"), aln_path, work / f"{key}.trimmed.fna"):
                shutil.copyfile(aln_path, outdir / f"{key}.untrimmed.aln.fna")
                shutil.copyfile(work / f"{key}.trimmed.fna", aln_path)
        else:
            kind = "codon" if is_prot and ac.get("codon_aware", True) else "nt"

        aln = read_fasta(aln_path)
        cols, med, pct_gaps = stats(aln)

        # ---- tree ------------------------------------------------------------
        if args.force or not nwk_path.exists() or stale or \
                nwk_path.stat().st_mtime < aln_path.stat().st_mtime:
            r = run([fasttree] + trc["fasttree_nt_opts"].split() + [str(aln_path)])
            nwk_path.write_text(r.stdout)
        n_trees += 1
        print(f"{tag}: {len(nt):>4} tips ({n_samples} samples + {n_refs} refs, "
              f"{n_org} organisms)  {cols:>5} cols  {pct_gaps:>4.1f}% gaps  ({kind})",
              file=sys.stderr)

        rows.append({"target": target, "marker": marker, "marker_type": row["marker_type"],
                     "n_tips": len(aln), "n_samples": n_samples, "n_refs": n_refs,
                     "n_outgroups": n_out, "n_organisms": n_org, "aln_cols": cols,
                     "median_ungapped_nt": med, "pct_gaps": pct_gaps,
                     "codon_aware": "yes" if kind == "codon" else "no",
                     "trimmed": ac.get("trim", "none"), "tree_built": "yes", "note": note,
                     "alignment": str(aln_path.relative_to(ROOT)),
                     "tree": str(nwk_path.relative_to(ROOT))})

    sum_path = root / "tree_summary.tsv"
    if (args.target or args.marker) and sum_path.exists():
        import csv
        done = {(r["target"], r["marker"]) for r in rows}
        with open(sum_path) as fh:
            others = [r for r in csv.DictReader(fh, delimiter="\t")
                      if (r["target"], r["marker"]) not in done]
        rows = others + rows
    rows.sort(key=lambda r: (-int(r["n_tips"]), r["target"], r["marker"]))
    write_tsv(sum_path, rows, SUM_COLS)
    print(f"\n{n_trees} trees built -> {root.relative_to(ROOT)}/<target>/", file=sys.stderr)
    print("next: step 19b roots each tree and draws the figures", file=sys.stderr)


if __name__ == "__main__":
    main()
