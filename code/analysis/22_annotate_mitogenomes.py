#!/usr/bin/env python3
"""Step 22 -- turn the step-21 gene loci into a final annotation per sample.

Two annotation sources, in preference order:

  1. Real MITOS2 (https://gitlab.com/Bernt/MITOS), if `mito.annotate.mitos_refdata_dir`
     in config/analysis_params.yaml points at an existing directory AND the
     `runmitos` command is on PATH. Run per sample on its selected mitogenome
     contig, output parsed from `result.bed`.
  2. Otherwise (the path actually used in this project -- see below), the
     homology loci from step 21 directly: each kept, selected-contig locus is
     extracted, oriented, and -- for the 13 protein-coding genes -- placed into
     its correct reading frame with the INVERTEBRATE MITOCHONDRIAL genetic code
     (NCBI transl_table 5; standard-code translation would misread every
     in-frame TGA, a normal Trp codon here, as a premature stop).

WHY THE FALLBACK IS WHAT ACTUALLY RUNS HERE: MITOS2's conda package installs
cleanly, but its reference data is distributed only from Zenodo
(zenodo.org/record/4284483), which this compute node could not reach -- every
attempt (see README_mito.md) timed out at the TLS layer. Point
`mitos_refdata_dir` at a copy fetched elsewhere and this step switches to real
MITOS2 automatically, no code change needed.

WHAT THE FALLBACK DOES NOT GIVE YOU: tRNA and gene-order annotation, exact
start/stop-codon boundaries (gene edges are the homology-supported span from
step 21, not a codon-precise call), and no circularity check. It is sufficient
for the per-gene alignments and partitioned tree in steps 23-25, which is the
explicit target of this workflow; it is not a substitute for a full MITOS2/
GenBank-quality annotation.

Outputs
  results/mito/genes/<GENE>.fna / .faa      one FASTA per gene, all samples
  results/mito/gff/<sample>.gff3            per-sample annotation, either source
  results/mito/annotation_summary.tsv       one row per sample
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (INVERT_MITO_TABLE, cfg_mito, frame_mito_pcg,  # noqa: E402
                         marker_name, mito_root, read_fasta, require, revcomp, run,
                         write_fasta)

SUMMARY_COLS = ["sample_id", "tick_species", "annotation_source", "contig_id",
                "contig_len", "n_genes_annotated", "genes_annotated",
                "n_flagged_internal_stops"]


def load_selection(root: Path) -> list[dict]:
    with open(root / "mito_selection.tsv") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["status"] == "ok"]


def load_loci(root: Path) -> dict[str, list[dict]]:
    """sample_id -> its kept, selected-contig loci (the homology annotation)."""
    out: dict[str, list[dict]] = defaultdict(list)
    with open(root / "mito_loci.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] == "kept" and r["is_selected_contig"] == "1":
                out[r["sample_id"]].append(r)
    return out


def mitos_available(ac: dict) -> str | None:
    refdata = Path(ac["mitos_refdata_dir"]) if ac.get("mitos_refdata_dir") else None
    if not refdata or not refdata.is_dir():
        return None
    from shutil import which
    return which("runmitos")


def run_mitos(exe: str, genome_fa: Path, outdir: Path, code: int, refdata: Path) -> list[dict]:
    """Run real MITOS2 and parse its result.bed into the same schema as the
    homology loci. UNTESTED beyond argument construction: this project's compute
    node cannot reach the Zenodo reference data (see module docstring), so this
    path has never been exercised end-to-end. If your MITOS2 version's
    result.bed columns differ, fix the parsing below -- the invocation itself
    is per MITOS2's own --help.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    # -r/--refseqver takes just the reference-data directory's OWN name, taken
    # from `refdata` itself so it can never drift out of sync with -R/--refdir
    # (its parent) -- these two must name the same directory between them.
    run([exe, "-i", str(genome_fa), "-o", str(outdir), "-c", str(code),
         "--linear", "-r", refdata.name, "--refdir", str(refdata.parent)])
    bed = outdir / "result.bed"
    if not bed.exists():
        raise FileNotFoundError(f"MITOS2 did not write {bed}")
    loci = []
    with open(bed) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            contig, start, end, name, score, strand = p[:6]
            loci.append({"gene": name, "contig_id": contig, "locus_start": int(start) + 1,
                        "locus_end": int(end), "strand": strand})
    return loci


def write_gff3(path: Path, sample: str, contig: str, loci: list[dict], source: str) -> None:
    lines = ["##gff-version 3"]
    for l in sorted(loci, key=lambda x: int(x["locus_start"])):
        typ = "CDS" if l.get("gene_type", "pcg") == "pcg" or l["gene"] not in ("rrnL", "rrnS") else "rRNA"
        lines.append("\t".join([contig, source, typ, str(l["locus_start"]), str(l["locus_end"]),
                                ".", l["strand"], "0" if typ == "CDS" else ".",
                                f"ID={l['gene']};Name={l['gene']}"]))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append", help="limit to this sample (repeatable)")
    args = ap.parse_args()

    mc = cfg_mito()
    ac = mc["annotate"]
    code = int(mc["genetic_code"])
    table = INVERT_MITO_TABLE if code == 5 else None
    if table is None:
        sys.exit(f"ERROR: mito.genetic_code={code} has no translation table wired up "
                 f"in markers_lib.py -- add one before using a non-invertebrate code.")

    root = mito_root()
    gene_dir, gff_dir = root / "genes", root / "gff"
    gene_dir.mkdir(parents=True, exist_ok=True)
    gff_dir.mkdir(parents=True, exist_ok=True)

    selection = load_selection(root)
    if args.sample:
        selection = [s for s in selection if s["sample_id"] in set(args.sample)]
    homology_loci = load_loci(root)

    mitos_exe = mitos_available(ac)
    refdata = Path(ac["mitos_refdata_dir"]) if mitos_exe else None
    print(f"MITOS2: {'using ' + mitos_exe if mitos_exe else 'refdata not configured/found -- using homology fallback'}",
          file=sys.stderr)

    by_gene_nt: dict[str, dict[str, str]] = defaultdict(dict)
    by_gene_aa: dict[str, dict[str, str]] = defaultdict(dict)
    # A --sample run must not drop the other samples from the per-gene FASTAs,
    # which aggregate across the whole cohort: seed each with what is already on
    # disk, keeping only entries for samples OUTSIDE this run's scope.
    scope = set(args.sample) if args.sample else None
    if not scope:
        # A full run rebuilds every gene's FASTA from scratch, authoritatively.
        # Without this, a gene/sample whose classification changes between runs
        # (e.g. flagged -> unflagged, as happened when the stop-codon-trim logic
        # below was added) leaves its STALE file from the previous run on disk
        # untouched, because the write loop only overwrites keys that have
        # entries in THIS run -- a gene with zero flagged calls this time never
        # gets its old, now-wrong *.flagged.fna file removed otherwise.
        for f in list(gene_dir.glob("*.fna")) + list(gene_dir.glob("*.faa")):
            f.unlink()
    if scope:
        for f in gene_dir.glob("*.fna"):
            seqs = {n: s for n, s in read_fasta(f).items() if n.split("|", 1)[0] not in scope}
            if seqs:
                by_gene_nt[f.stem] = seqs
        for f in gene_dir.glob("*.faa"):
            seqs = {n: s for n, s in read_fasta(f).items() if n.split("|", 1)[0] not in scope}
            if seqs:
                by_gene_aa[f.stem] = seqs
    summary: list[dict] = []
    n_flagged_total = 0

    for i, s in enumerate(selection, 1):
        sid, sp, contig = s["sample_id"], s["tick_species"], s["selected_contig"]
        genome_fa = ROOT / s["fasta"]
        genome = read_fasta(genome_fa)
        header = next(iter(genome))
        full_seq = genome[header]
        tag = f"[{i:>3}/{len(selection)}] {sid}"

        source = "homology_blast"
        loci = homology_loci.get(sid, [])
        if mitos_exe:
            try:
                mitos_loci = run_mitos(mitos_exe, genome_fa, root / "mitos" / sid, code,
                                       refdata)
                if mitos_loci:
                    loci, source = mitos_loci, "MITOS2"
            except Exception as e:  # noqa: BLE001 -- never let an unverified path abort the run
                print(f"{tag}: MITOS2 failed ({e}); falling back to homology calls",
                      file=sys.stderr)

        n_flagged = 0
        genes_here = []
        for l in loci:
            gene = l["gene"]
            lo, hi = int(l["locus_start"]), int(l["locus_end"])
            if hi - lo + 1 < ac["min_locus_len_nt"]:
                continue
            seq = full_seq[lo - 1:hi]
            if l["strand"] == "-":
                seq = revcomp(seq)

            is_pcg = gene not in ("rrnL", "rrnS")
            flagged = False
            if is_pcg:
                # frame, trim a trailing natural stop, flag real internal stops
                # -- see frame_mito_pcg; shared with step 20b's references
                seq, aa, flagged = frame_mito_pcg(seq, table)
                n_flagged += int(flagged)

            name = marker_name(sid, sp, "mitochondrion", gene, contig, len(seq))
            fkey = gene if not flagged else f"{gene}.flagged"
            by_gene_nt[fkey][name] = seq
            if is_pcg:
                by_gene_aa[fkey][name] = aa
            genes_here.append(gene)

        write_gff3(gff_dir / f"{sid}.gff3", sid, contig, loci, source)
        n_flagged_total += n_flagged
        summary.append({"sample_id": sid, "tick_species": sp, "annotation_source": source,
                        "contig_id": contig, "contig_len": s["contig_len"],
                        "n_genes_annotated": len(genes_here),
                        "genes_annotated": ",".join(sorted(genes_here)),
                        "n_flagged_internal_stops": n_flagged})
        if i % 25 == 0 or i == len(selection):
            print(f"{tag}: {len(genes_here)} genes annotated ({source})", file=sys.stderr)

    for gene, seqs in sorted(by_gene_nt.items()):
        write_fasta(gene_dir / f"{gene}.fna", seqs)
    for gene, seqs in sorted(by_gene_aa.items()):
        write_fasta(gene_dir / f"{gene}.faa", seqs)

    summary_path = root / "annotation_summary.tsv"
    if scope and summary_path.exists():
        with open(summary_path) as fh:
            others = [r for r in csv.DictReader(fh, delimiter="\t")
                     if r["sample_id"] not in scope]
        summary = others + summary
    summary.sort(key=lambda r: r["sample_id"])
    write_tsv(summary_path, summary, SUMMARY_COLS)

    print(f"\n{len(summary)} samples annotated, {n_flagged_total} gene calls flagged "
          f"for internal stop codons (held in *.flagged.fna, out of the alignments)",
          file=sys.stderr)
    for gene in ["COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L",
                "ND5", "ND6", "ATP6", "ATP8", "rrnL", "rrnS"]:
        n = len(by_gene_nt.get(gene, {}))
        print(f"  {gene:<8} {n:>4} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
