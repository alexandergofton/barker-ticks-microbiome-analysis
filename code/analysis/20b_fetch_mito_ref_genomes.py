#!/usr/bin/env python3
"""Step 20b -- download the published reference mitogenomes and annotate them.

The panel is metadata/mito_ref_genomes.tsv (curated from the raw list in
results/mito/mito_ref_seqs.txt; rows with use=no are skipped and say why).
These are complete, named-species tick mitogenomes. In steps 23-26 they are
aligned and treed WITH the samples, so every sample tip sits among species of
known identity and the tree can be rooted between the tick families.

Each reference is annotated the SAME way as the samples, not from its GenBank
features: tblastn/blastn with the step-20 probes, HSPs merged into loci by step
21's own functions, PCGs framed and stop-trimmed by the same
`frame_mito_pcg` as step 22. Two reasons:
  - many listed records (the Kelava et al. MT3718xx series) are GenBank
    "UNVERIFIED" sequences with no CDS features at all;
  - identical gene-boundary rules for references and samples keep the codon
    alignments free of an artificial references-vs-samples edge pattern.

Identity check: each downloaded record's GenBank title must contain the
expected species epithet (or the row's `genbank_check` text). A mismatch is
refused, not used -- one listed accession (NC_028521) turned out to be a frog.

Outputs
  results/mito/ref_genomes/fasta/<acc>.fna        downloaded record (cached)
  results/mito/ref_genomes/hits/<acc>.{pcg,rrna}.blast.tsv
  results/mito/ref_genomes/genes/<GENE>.fna/.faa  one FASTA per gene, all refs
  results/mito/ref_genomes/ref_loci.tsv           every locus called
  results/mito/ref_genomes/ref_genome_summary.tsv one row per panel row
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (INVERT_MITO_TABLE, Entrez, cfg_mito, frame_mito_pcg,  # noqa: E402
                         load_mito_panel, marker_name, mito_root, read_fasta, require,
                         revcomp, run, write_fasta)

# step 21's search functions, imported rather than copied so references and
# samples can never drift apart in how a gene locus is called
_spec = importlib.util.spec_from_file_location(
    "find_mito", Path(__file__).resolve().parent / "21_find_mito_candidates.py")
fm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fm)

PCG_GENES = ["COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L",
             "ND5", "ND6", "ATP6", "ATP8"]
SUMMARY_COLS = ["ref_id", "accession", "use", "status", "family", "genus", "subgenus",
                "species", "genbank_title", "length", "n_pcg", "n_genes",
                "genes_annotated", "flagged_genes", "note"]


def load_panel() -> list[dict]:
    path = ROOT / "metadata" / "mito_ref_genomes.tsv"
    with open(path) as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    for r in rows:
        for k in list(r):
            r[k] = (r[k] or "").strip()
    return rows


def ref_id(acc: str) -> str:
    """Tip/sample id for a reference. The REF_ prefix is how steps 24-26 tell
    references from samples."""
    return f"REF_{acc}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-download cached records")
    ap.add_argument("--force", action="store_true", help="re-run BLAST on cached records")
    args = ap.parse_args()

    mc = cfg_mito()
    cc, rc = mc["candidates"], mc["refs"]
    table = INVERT_MITO_TABLE
    min_len = int(mc["annotate"]["min_locus_len_nt"])
    root = mito_root()
    out = root / "ref_genomes"
    fa_dir, hit_dir, db_dir, gene_dir = out / "fasta", out / "hits", out / "blastdb", out / "genes"
    for d in (fa_dir, hit_dir, db_dir, gene_dir):
        d.mkdir(parents=True, exist_ok=True)

    makeblastdb, blastn, tblastn = require("makeblastdb"), require("blastn"), require("tblastn")
    panel_genes = load_mito_panel()
    thresholds = {r["gene"]: (r["min_aln_len"], r["min_pident"]) for r in panel_genes}
    pcg_probes, rrna_probes = fm.build_combined_probes(root / "refs", panel_genes)

    api_key = os.environ.get(rc.get("api_key_env", "NCBI_API_KEY")) or None
    ez = Entrez(rc["email"], api_key, tool="barker_ticks_mito")

    panel = load_panel()
    summary, all_loci = [], []
    by_gene_nt: dict[str, dict[str, str]] = defaultdict(dict)
    by_gene_aa: dict[str, dict[str, str]] = defaultdict(dict)
    use = [r for r in panel if r["use"] == "yes"]

    for i, r in enumerate(use, 1):
        acc, sp = r["accession"], r["species"]
        rid = ref_id(acc)
        tag = f"[{i:>3}/{len(use)}] {acc} {sp}"
        row = {k: r.get(k, "") for k in ("accession", "use", "family", "genus", "subgenus",
                                         "species", "note")}
        row.update(ref_id=rid, genbank_title="", length="", n_pcg=0, n_genes=0,
                   genes_annotated="", flagged_genes="")

        fa = fa_dir / f"{acc}.fna"
        if args.refresh or not fa.exists():
            txt = ez.efetch_fasta("nuccore", [acc])
            if not txt.startswith(">"):
                print(f"{tag}: download FAILED", file=sys.stderr)
                summary.append({**row, "status": "download_failed"})
                continue
            fa.write_text(txt)
        recs = read_fasta(fa)
        title, seq = next(iter(recs.items()))
        row["genbank_title"], row["length"] = title, len(seq)

        expect = r["genbank_check"] or sp.split()[-1]
        if expect.lower() not in title.lower():
            print(f"{tag}: REFUSED -- GenBank title does not mention '{expect}': {title}",
                  file=sys.stderr)
            summary.append({**row, "status": "identity_mismatch"})
            continue

        # one-record BLAST db, same search settings as the samples (step 21)
        db = db_dir / acc
        if args.force or args.refresh or not Path(f"{db}.nin").exists():
            run([makeblastdb, "-in", str(fa), "-dbtype", "nucl", "-out", str(db)])
        hits = {}
        for kind, exe, probes, ev, extra in (
                ("pcg", tblastn, pcg_probes, cc["evalue_pcg"], []),
                ("rrna", blastn, rrna_probes, cc["evalue_rrna"], ["-task", "blastn"])):
            hf = hit_dir / f"{acc}.{kind}.blast.tsv"
            if args.force or args.refresh or not hf.exists():
                run([exe, *extra, "-query", str(probes), "-db", str(db), "-out", str(hf),
                     "-outfmt", fm.FMT, "-evalue", str(ev),
                     "-num_threads", str(cc["threads"]), "-max_target_seqs", "100000"])
            hits[kind] = fm.parse_hits(hf, thresholds)
        loci = (fm.assemble_loci(hits["pcg"], cc["merge_gap_pcg"], "pcg")
                + fm.assemble_loci(hits["rrna"], cc["merge_gap_rrna"], "rrna"))
        loci = fm.cap_one_locus_per_gene(loci)
        for l in loci:
            l["sample_id"], l["is_selected_contig"] = rid, 1
        all_loci.extend(loci)

        genes, flagged = [], []
        for l in (l for l in loci if l["status"] == "kept"):
            gene = l["gene"]
            lo, hi = int(l["locus_start"]), int(l["locus_end"])
            if hi - lo + 1 < min_len:
                continue
            s = seq[lo - 1:hi]
            if l["strand"] == "-":
                s = revcomp(s)
            aa, bad = None, False
            if gene in PCG_GENES:
                s, aa, bad = frame_mito_pcg(s, table)
            if bad:
                # a published, complete genome with a broken gene here means a
                # mis-called locus, not a pseudogene -- keep it out of the tree
                flagged.append(gene)
                continue
            name = marker_name(rid, sp, "mitochondrion", gene, acc, len(s))
            by_gene_nt[gene][name] = s
            if aa is not None:
                by_gene_aa[gene][name] = aa
            genes.append(gene)

        row.update(status="ok", n_pcg=sum(g in PCG_GENES for g in genes),
                   n_genes=len(genes), genes_annotated=",".join(sorted(genes)),
                   flagged_genes=",".join(sorted(flagged)))
        summary.append(row)
        print(f"{tag}: {row['n_pcg']}/13 PCGs, {len(genes) - row['n_pcg']}/2 rRNAs"
              f"{'  flagged: ' + row['flagged_genes'] if flagged else ''}", file=sys.stderr)

    for r in panel:
        if r["use"] != "yes":
            summary.append({"ref_id": ref_id(r["accession"]), "accession": r["accession"],
                            "use": "no", "status": "not_used", "family": r["family"],
                            "genus": r["genus"], "subgenus": r["subgenus"],
                            "species": r["species"], "note": r["note"]})

    for f in list(gene_dir.glob("*.fna")) + list(gene_dir.glob("*.faa")):
        f.unlink()
    for gene, seqs in by_gene_nt.items():
        write_fasta(gene_dir / f"{gene}.fna", seqs)
    for gene, seqs in by_gene_aa.items():
        write_fasta(gene_dir / f"{gene}.faa", seqs)
    write_tsv(out / "ref_loci.tsv", all_loci, fm.LOCI_COLS)
    write_tsv(out / "ref_genome_summary.tsv", summary, SUMMARY_COLS)

    ok = [s for s in summary if s["status"] == "ok"]
    print(f"\n{len(ok)}/{len(use)} reference mitogenomes downloaded and annotated "
          f"({len(panel) - len(use)} panel rows not used, see the table)", file=sys.stderr)
    bad = [s for s in summary if s["status"] not in ("ok", "not_used")]
    for s in bad:
        print(f"  {s['accession']}: {s['status']}", file=sys.stderr)
    for gene in PCG_GENES + ["rrnL", "rrnS"]:
        print(f"  {gene:<8} {len(by_gene_nt.get(gene, {})):>4} refs", file=sys.stderr)


if __name__ == "__main__":
    main()
