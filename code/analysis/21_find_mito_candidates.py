#!/usr/bin/env python3
"""Step 21 -- pick ONE mitogenome-candidate contig per sample, by homology.

For each sample: keep contigs in the plausible mitogenome length window
(`mito.candidates.min/max_contig_len`), tblastn the 13 protein-coding gene
probes and blastn the 2 rRNA probes (metadata/mito_genes.tsv, fetched by step
20) against them, and merge nearby HSPs of the same gene into a locus the same
way the marker workflow does (step 13). Unlike step 13, loci are NOT resolved
across different genes by overlap: mitochondrial PCGs are single-copy and
share no homology with each other or the rRNAs, so an HSP-boundary overlap
between two different genes is normal adjacent-gene imprecision, not a
competing paralog call -- see the comment at the resolve-overlaps call site
below for the data-loss bug that assumption caused when first ported over.

Unlike the marker workflow, the target here is not "every contig for this
taxon" but ONE winning contig per sample: the contig carrying the most distinct
mitochondrial genes, ties broken by contig coverage, then summed bitscore,
then length. Before
that ranking, each candidate contig is checked against its own nt-BLAST LCA
from the main pipeline (lca/<sample>.blast.lca, classified via
results/taxonomy/taxon_lineage.tsv): a contig whose LCA is anything other than
a tick (order `mito.candidates.lca_filter.tick_order`) is not eligible -- see
classify_contigs() for why. That
contig is written out as this sample's mitogenome candidate; every gene locus
found on it (this step's homology calls) is also this sample's provisional
annotation, used directly by step 22 unless real MITOS2 is available.

Why homology across the whole order, not blastn alone: the sample set spans
Ixodidae and Argasidae, tick families separated by >200 My. Nucleotide identity
between them is often too low for blastn to find the gene at all; tblastn on
the translated, well-conserved protein carries the search across that gap.

Outputs
  results/mito/candidates/filtered_contigs/<sample>.fa   length-windowed contigs
  results/mito/candidates/hits/<sample>.{pcg,rrna}.blast.tsv
  results/mito/mito_loci.tsv         every kept gene locus, all candidate contigs,
                                     with is_selected_contig marking the winner
  results/mito/mito_selection.tsv    one row per sample: winning contig + evidence
  results/mito/genome/<sample>.mito.fna   the winning contig, relabelled
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ROOT, canonical, load_aliases, load_excluded, resolve_artefact,  # noqa: E402
                    write_tsv)
from markers_lib import (cfg_mito, load_mito_panel, mito_root, read_fasta, require,  # noqa: E402
                         run, tick_species_map, write_fasta)

FMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
       "evalue bitscore qlen slen")
COLS = FMT.split()[1:]

LOCI_COLS = ["sample_id", "gene", "gene_type", "contig_id", "contig_len",
             "locus_start", "locus_end", "strand", "locus_len", "n_hsps",
             "sum_bitscore", "max_bitscore", "best_pident", "best_ref",
             "ref_len", "ref_cov_pct", "status", "suppressed_by",
             "contig_lca", "contig_lca_class", "is_selected_contig"]

SEL_COLS = ["sample_id", "tick_species", "tick_species_source", "status",
            "n_candidate_contigs", "n_gene_contigs", "n_non_tick_excluded",
            "non_tick_excluded", "selected_contig", "selected_contig_lca",
            "selected_contig_lca_class", "contig_len", "n_genes_found",
            "genes_found", "sum_bitscore", "fasta"]

# Eligibility tier of a contig by its nt-BLAST LCA (lower = preferred).
# Contigs with no nt-BLAST hit at all are allowed but rank below confirmed
# tick contigs: several complete Ornithodoros mitogenomes on this cohort have
# no nt hit (no close relative in nt), so "no hit" cannot mean "not tick".
LCA_TIER = {"tick": 0, "no_blast_hit": 1, "no_lca_file": 1, "lca_filter_off": 1}


def build_combined_probes(refs_dir: Path, panel: list[dict]) -> tuple[Path, Path]:
    """One query FASTA per BLAST program, headers prefixed `GENE::id` so a hit's
    qseqid recovers which gene it came from without 15 separate blast calls."""
    pcg_out, rrna_out = refs_dir / "_combined_pcg.faa", refs_dir / "_combined_rrna.fna"
    pcg, rrna = {}, {}
    for row in panel:
        ext = "faa" if row["type"] == "pcg" else "fna"
        f = refs_dir / f"{row['key']}.{ext}"
        if not f.exists():
            sys.exit(f"ERROR: {f.relative_to(ROOT)} missing -- run step 20 first.")
        target = pcg if row["type"] == "pcg" else rrna
        for name, seq in read_fasta(f).items():
            target[f"{row['gene']}::{name}"] = seq
    write_fasta(pcg_out, pcg)
    write_fasta(rrna_out, rrna)
    return pcg_out, rrna_out


def parse_hits(path: Path, thresholds: dict[str, tuple[int, float]]) -> list[dict]:
    """Read BLAST tabular output, keeping HSPs that clear THAT gene's own
    min_aln_len/min_pident (looked up from the `GENE::` prefix on qseqid)."""
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(COLS):
                continue
            h = dict(zip(COLS, p))
            gene = h["qseqid"].split("::", 1)[0]
            min_len, min_pid = thresholds.get(gene, (0, 0.0))
            h["pident"] = float(h["pident"]); h["length"] = int(h["length"])
            h["bitscore"] = float(h["bitscore"]); h["evalue"] = float(h["evalue"])
            h["sstart"] = int(h["sstart"]); h["send"] = int(h["send"])
            h["qlen"] = int(h["qlen"]); h["slen"] = int(h["slen"])
            h["gene"] = gene
            if h["length"] < min_len or h["pident"] < min_pid:
                continue
            out.append(h)
    return out


def assemble_loci(hits: list[dict], merge_gap: int, gene_type: str) -> list[dict]:
    """Merge HSPs of the SAME gene into a locus, per contig and strand."""
    groups: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for h in hits:
        strand = "+" if h["sstart"] <= h["send"] else "-"
        groups[(h["sseqid"], strand, h["gene"])].append(h)

    out = []
    for (contig, strand, gene), hs in groups.items():
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
            unit = 3 if gene_type == "pcg" else 1
            cov = 100.0 * (hi - lo + 1) / max(1, top["qlen"] * unit)
            out.append({
                "gene": gene, "gene_type": gene_type, "contig_id": contig,
                "contig_len": top["slen"], "locus_start": lo, "locus_end": hi,
                "strand": strand, "locus_len": hi - lo + 1, "n_hsps": len(hsps),
                "sum_bitscore": round(sum(h["bitscore"] for h in hsps), 1),
                "max_bitscore": round(top["bitscore"], 1),
                "best_pident": round(max(h["pident"] for h in hsps), 1),
                "best_ref": top["qseqid"].split("::", 1)[1], "ref_len": top["qlen"],
                "ref_cov_pct": round(min(cov, 999.9), 1),
                "status": "kept", "suppressed_by": "",
            })
    return out


def cap_one_locus_per_gene(loci: list[dict]) -> list[dict]:
    """At most one KEPT locus per (contig, gene).

    A gene is single-copy in a mitochondrial genome; two disjoint loci for the
    same gene on the same contig (rare -- a weak secondary hit elsewhere, or a
    genuine internal assembly break) would otherwise both survive
    resolve_overlaps, since it only catches OVERLAPPING competing claims. Two
    surviving loci for one gene would collide under the same output name in
    step 22 and silently overwrite each other; keep only the stronger one.
    """
    by_cg: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for l in loci:
        if l["status"] == "kept":
            by_cg[(l["contig_id"], l["gene"])].append(l)
    for group in by_cg.values():
        if len(group) < 2:
            continue
        group.sort(key=lambda x: -x["sum_bitscore"])
        for extra in group[1:]:
            extra["status"] = "duplicate_suppressed"
            extra["suppressed_by"] = f"{group[0]['gene']}@bits{group[0]['max_bitscore']}"
    return loci


def load_lineage(path: Path) -> dict[str, dict]:
    with open(path) as fh:
        return {r["taxid"]: r for r in csv.DictReader(fh, delimiter="\t")}


def classify_contigs(sid: str, contigs: set[str], lineage: dict[str, dict] | None,
                     tick_order: str, aliases: dict[str, str]) -> dict[str, tuple[str, str]]:
    """contig -> (lca_name, class) from the sample's own whole-assembly nt BLAST.

    class is `tick` (LCA inside `tick_order`), `non_tick` (any other LCA,
    including one above the tick order such as Bilateria), `no_blast_hit`, or
    `no_lca_file`.

    Why this gate exists: the probes are tick proteins, but a vertebrate
    mitogenome carries the same 13 genes and clears the E-value gate at
    40-60% protein identity, so a blood-meal host mitogenome (Varanus,
    Petrogale, Isoodon on this cohort) can out-score a tick fragment on gene
    count alone. Bacterial contigs with 2-3 weak hits win the same way in
    samples with no real tick mitogenome contig. The tblastn evidence cannot
    tell these apart; the contig's nt-BLAST LCA can.
    """
    if lineage is None:
        return {c: ("", "lca_filter_off") for c in contigs}
    f = resolve_artefact(sid, ROOT / "lca", ".blast.lca", aliases)
    if f is None:
        return {c: ("", "no_lca_file") for c in contigs}
    out = {c: ("", "no_blast_hit") for c in contigs}
    with open(f) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["query"] in out:
                t = lineage.get(r["lca_taxid"], {})
                cls = "tick" if t.get("order") == tick_order else "non_tick"
                out[r["query"]] = (r["lca_name"], cls)
    return out


def contig_cov(contig_id: str) -> float:
    """SPAdes k-mer coverage from the NODE_..._cov_X name; 0 if absent."""
    m = re.search(r"_cov_([0-9.]+)$", contig_id)
    return float(m.group(1)) if m else 0.0


def choose_contig(loci: list[dict]) -> str | None:
    """Winning contig = best LCA tier (confirmed tick before no-nt-hit), then
    most distinct genes, then highest coverage, then highest summed bitscore,
    then longest. `loci` here is already filtered to status == kept and to
    LCA-eligible contigs.

    Coverage before bitscore: when two complete tick mitogenomes of different
    genera are in one assembly (a mixed or contaminated tube), bitscore only
    says which is closer to the probe set, while coverage says which tick
    contributed most of the DNA. 8164-M1-Amorelia picked a 5x Bothriocroton
    mitogenome over a 10x Amblyomma one on bitscore alone."""
    per_contig: dict[str, list[dict]] = defaultdict(list)
    for l in loci:
        per_contig[l["contig_id"]].append(l)
    if not per_contig:
        return None
    ranked = sorted(per_contig.items(), key=lambda kv: (
        LCA_TIER[kv[1][0]["contig_lca_class"]],
        -len({l["gene"] for l in kv[1]}),
        -contig_cov(kv[0]),
        -sum(l["sum_bitscore"] for l in kv[1]),
        -kv[1][0]["contig_len"],
    ))
    return ranked[0][0]


def find_sample_fastas() -> dict[str, Path]:
    """Every assembled sample except those in metadata/excluded_samples.tsv."""
    aliases, excluded = load_aliases(), load_excluded()
    out: dict[str, Path] = {}
    for f in sorted((ROOT / "contigs").glob("*.contigs.fa")):
        sid = canonical(f.name[: -len(".contigs.fa")], aliases)
        if sid not in excluded:
            out.setdefault(sid, f)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append", help="limit to this sample (repeatable)")
    ap.add_argument("--force", action="store_true", help="re-run filtering/BLAST")
    args = ap.parse_args()

    mc = cfg_mito()
    cc = mc["candidates"]
    root = mito_root()
    refs_dir, cand_dir = root / "refs", root / "candidates"
    filt_dir, hit_dir, db_dir = cand_dir / "filtered_contigs", cand_dir / "hits", cand_dir / "blastdb"
    genome_dir = root / "genome"
    for d in (filt_dir, hit_dir, db_dir, genome_dir):
        d.mkdir(parents=True, exist_ok=True)

    seqkit, makeblastdb = require("seqkit"), require("makeblastdb")
    blastn, tblastn = require("blastn"), require("tblastn")

    panel = load_mito_panel()
    thresholds = {r["gene"]: (r["min_aln_len"], r["min_pident"]) for r in panel}
    pcg_probes, rrna_probes = build_combined_probes(refs_dir, panel)

    samples = find_sample_fastas()
    if args.sample:
        samples = {s: p for s, p in samples.items() if s in set(args.sample)}
    species = tick_species_map()
    aliases = load_aliases()
    lf = cc.get("lca_filter", {}) or {}
    lineage = None
    if lf.get("enabled", True):
        lin_path = ROOT / "results" / "taxonomy" / "taxon_lineage.tsv"
        if not lin_path.exists():
            sys.exit(f"ERROR: {lin_path.relative_to(ROOT)} missing -- run step 03 first, "
                     "or set mito.candidates.lca_filter.enabled: false.")
        lineage = load_lineage(lin_path)
    tick_order = lf.get("tick_order", "Ixodida")
    print(f"{len(samples)} samples with assemblies", file=sys.stderr)

    selection: list[dict] = []
    all_loci: list[dict] = []

    for i, (sid, fasta) in enumerate(sorted(samples.items()), 1):
        sp, sp_src = species.get(sid, ("unknown_species", "not_in_metadata"))
        filt_fa = filt_dir / f"{sid}.fa"
        if args.force or not filt_fa.exists():
            run([seqkit, "seq", "-m", str(cc["min_contig_len"]), "-M", str(cc["max_contig_len"]),
                 str(fasta), "-o", str(filt_fa)])
        n_candidates = len(read_fasta(filt_fa))
        tag = f"[{i:>3}/{len(samples)}] {sid}"

        if n_candidates == 0:
            print(f"{tag}: no contigs in {cc['min_contig_len']}-{cc['max_contig_len']} bp window",
                  file=sys.stderr)
            selection.append({"sample_id": sid, "tick_species": sp,
                              "tick_species_source": sp_src, "status": "no_candidate_contigs",
                              "n_candidate_contigs": 0, "n_gene_contigs": 0,
                              "n_non_tick_excluded": 0, "non_tick_excluded": "",
                              "selected_contig": "", "selected_contig_lca": "",
                              "selected_contig_lca_class": "", "contig_len": "",
                              "n_genes_found": 0, "genes_found": "", "sum_bitscore": "", "fasta": ""})
            continue

        db = db_dir / sid
        if args.force or not Path(f"{db}.nin").exists():
            run([makeblastdb, "-in", str(filt_fa), "-dbtype", "nucl", "-out", str(db)])

        pcg_hits_f = hit_dir / f"{sid}.pcg.blast.tsv"
        if args.force or not pcg_hits_f.exists():
            run([tblastn, "-query", str(pcg_probes), "-db", str(db), "-out", str(pcg_hits_f),
                 "-outfmt", FMT, "-evalue", str(cc["evalue_pcg"]),
                 "-num_threads", str(cc["threads"]), "-max_target_seqs", "100000"])
        rrna_hits_f = hit_dir / f"{sid}.rrna.blast.tsv"
        if args.force or not rrna_hits_f.exists():
            run([blastn, "-task", "blastn", "-query", str(rrna_probes), "-db", str(db),
                 "-out", str(rrna_hits_f), "-outfmt", FMT, "-evalue", str(cc["evalue_rrna"]),
                 "-num_threads", str(cc["threads"]), "-max_target_seqs", "100000"])

        pcg_hits = parse_hits(pcg_hits_f, thresholds)
        rrna_hits = parse_hits(rrna_hits_f, thresholds)
        loci = (assemble_loci(pcg_hits, cc["merge_gap_pcg"], "pcg")
                + assemble_loci(rrna_hits, cc["merge_gap_rrna"], "rrna"))
        # NOTE: unlike the pathogen marker workflow (Rickettsia ompA/ompB, Anaplasma
        # msp2/p44), there is no cross-gene-family paralogy to resolve here --
        # mitochondrial PCGs are single-copy and not homologous to each other or
        # to the rRNAs. An earlier version applied the marker workflow's
        # highest-bitscore-wins overlap resolution here too, and it did real
        # damage: 146/146 "conflicts" it found were adjacent genes' HSP-approximate
        # boundaries overlapping by a few dozen bp (ND1 losing to its neighbour
        # rrnL in 87 samples, at 84% identity / E=1.7e-123 -- not a false hit,
        # just adjacent), silently deleting a strong, correct ND1 call in most
        # samples that had one. Overlap between genuinely different, adjacent
        # mitochondrial genes is normal and expected; each gene's own
        # independently-called span is kept as-is. cap_one_locus_per_gene still
        # protects against the one thing that IS a real error: the SAME gene
        # claiming two disjoint regions on one contig.
        loci = cap_one_locus_per_gene(loci)
        for l in loci:
            l["sample_id"] = sid

        lca = classify_contigs(sid, {l["contig_id"] for l in loci}, lineage, tick_order, aliases)
        for l in loci:
            l["contig_lca"], l["contig_lca_class"] = lca[l["contig_id"]]

        kept = [l for l in loci if l["status"] == "kept"]
        gene_contigs = {l["contig_id"] for l in kept}
        excluded = sorted(c for c in gene_contigs if lca[c][1] == "non_tick")
        winner = choose_contig([l for l in kept if l["contig_lca_class"] != "non_tick"])
        for l in loci:
            l["is_selected_contig"] = int(l["contig_id"] == winner)
        all_loci.extend(loci)
        excl_txt = ";".join(f"{c}={lca[c][0]}" for c in excluded)

        if winner is None:
            status = "no_tick_contig" if excluded else "no_genes_recovered"
            print(f"{tag}: {n_candidates} candidate contigs, "
                  + (f"genes only on {len(excluded)} non-tick contig(s): {excl_txt}"
                     if excluded else "no gene recovered on any"), file=sys.stderr)
            selection.append({"sample_id": sid, "tick_species": sp,
                              "tick_species_source": sp_src, "status": status,
                              "n_candidate_contigs": n_candidates,
                              "n_gene_contigs": len(gene_contigs),
                              "n_non_tick_excluded": len(excluded), "non_tick_excluded": excl_txt,
                              "selected_contig": "", "selected_contig_lca": "",
                              "selected_contig_lca_class": "",
                              "contig_len": "", "n_genes_found": 0, "genes_found": "",
                              "sum_bitscore": "", "fasta": ""})
            continue

        win_loci = [l for l in kept if l["contig_id"] == winner]
        genes_found = sorted({l["gene"] for l in win_loci})
        seq = read_fasta(filt_fa)[winner]
        out_fa = genome_dir / f"{sid}.mito.fna"
        write_fasta(out_fa, {f"{sid}|{sp.replace(' ', '_')}|{winner}|{len(seq)}bp": seq})

        print(f"{tag}: {n_candidates} candidates -> {winner} ({len(seq):,} bp, "
              f"{len(genes_found)}/13+2 genes, LCA {lca[winner][0] or lca[winner][1]})"
              + (f"; {len(excluded)} non-tick excluded" if excluded else ""), file=sys.stderr)
        selection.append({"sample_id": sid, "tick_species": sp,
                          "tick_species_source": sp_src, "status": "ok",
                          "n_candidate_contigs": n_candidates,
                          "n_gene_contigs": len(gene_contigs),
                          "n_non_tick_excluded": len(excluded), "non_tick_excluded": excl_txt,
                          "selected_contig": winner, "selected_contig_lca": lca[winner][0],
                          "selected_contig_lca_class": lca[winner][1],
                          "contig_len": len(seq), "n_genes_found": len(genes_found),
                          "genes_found": ",".join(genes_found),
                          "sum_bitscore": round(sum(l["sum_bitscore"] for l in win_loci), 1),
                          "fasta": str(out_fa.relative_to(ROOT))})

    # A --sample run must leave the other samples' rows in place -- both tables
    # are read whole downstream (step 22 onward).
    sel_path, loci_path = root / "mito_selection.tsv", root / "mito_loci.tsv"
    if args.sample:
        scope = set(args.sample)
        if sel_path.exists():
            with open(sel_path) as fh:
                others = [r for r in csv.DictReader(fh, delimiter="\t")
                         if r["sample_id"] not in scope]
            for r in others:
                r["n_genes_found"] = int(r["n_genes_found"] or 0)
            selection = others + selection
        if loci_path.exists():
            with open(loci_path) as fh:
                others = [r for r in csv.DictReader(fh, delimiter="\t")
                         if r["sample_id"] not in scope]
            for r in others:
                r["is_selected_contig"] = int(r["is_selected_contig"])
            all_loci = others + all_loci

    selection.sort(key=lambda r: r["sample_id"])
    all_loci.sort(key=lambda l: (l["sample_id"], l["gene"], l["contig_id"]))
    write_tsv(sel_path, selection, SEL_COLS)
    write_tsv(loci_path, all_loci, LOCI_COLS)

    ok = [s for s in selection if s["status"] == "ok"]
    print(f"\n{len(ok)}/{len(selection)} samples yielded a candidate mitogenome",
          file=sys.stderr)
    nt = [s for s in selection if s["status"] == "no_tick_contig"]
    if nt:
        print(f"{len(nt)} samples had mitochondrial genes only on non-tick contigs: "
              + ", ".join(s["sample_id"] for s in nt), file=sys.stderr)
    if ok:
        import statistics
        ng = [s["n_genes_found"] for s in ok]
        print(f"genes recovered per sample: median {statistics.median(ng)}, "
              f"min {min(ng)}, max {max(ng)} (of 15 total)", file=sys.stderr)


if __name__ == "__main__":
    main()
