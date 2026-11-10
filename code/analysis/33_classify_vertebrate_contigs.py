#!/usr/bin/env python3
"""Step 33 -- name every contig that carries vertebrate mitochondrial or
nuclear-rRNA sequence, and decide which evidence tier it belongs to.

Input per sample: the step-30 cache of vertebrate HSPs, the contig's nt LCA
(lca/<s>.blast.lca) and its mapped-read count (mapping/<s>.coverage.tsv).

Per contig, separately for mtDNA and nuclear-rRNA subjects
(bloodmeal_lib.locus_class -- strict, so nuclear genes for mitochondrial
proteins and PREDICTED transcripts are not counted):

  1. Competitive check. The best vertebrate HSP must fall inside the window the
     contig's own nt LCA used: bitscore >= (1 - lca_run.top_percent/100) x the
     contig's top nt bitscore over ALL taxa. A tick mitogenome contig that also
     hits a vertebrate COX1 at 75% is far outside that window and is rejected
     (`not_best_hit`), as is a bacterial contig with a short vertebrate hit.
  2. Name. LCA of the vertebrate subjects within the same window of the best
     vertebrate HSP (lca.py's rule, applied to this locus' hits only).
  3. Identity cap. The LCA is truncated to the deepest rank the best HSP's
     identity supports (`bloodmeal.contig.identity_caps`): nt has a mitogenome
     for few Australian squamates, so a 91% hit to the one Varanus that has
     one must read "Varanidae", not "Varanus komodoensis". 18S uses its own,
     much stricter caps (`identity_caps_18s`).
  4. Tier. mtDNA with best HSP >= tier1_min_aln -> tier1_mito; else COX1
     >= tier2_min_aln -> tier2_COI; 18S >= tier2_min_aln -> tier2_18S;
     anything shorter -> below_floor (listed, never used).

Every vertebrate-mito contig that passes the competitive check is then
blasted against the RefSeq vertebrate mitogenomes (step 31) to list the genes
it spans and give a curated second-opinion best hit.

Outputs (results/bloodmeal/)
  tables/vertebrate_contigs.tsv          one row per contig x locus kind, incl. rejected
  tables/nuclear_vertebrate_background.tsv  per sample x class: reads on contigs
                                         whose nt LCA is vertebrate (any locus) --
                                         the nuclear signal, context only
  contigs/vertebrate_mito_contigs.fa     the mtDNA contigs that passed
  contigs/refseq_hits.tsv                raw blastn vs RefSeq (cached)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_aliases, load_config, resolve_artefact, write_tsv  # noqa: E402
from bloodmeal_lib import (RANKS, Taxonomy, bloodmeal_root, cfg_bloodmeal, contig_len,  # noqa: E402
                           identity_rank, locus_class, read_tsv)
from markers_lib import require, run  # noqa: E402

COLS = ["sample_id", "contig", "contig_len", "mapped_reads", "kind", "genes_nt", "tier",
        "status", "best_pident", "best_aln_len", "best_bitscore", "contig_top_bitscore",
        "best_subject", "best_subject_title", "n_hits_in_window", "lca_rank_raw", "lca_name_raw",
        "identity_rank", "call_taxid", "call_rank", "call_name",
        *[f"call_{r}" for r in RANKS],
        "contig_nt_lca", "contig_nt_lca_rank",
        "refseq_best", "refseq_organism", "refseq_pident", "refseq_aln_len", "refseq_genes"]

BG_COLS = ["sample_id", "class", "n_contigs", "mapped_reads", "share_of_vertebrate_reads"]


def load_lca(path: Path | None) -> dict[str, tuple[str, str, str, float]]:
    out = {}
    if path is None:
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["query"]] = (r["lca_taxid"], r["lca_name"], r["lca_rank"],
                               float(r["top_bitscore"] or 0))
    return out


def load_reads(path: Path | None, wanted: set[str] | None) -> dict[str, int]:
    out = {}
    if path is None:
        return out
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.split("\t", 3)
            if wanted is None or p[0] in wanted:
                out[p[0]] = int(float(p[2]))
    return out


def classify_sample(sid: str, cache: Path, lca: dict, reads: dict, tax: Taxonomy,
                    cc: dict, top_pct: float) -> list[dict]:
    hits: dict[tuple[str, str], list[dict]] = defaultdict(list)
    with gzip.open(cache, "rt") as fh:
        for h in csv.DictReader(fh, delimiter="\t"):
            kind, gene = locus_class(h["saccver"], h["stitle"])
            if kind == "other" or float(h["pident"]) < float(cc["min_pident"]):
                continue
            if kind == "nuc_rRNA" and gene != "18S":
                kind = "nuc_rRNA_other"
            h["gene"] = gene
            for k in ("pident", "bitscore"):
                h[k] = float(h[k])
            h["length"] = int(h["length"])
            hits[(h["contig"], kind)].append(h)

    caps = cc["identity_caps"]
    rows = []
    for (contig, kind), hs in hits.items():
        best = max(hs, key=lambda h: h["bitscore"])
        lt = lca.get(contig, ("", "", "", best["bitscore"]))
        top = max(lt[3], best["bitscore"])
        window = [h for h in hs if h["bitscore"] >= best["bitscore"] * (1 - top_pct / 100)]
        raw = tax.lca({h["staxid"] for h in window})
        id_rank = identity_rank(best["pident"],
                                cc["identity_caps_18s"] if kind == "nuc_rRNA" else caps)
        called = tax.cap(raw, id_rank)
        lin, raw_lin = tax.lineage(called), tax.lineage(raw)

        genes = sorted({h["gene"] for h in hs if h["gene"]})
        if best["bitscore"] < top * (1 - top_pct / 100):
            status = "not_best_hit"
        elif kind == "nuc_rRNA_other":
            status = "not_used"
        else:
            status = "ok"
        if kind == "mtDNA":
            tier = ("tier1_mito" if best["length"] >= cc["tier1_min_aln"] else
                    "tier2_COI" if "COX1" in genes and best["length"] >= cc["tier2_min_aln"]
                    else "below_floor")
        elif kind == "nuc_rRNA":
            tier = "tier2_18S" if best["length"] >= cc["tier2_min_aln"] else "below_floor"
        else:
            tier = ""
        rows.append({
            "sample_id": sid, "contig": contig, "contig_len": contig_len(contig),
            "mapped_reads": reads.get(contig, 0), "kind": kind, "genes_nt": ",".join(genes),
            "tier": tier, "status": status, "best_pident": best["pident"],
            "best_aln_len": best["length"], "best_bitscore": best["bitscore"],
            "contig_top_bitscore": lt[3], "best_subject": best["saccver"],
            "best_subject_title": best["stitle"][:150], "n_hits_in_window": len(window),
            "lca_rank_raw": raw_lin["rank"], "lca_name_raw": raw_lin["name"],
            "identity_rank": id_rank, "call_taxid": lin["taxid"], "call_rank": lin["rank"],
            "call_name": lin["name"], **{f"call_{r}": lin[r] for r in RANKS},
            "contig_nt_lca": lt[1], "contig_nt_lca_rank": lt[2]})
    return rows


def nuclear_background(sid: str, lca_path: Path | None, reads: dict,
                       lineage: dict[str, str]) -> list[dict]:
    """Reads on contigs whose whole-contig nt LCA is vertebrate, by class."""
    if lca_path is None:
        return []
    by: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    with open(lca_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cls = lineage.get(r["lca_taxid"])
            if cls is None:
                continue
            by[cls or "Vertebrata (unresolved)"][0] += 1
            by[cls or "Vertebrata (unresolved)"][1] += reads.get(r["query"], 0)
    total = sum(v[1] for v in by.values()) or 1
    return [{"sample_id": sid, "class": c, "n_contigs": v[0], "mapped_reads": v[1],
             "share_of_vertebrate_reads": round(v[1] / total, 4)}
            for c, v in sorted(by.items(), key=lambda x: -x[1][1])]


def vertebrate_classes() -> dict[str, str]:
    """LCA taxid -> class name for every vertebrate taxid in step 03's lineage table."""
    out = {}
    for r in read_tsv(ROOT / "results" / "taxonomy" / "taxon_lineage.tsv"):
        if r["phylum"] == "Chordata" and r.get("is_bloodmeal") == "1":
            out[r["taxid"]] = r["class"] if r["class"] not in ("Ascidiacea", "Leptocardii",
                                                               "Appendicularia") else None
    return {k: v for k, v in out.items() if v is not None}


# ---------------------------------------------------------------------------
# RefSeq second opinion + gene span
# ---------------------------------------------------------------------------
def extract_contigs(want: dict[str, set[str]], out_fa: Path) -> None:
    aliases = load_aliases()
    seqkit = require("seqkit")
    with open(out_fa, "w") as fo:
        for sid, ids in sorted(want.items()):
            fa = resolve_artefact(sid, ROOT / "contigs", ".contigs.fa", aliases)
            if fa is None:
                print(f"  {sid}: no contigs file, {len(ids)} contigs not extracted",
                      file=sys.stderr)
                continue
            idf = out_fa.with_suffix(".ids")
            idf.write_text("\n".join(sorted(ids)) + "\n")
            r = subprocess.run([seqkit, "grep", "-f", str(idf), str(fa)],
                               capture_output=True, text=True, check=True)
            # prefix the sample so contig ids are unique across samples
            fo.write(r.stdout.replace(">", f">{sid}|"))
            idf.unlink()


def refseq_blast(fa: Path, out: Path, cc: dict) -> None:
    db = bloodmeal_root() / "refs" / "blastdb" / "vert_mito_refseq"
    run([require("blastn"), "-task", cc.get("refseq_task", "dc-megablast"),
         "-query", str(fa), "-db", str(db), "-out", str(out),
         "-evalue", str(cc.get("refseq_evalue", 1e-10)), "-max_target_seqs", "5",
         "-num_threads", str(cc.get("threads", 2)),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"])


def refseq_annotation(hits_path: Path, cc: dict) -> dict[tuple[str, str], dict]:
    root = bloodmeal_root() / "refs"
    org = {r["accession"]: r["organism"] for r in read_tsv(root / "vert_mito_refseq.tsv")}
    by_q: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    with open(hits_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            by_q[p[0]][p[1]].append((float(p[2]), int(p[3]), int(p[6]), int(p[7]), float(p[9])))
    best_subj = {}
    for q, subj in by_q.items():
        s, hs = max(subj.items(), key=lambda x: sum(h[4] for h in x[1]))
        best_subj[q] = (s, hs)

    feats: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    need = {s for s, _ in best_subj.values()}
    for r in read_tsv(root / "vert_mito_features.tsv"):
        if r["accession"] in need and r["gene"] != "tRNA":
            feats[r["accession"]].append((r["gene"], int(r["start"]), int(r["end"])))

    min_ov = float(cc.get("min_gene_overlap", 0.3))
    out = {}
    for q, (s, hs) in best_subj.items():
        spans = [(min(a, b), max(a, b)) for _, _, a, b, _ in hs]
        genes = []
        for g, gs, ge in feats.get(s, []):
            cov = sum(max(0, min(ge, e) - max(gs, b) + 1) for b, e in spans)
            if cov / (ge - gs + 1) >= min_ov:
                genes.append(g)
        aln = sum(h[1] for h in hs)
        pid = sum(h[0] * h[1] for h in hs) / aln
        sid, contig = q.split("|", 1)
        out[(sid, contig)] = {"refseq_best": s, "refseq_organism": org.get(s, ""),
                              "refseq_pident": round(pid, 2), "refseq_aln_len": aln,
                              "refseq_genes": ",".join(dict.fromkeys(genes))}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", action="append", help="limit to this sample (repeatable)")
    ap.add_argument("--force-refseq", action="store_true", help="re-run the RefSeq blastn")
    args = ap.parse_args()

    cfg = load_config()
    cc = cfg_bloodmeal()["contig"]
    top_pct = float(cfg["lca_run"]["top_percent"])
    root = bloodmeal_root()
    aliases = load_aliases()
    caches = sorted((root / "vert_hits").glob("*.tsv.gz"))
    if args.sample:
        caches = [c for c in caches if c.name[:-7] in set(args.sample)]
    if not caches:
        sys.exit("ERROR: no step-30 caches under results/bloodmeal/vert_hits/ -- run step 30")

    print("loading taxonomy", file=sys.stderr)
    tax = Taxonomy()
    vclass = vertebrate_classes()
    rows, bg = [], []
    for i, cache in enumerate(caches, 1):
        sid = cache.name[:-len(".tsv.gz")]
        lca_p = resolve_artefact(sid, ROOT / "lca", ".blast.lca", aliases)
        cov_p = resolve_artefact(sid, ROOT / "mapping", ".coverage.tsv", aliases)
        lca = load_lca(lca_p)
        reads = load_reads(cov_p, None)
        r = classify_sample(sid, cache, lca, reads, tax, cc, top_pct)
        rows += r
        bg += nuclear_background(sid, lca_p, reads, vclass)
        ok = sum(x["status"] == "ok" and x["kind"] == "mtDNA" for x in r)
        print(f"  [{i}/{len(caches)}] {sid}: {len(r)} vertebrate loci, {ok} mtDNA ok",
              file=sys.stderr, flush=True)
        del reads, lca

    # ---- RefSeq gene span / second opinion ---------------------------------
    cdir = root / "contigs"
    cdir.mkdir(parents=True, exist_ok=True)
    fa, hits = cdir / "vertebrate_mito_contigs.fa", cdir / "refseq_hits.tsv"
    if args.force_refseq:
        fa.unlink(missing_ok=True)
        hits.unlink(missing_ok=True)
    # Incremental: only contigs not already extracted + blasted are processed,
    # so adding samples does not re-read 70 GB of assemblies.
    have = set()
    if fa.exists() and hits.exists():
        with open(fa) as fh:
            have = {line[1:].split()[0] for line in fh if line.startswith(">")}
    want: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        if (r["kind"] == "mtDNA" and r["status"] == "ok"
                and f"{r['sample_id']}|{r['contig']}" not in have):
            want[r["sample_id"]].add(r["contig"])
    if want:
        print(f"extracting {sum(map(len, want.values()))} new vertebrate mito contigs "
              f"from {len(want)} samples", file=sys.stderr)
        new_fa, new_hits = fa.with_suffix(".new.fa"), hits.with_suffix(".new.tsv")
        extract_contigs(want, new_fa)
        print("blastn vs RefSeq vertebrate mitogenomes", file=sys.stderr)
        refseq_blast(new_fa, new_hits, cc)
        for new, old in ((new_fa, fa), (new_hits, hits)):
            with open(old, "a") as fo, open(new) as fi:
                fo.write(fi.read())
            new.unlink()
    ann = refseq_annotation(hits, cc) if hits.exists() else {}
    for r in rows:
        r.update(ann.get((r["sample_id"], r["contig"]), {}))

    rows.sort(key=lambda r: (r["sample_id"], r["kind"], -r["mapped_reads"]))
    tables = root / "tables"
    if args.sample:  # partial run: do not overwrite the cohort tables
        tables = root / "tables" / "partial"
    write_tsv(tables / "vertebrate_contigs.tsv", rows, COLS)
    write_tsv(tables / "nuclear_vertebrate_background.tsv", bg, BG_COLS)


if __name__ == "__main__":
    main()
