#!/usr/bin/env python3
"""Step 33c -- re-name each host mtDNA contig against the candidate families'
full mitochondrial record sets (step 33b) plus the RefSeq vertebrate
mitogenomes (step 31).

Naming is by identity, compared like for like:

  1. The contig is cut into windows of `refine.window` bp. Every HSP of
     >= `min_aln` bp gives its subject's taxon an identity in each window it
     covers (>= `window_min_overlap` of the window, or of the HSP if that is
     shorter); a taxon keeps its best value per window.
  2. Candidates = taxa within `delta_pident` of the best identity in at least
     one window.
  3. Pairwise, over SHARED windows only: A eliminates B when A's mean identity
     there exceeds B's by more than `delta_pident`. This is what makes the
     comparison fair between records of different extent: a congener's complete
     mitogenome at 95% is eliminated by the true species' 99% CYTB record on
     the CYTB windows, and cannot win by covering 50 other windows where the
     true species simply has no sequence.
  4. The surviving taxa name the contig (their LCA), capped at the rank their
     identity supports (`contig.identity_caps`). Two survivors that never
     overlap (one known only from 12S, the other only from CYTB) cannot be
     separated and are named at their common ancestor.

Outputs (results/bloodmeal/)
  refs/blastdb/host_mito_refine.*      family records + RefSeq, one database
  contigs/refine_hits.tsv              raw blastn (cached)
  tables/refined_contig_calls.tsv      one row per contig
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv  # noqa: E402
from bloodmeal_lib import (RANKS, Taxonomy, bloodmeal_root, cfg_bloodmeal,  # noqa: E402
                           identity_rank, read_tsv)
from markers_lib import read_fasta, require, run, write_fasta  # noqa: E402

COLS = ["sample_id", "contig", "contig_len", "n_windows", "n_windows_hit", "n_candidates",
        "survivors", "refined_pident", "refined_identity_rank", "refined_taxid",
        "refined_rank", "refined_name", *[f"refined_{r}" for r in RANKS],
        "best_single_hit", "best_single_hit_species", "best_single_hit_pident",
        "best_single_hit_len"]


def build_db(root: Path) -> tuple[Path, dict[str, str]]:
    refs = root / "refs"
    taxid: dict[str, str] = {r["accession"]: r["taxid"]
                             for r in read_tsv(refs / "vert_mito_refseq.tsv")}
    inputs = [refs / "vert_mito_refseq.fa"]
    for tsv in sorted((refs / "families").glob("*.tsv")):
        if tsv.name.startswith("_"):
            continue
        for r in read_tsv(tsv):
            taxid.setdefault(r["accession"], r["taxid"])
        inputs.append(tsv.with_suffix(".fa"))
    fa = refs / "host_mito_refine.fa"
    db = refs / "blastdb" / "host_mito_refine"
    newest = max(p.stat().st_mtime for p in inputs if p.exists())
    if not Path(f"{db}.nsq").exists() or Path(f"{db}.nsq").stat().st_mtime < newest:
        # streamed: the inputs total a few hundred MB and this node has ~1 GB free
        seen: set[str] = set()
        with open(fa, "w") as fo:
            for p in inputs:
                if not p.exists():
                    continue
                keep = False
                with open(p) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            acc = line[1:].split()[0]
                            keep = acc not in seen
                            seen.add(acc)
                            line = f">{acc}\n"
                        if keep:
                            fo.write(line)
        run([require("makeblastdb"), "-in", str(fa), "-dbtype", "nucl", "-out", str(db),
             "-title", "candidate host family mito records + RefSeq vertebrate mitogenomes"])
        print(f"  refine database: {len(seen)} records", file=sys.stderr)
    return db, taxid


def vote(hsps: list[tuple], qlen: int, rc: dict, taxid_of: dict[str, str],
         tax: Taxonomy) -> dict:
    W = int(rc["window"])
    n = max(1, round(qlen / W))
    edges = [round(i * qlen / n) for i in range(n + 1)]
    min_aln = min(int(rc["min_aln"]), int(0.8 * qlen))
    ov = float(rc["window_min_overlap"])
    delta = float(rc["delta_pident"])

    ident: dict[str, dict[int, float]] = defaultdict(dict)   # taxid -> window -> pident
    best_single = None
    for sacc, pid, length, qs, qe, bits in hsps:
        if length < min_aln:
            continue
        t = taxid_of.get(sacc)
        if not t:
            continue
        if best_single is None or bits > best_single[5]:
            best_single = (sacc, pid, length, qs, qe, bits)
        lo, hi = min(qs, qe), max(qs, qe)
        for w in range(n):
            a, b = edges[w] + 1, edges[w + 1]
            # a short HSP lying inside one window still counts there
            if min(hi, b) - max(lo, a) + 1 >= ov * min(b - a + 1, hi - lo + 1):
                if pid > ident[t].get(w, -1):
                    ident[t][w] = pid
    if not ident:
        return {}
    windows = {w for d in ident.values() for w in d}
    cands = set()
    for w in windows:
        m = max(d[w] for d in ident.values() if w in d)
        cands |= {t for t, d in ident.items() if w in d and d[w] >= m - delta}

    survivors = []
    for b in cands:
        beaten = False
        for a in cands:
            if a == b:
                continue
            shared = ident[a].keys() & ident[b].keys()
            if shared and sum(ident[a][w] - ident[b][w] for w in shared) / len(shared) > delta:
                beaten = True
                break
        if not beaten:
            survivors.append(b)
    mean_id = {t: sum(ident[t].values()) / len(ident[t]) for t in survivors}
    pid = max(mean_id.values())
    raw = tax.lca(survivors)
    id_rank = identity_rank(pid, cfg_bloodmeal()["contig"]["identity_caps"])
    called = tax.cap(raw, id_rank)
    lin = tax.lineage(called)
    names = sorted(survivors, key=lambda t: -mean_id[t])
    bs_t = taxid_of.get(best_single[0], "")
    return {"n_windows": n, "n_windows_hit": len(windows), "n_candidates": len(cands),
            "survivors": "; ".join(f"{tax.names.get(int(t), t)} {mean_id[t]:.1f}% "
                                   f"({len(ident[t])}w)" for t in names[:6])
            + (f"; +{len(names) - 6} more" if len(names) > 6 else ""),
            "refined_pident": round(pid, 2), "refined_identity_rank": id_rank,
            "refined_taxid": lin["taxid"], "refined_rank": lin["rank"],
            "refined_name": lin["name"], **{f"refined_{r}": lin[r] for r in RANKS},
            "best_single_hit": best_single[0],
            "best_single_hit_species": tax.lineage(bs_t)["name"] if bs_t else "",
            "best_single_hit_pident": best_single[1], "best_single_hit_len": best_single[2]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="re-run blastn")
    args = ap.parse_args()

    rc = cfg_bloodmeal()["refine"]
    threads = str(cfg_bloodmeal()["contig"].get("threads", 2))
    root = bloodmeal_root()
    want = {(r["sample_id"], r["contig"]) for r in read_tsv(root / "tables" / "vertebrate_contigs.tsv")
            if r["kind"] == "mtDNA" and r["status"] == "ok"
            and r["tier"] in ("tier1_mito", "tier2_COI")}
    allfa = read_fasta(root / "contigs" / "vertebrate_mito_contigs.fa")
    seqs = {k.split()[0]: v for k, v in allfa.items()
            if tuple(k.split()[0].split("|", 1)) in want}
    missing = len(want) - len(seqs)
    if missing:
        print(f"WARNING: {missing} contigs not in vertebrate_mito_contigs.fa -- re-run step 33",
              file=sys.stderr)
    q = root / "contigs" / "refine_queries.fa"
    write_fasta(q, seqs)

    db, taxid_of = build_db(root)
    hits = root / "contigs" / "refine_hits.tsv"
    if args.force or not hits.exists() or hits.stat().st_mtime < q.stat().st_mtime:
        print(f"blastn {len(seqs)} host mito contigs vs refine database", file=sys.stderr)
        run([require("blastn"), "-task", rc["task"], "-query", str(q), "-db", str(db),
             "-out", str(hits), "-evalue", str(rc["evalue"]),
             "-max_target_seqs", str(rc["max_target_seqs"]), "-num_threads", threads,
             "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"])

    by_q: dict[str, list[tuple]] = defaultdict(list)
    with open(hits) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            by_q[p[0]].append((p[1], float(p[2]), int(p[3]), int(p[4]), int(p[5]), float(p[9])))

    tax = Taxonomy()
    rows = []
    for qid, seq in sorted(seqs.items()):
        sid, contig = qid.split("|", 1)
        row = {"sample_id": sid, "contig": contig, "contig_len": len(seq)}
        row.update(vote(by_q.get(qid, []), len(seq), rc, taxid_of, tax))
        rows.append(row)
    write_tsv(root / "tables" / "refined_contig_calls.tsv", rows, COLS)


if __name__ == "__main__":
    main()
