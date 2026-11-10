#!/usr/bin/env python3
"""Step 18 -- cut the marker locus out of every downloaded reference.

Step 17 downloads whatever the panel queries returned: gene records, whole
genomes, and an outgroup. None of it is trustworthy as-is -- a record matched by
title may not carry the locus at all, and a 1.3 Mb genome is 99.9% not the gene.
This step blasts each record and keeps only the part that is the marker:

  1. The queries are the step-11 search probes for this locus plus the longest
     few sample sequences (`curate.max_sample_queries`). blastn for rRNA loci,
     tblastn for protein-coding ones -- the same programs step 13 used against
     the contigs, against the same panel floors.
  2. HSPs on one record and one strand are merged into a locus, and the
     best-scoring locus per record is kept. A genome holds one groEL and twenty
     IS1111 copies; in both cases the best hit is the answer, because a second
     copy of the same locus is not extra phylogenetic information.
  3. The locus is cut out, reverse-complemented if it sat on the minus strand,
     and -- for a protein-coding marker -- trimmed to the forward frame with the
     fewest internal stops and to a whole number of codons, exactly as step 14
     treats a sample sequence. That is what lets references and samples share
     one codon alignment in step 19.
  4. A reference whose locus still carries an internal stop codon is dropped
     (`curate.drop_internal_stops`): it is a pseudogene or a mis-annotation, and
     it would distort the tree it was added to.
  5. What survives is capped per organism and filled round-robin across
     organisms, so the reference set is broad rather than forty strains of the
     best-studied species.

A record that fails every test is simply not a reference for this locus. That is
a normal outcome of a loose query and is reported, never silently dropped:
`tree_reference_sequences.tsv` lists every reference that made it, with the
identity and coverage that got it there.

Outputs
  results/markers/trees/refs/sets/<target>__<marker>.refs.fna   (+ .faa)
  results/markers/trees/tree_reference_sequences.tsv
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, write_tsv  # noqa: E402
from markers_lib import (best_frame, internal_stops, read_fasta, require, revcomp,  # noqa: E402
                         run, slug, translate, write_fasta)
from trees_lib import (BLAST_FMT, OUTGROUP_TAG, REF_TAG, best_locus_per_subject,  # noqa: E402
                       cfg_trees, organism_from_title, parse_blast, ref_name,
                       sample_fasta, loci_with_sequences, trees_root)

SEQ_COLS = ["target", "marker", "marker_type", "seq_id", "accession", "organism",
            "tier", "role", "record_len", "locus_start", "locus_end", "strand",
            "seq_len_nt", "aa_len", "pident", "ref_cov_pct", "max_bitscore",
            "best_query", "defline"]


def extract(seq_full: str, locus: dict, is_prot: bool, flank: int) -> tuple[str, str, int, int]:
    """Cut the locus out of a record; returns (nt, aa, start, end), nt in frame."""
    start = max(1, locus["start"] - flank)
    end = min(len(seq_full), locus["end"] + flank)
    sub = seq_full[start - 1:end]
    if locus["strand"] == "-":
        sub = revcomp(sub)
    aa = ""
    if is_prot:
        frame, _ = best_frame(sub)
        sub = sub[frame:]
        sub = sub[: len(sub) - len(sub) % 3]
        aa = translate(sub, 0)
    return sub.upper(), aa, start, end


def spread(cands: list[dict], max_per_org: int, want: int) -> list[dict]:
    """Round-robin across organisms, best candidate first within each."""
    by_org: dict[str, list[dict]] = defaultdict(list)
    for c in sorted(cands, key=lambda c: (-c["max_bitscore"], -c["seq_len_nt"])):
        by_org[c["organism"]].append(c)
    order = sorted(by_org, key=lambda o: -by_org[o][0]["max_bitscore"])
    picked: list[dict] = []
    for rnd in range(max_per_org):
        for org in order:
            if len(picked) >= want:
                return picked
            if rnd < len(by_org[org]):
                picked.append(by_org[org][rnd])
    return picked


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", action="append", help="limit to this target (repeatable)")
    ap.add_argument("--marker", action="append", help="limit to this marker (repeatable)")
    ap.add_argument("--force", action="store_true", help="re-run BLAST for cached hits")
    args = ap.parse_args()

    tc = cfg_trees()
    cc = tc["curate"]
    root = trees_root()
    refs, hits_dir = root / "refs", root / "hits"
    sets_dir, db_dir, q_dir = refs / "sets", refs / "blastdb", refs / "queries"
    for d in (hits_dir, sets_dir, db_dir, q_dir):
        d.mkdir(parents=True, exist_ok=True)

    makeblastdb, blastn, tblastn = require("makeblastdb"), require("blastn"), require("tblastn")

    loci = loci_with_sequences(args.target, args.marker)
    if not loci:
        sys.exit("ERROR: no target x marker has sample sequences -- run run_markers.sh first.")

    record_cache: dict[Path, dict[str, str]] = {}

    def records(path: Path) -> dict[str, str]:
        if path not in record_cache:
            record_cache[path] = read_fasta(path) if path.exists() else {}
        return record_cache[path]

    rows: list[dict] = []
    summary: list[tuple[str, str, int, int, int]] = []

    for i, row in enumerate(loci, 1):
        target, marker = row["target"], row["marker"]
        key = f"{slug(target)}__{slug(marker)}"
        is_prot = row["marker_type"] == "protein"
        ext = "faa" if is_prot else "fna"
        unit = 3 if is_prot else 1
        tag = f"[{i:>3}/{len(loci)}] {target} {marker}"

        # ---- queries: the step-11 probes plus the longest sample sequences ---
        queries: dict[str, str] = {}
        probe = ROOT / "results" / "markers" / "refs" / f"{key}.{ext}"
        queries.update({f"probe_{n.split()[0]}": s for n, s in records(probe).items()})
        samples = records(sample_fasta(target, marker, ext))
        for n, s in sorted(samples.items(), key=lambda kv: -len(kv[1]))[:int(cc["max_sample_queries"])]:
            queries[n] = s
        if not queries:
            print(f"{tag}: no probes and no sample sequences -- skipped", file=sys.stderr)
            continue
        qfile = q_dir / f"{key}.{ext}"
        write_fasta(qfile, queries)

        # ---- one pass per reference tier -------------------------------------
        sources = [("gene", refs / "gene" / f"{key}.fna", REF_TAG),
                   ("genome", refs / "genome" / f"{slug(target)}.fna", REF_TAG),
                   ("outgroup", refs / "outgroup" / f"{key}.fna", OUTGROUP_TAG)]
        cands: list[dict] = []
        n_dropped_stop = n_dropped_len = 0

        for tier, fa, tag_name in sources:
            if not fa.exists() or fa.stat().st_size == 0:
                continue
            db = db_dir / (f"{key}.{tier}" if tier != "genome" else f"{slug(target)}.genome")
            if args.force or not Path(f"{db}.nin").exists():
                run([makeblastdb, "-in", str(fa), "-dbtype", "nucl", "-out", str(db)])
            hit_file = hits_dir / f"{key}.{tier}.blast.tsv"
            if args.force or not hit_file.exists():
                prog = tblastn if is_prot else blastn
                cmd = [prog, "-query", str(qfile), "-db", str(db), "-out", str(hit_file),
                       "-outfmt", BLAST_FMT, "-num_threads", str(cc["threads"]),
                       "-evalue", str(cc["evalue_protein"] if is_prot else cc["evalue_rrna"]),
                       "-max_target_seqs", "10000"]
                if not is_prot:
                    cmd += ["-task", "blastn"]   # not megablast: cross-genus rRNA
                run(cmd)

            min_pid = (float(cc["outgroup_min_pident"]) if tier == "outgroup"
                       else float(row["min_pident"]))
            loci_by_rec = best_locus_per_subject(
                parse_blast(hit_file, int(row["min_aln_len"]), min_pid),
                int(cc["merge_gap_protein"] if is_prot else cc["merge_gap_rrna"]), unit)

            seqs = records(fa)
            by_id = {n.split()[0]: n for n in seqs}
            for rec_id, locus in loci_by_rec.items():
                defline = by_id.get(rec_id)
                if defline is None:
                    continue
                nt, aa, start, end = extract(
                    seqs[defline], locus, is_prot,
                    0 if is_prot else int(cc["flank_rrna"]))
                if not (int(cc["min_ref_len_nt"]) <= len(nt) <= int(cc["max_ref_len_nt"])):
                    n_dropped_len += 1
                    continue
                if is_prot and cc["drop_internal_stops"] and internal_stops(aa):
                    n_dropped_stop += 1
                    continue
                organism = organism_from_title(defline)
                cands.append({
                    "target": target, "marker": marker, "marker_type": row["marker_type"],
                    "accession": rec_id, "organism": organism, "tier": tier,
                    "role": "outgroup" if tag_name == OUTGROUP_TAG else "reference",
                    "tag": tag_name, "record_len": locus["subject_len"],
                    "locus_start": start, "locus_end": end, "strand": locus["strand"],
                    "seq_len_nt": len(nt), "aa_len": len(aa) if is_prot else "",
                    "pident": locus["best_pident"], "ref_cov_pct": locus["ref_cov_pct"],
                    "max_bitscore": locus["max_bitscore"], "best_query": locus["best_query"],
                    "defline": defline, "nt": nt, "aa": aa,
                })

        # ---- dedupe, cap, name ----------------------------------------------
        seen: set[str] = set()
        uniq: list[dict] = []
        for c in sorted(cands, key=lambda c: (-c["max_bitscore"], -c["seq_len_nt"])):
            if c["nt"] in seen:
                continue
            seen.add(c["nt"])
            uniq.append(c)

        outgroups = [c for c in uniq if c["role"] == "outgroup"][:int(cc["max_outgroups"])]
        keep = spread([c for c in uniq if c["role"] == "reference"],
                      int(cc["max_per_organism"]), int(cc["max_refs_per_locus"])) + outgroups

        nt_out: dict[str, str] = {}
        aa_out: dict[str, str] = {}
        for c in keep:
            name = ref_name(c["tag"], c["accession"], c["organism"], target, marker,
                            c["seq_len_nt"])
            # Two records of the same organism and identical trimmed length would
            # collide; the accession is already in the name, so a collision means
            # the same record twice.
            if name in nt_out:
                continue
            nt_out[name] = c["nt"]
            if c["aa"]:
                aa_out[name] = c["aa"]
            rows.append({**{k: v for k, v in c.items() if k in SEQ_COLS}, "seq_id": name})

        if nt_out:
            write_fasta(sets_dir / f"{key}.refs.fna", nt_out)
            if aa_out:
                write_fasta(sets_dir / f"{key}.refs.faa", aa_out)
        n_org = len({c["organism"] for c in keep if c["role"] == "reference"})
        summary.append((target, marker, len(nt_out), n_org, len(outgroups)))
        print(f"{tag}: {len(cands):>4} candidate -> {len(nt_out):>3} refs "
              f"({n_org} organisms, {len(outgroups)} outgroup)"
              + (f", dropped {n_dropped_stop} internal-stop" if n_dropped_stop else "")
              + (f", {n_dropped_len} out-of-range" if n_dropped_len else ""),
              file=sys.stderr)

    seq_path = root / "tree_reference_sequences.tsv"
    if (args.target or args.marker) and seq_path.exists():
        import csv
        done = {(r["target"], r["marker"]) for r in rows}
        with open(seq_path) as fh:
            others = [r for r in csv.DictReader(fh, delimiter="\t")
                      if (r["target"], r["marker"]) not in done]
        print(f"  merging with {len(others):,} references from loci not in this run",
              file=sys.stderr)
        rows = sorted(others + rows, key=lambda r: (r["target"], r["marker"], r["seq_id"]))
    write_tsv(seq_path, rows, SEQ_COLS)

    bare = [s for s in summary if s[2] == 0]
    if bare:
        print(f"\n{len(bare)} loci ended with NO references: "
              + ", ".join(f"{t}/{m}" for t, m, *_ in bare[:12])
              + (" ..." if len(bare) > 12 else "")
              + "\n  Either the query returned nothing usable or nothing in it "
                "carries the locus at the panel's min_pident/min_aln_len. Those "
                "trees will hold sample sequences only.", file=sys.stderr)
    print(f"\n{sum(s[2] for s in summary):,} references over {len(summary)} loci "
          f"-> {(sets_dir).relative_to(ROOT)}/", file=sys.stderr)


if __name__ == "__main__":
    main()
