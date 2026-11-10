#!/usr/bin/env python3
"""Step 26c -- is a genus-discordant mitogenome a contaminant, a mixed tube, or
the wrong tick altogether?

Step 26b flags samples whose mitogenome is a tick of another genus than the
curated one. The mitogenome alone cannot say why. This step asks the sample's
own whole-assembly nt BLAST (blast/<sample>.blast.txt) a second, independent
question: which tick genus do its NUCLEAR rRNA contigs (18S/28S/5.8S/ITS)
belong to?

Why nuclear rRNA and not "percentage of tick hits per genus": nt holds whole
nuclear genomes for only a few tick genera (Dermacentor, Rhipicephalus,
Ixodes, Haemaphysalis, Hyalomma). On this cohort ~70% of all tick HSPs in a
correctly labelled Amblyomma or Bothriocroton sample go to Dermacentor, and
1% to the curated genus -- the raw percentage measures nt's composition, not the
sample's. Restricting to high-identity hits (`min_pident`, default 95%) and
splitting by subject type removes that bias: tick nuclear rRNA and
mitochondrial genes are in nt for nearly every genus, and at >=95% identity a
hit names the genus actually present.

Per sample, three signals are compared, each as a share of
contig length x k-mer coverage (a proxy for read share) among contigs whose
best high-identity tick hit is that genus:

  rrna_*   nuclear rRNA contigs -- which tick's nuclear DNA is in the tube
  mito_*   mitochondrial-subject hits on contigs of any length (incl. the
           fragments < 3 kb that step 21 never sees)
  nearest sample mitogenome (p-distance, step 26): a foreign mitogenome
           0.00% from another sample on the SAME sequencing lane points at that
           sample as the source

Calls (thresholds: `mito.audit` in config/analysis_params.yaml):

  CONTAMINANT_MITOGENOME  nuclear rRNA is the curated genus (>= rrna_major);
                          the mitogenome is a foreign tick's. The tick ID
                          stands; the selected mitogenome does not represent it.
  MIXED_SAMPLE            both genera carry real nuclear rRNA share. Which tick
                          is "the sample" cannot be decided from sequence.
  WRONG_TICK              curated genus absent (< rrna_minor) from both nuclear
                          rRNA and mitochondrial hits: tube swap, mislabel or
                          misidentification.
  UNRESOLVED              no nuclear rRNA tick hit (or no BLAST file) to decide.
  ID_SUPPORTED[_NO_MITOGENOME]  nuclear rRNA is the curated genus and so is the
                          selected mitogenome (or none was selected).
  UNTESTABLE_GENUS        the curated genus has no nuclear rRNA in nt at all.

CONTAMINANT_MITOGENOME, WRONG_TICK, MIXED_SAMPLE and UNRESOLVED mean
`exclude_from_mito_tree = yes` (when the sample has a mitogenome); the last three
also mean `exclude_from_species_analyses = yes` (the curated
species label cannot be trusted for that tube). CONTAMINANT_MITOGENOME samples
carry another sample's DNA, so their microbiome profile may carry that
sample's microbes too (`microbiome_caution`).

Outputs
  results/mito/audit/tick_hits/<sample>.tsv   cache: best high-identity tick hit
                                              per (contig, subject type)
  results/mito/tables/tick_genus_signal.tsv   sample x type x genus shares, every sample
  results/mito/tables/mito_discordance_audit.tsv   one row per audited sample
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_aliases, resolve_artefact, write_tsv  # noqa: E402
from markers_lib import cfg_mito, mito_root  # noqa: E402

TICK_TAXID = "6935"  # Ixodida
TYPES = ("mito", "nuc_rRNA", "nuclear")

SIGNAL_COLS = ["sample_id", "curated_genus", "subject_type", "genus", "n_contigs",
               "weight", "share"]
AUDIT_COLS = ["sample_id", "tick_species", "curated_genus", "life_stage", "n_ticks_pooled",
              "host_species", "seq_run", "verdict_26b", "mito_status", "mito_genus",
              "mito_contig", "mito_cov", "n_genes", "rrna_top", "rrna_curated_share",
              "rrna_n_contigs", "mito_hits_top", "mito_hits_curated_share",
              "best_curated_mito_candidate", "nearest_sample", "nearest_sample_species",
              "nearest_sample_pdist", "nearest_same_run", "call",
              "exclude_from_mito_tree", "exclude_from_species_analyses",
              "microbiome_caution", "reason"]


def subject_type(title: str) -> str:
    t = title.lower()
    if "mitochond" in t or re.search(r"cytochrome|\bcox[123]\b|\bcoi\b|12s|16s ribosomal", t):
        return "mito"
    if re.search(r"18s|28s|5\.8s|internal transcribed|\bits[12]?\b|ribosomal rna|rrna", t):
        return "nuc_rRNA"
    return "nuclear"


def tick_genus_map() -> dict[str, str]:
    """Every taxid under Ixodida -> its genus name (merged taxids included)."""
    parent_kids: dict[str, list[str]] = defaultdict(list)
    rank = {}
    with open(ROOT / "taxdump" / "nodes.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t")
            parent_kids[p[1]].append(p[0]); rank[p[0]] = p[2]
    names = {}
    with open(ROOT / "taxdump" / "names.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t")
            if p[3].startswith("scientific name"):
                names[p[0]] = p[1]
    out, stack = {}, [(TICK_TAXID, "")]
    while stack:
        t, g = stack.pop()
        if rank.get(t) == "genus":
            g = names.get(t, "")
        out[t] = g or "Ixodida_unassigned"
        stack.extend((k, g) for k in parent_kids.get(t, []) if k != t)
    with open(ROOT / "taxdump" / "merged.dmp") as fh:
        for line in fh:
            p = line.split("\t|\t")
            new = p[1].split("\t")[0]
            if new in out:
                out[p[0]] = out[new]
    return out


def extract(sid: str, blast: Path, cache: Path, genus_tsv: Path, min_pid: float,
            min_len: int) -> None:
    """Stream the raw BLAST once with awk (files are up to ~60 GB), keep tick
    HSPs above the identity/length floor, reduce to the best hit per
    (contig, subject type)."""
    awk = ("NR==FNR{g[$1]=$2; next} $8>=P && $9>=L {n=split($7,t,\";\");"
           " for(i=1;i<=n;i++) if(t[i] in g){print $1\"\\t\"g[t[i]]\"\\t\"$17\"\\t\"$3; break}}")
    proc = subprocess.Popen(["awk", "-F", "\t", "-v", f"P={min_pid}", "-v", f"L={min_len}",
                             awk, str(genus_tsv), str(blast)],
                            stdout=subprocess.PIPE, text=True, env={"LC_ALL": "C"})
    best: dict[tuple[str, str], tuple[str, float]] = {}
    for line in proc.stdout:
        contig, genus, bits, title = line.rstrip("\n").split("\t", 3)
        key, b = (contig, subject_type(title)), float(bits)
        if key not in best or b > best[key][1]:
            best[key] = (genus, b)
    if proc.wait() != 0:
        raise RuntimeError(f"awk failed on {blast}")
    tmp = cache.with_suffix(".tmp")
    with open(tmp, "w") as fh:
        fh.write("contig\tsubject_type\tgenus\tbitscore\n")
        for (c, k), (g, b) in sorted(best.items()):
            fh.write(f"{c}\t{k}\t{g}\t{b}\n")
    tmp.rename(cache)


def contig_weight(contig: str) -> float:
    m = re.search(r"length_(\d+)_cov_([0-9.]+)", contig)
    return int(m.group(1)) * float(m.group(2)) if m else 0.0


def contig_cov(contig: str) -> float:
    m = re.search(r"_cov_([0-9.]+)$", contig or "")
    return float(m.group(1)) if m else 0.0


def genus_signal(cache: Path) -> dict[str, tuple[Counter, Counter]]:
    """subject type -> (genus -> weight, genus -> n contigs)."""
    out = {k: (Counter(), Counter()) for k in TYPES}
    with open(cache) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            w, n = out[r["subject_type"]]
            w[r["genus"]] += contig_weight(r["contig"]); n[r["genus"]] += 1
    return out


def share(c: Counter, g: str) -> float | None:
    t = sum(c.values())
    return c.get(g, 0) / t if t else None


def top(c: Counter, k: int = 3) -> str:
    t = sum(c.values())
    return ", ".join(f"{g} {100 * v / t:.0f}%" for g, v in c.most_common(k)) if t else ""


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def nearest_samples(pdist: Path, wanted: set[str]) -> dict[str, tuple[str, float]]:
    out = {}
    with open(pdist) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = [h.split("|")[0] for h in next(rd)[1:]]
        for row in rd:
            s = row[0].split("|")[0]
            if s not in wanted:
                continue
            d = [(float(x), h) for x, h in zip(row[1:], hdr)
                 if h != s and not h.startswith("REF_") and x not in ("", "NA")]
            if d:
                v, h = min(d)
                out[s] = (h, v)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append", help="limit extraction to this sample (repeatable)")
    ap.add_argument("--jobs", type=int, default=4, help="parallel BLAST-file scans (I/O bound)")
    ap.add_argument("--force", action="store_true", help="re-scan raw BLAST files")
    args = ap.parse_args()

    ac = cfg_mito().get("audit", {}) or {}
    min_pid, min_len = float(ac.get("min_pident", 95)), int(ac.get("min_aln_len", 150))
    rrna_major, rrna_minor = float(ac.get("rrna_major", 0.8)), float(ac.get("rrna_minor", 0.05))

    root = mito_root()
    cache_dir = root / "audit" / "tick_hits"
    cache_dir.mkdir(parents=True, exist_ok=True)
    aliases = load_aliases()

    meta = {r["sample_id"]: r for r in read_tsv(ROOT / "metadata" / "sample_metadata.tsv")}
    sel = {r["sample_id"]: r for r in read_tsv(root / "mito_selection.tsv")}
    vm = {r["sample_id"]: r for r in read_tsv(root / "tables" / "mito_vs_metadata.tsv")}
    lca_genus = {r["lca_name"]: (r["genus"] or r["lca_name"])
                 for r in read_tsv(ROOT / "results" / "taxonomy" / "taxon_lineage.tsv")
                 if r["order"] == "Ixodida"}

    samples = sorted(sel)
    if args.sample:
        samples = [s for s in samples if s in set(args.sample)]

    # ---- 1. cached per-sample scan of the raw BLAST ------------------------
    todo = []
    for s in samples:
        cache = cache_dir / f"{s}.tsv"
        blast = resolve_artefact(s, ROOT / "blast", ".blast.txt", aliases)
        if blast is not None and (args.force or not cache.exists()):
            todo.append((s, blast, cache))
    if todo:
        genus_tsv = cache_dir / "_ixodida_taxid_genus.tsv"
        with open(genus_tsv, "w") as fh:
            fh.writelines(f"{t}\t{g}\n" for t, g in tick_genus_map().items())
        print(f"scanning {len(todo)} raw BLAST files ({args.jobs} at a time)", file=sys.stderr)
        with ThreadPoolExecutor(args.jobs) as ex:
            futs = {ex.submit(extract, s, b, c, genus_tsv, min_pid, min_len): s for s, b, c in todo}
            for i, f in enumerate(futs, 1):
                f.result()
                print(f"  [{i}/{len(todo)}] {futs[f]}", file=sys.stderr)

    # ---- 2. genus signal, every sample -------------------------------------
    signal, sig_rows = {}, []
    for s in samples:
        cache = cache_dir / f"{s}.tsv"
        if not cache.exists():
            continue
        signal[s] = genus_signal(cache)
        cg = (meta.get(s, {}).get("tick_genus") or "")
        for k in TYPES:
            w, n = signal[s][k]
            tot = sum(w.values())
            for g, v in w.most_common():
                sig_rows.append({"sample_id": s, "curated_genus": cg, "subject_type": k,
                                 "genus": g, "n_contigs": n[g], "weight": round(v),
                                 "share": round(v / tot, 4) if tot else ""})
    write_tsv(root / "tables" / "tick_genus_signal.tsv", sig_rows, SIGNAL_COLS)

    # ---- 3. audit: 26b-discordant samples + any sample whose nuclear rRNA
    #         genus disagrees with the curated genus ------------------------
    def mito_genus(s: str) -> str:
        r = sel.get(s, {})
        if r.get("status") != "ok":
            return ""
        if r.get("selected_contig_lca_class") == "tick":
            return lca_genus.get(r["selected_contig_lca"], r["selected_contig_lca"])
        ref = vm.get(s, {}).get("nearest_ref_species", "")
        return f"{ref.split(' ')[0]} (nearest ref)" if ref else ""

    audit_set = []
    for s in samples:
        cg = meta.get(s, {}).get("tick_genus") or ""
        v = vm.get(s, {}).get("verdict", "")
        rr = signal.get(s, {}).get("nuc_rRNA", (Counter(), Counter()))[0]
        rr_disagree = cg and sum(rr.values()) and (share(rr, cg) or 0) < rrna_major
        if (v in ("GENUS_DISCORDANT", "NOT_TICK_CONTIG") or rr_disagree
                or sel[s]["status"] == "no_tick_contig"):
            audit_set.append(s)
    # A genus with no >=95% nuclear-rRNA hit in ANY sample is not in nt (e.g.
    # Cryptocroton): its ticks' rRNA lands on the nearest represented genus,
    # so absence of the curated genus there proves nothing.
    rrna_genera = {g for sg in signal.values() for g in sg["nuc_rRNA"][0]}

    pd_path = root / "tables" / "mito_distance_pdistance.tsv"
    nearest = nearest_samples(pd_path, set(audit_set)) if pd_path.exists() else {}

    loci_by_contig: dict[tuple[str, str], set] = defaultdict(set)
    with open(root / "mito_loci.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] == "kept" and r["contig_lca_class"] == "tick":
                loci_by_contig[(r["sample_id"], r["contig_id"] + "\t" + r["contig_lca"])].add(r["gene"])

    rows = []
    for s in audit_set:
        m, r = meta.get(s, {}), sel.get(s, {})
        cg = m.get("tick_genus") or ""
        sg = signal.get(s)
        rr, rn = sg["nuc_rRNA"] if sg else (Counter(), Counter())
        mh = sg["mito"][0] if sg else Counter()
        rs, ms = share(rr, cg), share(mh, cg)
        mg = mito_genus(s)

        # best tick-LCA candidate contig of the curated genus (step 21 loci)
        cands = []
        for (sid, key), genes in loci_by_contig.items():
            if sid != s:
                continue
            contig, lca = key.split("\t")
            if cg and lca_genus.get(lca) == cg:
                cands.append((len(genes), contig_cov(contig), contig))
        best_cur = max(cands) if cands else None

        nn = nearest.get(s)
        nn_run = meta.get(nn[0], {}).get("seq_run", "") if nn else ""
        same_run = bool(nn and m.get("seq_run") and nn_run == m.get("seq_run"))

        has_mito = r.get("status") == "ok"
        mito_is_curated = bool(cg) and mg.split(" ")[0] == cg
        if not cg:
            call, reason = "NO_CURATED_ID", f"no curated genus; nuclear rRNA {top(rr)}"
        elif cg not in rrna_genera:
            call, reason = "UNTESTABLE_GENUS", f"no {cg} nuclear rRNA in nt; nearest hits {top(rr)}"
        elif rs is None:
            call, reason = "UNRESOLVED", "no >=95% tick nuclear-rRNA hit (or no BLAST file)"
        elif rs >= rrna_major:
            call = ("ID_SUPPORTED_NO_MITOGENOME" if not has_mito else
                    "ID_SUPPORTED" if mito_is_curated else "CONTAMINANT_MITOGENOME")
            reason = f"nuclear rRNA {100 * rs:.0f}% {cg}"
        elif rs < rrna_minor and (ms or 0) < rrna_minor:
            call, reason = "WRONG_TICK", f"nuclear rRNA {100 * rs:.0f}% {cg}; {top(rr, 1)}"
        else:
            call = "MIXED_SAMPLE"
            reason = f"nuclear rRNA {top(rr)}; mito hits {top(mh)}"
        if has_mito and nn and nn[1] == 0 and same_run and not mito_is_curated:
            reason += f"; mitogenome identical to {nn[0]} on the same lane"

        excl_tree = has_mito and call in ("CONTAMINANT_MITOGENOME", "WRONG_TICK",
                                          "MIXED_SAMPLE", "UNRESOLVED")
        excl_species = call in ("WRONG_TICK", "MIXED_SAMPLE", "UNRESOLVED")
        rows.append({
            "sample_id": s, "tick_species": m.get("tick_species", ""), "curated_genus": cg,
            "life_stage": m.get("life_stage", ""), "n_ticks_pooled": m.get("n_ticks_pooled", ""),
            "host_species": m.get("host_species", ""), "seq_run": m.get("seq_run", ""),
            "verdict_26b": vm.get(s, {}).get("verdict", ""), "mito_status": r.get("status", ""),
            "mito_genus": mg, "mito_contig": r.get("selected_contig", ""),
            "mito_cov": contig_cov(r.get("selected_contig", "")) or "",
            "n_genes": r.get("n_genes_found", ""),
            "rrna_top": top(rr), "rrna_curated_share": "" if rs is None else round(rs, 3),
            "rrna_n_contigs": sum(rn.values()),
            "mito_hits_top": top(mh), "mito_hits_curated_share": "" if ms is None else round(ms, 3),
            "best_curated_mito_candidate": (f"{best_cur[2]} ({best_cur[0]} genes)" if best_cur else ""),
            "nearest_sample": nn[0] if nn else "",
            "nearest_sample_species": meta.get(nn[0], {}).get("tick_species", "") if nn else "",
            "nearest_sample_pdist": nn[1] if nn else "", "nearest_same_run": int(same_run),
            "call": call, "exclude_from_mito_tree": "yes" if excl_tree else "no",
            "exclude_from_species_analyses": "yes" if excl_species else "no",
            "microbiome_caution": "yes" if call in ("CONTAMINANT_MITOGENOME", "MIXED_SAMPLE") else "no",
            "reason": reason,
        })
    order = {"WRONG_TICK": 0, "MIXED_SAMPLE": 1, "UNRESOLVED": 2, "CONTAMINANT_MITOGENOME": 3,
             "ID_SUPPORTED_NO_MITOGENOME": 4, "ID_SUPPORTED": 5, "UNTESTABLE_GENUS": 6,
             "NO_CURATED_ID": 7}
    rows.sort(key=lambda x: (order[x["call"]], x["curated_genus"], x["sample_id"]))
    write_tsv(root / "tables" / "mito_discordance_audit.tsv", rows, AUDIT_COLS)
    for c, n in Counter(r["call"] for r in rows).most_common():
        print(f"  {c:<24} {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
