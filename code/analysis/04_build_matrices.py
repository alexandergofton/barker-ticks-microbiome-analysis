#!/usr/bin/env python3
"""Step 4 -- join LCA taxonomy to read mapping and build the abundance matrices.

This is the step the project was missing. Until now taxonomy (lca/) and
abundance (mapping/) were never joined, so the Pavian kreports counted CONTIGS
-- one contig = one "read" -- and a 300 bp contig at 2x coverage scored the same
as a 500 kb contig at 300x.

PER SAMPLE
  1. inner join lca/<s>.blast.lca  x  mapping/<s>.coverage.tsv  on contig id
  2. apply the contig evidence filter (config/analysis_params.yaml)
  3. renormalise TPM twice:
       tpm_total     -- over all retained contigs incl. host, sums to 1e6
       tpm_microbial -- over microbial contigs only, sums to 1e6  <- primary assay
  4. aggregate to taxon, apply the taxon detection thresholds
  5. aggregate to each export rank

WHY THE EVIDENCE FILTER EXISTS
  TPM divides by contig length, so short thinly-covered contigs are handed a
  disproportionate share of the budget. Measured on 0143-F1-Am-sp, contigs with
  an LCA call hold 95.4% of mapped READS but only 60.3% of TPM mass -- the
  remainder is short-contig assembly noise. Filtering before renormalising is
  what stops that noise from setting the denominator.

WHY READ COUNTS ARE STILL CARRIED
  TPM is the agreed abundance unit, but TPM alone cannot support a presence
  call: a 400 bp contig with 6 reads can outrank a real genome. mapped_reads,
  contig_len, n_contigs and bitscore travel with every row as DETECTION
  EVIDENCE and gate the presence call. They are not an alternative abundance.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_aliases, load_config, resolve_artefact, write_tsv  # noqa: E402

EXPORT_RANKS = ["species", "genus", "family", "phylum"]
BOWTIE_RE = re.compile(r"([\d.]+)% overall alignment rate")


ALIASES = load_aliases()


def read_alignment_rate(sample: str) -> float | None:
    p = resolve_artefact(sample, ROOT / "mapping", ".bowtie2.log", ALIASES)
    if p is None:
        return None
    m = BOWTIE_RE.search(p.read_text())
    return float(m.group(1)) if m else None


def load_taxonomy() -> pd.DataFrame:
    tax = pd.read_csv(ROOT / "results" / "taxonomy" / "taxon_lineage.tsv",
                      sep="\t", dtype={"taxid": int}, keep_default_na=False)
    return tax.set_index("taxid")


def rank_label(df: pd.DataFrame, rank: str, order: list[str]) -> pd.Series:
    """Name at `rank`, or the nearest resolved ancestor marked unresolved.

    An LCA that stopped at family cannot be pushed down to species. Binning it
    as '<family> (unresolved)' keeps it in the denominator so the column still
    sums to 1e6 and the unresolved share stays visible, instead of silently
    dropping it and inflating everything else.

    Vectorised: ~25% of LCA calls sit above species, so this runs over the whole
    frame on every rank and a per-row apply is far too slow.
    """
    out = df[rank].fillna("").astype(str)
    missing = out == ""
    for up in order[order.index(rank) + 1:]:
        if not missing.any():
            break
        anc = df[up].fillna("").astype(str)
        fill = missing & (anc != "")
        out = out.mask(fill, anc + " (unresolved)")
        missing = missing & ~fill
    return out.mask(missing, "unresolved")


def diagnose_skip(sample: str, cfg_filter: dict) -> str:
    """Explain why a sample produced no usable contigs."""
    lca_p = resolve_artefact(sample, ROOT / "lca", ".blast.lca", ALIASES)
    cov_p = resolve_artefact(sample, ROOT / "mapping", ".coverage.tsv", ALIASES)
    if lca_p is None:
        return "no LCA file found (checked canonical id and aliases)"
    if cov_p is None:
        return "no coverage file found (checked canonical id and aliases)"
    lca = pd.read_csv(lca_p, sep="\t", usecols=["query", "top_bitscore"])
    cov = pd.read_csv(cov_p, sep="\t", usecols=["contig_id", "contig_len", "mapped_reads"])
    j = cov.merge(lca, left_on="contig_id", right_on="query", how="inner")
    if len(j) == 0:
        return f"LCA and coverage share no contig ids ({len(lca)} lca, {len(cov)} cov)"
    return (f"no contig passes the evidence filter: {len(j)} joined contigs, "
            f"max_len={int(j.contig_len.max())} (need >={cfg_filter['min_contig_len']}), "
            f"max_reads={int(j.mapped_reads.max())} (need >={cfg_filter['min_mapped_reads']}), "
            f"total_mapped_reads={int(cov.mapped_reads.sum())} "
            f"-- library too shallow to analyse")


def process_sample(sample: str, tax: pd.DataFrame, cfg_filter: dict) -> tuple | None:
    # Files may be stored under an alias stem rather than the canonical id.
    lca_p = resolve_artefact(sample, ROOT / "lca", ".blast.lca", ALIASES)
    cov_p = resolve_artefact(sample, ROOT / "mapping", ".coverage.tsv", ALIASES)
    if lca_p is None or cov_p is None:
        return None

    lca = pd.read_csv(lca_p, sep="\t", usecols=["query", "lca_taxid", "n_hits", "top_bitscore"])
    cov = pd.read_csv(cov_p, sep="\t",
                      usecols=["contig_id", "contig_len", "mapped_reads", "tpm"])

    cov_tpm_all = cov["tpm"].sum()
    cov_reads_all = cov["mapped_reads"].sum()

    df = cov.merge(lca, left_on="contig_id", right_on="query", how="inner")
    n_lca, n_joined = len(lca), len(df)

    # ---- contig evidence filter -------------------------------------------
    keep = (
        (df["contig_len"] >= cfg_filter["min_contig_len"])
        & (df["mapped_reads"] >= cfg_filter["min_mapped_reads"])
        & (df["top_bitscore"] >= cfg_filter["min_top_bitscore"])
    )
    df["passed_filter"] = keep.astype(int)
    kept = df[keep].copy()
    if kept.empty:
        return None

    # ---- attach scope ------------------------------------------------------
    kept["lca_taxid"] = kept["lca_taxid"].astype(int)
    kept = kept.join(tax[["scope_class", "superkingdom_group"]], on="lca_taxid")
    kept["scope_class"] = kept["scope_class"].fillna("ambiguous")
    kept["superkingdom_group"] = kept["superkingdom_group"].fillna("unclassified")

    # ---- renormalise -------------------------------------------------------
    tot = kept["tpm"].sum()
    kept["tpm_total"] = kept["tpm"] / tot * 1e6 if tot > 0 else 0.0
    micro = kept["scope_class"] == "microbial"
    mtot = kept.loc[micro, "tpm"].sum()
    kept["tpm_microbial"] = 0.0
    if mtot > 0:
        kept.loc[micro, "tpm_microbial"] = kept.loc[micro, "tpm"] / mtot * 1e6

    # ---- aggregate to taxon ------------------------------------------------
    agg = (kept.groupby("lca_taxid")
              .agg(tpm_total=("tpm_total", "sum"),
                   tpm_microbial=("tpm_microbial", "sum"),
                   n_contigs=("contig_id", "size"),
                   sum_mapped_reads=("mapped_reads", "sum"),
                   max_contig_len=("contig_len", "max"),
                   max_bitscore=("top_bitscore", "max"))
              .reset_index())
    agg.insert(0, "sample_id", sample)

    # ---- fraction summary --------------------------------------------------
    by_scope = kept.groupby("scope_class")[["tpm", "mapped_reads"]].sum()
    by_grp = kept.groupby("superkingdom_group")[["tpm", "mapped_reads"]].sum()
    frac = {
        "sample_id": sample,
        "contigs_total": len(cov),
        "contigs_with_lca": n_lca,
        "contigs_joined": n_joined,
        "contigs_passing_filter": int(keep.sum()),
        "pct_contigs_passing": round(100 * keep.sum() / max(n_joined, 1), 2),
        "reads_mapped_total": int(cov_reads_all),
        "alignment_rate_pct": read_alignment_rate(sample),
        # TPM shares are of the WHOLE coverage table, so unassigned is visible
        "tpm_pct_no_lca": round(100 * (cov_tpm_all - df["tpm"].sum()) / cov_tpm_all, 3),
        "tpm_pct_filtered_out": round(100 * df.loc[~keep, "tpm"].sum() / cov_tpm_all, 3),
        "read_pct_no_lca": round(100 * (cov_reads_all - df["mapped_reads"].sum()) / max(cov_reads_all, 1), 3),
    }
    for cls in ["host", "bloodmeal", "microbial", "other_eukaryote", "ambiguous"]:
        frac[f"tpm_pct_{cls}"] = round(100 * by_scope["tpm"].get(cls, 0) / cov_tpm_all, 3)
        frac[f"read_pct_{cls}"] = round(100 * by_scope["mapped_reads"].get(cls, 0) / max(cov_reads_all, 1), 3)
    for g in ["Bacteria", "Archaea", "Viruses", "Eukaryota"]:
        frac[f"tpm_pct_{g.lower()}"] = round(100 * by_grp["tpm"].get(g, 0) / cov_tpm_all, 3)

    kept_out = kept[["contig_id", "contig_len", "mapped_reads", "tpm",
                     "tpm_total", "tpm_microbial", "lca_taxid", "scope_class",
                     "n_hits", "top_bitscore"]].copy()
    kept_out.insert(0, "sample_id", sample)
    return agg, frac, kept_out


def main() -> None:
    cfg = load_config()
    cfg_filter = cfg["contig_filter"]
    tax_cfg = cfg["taxon_filter"]
    tax = load_taxonomy()

    audit = pd.read_csv(ROOT / "results" / "audit" / "sample_completeness.tsv", sep="\t")
    samples = sorted(audit.loc[audit["status"] == "complete", "sample_id"])
    print(f"{len(samples)} complete samples", file=sys.stderr)

    tbl = ROOT / "results" / "tables"
    mat = ROOT / "results" / "matrices"
    tbl.mkdir(parents=True, exist_ok=True)
    mat.mkdir(parents=True, exist_ok=True)

    aggs, fracs, excluded = [], [], []
    contig_out = gzip.open(tbl / "contig_taxon_long.tsv.gz", "wt")
    header_written = False

    for i, s in enumerate(samples, 1):
        res = process_sample(s, tax, cfg_filter)
        if res is None:
            # Record rather than merely log: a sample vanishing from the matrix
            # must be traceable to a reason, not inferred from a row count.
            reason = diagnose_skip(s, cfg_filter)
            excluded.append({"sample_id": s, "stage": "matrix_build", "reason": reason})
            print(f"  [{i}/{len(samples)}] {s}: EXCLUDED -- {reason}", file=sys.stderr)
            continue
        agg, frac, contigs = res
        aggs.append(agg)
        fracs.append(frac)
        contigs.to_csv(contig_out, sep="\t", index=False, header=not header_written)
        header_written = True
        if i % 20 == 0 or i == len(samples):
            print(f"  [{i}/{len(samples)}] {s}", file=sys.stderr)
    contig_out.close()

    taxon = pd.concat(aggs, ignore_index=True)

    # ---- taxon-level detection threshold -----------------------------------
    detected = (
        (taxon["n_contigs"] >= tax_cfg["min_contigs_per_taxon"])
        & (taxon["sum_mapped_reads"] >= tax_cfg["min_reads_per_taxon"])
    )
    taxon["detected"] = detected.astype(int)
    print(f"\ntaxon x sample observations: {len(taxon):,} "
          f"({detected.sum():,} pass detection thresholds)", file=sys.stderr)

    taxon = taxon.join(tax, on="lca_taxid", rsuffix="_tax")
    taxon.to_csv(tbl / "taxon_abundance_long.tsv.gz", sep="\t", index=False,
                 compression="gzip")
    print(f"  wrote results/tables/taxon_abundance_long.tsv.gz ({len(taxon):,} rows)",
          file=sys.stderr)

    write_tsv(tbl / "fraction_summary.tsv", fracs, list(fracs[0].keys()))
    write_tsv(ROOT / "results" / "audit" / "excluded_samples.tsv",
              excluded or [{"sample_id": "", "stage": "", "reason": "none"}],
              ["sample_id", "stage", "reason"])

    # ---- wide matrices per rank -------------------------------------------
    det = taxon[taxon["detected"] == 1].copy()
    micro = det[det["is_microbial"] == 1].copy()
    order = EXPORT_RANKS + ["domain"]

    for rank in EXPORT_RANKS:
        # microbial TPM, renormalised again after detection filtering
        m = micro.assign(label=rank_label(micro, rank, order))
        w = m.pivot_table(index="label", columns="sample_id",
                          values="tpm_microbial", aggfunc="sum", fill_value=0.0)
        w = w.div(w.sum(axis=0).replace(0, pd.NA), axis=1).fillna(0.0) * 1e6
        w.to_csv(mat / f"tpm_microbial_{rank}.tsv", sep="\t",
                 float_format="%.6f")
        print(f"  wrote results/matrices/tpm_microbial_{rank}.tsv  "
              f"({w.shape[0]} taxa x {w.shape[1]} samples)", file=sys.stderr)

        # matching read-evidence matrix (NOT an abundance -- gating/QC only)
        e = m.pivot_table(index="label", columns="sample_id",
                          values="sum_mapped_reads", aggfunc="sum", fill_value=0)
        e = e.reindex(index=w.index, columns=w.columns, fill_value=0).astype(int)
        e.to_csv(mat / f"evidence_reads_{rank}.tsv", sep="\t")

        n = m.pivot_table(index="label", columns="sample_id",
                          values="n_contigs", aggfunc="sum", fill_value=0)
        n.reindex(index=w.index, columns=w.columns, fill_value=0).astype(int) \
         .to_csv(mat / f"evidence_ncontigs_{rank}.tsv", sep="\t")

    print("\ndone", file=sys.stderr)


if __name__ == "__main__":
    main()
