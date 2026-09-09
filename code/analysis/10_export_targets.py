#!/usr/bin/env python3
"""Step 10 -- export pathogen / endosymbiont target tables.

This is the handoff into the targeted follow-up analysis. extracting_pathogens/
already holds per-genus contig pulls for Anaplasma, Coxiella, Ehrlichia,
Francisella, Midichloria, Neoehrlichia, Rickettsia and Spiroplasma; this step
puts those genera on the same quantitative footing as everything else and adds
the piroplasms, Bartonella and the tick-borne virus families.

Outputs
  results/targets/target_prevalence.tsv    target x sample: call + full evidence
  results/targets/target_summary.tsv       one row per target across the cohort
  results/targets/target_contig_index.tsv  every contig assigned to a target,
                                           with its path -- the direct input to
                                           per-taxon phylogenetics or binning

A PRESENCE CALL IS NEVER MADE ON TPM ALONE. It requires the contig and read
evidence thresholds in config/analysis_params.yaml to be met, because a 400 bp
contig with 6 reads can outrank a real genome on a length-normalised metric.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_aliases, load_config, resolve_artefact, write_tsv  # noqa: E402

TARGETS_PATH = ROOT / "metadata" / "target_taxa.tsv"
LINEAGE_COLS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_targets() -> pd.DataFrame:
    rows = []
    with open(TARGETS_PATH) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] == "taxon":
                continue
            rows.append({"target": p[0], "rank_hint": p[1], "category": p[2],
                         "note": p[3] if len(p) > 3 else ""})
    return pd.DataFrame(rows)


def main() -> None:
    cfg = load_config()
    tf = cfg["taxon_filter"]
    aliases = load_aliases()
    out = ROOT / "results" / "targets"
    out.mkdir(parents=True, exist_ok=True)

    targets = load_targets()
    print(f"{len(targets)} target taxa", file=sys.stderr)

    taxon = pd.read_csv(ROOT / "results" / "tables" / "taxon_abundance_long.tsv.gz", sep="\t")
    for c in LINEAGE_COLS:
        if c in taxon:
            taxon[c] = taxon[c].fillna("")

    # ---- map every taxon row to a target ----------------------------------
    # Match at the hinted rank first, then anywhere in the lineage, so a
    # family-level target (Nairoviridae) picks up all genera beneath it.
    assign: dict[int, str] = {}
    for _, t in targets.iterrows():
        name, hint = t["target"], t["rank_hint"]
        cols = ([hint] if hint in LINEAGE_COLS else []) + \
               [c for c in LINEAGE_COLS if c != hint]
        # NCBI prefixes provisional genera with "Candidatus", so the genus of
        # Midichloria mitochondrii is literally "Candidatus Midichloria". An
        # exact match on "Midichloria" silently finds nothing -- which is
        # precisely the wrong failure mode for a pathogen survey, since it
        # reports absence rather than erroring. Accept both spellings.
        variants = {name}
        if name.startswith("Candidatus "):
            variants.add(name[len("Candidatus "):])
        else:
            variants.add(f"Candidatus {name}")
        hit = pd.Series(False, index=taxon.index)
        for c in cols:
            if c not in taxon:
                continue
            hit |= taxon[c].isin(variants)
        n_hit = int(hit.sum())
        if n_hit == 0:
            print(f"  no signal anywhere for target: {name}", file=sys.stderr)
        for i in taxon.index[hit]:
            assign.setdefault(i, name)
    taxon["target"] = pd.Series(assign)
    tt = taxon[taxon["target"].notna()].copy()
    print(f"  {len(tt):,} taxon x sample rows map to a target "
          f"({tt['target'].nunique()} targets seen)", file=sys.stderr)

    # ---- per target x sample ----------------------------------------------
    per = (tt.groupby(["target", "sample_id"])
             .agg(tpm_microbial=("tpm_microbial", "sum"),
                  tpm_total=("tpm_total", "sum"),
                  n_contigs=("n_contigs", "sum"),
                  total_mapped_reads=("sum_mapped_reads", "sum"),
                  longest_contig=("max_contig_len", "max"),
                  max_bitscore=("max_bitscore", "max"),
                  n_taxa=("lca_taxid", "nunique"))
             .reset_index())
    per["present"] = (
        (per["n_contigs"] >= tf["min_contigs_per_taxon"])
        & (per["total_mapped_reads"] >= tf["min_reads_per_taxon"])
    ).astype(int)
    # A blunt confidence tier, so a marginal call is never read as a firm one.
    per["evidence"] = pd.cut(
        per["total_mapped_reads"],
        bins=[-1, tf["min_reads_per_taxon"] - 1, 500, 5000, float("inf")],
        labels=["below_threshold", "low", "moderate", "high"])
    per = per.merge(targets[["target", "category"]], on="target", how="left")
    per = per.sort_values(["target", "tpm_microbial"], ascending=[True, False])
    per.to_csv(out / "target_prevalence.tsv", sep="\t", index=False,
               float_format="%.4f")
    print(f"  wrote results/targets/target_prevalence.tsv ({len(per):,} rows)", file=sys.stderr)

    # ---- cohort summary ----------------------------------------------------
    n_samples = taxon["sample_id"].nunique()
    summ = (per.groupby(["target", "category"])
              .agg(n_samples_present=("present", "sum"),
                   n_samples_any_signal=("sample_id", "nunique"),
                   median_tpm_when_present=("tpm_microbial",
                                            lambda s: s[s > 0].median() if (s > 0).any() else 0),
                   max_tpm=("tpm_microbial", "max"),
                   total_reads=("total_mapped_reads", "sum"),
                   longest_contig=("longest_contig", "max"))
              .reset_index())
    summ["prevalence_pct"] = (100 * summ["n_samples_present"] / n_samples).round(2)
    summ = summ.sort_values("n_samples_present", ascending=False)
    summ.to_csv(out / "target_summary.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"  wrote results/targets/target_summary.tsv", file=sys.stderr)

    # ---- contig index ------------------------------------------------------
    tax_of_target = tt[["lca_taxid", "target"]].drop_duplicates()
    contigs = pd.read_csv(ROOT / "results" / "tables" / "contig_taxon_long.tsv.gz",
                          sep="\t", dtype={"sample_id": str}, low_memory=False)
    ci = contigs.merge(tax_of_target, on="lca_taxid", how="inner")
    lin = pd.read_csv(ROOT / "results" / "taxonomy" / "taxon_lineage.tsv", sep="\t",
                      usecols=["taxid", "lca_name", "genus", "species"])
    ci = ci.merge(lin, left_on="lca_taxid", right_on="taxid", how="left")

    def contig_path(s: str) -> str:
        p = resolve_artefact(s, ROOT / "contigs", ".contigs.fa", aliases)
        return str(p) if p else ""

    paths = {s: contig_path(s) for s in ci["sample_id"].unique()}
    ci["contigs_fasta"] = ci["sample_id"].map(paths)
    ci = ci[["sample_id", "target", "contig_id", "contig_len", "mapped_reads",
             "tpm_microbial", "lca_taxid", "lca_name", "genus", "species",
             "top_bitscore", "contigs_fasta"]]
    ci = ci.sort_values(["target", "sample_id", "contig_len"],
                        ascending=[True, True, False])
    ci.to_csv(out / "target_contig_index.tsv", sep="\t", index=False,
              float_format="%.4f")
    print(f"  wrote results/targets/target_contig_index.tsv ({len(ci):,} contigs)",
          file=sys.stderr)

    # per-target seqid lists, ready for `seqtk subseq` / `samtools faidx`
    sd = out / "seqids"; sd.mkdir(exist_ok=True)
    for tgt, g in ci.groupby("target"):
        with open(sd / f"{tgt.replace(' ', '_')}.contigs.tsv", "w") as fh:
            fh.write("sample_id\tcontig_id\tcontigs_fasta\n")
            for _, r in g.iterrows():
                fh.write(f"{r['sample_id']}\t{r['contig_id']}\t{r['contigs_fasta']}\n")
    print(f"  wrote {len(ci.groupby('target'))} per-target seqid lists to "
          f"results/targets/seqids/", file=sys.stderr)

    # ---- console summary ---------------------------------------------------
    print(f"\nprevalence across {n_samples} samples", file=sys.stderr)
    print(f"  {'target':<22} {'category':<13} {'n':>4} {'%':>6}  {'max TPM':>10} "
          f"{'longest contig':>14}", file=sys.stderr)
    for _, r in summ.iterrows():
        if r["n_samples_present"] == 0:
            continue
        print(f"  {r['target']:<22} {r['category']:<13} {int(r['n_samples_present']):>4} "
              f"{r['prevalence_pct']:>6.1f}  {r['max_tpm']:>10.0f} "
              f"{int(r['longest_contig']):>14,}", file=sys.stderr)


if __name__ == "__main__":
    main()
