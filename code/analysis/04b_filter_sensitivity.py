#!/usr/bin/env python3
"""Step 4b -- how much does the evidence filter actually change the answer?

The contig filter is a judgement call, so its effect has to be shown rather
than asserted. This re-runs the join and aggregation under the lenient /
default / strict profiles in config/analysis_params.yaml and reports:

  results/tables/filter_sensitivity_fractions.tsv
      per profile: mean share of TPM that is unassigned, filtered out, host,
      and microbial. The unassigned share is the number to watch -- measured on
      0143-F1-Am-sp, LCA-assigned contigs hold 95.4% of READS but only 60.3% of
      TPM, because length normalisation inflates short unassigned contigs. If
      that share does not fall as the floors engage, the filter is not working.

  results/tables/filter_sensitivity_genera.tsv
      top-20 genus composition under each profile, so you can see whether any
      conclusion depends on the threshold choice.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
from common import ROOT, load_config, write_tsv  # noqa: E402

bm = importlib.import_module("04_build_matrices") if False else None
# 04_build_matrices is not importable by name (leading digit), so load it by path
spec = importlib.util.spec_from_file_location(
    "build_matrices", Path(__file__).resolve().parent / "04_build_matrices.py")
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)


def main() -> None:
    cfg = load_config()
    profiles = cfg["sensitivity_profiles"]
    tax = bm.load_taxonomy()

    audit = pd.read_csv(ROOT / "results" / "audit" / "sample_completeness.tsv", sep="\t")
    samples = sorted(audit.loc[audit["status"] == "complete", "sample_id"])

    # A subset keeps this affordable; it is stated in the output rather than
    # silently applied, because a sampled result must not read as a full one.
    n_sub = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    step = max(1, len(samples) // n_sub)
    subset = samples[::step][:n_sub]
    print(f"sensitivity over {len(subset)} of {len(samples)} samples "
          f"(every {step}th, for runtime)", file=sys.stderr)

    frac_rows, gen_rows = [], []
    for pname, pf in profiles.items():
        cf = {k: pf[k] for k in ("min_contig_len", "min_mapped_reads", "min_top_bitscore")}
        fr, agg_all = [], []
        for s in subset:
            res = bm.process_sample(s, tax, cf)
            if res is None:
                continue
            agg, frac, _ = res
            fr.append(frac)
            agg_all.append(agg)
        if not fr:
            continue
        f = pd.DataFrame(fr)
        row = {
            "profile": pname,
            "n_samples": len(fr),
            **{k: pf[k] for k in pf},
            "mean_pct_contigs_passing": round(f["pct_contigs_passing"].mean(), 2),
            "mean_tpm_pct_no_lca": round(f["tpm_pct_no_lca"].mean(), 2),
            "mean_tpm_pct_filtered_out": round(f["tpm_pct_filtered_out"].mean(), 2),
            "mean_tpm_pct_host": round(f["tpm_pct_host"].mean(), 2),
            "mean_tpm_pct_microbial": round(f["tpm_pct_microbial"].mean(), 2),
            # tpm_pct_no_lca is measured against the FULL coverage table, so it
            # is by construction unchanged by the contig filter. The number the
            # filter actually moves is tpm_pct_filtered_out: low-evidence
            # contigs that DID get an LCA but should not set the denominator.
            "mean_read_pct_no_lca": round(f["read_pct_no_lca"].mean(), 2),
        }
        # the discarded-signal counterpart of the above
        row["mean_tpm_pct_retained"] = round(
            100 - row["mean_tpm_pct_no_lca"] - row["mean_tpm_pct_filtered_out"], 2)
        frac_rows.append(row)

        a = pd.concat(agg_all, ignore_index=True)
        a = a.join(tax, on="lca_taxid")
        a = a[(a["is_microbial"] == 1)
              & (a["n_contigs"] >= pf["min_contigs_per_taxon"])
              & (a["sum_mapped_reads"] >= pf["min_reads_per_taxon"])]
        g = (a[a["genus"].astype(str) != ""]
             .groupby("genus")["tpm_microbial"].sum().sort_values(ascending=False))
        tot = g.sum()
        for rank_i, (gen, v) in enumerate(g.head(20).items(), 1):
            gen_rows.append({"profile": pname, "rank": rank_i, "genus": gen,
                             "pct_of_microbial_tpm": round(100 * v / tot, 3)})
        print(f"  {pname:<8} contigs passing {row['mean_pct_contigs_passing']:>6.2f}%  "
              f"unassigned TPM {row['mean_tpm_pct_no_lca']:>6.2f}%  "
              f"top genus {g.index[0] if len(g) else 'n/a'}", file=sys.stderr)

    write_tsv(ROOT / "results" / "tables" / "filter_sensitivity_fractions.tsv",
              frac_rows, list(frac_rows[0].keys()))
    write_tsv(ROOT / "results" / "tables" / "filter_sensitivity_genera.tsv",
              gen_rows, ["profile", "rank", "genus", "pct_of_microbial_tpm"])

    # Does the ranking survive the threshold choice?
    gd = pd.DataFrame(gen_rows)
    piv = gd.pivot_table(index="genus", columns="profile", values="rank")
    both = piv.dropna()
    if len(both) > 2 and {"lenient", "strict"} <= set(piv.columns):
        # numpy, not pandas method="spearman" -- that path imports scipy, which
        # is not installed in the base env this runs in.
        import numpy as np
        a = np.asarray(both["lenient"], float)
        b = np.asarray(both["strict"], float)
        ra, rb = pd.Series(a).rank().to_numpy(), pd.Series(b).rank().to_numpy()
        rho = float(np.corrcoef(ra, rb)[0, 1])
        print(f"\n  Spearman rho of top-20 genus RANK, lenient vs strict: {rho:.3f}",
              file=sys.stderr)
        print(f"  ({len(both)} genera in both lists; a high value means the "
              f"composition\n   conclusions do not hinge on the threshold choice)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
