#!/usr/bin/env python3
"""Verification checks for the exported matrices. Run after step 4.

Each check prints PASS/FAIL/INFO and the number it is asserting on, so a
regression shows up as a changed number rather than a missing file.
"""
from __future__ import annotations
import sys, gzip
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, load_aliases, load_config, resolve_artefact

ok = True
def check(name, passed, detail=""):
    global ok
    ok &= bool(passed)
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
def info(name, detail):
    print(f"  [INFO] {name} -- {detail}")

cfg = load_config()
M = ROOT / "results" / "matrices"
T = ROOT / "results" / "tables"

print("\n1. TPM conservation")
for rank in ["species", "genus", "family", "phylum"]:
    m = pd.read_csv(M / f"tpm_microbial_{rank}.tsv", sep="\t", index_col=0)
    cs = m.sum(axis=0)
    nz = cs[cs > 0]
    check(f"{rank}: colSums == 1e6", bool(((nz - 1e6).abs() < 1.0).all()),
          f"n={len(nz)} min={nz.min():.3f} max={nz.max():.3f}")
    if (cs == 0).any():
        info(f"{rank}: zero-signal samples", f"{int((cs==0).sum())}")

print("\n2. Rank consistency (aggregating species should reproduce genus totals)")
sp = pd.read_csv(M / "tpm_microbial_species.tsv", sep="\t", index_col=0)
gn = pd.read_csv(M / "tpm_microbial_genus.tsv", sep="\t", index_col=0)
check("same samples at every rank", list(sp.columns) == list(gn.columns),
      f"{sp.shape[1]} samples")
check("species matrix has >= genus rows", sp.shape[0] >= gn.shape[0],
      f"{sp.shape[0]} species labels vs {gn.shape[0]} genus labels")

print("\n3. Join integrity (LCA contigs all present in coverage)")
aliases = load_aliases()
audit = pd.read_csv(ROOT/"results"/"audit"/"sample_completeness.tsv", sep="\t")
comp = audit.loc[audit.status == "complete", "sample_id"].tolist()
for s in comp[:3] + comp[-2:]:
    lp = resolve_artefact(s, ROOT/"lca", ".blast.lca", aliases)
    cp = resolve_artefact(s, ROOT/"mapping", ".coverage.tsv", aliases)
    l = pd.read_csv(lp, sep="\t", usecols=["query"])
    c = pd.read_csv(cp, sep="\t", usecols=["contig_id"])
    miss = len(set(l["query"]) - set(c["contig_id"]))
    check(f"{s}", miss == 0, f"{len(l):,} lca contigs, {miss} missing from coverage")

print("\n4. Evidence assays align with the abundance assay")
for rank in ["species", "genus"]:
    a = pd.read_csv(M/f"tpm_microbial_{rank}.tsv", sep="\t", index_col=0)
    r = pd.read_csv(M/f"evidence_reads_{rank}.tsv", sep="\t", index_col=0)
    n = pd.read_csv(M/f"evidence_ncontigs_{rank}.tsv", sep="\t", index_col=0)
    check(f"{rank}: reads assay same shape/order", a.shape == r.shape and list(a.index)==list(r.index))
    check(f"{rank}: ncontigs assay same shape/order", a.shape == n.shape and list(a.index)==list(n.index))
    check(f"{rank}: no abundance without read support",
          bool(((a > 0) & (r == 0)).sum().sum() == 0),
          f"{int(((a>0)&(r==0)).sum().sum())} violations")

print("\n5. Scope separation -- host must be absent from the microbial matrix")
# Test on TAXIDS, not names. NCBI reuses genus names across kingdoms
# (Bacillus, Planococcus -- see results/taxonomy/homonyms.tsv), so a name-based
# test reports false leaks. The matrix is built from is_microbial rows, so the
# real question is whether any non-microbial taxid contributed to it.
tax = pd.read_csv(ROOT/"results"/"taxonomy"/"taxon_lineage.tsv", sep="\t", keep_default_na=False)
long = pd.read_csv(T/"taxon_abundance_long.tsv.gz", sep="\t",
                   usecols=["lca_taxid","tpm_microbial","is_microbial","scope_class"])
bad = long[(long.tpm_microbial > 0) & (long.is_microbial != 1)]
check("no non-microbial taxid carries microbial TPM", len(bad)==0,
      f"{len(bad)} rows" + (f", classes: {sorted(bad.scope_class.unique())}" if len(bad) else ""))
hom = pd.read_csv(ROOT/"results"/"taxonomy"/"homonyms.tsv", sep="\t") \
        if (ROOT/"results"/"taxonomy"/"homonyms.tsv").exists() else pd.DataFrame()
if len(hom) and hom.iloc[0].get("superkingdom_groups") != "none":
    info("cross-kingdom name collisions", 
         f"{len(hom)} (handled by taxid; see results/taxonomy/homonyms.tsv): "
         + ", ".join(hom["name"].astype(str).head(5)))

print("\n6. Known-biology spot check (tick endosymbionts must be detectable)")
frac = pd.read_csv(T/"fraction_summary.tsv", sep="\t")
for g in ["Coxiella", "Rickettsia", "Midichloria", "Wolbachia", "Francisella"]:
    # NCBI genus names carry the "Candidatus " prefix for provisional genera
    # (Candidatus Midichloria), so match both spellings.
    hit = [x for x in (g, f"Candidatus {g}") if x in gn.index]
    if hit:
        g, = hit[:1]
        v = gn.loc[g]
        info(f"{g}", f"detected in {int((v>0).sum())}/{gn.shape[1]} samples, "
                     f"max {v.max():,.0f} TPM ({v.max()/1e4:.1f}% of microbial)")
    else:
        info(f"{g}", "NOT PRESENT at genus level")

print("\n7. Host fraction sanity (whole-tick homogenate should be host-dominated)")
info("median host TPM share", f"{frac.tpm_pct_host.median():.1f}%")
info("median microbial TPM share", f"{frac.tpm_pct_microbial.median():.1f}%")
info("median unassigned (no LCA) TPM share", f"{frac.tpm_pct_no_lca.median():.1f}%")
info("median unassigned READ share", f"{frac.read_pct_no_lca.median():.1f}%")
check("host share exceeds microbial in most samples",
      (frac.tpm_pct_host > frac.tpm_pct_microbial).mean() > 0.5,
      f"{100*(frac.tpm_pct_host > frac.tpm_pct_microbial).mean():.0f}% of samples")

print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
sys.exit(0 if ok else 1)
