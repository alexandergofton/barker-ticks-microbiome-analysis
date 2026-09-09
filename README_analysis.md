# Tick microbiome / pathogen composition — export and analysis

Joins the LCA taxonomy in `lca/` to the per-contig read mapping in `mapping/`,
producing a TPM abundance matrix and a `TreeSummarizedExperiment` for
descriptive compositional analysis, plus a handoff for the targeted
pathogen/endosymbiont work.

---

## The problem this solves

`code/lca_to_kreport.py` builds each Pavian kreport by counting **contigs** —
one contig counts as one "read". A taxon on a single 300 bp contig at 2×
coverage therefore scores the same as a taxon on a 500 kb contig at 300×. The
abundance information lives in `mapping/*.coverage.tsv` and was never used.

Everything here joins the two on contig id and weights taxonomy by TPM.

> The kreports are still regenerated and still work in Pavian. Treat them as a
> **richness** view only. Abundance comes from `results/matrices/`.

---

## Three things to know before reading any result

**1. ~88% of reads have no BLAST hit, and that is not a bug.**

Median across 118 samples: **94.2% of TPM and 87.9% of mapped reads sit on
contigs with no BLAST hit at all.** This was cross-checked against the
independent read-level Kraken2 track — for `0143-F1-Am-sp`, Kraken2 leaves
**96.9% of reads unclassified**, having never seen the assembly or the BLAST
output. Two methodologically independent classifiers fail on the same bulk of
the data.

The explanation is that Australian tick genomes are essentially absent from
`nt`. The tick's own DNA — most of a whole-tick homogenate — yields no hit and
lands in "unassigned".

Consequences, which are not optional caveats:
- `tpm_pct_host` (median 1.4%) **massively understates** host content. It
  counts only tick contigs similar enough to something already in `nt`.
- **No statement of the form "X% of this sample is microbial" is supportable.**
- The primary assay is `tpm_microbial`, renormalised *within* the microbial
  fraction. That is the object to analyse, and it is unaffected by the above,
  because it never uses the unassigned fraction as a denominator.

**2. There are no controls.** No extraction blanks, no no-template controls.
The experiment was not designed for microbiome quantification. Contaminant taxa
are **annotated, never removed** (`metadata/reagent_contaminants.tsv`).
Composition is descriptive.

**3. TPM alone never makes a presence call.** TPM divides by contig length, so a
400 bp contig with 6 reads can outrank a real genome. `mapped_reads`,
`n_contigs`, `contig_len` and `bitscore` travel with every row as **detection
evidence** and gate every presence call. They are not an alternative abundance
and must not be used for count-based differential abundance — they are
assembly-derived contig read counts, not a designed count matrix.

---

## Sample accounting

From `results/audit/audit_summary.txt`:

| status | n | meaning |
|---|---|---|
| complete | 120 | has both LCA and coverage |
| needs_lca | 12 | has BLAST + contigs + coverage, LCA never ran |
| needs_mapping | 24 | has LCA, read mapping never ran |
| excluded_no_reads | 9 | contigs exist, reads absent — cannot be mapped |
| excluded_no_blast | 2 | contigs exist, BLAST never produced output |

118 samples reach the abundance matrix (120 complete, minus `B5093_F2_S8`
— 10 contigs and 61 reads total, too shallow to analyse — plus the two alias
recoveries). Every exclusion is recorded with a diagnosed reason in
`results/audit/excluded_samples.tsv`.

**Why the 12 LCA backfills never ran:** `code/lca.py:31` does
`pd.read_csv(blast_file)` on the whole BLAST output. Seven of those files are
10–150 GB (largest: `4480-F1-Isimp` at 149 GB), so they OOM.
`code/analysis/lca_stream.py` streams by query block instead — measured **flat
2.34 GB RSS at 225k lines/s on the 149 GB file**, independent of input size.

Two filename aliases are resolved without renaming anything on disk
(`metadata/sample_aliases.tsv`):
- `8334-AMLoculosum-BriddledTurn7` (in `lca/`) ≡ `8334-AMLoculosum-BriddledTurn`
- `B6852_F1_38_1` (in `contigs/`) ≡ `B6852_F1_38_S1`

---

## Running it

### One-time setup

```bash
mamba env create -f env/tick_micro.yml      # R 4.3 + mia + vegan + decontam
python3 code/analysis/make_merged_filtered.py   # builds taxdump/merged.filtered.dmp
```

`merged.filtered.dmp` rescues staxids NCBI has merged since the local April 2026
dump; without it they silently became "unclassified". `nodes.dmp`/`names.dmp` are
deliberately **not** updated — the 120 existing LCA calls were made against them
and replacing them would make backfilled samples incomparable. The current
`merged.dmp` is filtered to targets that exist locally (2,024 of 100,581 rows
dropped), because taxopy raises `KeyError` on the rest.

### The two cluster steps

These are the only compute-heavy steps. Both are POSIX `sh` — this site's Slurm
runs batch scripts with `/bin/sh` regardless of the shebang, so they contain no
`[[ ]]`, `(( ))`, `mapfile`, `<<<`, arrays or `bc`.

```bash
bash code/analysis/submit_backfills.sh          # prints the sbatch lines
bash code/analysis/submit_backfills.sh --go     # submits them
```

The helper sizes `--array` to the current queue length and passes `TICK_ROOT`,
which matters because under Slurm `$0` is the spool copy
(`/cm/local/.../slurm_script`), not the script's path in the repo. Queue files
hold paths **relative** to the project root, so they work whether the filesystem
is mounted at `/data/hb-austick/...` or `/datasets/work/hb-austicks/...`.

Both scripts skip samples whose output already exists, so they are safe to
re-run. If `python3`/`bowtie2` need a module on the compute nodes, set
`PYTHON=` / `BOWTIE2_BIN=` in the environment.

### Everything else

```bash
bash run_all.sh
```

Re-run this after each backfill finishes to pull the new samples in.

---

## Pipeline steps

| script | does |
|---|---|
| `00_audit_samples.py` | reconciles every sample across every artefact; writes the two backfill queues |
| `make_merged_filtered.py` | builds `taxdump/merged.filtered.dmp` |
| `lca_stream.py` | streaming replacement for `code/lca.py`; constant memory |
| `01_backfill_lca.sh` | Slurm array over `lca_backfill_queue.tsv` |
| `02_backfill_mapping.sh` | Slurm array over `read_dir_priority.tsv`; wraps the existing `code/map_reads_to_contigs.sh` unchanged |
| `03_build_taxonomy.py` | lineage + scope + contaminant annotation per taxid |
| `04_build_matrices.py` | **the core join**; contig filter, TPM renormalisation, wide matrices |
| `04b_filter_sensitivity.py` | re-runs under lenient/default/strict to show the filter's effect |
| `05_make_metadata_template.py` | pre-fills the metadata template |
| `06_validate_metadata.R` | validates your filled-in metadata |
| `07_build_tse.R` | assembles the TreeSummarizedExperiment |
| `08_composition.R` | composition, contamination diagnostics, Bracken cross-check |
| `09_diversity.R` | alpha/beta diversity, PERMANOVA incl. batch terms |
| `10_export_targets.py` | pathogen/endosymbiont tables + contig index |
| `verify.py` | assertions over the exported matrices |

### Marker genes for the targets — `bash run_markers.sh`

Steps 11–16 are a separate, second-stage workflow that takes over where step 10
stops. Step 10 says *this genus is present, at this abundance*; it cannot say
which **species**, which for `Coxiella`, `Rickettsia` and `Francisella` — each
holding both endosymbionts and genuine pathogens — is the whole question.

Steps 11–16 find, extract, label and align the phylogenetic marker genes that
do resolve them (`gltA`, `ompA`/`ompB`/`sca4`, `flaB`, `msp4`, `rpoB`, `cox1`,
`cytb`, RdRp, and 16S/18S), per target, and hand you alignments to add your own
reference sequences to. **See `README_markers.md`.**

| script | does |
|---|---|
| `11_fetch_marker_refs.py` | reference probes per locus from NCBI Entrez; cached |
| `12_extract_target_contigs.py` | pulls target contigs from the assemblies, relabelled |
| `13_find_markers.py` | blastn / tblastn; assembles loci; resolves Sca-family paralogues |
| `14_extract_markers.py` | cuts out the loci, fixes the reading frame, names everything |
| `15_align_markers.py` | MAFFT per target × marker, codon-aware for protein-coding loci |
| `16_marker_report.py` | what was recovered, and what was not |

The panel of loci is `metadata/target_markers.tsv` (120 loci over 32 targets) and
the settings are the `markers:` block of `config/analysis_params.yaml`. Needs
`env/tick_markers.yml` for MAFFT.

---

## Key outputs

```
results/
  audit/         sample_completeness.tsv, excluded_samples.tsv, audit_summary.txt
  taxonomy/      taxon_lineage.tsv, homonyms.tsv
  tables/        contig_taxon_long.tsv.gz   <- source of truth, one row per contig
                 taxon_abundance_long.tsv.gz
                 fraction_summary.tsv       <- per-sample host/microbial/unassigned
                 filter_sensitivity_*.tsv
  matrices/      tpm_microbial_{species,genus,family,phylum}.tsv   <- PRIMARY
                 evidence_{reads,ncontigs}_*.tsv                   <- gating only
  tse/           tick_microbiome_tse.rds  (+ _genus/_family/_phylum)
  targets/       target_prevalence.tsv, target_summary.tsv,
                 target_contig_index.tsv, seqids/<target>.contigs.tsv
  markers/       sequences/<target>/<target>__<marker>.fna     <- extracted loci
                 alignments/<target>/<target>__<marker>.aln.fasta
                 marker_summary.tsv, marker_by_sample.tsv, marker_gaps.tsv
  figures/       01..09 png
```

### Using the TSE

```r
library(mia)
tse <- readRDS("results/tse/tick_microbiome_tse.rds")

assay(tse, "tpm_microbial")   # PRIMARY abundance, sums to 1e6 per sample
assay(tse, "reads")           # EVIDENCE ONLY — gating and confidence
assay(tse, "n_contigs")       # EVIDENCE ONLY

rowData(tse)$contaminant_flag # annotation; nothing has been removed
metadata(tse)$caveats
```

Agglomerated rank objects are relabelled to match `results/matrices/` and
cross-checked: **max per-sample TPM difference 0.0000 at every rank**. Row
counts differ slightly (1246 vs 1226 at genus) because `mia` keys rows on the
full lineage, which correctly keeps cross-kingdom homonyms apart — `Bacillus`
and `Planococcus` are each both a bacterial genus and an arthropod genus in
NCBI (`results/taxonomy/homonyms.tsv`).

---

## What you need to supply

**`metadata/sample_metadata.tsv`** — fill in
`metadata/sample_metadata_TEMPLATE.tsv` and save under that name.

Pre-filled for you: `sample_id`, `status`, `flowcell`, `index_i7`, `index_i5`,
`lane`, `seq_run`, `sample_sheet_index`, `seq_submission`, plus **guesses** —
`tick_species_guess` (84/167, from `metadata/tick_species_abbrev.tsv`),
`life_stage_guess` and `n_ticks_guess` (125/167, parsed from the `F1`/`2F`/`N3`
tokens). The guesses are provisional and are superseded by your collaborator's
mitogenome identifications.

To fill by hand: `tick_species`, `life_stage`, `n_ticks_pooled`, `host_species`,
`collection_site`, `state`, `latitude`, `longitude`, `collection_date`,
`dna_conc_ng_ul`, `extraction_batch`.

**`dna_conc_ng_ul` is worth more than it looks.** It enables
`decontam::isContaminant(method="frequency")`, which needs concentrations rather
than blanks — the only contamination test this design can actually support.
Without it, step 8 skips that analysis.

---

## Current findings (118 samples, provisional metadata)

**Batch effect is real but small.** PERMANOVA on Bray–Curtis: flowcell/lane
p=0.001–0.005 but R²≈0.029. Report it; do not adjust it away.

**Tick species dominates.** On the provisional id-string guess (70/118 samples,
21 levels): R²=0.38–0.40, p=0.007 — an order of magnitude larger than the
technical terms. Worth confirming once real identifications are in.

**Richness is depth-confounded.** Observed genera vs contigs-passing-filter:
Spearman ρ=0.73. Figure 06 plots this explicitly rather than hiding it. Do not
interpret richness differences as biology without accounting for it.

**Assembly-LCA vs read-level Bracken:** median Spearman ρ=0.391 across 85
shared samples. Moderate, not strong. The two methods measure different things
(assembled contigs vs raw reads) and disagree most where assembly is poor —
this is itself a result, and a reason to treat low-abundance calls cautiously.

**Targets** (`results/targets/target_summary.tsv`), prevalence across 119 samples:

| target | category | n | % | longest contig |
|---|---|---|---|---|
| *Rickettsia* | both | 65 | 54.6 | 1,157,618 |
| *Francisella* | both | 48 | 40.3 | 682,921 |
| Flaviviridae | pathogen | 29 | 24.4 | 25,116 |
| Phenuiviridae | pathogen | 23 | 19.3 | 11,979 |
| *Coxiella* | both | 19 | 16.0 | 606,599 |
| *Arsenophonus* | endosymbiont | 16 | 13.4 | 7,436 |
| *Rickettsiella* | endosymbiont | 14 | 11.8 | 228,266 |
| *Babesia* | pathogen | 11 | 9.2 | 4,338 |
| Nairoviridae | pathogen | 11 | 9.2 | 8,566 |
| *Hepatozoon* | pathogen | 8 | 6.7 | 9,430 |
| *Borrelia* | pathogen | 6 | 5.0 | 12,332 |
| *Candidatus* Midichloria | endosymbiont | 6 | 5.0 | 72,486 |

The megabase-scale *Rickettsia*, *Francisella* and *Coxiella* contigs indicate
near-complete genome recovery in the samples that carry them.

> **A genus-level hit in *Coxiella*, *Rickettsia* or *Francisella* does not
> establish pathogen presence.** Each genus contains obligate tick
> endosymbionts *and* genuine pathogens (*C. burnetii*, *F. tularensis*,
> spotted-fever group) that are near-identical over short contigs. Resolving
> which is present is exactly what the follow-up analysis is for —
> `results/targets/target_contig_index.tsv` and
> `results/targets/seqids/<target>.contigs.tsv` give you the contigs and their
> FASTA paths to start from.
>
> The same caution applies to *Borrelia* in an Australian context.

---

## Tuning

`config/analysis_params.yaml` holds every threshold. The contig evidence filter
defaults to `min_contig_len: 300`, `min_mapped_reads: 10`,
`min_top_bitscore: 100`; a taxon is called present at `>= 1` contig and `>= 50`
reads.

### How much does the filter change the answer?

From `results/tables/filter_sensitivity_*.tsv` (30-sample subset):

| profile | contigs passing | TPM filtered out | top genera |
|---|---|---|---|
| lenient | 100.0% | 0.0% | *Rickettsia* 17%, *Klebsiella* 9%, *Francisella* 8%, *Iltovirus* 8% |
| **default** | **28.3%** | **10.2%** | ***Rickettsia* 22%, *Francisella* 14%, *Iltovirus* 10%, *Paenibacillus* 7%** |
| strict | 3.3% | 16.5% | *Rickettsia* 28%, *Francisella* 23%, *Paenibacillus* 6%, *Enterobacter* 5% |

Top-20 genus rank correlation: lenient↔default ρ=0.78, default↔strict ρ=0.69.
**The dominant taxa are robust to the threshold choice.** *Rickettsia* and
*Francisella* *grow* as the filter tightens (17%→28%, 8%→23%) — the expected
behaviour for real genomes with deep, long-contig support, while thin
assemblies fall away. *Klebsiella* and *Iltovirus* drop out entirely under
strict and should be treated with suspicion.

Note that `tpm_pct_no_lca` (~80% here) is **unchanged** across profiles by
construction — it is measured against the full coverage table, and the contig
filter cannot affect contigs that never had an LCA. What the filter moves is
`tpm_pct_filtered_out`: assigned-but-low-evidence contigs that should not be
setting the denominator.
