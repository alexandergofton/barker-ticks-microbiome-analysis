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

### Samples removed by decision — `metadata/excluded_samples.tsv`

Ten samples are removed from **every** analysis (microbiome, markers, trees,
mitogenome) after the step-26c mitogenome audit (`README_mito.md`): wrong tick
in the tube, a mixed tube, or a contaminating tick's mitogenome. The file gives
each one with its reason. Step 00 gives them the status `excluded_manual`,
which every later step honours (matrices, metadata, targets, maps). Step 21
reads `contigs/` directly, so it skips them too (`common.py:load_excluded`).
Nothing on disk is touched. To restore a sample, delete its row and re-run.

### Contigs re-assigned by decision — `metadata/contig_reassignments.tsv`

Some contigs carry a wrong nt LCA because the **viral reference records** contain
tick genomic repeats (see `extracting_pathogens/*/are_these_real.html`). They
are re-assigned contig by contig. Steps 03 and 04 apply the table on top of
`lca/<s>.blast.lca` (`common.py:load_contig_reassignments`), so a re-assigned
contig is filtered exactly like any other contig of its new taxon. The LCA files
themselves are not edited. Current contents (27,080 contigs, all → Ixodida, so
scope `host`):

- every *Cheeloo Jingmen-like virus* contig (16,699);
- every *North shore virus* contig (10,381) except the 9 with RdRp protein
  support (8143, 8337-Iholo-Xbtp-Awzh, B5119_F1_S15, B6836_N1_S13).

The table is regenerated with
`python3 extracting_pathogens/curate_tick_repeat_contigs.py`. To undo a
re-assignment, remove its rows and re-run from step 03.

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
| `05b_import_metadata.py` | builds `metadata/sample_metadata.tsv` from the curated CSV; logs issues |
| `06_validate_metadata.R` | validates your filled-in metadata |
| `07_build_tse.R` | assembles the TreeSummarizedExperiment |
| `08_composition.R` | composition, contamination diagnostics, Bracken cross-check |
| `09_diversity.R` | alpha/beta diversity, PERMANOVA incl. batch terms |
| `10_export_targets.py` | pathogen/endosymbiont tables + contig index |
| `10b_tick_group_comparison.R` | composition + endosymbionts compared between tick genera/species |
| `10c_sampling_maps.R` | static and interactive sampling maps |
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

### Trees for those markers — `bash run_trees.sh`

Steps 17–19 close the loop steps 11–16 leave open. Step 15 hands you an
alignment of sample sequences with nothing to compare them to; these steps
download the named-species references each locus needs, cut the locus out of
them, align them with the samples and build a tree. **See `README_trees.md`.**

| script | does |
|---|---|
| `17_fetch_tree_refs.py` | downloads the reference sequences that go INTO the trees; cached |
| `18_curate_tree_refs.py` | blasts each downloaded record and cuts the marker locus out of it |
| `19_build_trees.py` | aligns samples + references, one FastTree per locus |
| `19b_tree_figures.R` | roots each tree and draws it |

Two reference panels drive it: `metadata/tree_references.tsv` (per target ×
marker) and `metadata/tree_ref_genomes.tsv` (per target — several loci exist
only inside complete genomes, and step 18 cuts them out). Settings are the
`trees:` block of `config/analysis_params.yaml`.

These are **screening trees** — FastTree, SH-like local support, seconds per
locus. They resolve endosymbiont from pathogen well enough to decide which loci
deserve a proper partitioned analysis; they are not the paper's figure.

### Tick mitogenome phylogeny — `bash run_mito.sh`

A separate, independent workflow: extracts one mitochondrial genome candidate
per sample, annotates its genes, aligns each of the 13 protein-coding genes
separately, and builds one partitioned ML tree in IQ-TREE (ModelFinder, one
substitution model per gene) -- a mitogenome phylogeny of the ticks
themselves. **See `README_mito.md`**, including why MITOS2 (the intended
annotator) is installed but running on a homology-based fallback here: its
reference data is Zenodo-hosted and unreachable from this compute node.

| script | does |
|---|---|
| `20_fetch_mito_refs.py` | reference probes for 13 PCGs + 2 rRNAs; cached |
| `21_find_mito_candidates.py` | picks ONE mitogenome contig per sample by homology |
| `22_annotate_mitogenomes.py` | final per-gene sequences; MITOS2 if available, else homology |
| `23_align_mito_genes.py` | MAFFT per gene, codon-aware for the 13 PCGs |
| `24_concat_partition.py` | concatenates the PCGs into one partitioned supermatrix |
| `25_run_iqtree.py` | ModelFinder + partitioned ML tree, UFBoot + SH-aLRT |
| `26_mito_figures.R` | tree figure + two distance tables (p-distance, patristic) |

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
                 trees/<target>/<target>__<marker>.{aln.fna,nwk,rooted.nwk,png}
                 trees/tree_summary.tsv, trees/tree_reference_sequences.tsv
  mito/          genome/<sample>.mito.fna, genes/<GENE>.fna, tree/mito.treefile
                 figures/mito_tree.png, tables/mito_distance_{pdistance,patristic}.tsv
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

## Sample metadata

**The source of truth is `metadata/barker_microbe_analysis_metadata.csv`**:
one row per sample, with `tick_genus`, `tick_species` (epithet), `life_stage`,
`n_ticks_pooled`, `host_type`, `host_species`, `location`, `latitude` and
`longitude`. Edit that file, not anything derived from it.

Step 5b (`05b_import_metadata.py`) joins it to the fields step 5 derives
automatically: flowcell, lane, index and `seq_submission`, which are the only
technical covariates this design has, plus the id-string guesses. It writes
`metadata/sample_metadata.tsv`, which is **regenerated on every run**. Step 5b
only fixes formatting: stray whitespace, and a `"lat, lon"` pair pasted into
the latitude column. It never changes content. Anything that looks wrong goes
to `results/audit/metadata_import_issues.tsv` with level `review`, for you to
check against the source. Examples are a species that disagrees with the
sample-id string, or a `host_type` that disagrees with the host named.

Derived columns: `tick_species` is the binomial (`Ixodes tasmani`, or
`Ixodes sp.`), `tick_epithet` keeps the bare epithet, `species_identified`
is `no` for `sp.`, and `tick_family` is set from the genus.

Not supplied, and therefore not analysed: collection date, state, DNA
concentration, extraction batch. Without `dna_conc_ng_ul`, step 8 cannot run
the decontam frequency method, which is the one contamination test this design
supports without blanks.

---

## Tick genus / species comparison — step 10b

`10b_tick_group_comparison.R` writes tables to `results/tables/tick_groups/`
and figures 10–15. The settings are in the `tick_groups:` block of the config.

| output | question |
|---|---|
| `10_samples_by_tick_and_host.png`, `sample_counts.tsv` | who was sampled; genus × host type confounding |
| `11a/11b_composition_by_tick_{genus,species}.png`, `mean_composition_by_*.tsv` | mean composition per group (heatmap, top 20 genera) |
| `12_pcoa_{bray,aitchison}_tick_genus.png`, `permanova.tsv`, `permanova_pairwise_genus.tsv` | does composition differ between genera, between species within a genus, and by host within a genus? |
| `13a/13b_endosymbiont_prevalence_by_*.png`, `endosymbiont_prevalence_by_*.tsv`, `endosymbiont_fisher_tests.tsv` | which endosymbionts each tick genus/species carries |
| `14_endosymbiont_abundance_by_genus.png`, `endosymbiont_abundance_kruskal.tsv` | how dominant an endosymbiont is in the ticks that carry it |
| `15_dominant_genus_by_tick_species.png`, `dominant_genus_per_sample.tsv` | the single most abundant microbial genus in each tick |
| `endosymbiont_cooccurrence.tsv` | symbiont pairs that co-occur or exclude each other, raw and stratified by genus |
| `indicator_genera_by_tick_genus.tsv` | microbial genera characteristic of a tick genus (IndVal.g) |
| `mantel_composition_vs_geography.tsv` | does composition track geographic distance, overall and within genus/species? |

How the tests are built:
- **Groups with fewer than 3 samples are drawn but never tested.**
  Unidentified (`sp.`) ticks count toward their genus but not toward any
  species.
- **Tick genus and host type are confounded.** *Amblyomma* is mostly on
  reptiles; *Ixodes* and *Haemaphysalis* are on mammals. Host type is
  therefore tested with permutations restricted to within-genus blocks
  (`host type | genus`). A marginal model (`~ genus + host + seq_submission`,
  `by = "margin"`) gives each term's unique share.
- **PERMDISP sits beside every single-factor PERMANOVA**
  (`dispersion_p`). When both are significant, part of the "difference" is
  spread, not location.
- **Endosymbiont presence uses the step-10 call** (reads and contigs), never
  TPM alone. *Rickettsia*, *Francisella* and *Coxiella* are genus-level
  calls (marked `*`). Each of these genera contains both endosymbionts and
  pathogens.
- **Co-occurrence is also tested stratified by tick genus** (exact
  Mantel–Haenszel). Two symbionts that merely favour the same tick genus are
  then not reported as associated.
- **"Keystone" here means the dominant maternally inherited symbionts.** This
  analysis does not infer keystone taxa from a co-occurrence network. With
  3–56 samples per genus and compositional, assembly-derived abundances, a
  network-centrality analysis would not be reliable.

## Sampling maps — step 10c

| output | |
|---|---|
| `results/maps/sampling_map.html` | interactive (Leaflet). Toggle tick genera and base layers; click a point for its metadata, dominant microbial genus and endosymbionts. Needs internet for the map tiles. |
| `16a_sampling_map.png` | all samples: colour = tick genus, shape = host type, size = samples at the site |
| `16b_sampling_map_by_genus.png` | one panel per tick genus |
| `16c_endosymbiont_map.png` | where each common endosymbiont was detected vs not |
| `results/maps/sampling_sites.tsv` | the mapped table |

Every sample with coordinates is mapped, including those never analysed (they
show faded in the HTML). Samples outside the static extent (the Falklands
*I. uriae*) are named in the caption. State borders are fetched from Natural
Earth once and cached in `results/maps/cache/`.

---

## Earlier findings (118 samples, provisional metadata)

> **Superseded.** These numbers predate the backfills (now 149 analysed
> samples) and the curated metadata. Current results, with every number
> computed from the outputs, are in **`report/tick_microbiome_report.html`**
> (source: `report/tick_microbiome_report.qmd`). This section is kept only as a
> record of the provisional analysis.

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

## Figure style and colour

Every figure is written as both `.png` and `.svg` with the same name.
`save_plot()` in `common.R` replaces `ggsave()`. `png_svg()` draws each
base-graphics figure (pheatmap, phylogenies) once per device. SVGs are
written with `svglite`, so their text stays editable. svglite lives in a
separate R library per env (`env/r_lib_svglite*`, see `env/README_svglite.md`),
not in the pinned conda envs. Without it, SVGs fall back to Cairo
(text as outlines).

Every ggplot uses `theme_tick()` from `code/analysis/common.R`, which is
`theme_classic()` plus light-grey major gridlines in both directions, set with
`theme_set()` when `common.R` is sourced. All colours come from the CSIRO brand
palette in `CSIRO_colours.R`, read into the named vector `CSIRO`:

- `CSIRO_CAT`: nine categorical hues, **chosen and ordered by computation**.
  They were simulated under protanopia and deuteranopia (Machado 2009) and scored
  by OKLab ΔE. Every pair stays ≥ 8.7 apart under both deficiencies and ≥ 15.8
  under normal vision, and each prefix is the best-separated set of that size.
  `csiro_cat(n)` returns the first n. A figure never uses more than nine: the
  rest fold into "other".
- `TICK_GENUS_PAL`: one fixed hue per tick genus, taken from `CSIRO_CAT`.
- `CSIRO_OTHER` (Mist) for other / unidentified / not detected. Steel was
  rejected because it collapses onto Mint and Fuchsia under deuteranopia.
- `CSIRO_SEQ` / `csiro_seq(n)`: a single-hue white → Midday → Ocean →
  Midnight Blue ramp for magnitudes.
- `MARK_OUTLINE` (Midnight Blue) outlines every filled mark. Gold and Light
  Mint have < 3:1 contrast on white, and brand values cannot be re-stepped.

Identity is never shown by colour alone: there is always a shape, a facet or a
direct label as well.

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
