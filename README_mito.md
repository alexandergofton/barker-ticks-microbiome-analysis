# Tick mitochondrial genome phylogenomics

Extracts one mitochondrial genome candidate per sample from the assemblies,
annotates its genes, aligns each protein-coding gene separately, and builds one
partitioned maximum-likelihood tree in IQ-TREE with a per-gene substitution
model -- a mitogenome-based phylogeny of the ticks themselves, independent of
the pathogen/endosymbiont work in `README_markers.md`.

```bash
bash run_mito.sh
```

---

## The steps

| script | does |
|---|---|
| `20_fetch_mito_refs.py` | reference probes for 13 PCGs + 2 rRNAs, from NCBI Entrez; cached |
| `20b_fetch_mito_ref_genomes.py` | downloads the published reference mitogenomes in `metadata/mito_ref_genomes.tsv` and annotates them exactly like the samples |
| `21_find_mito_candidates.py` | picks ONE mitogenome contig per sample by homology |
| `22_annotate_mitogenomes.py` | final per-gene sequences; MITOS2 if available, else homology fallback |
| `23_align_mito_genes.py` | MAFFT per gene, codon-aware for the 13 PCGs |
| `24_concat_partition.py` | concatenates the PCGs into one partitioned supermatrix |
| `25_run_iqtree.py` | ModelFinder + partitioned ML tree, ultrafast bootstrap + SH-aLRT |
| `26_mito_figures.R` | tree figure + two distance tables |
| `26b_mito_vs_metadata.R` | curated IDs vs the tree; nearest named reference per sample |
| `26c_mito_discordance_audit.py` | why a sample is discordant: contaminant mitogenome / mixed tube / wrong tick, from nuclear-rRNA BLAST hits |

The gene panel is `metadata/mito_genes.tsv`; settings are the `mito:` block of
`config/analysis_params.yaml`.

---

## 0. Published reference mitogenomes (step 20b)

The tree is built from the samples **plus 97 published tick mitogenomes** of
named species, so every sample tip sits among species of known identity and
the tree can be rooted between the tick families.

The panel is `metadata/mito_ref_genomes.tsv`, curated from the raw list in
`results/mito/mito_ref_seqs.txt` (121 lines). Every line is accounted for:
one row per distinct accession, `use=no` plus a reason for anything not
downloaded. What curation found in the raw list:

- **`NC_028521` (listed as *Dermacentor nuttalli*) is a frog**, *Rana
  huanrensis*. Replaced with `NC_028528`, the *D. nuttalli* RefSeq mitogenome
  (a one-digit typo).
- **`MT371800` (listed as *D. taiwanensis*) is *D. bellulus* in GenBank.**
  Used and labelled as listed, flagged in the table for a check.
- Truncated/garbled accessions resolved to the only matching record:
  `MT37180` → `MT371802` (*H. cornigera*), `MT3717` → `MT371797`
  (*A. geoemydae*), `OL74173` → `OL741735` (*A. latum*), `KC50326:` →
  `KC503261` (*R. microplus*), `NC_ 017757` → `NC_017757`,
  `CC 0326552` → `NC_026552` (*D. silvarum*). `KY45753` (*R. evertsi*) matches
  two clones from one study, `KY457537` and `KY457538`; `KY457537` is used.
- Duplicates dropped: `NC_018369` (= `JN248424`), `NC_088057` (= `OR416215`),
  `NC_023335` (= `KC503261`), plus ~15 accessions listed twice.
- **Two SRA read runs** (`SRR4301110` *A. aureolatum*, `SRR8074777`
  *A.* aff. *moreliae*) are not assembled genomes and are not used; each would
  need its mitogenome assembled from reads first.

Each record's GenBank title must contain the expected species epithet (or the
row's `genbank_check`), or step 20b refuses it -- the check that would have
caught the frog automatically.

References are annotated **by the same homology method as the samples**, not
from their GenBank features: step 21's own tblastn/blastn and locus-merging
functions, step 22's framing and stop-trim (`markers_lib.frame_mito_pcg`).
Many listed records (the Kelava et al. `MT3718xx` series) are GenBank
"UNVERIFIED" sequences with no CDS features at all, and identical boundary
rules keep the codon alignments free of an artificial references-vs-samples
edge pattern. A reference gene that fails the stop-codon check is left out
rather than flagged (on a published complete genome it means a mis-called
locus, not a pseudogene): ND1 in 6 references, ATP6 in 1.

Recovery on these complete genomes is itself informative: **ATP8 in 2/97 and
ND6 in 51/97**, the same two genes the samples lack. Those gaps are therefore a
limit of homology search on short/divergent genes, not missing sequence in
the samples.

Set `mito.ref_genomes.include_in_tree: false` to rebuild the samples-only tree.
The previous samples-only tree and its step-26b tables are kept in
`results/mito/samples_only/` for comparison.

---

## 1. Finding the mitogenome: one contig per sample, by homology

Every sample's contigs in a plausible mitogenome length window
(`mito.candidates.min/max_contig_len`, default 3,000-22,000 bp) are searched
with **tblastn** against 13 mitochondrial protein-coding gene probes and
**blastn** against 2 rRNA probes, all fetched at **order level** (Ixodida = every
tick). That matters: the sample set spans Ixodidae (hard ticks) and Argasidae
(soft ticks), families separated by upwards of 200 million years. Nucleotide
identity between them is often too low for blastn to find the gene at all --
tblastn on the translated, conserved protein is what actually carries the
search across that gap. A first pass at the marker workflow's within-genus
`min_pident` floors (40-45%) rejected real hits at this divergence -- a COX1
match at 41.5% identity over 508/512 aa, *E* = 5×10⁻¹⁰⁷, is not noise, it is
one of the strongest signals BLAST can report. The panel's `min_pident` was
lowered accordingly (30% protein / 65% rRNA); the *E*-value cutoffs
(`evalue_pcg`/`evalue_rrna`) do the real significance gating.

**Only tick contigs are eligible.** Each candidate contig's own nt-BLAST LCA
from the main pipeline (`lca/<sample>.blast.lca`, lineage from step 03) must
fall inside Ixodida (`mito.candidates.lca_filter`). The probes are tick
proteins, but a vertebrate mitogenome carries the same 13 genes and clears the
E-value gate, so before this gate blood-meal mitogenomes (*Varanus*,
*Petrogale*, *Isoodon*) and bacterial contigs with 2-3 weak hits won in 7
samples. A contig with no nt-BLAST hit at all stays eligible but ranks below
confirmed tick contigs: the complete *Ornithodoros* mitogenomes here
(`Oguryn1/2`, `B3648_M1_14_8_S19`) have no nt hit.

Among eligible contigs the winner is the one with the **most distinct genes**
recovered, ties broken by **coverage**, then summed bitscore, then length --
not simply "closest to 15 kb", so a contig with strong evidence for 12 genes at
13 kb beats a 14.8 kb contig that only hits one gene by chance. Coverage comes
before bitscore because in a mixed tube bitscore only says which tick is closer
to the probe set; `8164-M1-Amorelia` had picked a 5x *Bothriocroton* mitogenome
over its own 10x *Amblyomma* one on bitscore.

Unlike the marker workflow, overlapping HSPs from two *different* genes are
never treated as competing paralog calls here, and deliberately so. An earlier
version of this step ported the marker workflow's highest-bitscore-wins overlap
resolution over unchanged, and it did real damage: every one of the 146
"conflicts" it found on this cohort was two adjacent, unrelated mitochondrial
genes (mostly ND1 vs. its neighbour rrnL) whose independently blasted,
HSP-approximate boundaries happened to overlap by a few dozen bp -- not a
competing annotation of the same DNA, which is impossible between a
protein-coding gene and an rRNA. It was silently deleting a correct, strong
ND1 call (84% identity, *E* = 1.7×10⁻¹²³) in most samples that had one. Fixed by
removing cross-gene overlap resolution entirely: mitochondrial PCGs are
single-copy and non-homologous to each other, so there is no real paralog
conflict to resolve in the first place. The only same-gene safety check kept
is `cap_one_locus_per_gene`, for the genuinely possible case of one gene
producing two disjoint loci on the same contig.

### Result on this cohort

158 samples had assemblies. **140 yielded a mitogenome candidate** (3 had no
contig in the length window at all, 6 had contigs but no recognisable
mitochondrial gene on any of them, and 9 had genes only on non-tick contigs --
status `no_tick_contig`, with the rejected contigs and their LCAs in
`non_tick_excluded`); results/mito/mito_selection.tsv has the
per-sample outcome, including `n_candidate_contigs` and `n_genes_found` so a
"no mitogenome" result is never just a bare negative. Recovery is highly
sample-dependent -- a full 13-15 kb genome on one long, low-coverage contig
versus a handful of genes on a short fragment in a shallow assembly -- and that
variation is expected, not a bug.

Gene recovery across the 140 (of 15 possible genes; see `mito_selection.tsv`):
COX3/ATP6 highest (133 each), **ATP8 lowest by far (3/140)** -- it is the
shortest mitochondrial protein in animals (~50 aa) and routinely falls below
any BLAST significance floor regardless of threshold tuning, a known,
unavoidable limitation of homology search on this gene, not specific to this
pipeline. **ND6 is the next lowest (61/140)**, consistent with it being the
most divergent, hardest-to-detect NADH subunit across metazoans generally.

**137 of the 140 made it into the final partitioned tree**
(`tree.min_genes_for_tree=3` excludes 3 -- see `concat/excluded_samples.tsv`),
alongside the 97 published references: 234 taxa x 13 partitions, with
ModelFinder selecting a different model per partition
(`tree/mito.best_scheme.nex`). IQ-TREE took 2 h 56 min on 4 cores.

**The references validate the tree at depth**: Ixodidae is monophyletic,
*Nuttalliella* is sister to Ixodidae, and every reference genus is
monophyletic except *Ornithodoros* s.l. and *Haemaphysalis* (because of
*Alloceraea*), both known results.

**No selected contig is non-tick any more** (step 26b `NOT_TICK_CONTIG`: 0;
it was 7 before the LCA gate). Of those 7: 4904, 8177, 8267 and 5070 have no
tick mitogenome contig at all; 8141 and 8405 have only a 1-gene fragment
(excluded from the tree); `B4276_F1_28_7_S20` now carries a 9-gene
*Haemaphysalis* contig. `8186-M1-F1-Afimb` (LCA "Bilateria", a long branch)
and `B7115_N3_S40` (no nt hit, a long branch) also moved to tick contigs of
their curated genus. The pre-fix tree and tables are kept in
`results/mito/pre_lca_filter/`.

### Genus-discordant samples (step 26c)

**18 samples are genus-discordant** in the rebuilt tree (step 26b). The
mitogenome alone cannot say whether that is a contaminant, a mixed tube or the
wrong tick, so step 26c asks the sample's raw nt BLAST which tick genus its
**nuclear rRNA** (18S/28S/5.8S/ITS) contigs belong to.

The obvious version of this question -- what percentage of all tick hits go
to each genus -- does not work. nt holds whole nuclear genomes for only a few
tick genera, so in a correctly labelled *Amblyomma* or *Bothriocroton* sample
~70% of tick HSPs go to *Dermacentor* and ~1% to the curated genus. Only hits at
>= 95% identity are used, split by subject type. Tick nuclear rRNA and
mitochondrial genes are in nt for nearly every genus, and at that identity a
hit names the genus actually present. Shares are weighted by contig length x
coverage. `tables/tick_genus_signal.tsv` has every sample's breakdown.

`tables/mito_discordance_audit.tsv` (31 samples: the 26b-discordant ones, the
`no_tick_contig` ones, and any sample whose nuclear rRNA disagrees with its
curated genus):

| call | n | samples | meaning / action |
|---|---|---|---|
| `WRONG_TICK` | 12 | 0895-N1-7-Afimb, 8140-B-M1-Afimb, 8184-F1-N1-3-Alimb, 8188-M2-Afimb, 8316-B-F1-AmblyomaSP, 8336-Argasidae-4L-Boobie, B7024_M1_21_2_S11, B3015_F1_0_49_S28, B4560b_F1_S5, B4283_M1_S2, B5053a_N1_1_14_S10, B7002_1_F1_5_22_S9 | curated genus is 0% of both nuclear rRNA and mito hits: swap, mislabel or misidentification. **Remove** from species-keyed analyses and the mito tree |
| `MIXED_SAMPLE` | 4 | 8225-3L1, 8260-F1-Iconf, 5070-847-DECB-ParkhurstHbanc, B1201_M1_1_96_S24 | two genera with real nuclear-rRNA share; which tick is "the sample" cannot be decided. **Remove** |
| `UNRESOLVED` | 2 | B5056b_F1_S7, 8399-F-Benaraby-M-gig-14-4-24 | no raw BLAST file to test. **Remove** until checked |
| `CONTAMINANT_MITOGENOME` | 3 | 5121-M1-Alimb, 8187-F1-Afimb, 8261-F1-Atrig | nuclear rRNA 92-97% curated genus, mitogenome a foreign tick's. Tick ID stands; **drop from the mito tree**; microbiome may carry the source sample's microbes |
| `ID_SUPPORTED[_NO_MITOGENOME]` | 7 | 0143, 4904, 8130, 8177, 8256, 3263, 6587 | ID supported, no tick mitogenome contig |
| `UNTESTABLE_GENUS` | 2 | 8202-F1-Cpapua, B7103.F1.Cryptocroton | no *Cryptocroton* rRNA in nt (its hits land on *Amblyomma*) |
| `NO_CURATED_ID` | 1 | 8267 | nuclear rRNA is 97% *Haemaphysalis* -- usable to fill in its metadata |

What the evidence says about individual cases:

- **Same-lane identity.** A foreign mitogenome 0.00% from another sample on
  the same lane points at that sample as the source. `8187-F1-Afimb` carries
  `6700-10-CT-Hbanc`'s *Haemaphysalis* mitogenome, at 7x vs that sample's 102x,
  consistent with carry-over into a low-input library. `8184-F1-N1-3-Alimb` is
  identical to `8302-Argusrobertsi`, but at twice its coverage with 97% *Argas*
  rRNA: the tube holds an *Argas*, not a trace of one. In the pre-fix tree the
  *Varanus* blood-meal mitogenomes of `8177` and `5070-847-DECB` were identical
  too, and 32% of 5070's nuclear rRNA is *Amblyomma*, so 8177's DNA is in 5070.
- **`8225-3L1`** is 3 pooled larvae from a long-nosed bandicoot at PHRF, the
  same host and site as `8278-F1-Itasm`, whose *Ixodes* mitogenome it carries
  exactly. Most likely a real mixed pool (one *I. tasmani* larva), not lab
  contamination.
- **Mixtures where the curated genus is also present** (5121, 8187, 8261,
  8225): the foreign mitogenome has the *higher* coverage in every case, so the
  coverage tie-break cannot rescue them. It did rescue 8164 and 8140-A.
- **`8336-Argasidae-4L-Boobie`** comes out as *Alectorobius capensis*, a seabird
  soft tick, from a booby: likely a genus-level relabel within Argasidae, not an
  error. **`8316-B-F1-AmblyomaSP`** reads as *Dermacentor*/*Rhipicentor* and is
  >= 23% from every other sample: likely a misidentified tick from outside this
  cohort's fauna.

26c only flags. Nothing is removed automatically; the tree stays a record of
what the pipeline selected.

Two real bugs were caught and fixed while building this step, both from
porting logic out of the pathogen marker workflow without re-examining the
assumption behind it -- worth knowing about if you extend this workflow
further:

- **Cross-gene overlap "paralog" resolution doesn't belong here.** The marker
  workflow suppresses whichever of two competing gene calls has the lower
  bitscore when they overlap heavily, because Rickettsia's ompA/ompB/sca4
  probes genuinely cross-hit the same Sca-family DNA. Applied unchanged to the
  mitogenome, it suppressed a correct, strong ND1 call (84% identity,
  *E*=1.7×10⁻¹²³) in 87 of 149 samples, every time because ND1's HSP-approximate
  boundary happened to overlap its neighbour rrnL by a few dozen bp -- an
  artifact of independent probe searches, not a competing annotation, and
  something that cannot happen between a protein-coding gene and an rRNA in the
  first place. Removed entirely for this workflow; single-copy mitochondrial
  genes have no real paralogs to resolve.
- **A full annotation run must delete stale per-gene output files it no longer
  writes to.** When the stop-codon-trim fix below was added, some sequences
  moved from `GENE.flagged.fna` to `GENE.fna` between runs -- but the old
  `GENE.flagged.fna` was never touched by a run that produced zero flagged
  calls for that gene, so it kept reporting the previous run's (wrong) flagged
  set. Step 22 now clears every existing `genes/*.fna`/`*.faa` at the start of
  a full run before rebuilding them.

---

## 2. Annotation: MITOS2 vs the homology fallback

**MITOS2 (https://gitlab.com/Bernt/MITOS) is the intended tool and is
installed** (`mamba create -n mitos2 -c conda-forge -c bioconda "mitos>=2"`,
env `mitos2`). It could not be **run** here: its reference database is
distributed only from Zenodo (`zenodo.org/record/4284483`), and every attempt
to reach it from this compute node failed --

```
$ curl -sS --max-time 40 "https://zenodo.org/api/records/4284483"
# 504 Gateway Timeout after 31s, or no response at all following the doi.org redirect
```

(NCBI E-utilities, GitHub and GitLab raw content, and NCBI's FTP server were all
reachable in the same session, so this looks specific to zenodo.org from this
network, not a general egress block.)

**If you can reach Zenodo from another machine:** download the MITOS2
refdata (link above), copy it somewhere this project can see, and set

```yaml
# config/analysis_params.yaml, mito.annotate block
mitos_refdata_dir: "/path/to/refdata/refseq89_metazoa"
```

Re-run `code/analysis/22_annotate_mitogenomes.py` -- it detects the directory
and the `runmitos` command in the `mitos2` env automatically and switches to
real MITOS2 for every sample, no other change needed. **This path is
implemented but has never been exercised end-to-end** (there is no way to test
it without the reference data); if your installed MITOS2 version's
`result.bed` columns differ from what `run_mitos()` expects, fix the parsing
in `22_annotate_mitogenomes.py` -- the invocation itself follows MITOS2's own
`--help`.

**What actually ran, and what it does and does not give you.** Absent the
refdata, step 22 annotates by homology: the same tblastn/blastn gene calls from
step 21, re-oriented and placed into the correct reading frame using the
**invertebrate mitochondrial genetic code** (NCBI translation table 5 --
`markers_lib.INVERT_MITO_TABLE`). This code differs from the standard one at
four codons (ATA=Met, TGA=Trp, AGA/AGG=Ser); translating a tick mitochondrial
gene with the standard code misreads every in-frame TGA -- a normal Trp codon
here -- as a premature stop, which would flag most real genes as broken. Getting
this right is why `markers_lib.translate`/`best_frame` now take a `table`
argument instead of assuming one genetic code project-wide.

This gives you, per sample: the 13 PCG + 2 rRNA sequences, in frame, with
internal stop codons flagged (`*.flagged.fna`, held out of the alignments) --
enough to build the phylogeny in steps 23-25. It does **not** give you: tRNA
annotation, gene order, a circularity call, or codon-exact start/stop
boundaries -- gene edges are the homology-supported span from step 21, not a
verified start/stop codon. If you need those, MITOS2 (or tRNAscan-SE/ARWEN for
tRNAs specifically) is where to get them once the refdata is in hand.

That last point has one specific, common consequence worth knowing: because a
gene's boundary is a BLAST hit span rather than a verified stop codon, it
routinely runs one codon past the gene's own true, natural stop. The first
version of this logic flagged every such case as broken -- 63 of 149 samples'
worth, mostly full-length, correct genes with one extra spurious codon after a
real in-frame stop. Step 22 now trims there instead: a single stop codon in the
last 15% of the ORF is treated as the gene's own end and cut off (not a
frameshift), and only a sequence with multiple stops, or one stop well before
the expected end, is actually flagged. That distinction reduced the flagged
count on this cohort from 63 to **7** -- 6 ND1, 1 COX3 -- which are worth a
manual look as likely genuine frameshifts, pseudogenes, or mis-assemblies
rather than boundary artifacts.

---

## 3. Alignment and the partitioned tree

Each of the 13 protein-coding genes is aligned **separately**: MAFFT on the
translation, back-translated onto the in-frame nucleotides (the same
codon-aware machinery as the marker-gene workflow, `markers_lib.align` /
`backtranslate`), so every alignment column-triplet is one codon and no gap
falls mid-codon. The 2 rRNAs are aligned directly as nucleotide and kept
alongside for anyone who wants them, but are **not part of the tree** -- only
the 13 PCGs are concatenated, matching what was asked for.

The supermatrix (`results/mito/concat/mito_concat.fasta`) is genes side by
side, in a fixed order recorded in `mito_partitions.txt`. A sample missing a
gene gets an all-gap block of that gene's aligned length, not a dropped row --
dropping it would shrink every OTHER gene's taxon set to whatever that one
sample happened to have. A sample recovering fewer than
`tree.min_genes_for_tree` genes overall (default 3 of 13) is excluded from the
tree entirely instead (`results/mito/concat/excluded_samples.tsv`) -- a taxon
that is almost pure missing data is exactly the profile that produces a
rogue/long-branch artifact in ML tree search, not a taxon ML can place
reliably.

IQ-TREE runs **edge-linked** partitioning (`-p`): every partition shares the
same topology and branch lengths -- correct here because all 13 genes are
inherited together on one non-recombining mitochondrial genome and therefore
share one genealogy -- while ModelFinder (`-m MFP`) selects a substitution
model **independently per partition**. Branch support is ultrafast bootstrap +
SH-aLRT together. IQ-TREE's tree is unrooted; the R steps root it on the
**Argasidae reference mitogenomes** (`mito.tree.outgroup_family`), i.e. on the
Argasidae/Ixodidae split, the deepest well-established node among living
ticks. Only reference tips define the outgroup, so a mislabelled sample cannot
break the rooting; if those tips are not monophyletic the figure falls back to
midpoint rooting and says so. *Nuttalliella namaqua* is deliberately not part
of the outgroup: its position relative to the two families is the contested
question, so it is left for the data to place.

---

## 4. Figure and distance tables (R)

`26_mito_figures.R` (the `tick_micro` env; adds `r-ape` + `r-phangorn`) writes:

- `results/mito/figures/mito_tree.png` / `.pdf` -- the rooted tree, sample
  tips coloured by tick genus, references in black bold, SH-aLRT/UFBoot
  support at internal nodes.
- `results/mito/tables/mito_distance_pdistance.tsv` -- pairwise uncorrected
  p-distance (%) from the concatenated alignment, computed with **pairwise
  deletion**. This is not a cosmetic choice: with missing genes represented as
  all-gap blocks, comparing two samples over the full matrix without pairwise
  deletion would count every gene one of them lacks as a difference, inflating
  distances between exactly the pairs with the least data in common.
- `results/mito/tables/mito_distance_patristic.tsv` -- tree-based (patristic)
  distance between every pair of tips, from the fitted branch lengths.

---

## Sequence naming

Every gene sequence, contig and figure tip is named the same way as the marker
workflow (`markers_lib.marker_name`, target fixed to `mitochondrion`):

```
<sample_id>|<tick_species>|mitochondrion|<gene>|<contig_id>|<len>bp

2331-Apersicus|Argas_persicus|mitochondrion|COX1|NODE_908_length_14562_cov_12.2|1536bp
```

Tick species is the same guess-until-corrected value used throughout this
project (`tick_species_guess` until `metadata/sample_metadata.tsv` is filled
in) -- see `README_markers.md` for the caveat, which applies identically here.

---

## Outputs

```
results/mito/
  mito_ref_seqs.txt       the raw reference list (input; curated into metadata/mito_ref_genomes.tsv)
  refs/                   reference probes + reference_manifest.tsv
  ref_genomes/            step 20b: fasta/<acc>.fna, hits/, genes/<GENE>.fna/.faa,
                          ref_loci.tsv, ref_genome_summary.tsv
  samples_only/           the previous samples-only tree + its 26b tables
  candidates/              filtered_contigs/, hits/, blastdb/  (per-sample, step 21 audit trail)
  genome/<sample>.mito.fna  the ONE selected contig per sample
  mito_selection.tsv       one row per sample: winning contig + evidence + status
  mito_loci.tsv            every candidate locus, all contigs, incl. suppressed
  genes/<GENE>.fna/.faa    one FASTA per gene, all samples (+ *.flagged.fna)
  gff/<sample>.gff3         per-sample annotation (homology or MITOS2)
  annotation_summary.tsv    one row per sample: genes annotated, source, flags
  alignments/<GENE>.aln.fasta / .aa.aln.fasta  (+ <GENE>.input.* = samples + refs as aligned)
  concat/mito_concat.fasta, mito_partitions.txt, mito_partition_summary.tsv,
         excluded_samples.tsv
  tree/mito.treefile, mito.iqtree (per-partition models), mito.contree, ...
  figures/mito_tree.png, .pdf
  tables/mito_distance_pdistance.tsv, mito_distance_patristic.tsv
  tables/mito_vs_metadata.tsv, mito_species_monophyly.tsv   (step 26b)
  tables/tick_genus_signal.tsv, mito_discordance_audit.tsv  (step 26c)
  audit/tick_hits/<sample>.tsv  step 26c cache: best >=95% tick hit per contig x subject type
  pre_lca_filter/           tree + tables from before the step-21 LCA gate
```

## Setup

```bash
# alignment + tree (shared with the marker workflow)
mamba env create -f env/tick_markers.yml

# figures + distances
mamba install -n tick_micro -c conda-forge r-ape r-phangorn   # or: mamba env update -f env/tick_micro.yml

# MITOS2 -- installed, but see section 2 above before expecting it to run
mamba create --strict-channel-priority -c conda-forge -c bioconda -n mitos2 "mitos>=2"
```

Re-running is cheap: reference probes, BLAST hits, and length-filtered contig
sets are all cached per sample/gene. Use `--force` on step 21 or `--refresh` on
step 20 to redo them; `--sample` / `--gene` scope any step to a subset without
touching the other samples' rows in the shared tables.
