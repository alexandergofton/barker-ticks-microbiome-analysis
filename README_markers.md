# Phylogenetic marker genes for the target taxa

Finds, extracts, labels and aligns the marker genes that actually resolve each
target — 16S/18S plus the discriminatory protein-coding loci — from the contigs
already assigned to that target by `10_export_targets.py`. The output is a set
of alignments waiting for your reference sequences.

```bash
bash run_markers.sh                                    # all targets
bash run_markers.sh --target Rickettsia --target Coxiella
```

---

## Why this exists

`results/targets/` establishes that a genus is present and how much of it there
is. It cannot say **which species**, and for three of the targets that is the
entire question:

> *Coxiella*, *Rickettsia* and *Francisella* each contain obligate tick
> endosymbionts **and** genuine human pathogens that are near-identical over a
> short contig. A genus-level hit does not establish pathogen presence.
> — `metadata/target_taxa.tsv`

16S rRNA does not settle it. It separates *Rickettsia* from *Anaplasma* and
stops there — it barely distinguishes spotted-fever-group species from each
other, and Coxiella-like endosymbionts sit within a few bases of *C. burnetii*.
`gltA`, `ompA`, `ompB` and `sca4` do the work for *Rickettsia*; `rpoB` and
`icd` for *Coxiella*; `flaB` for *Borrelia*; `msp4` for *Anaplasma*; `cox1` and
`cytb` for the piroplasms. So the panel is **defined per target**, and each
locus is searched independently.

There is no guarantee any given locus is present. Whether `gltA` was assembled
in a given sample depends on coverage and on what SPAdes reconstructed. That is
the expected case, not a failure, and `marker_gaps.tsv` reports it explicitly.

---

## The panel

`metadata/target_markers.tsv` — one row per target × marker, 120 loci over 32
targets. **Edit this file to change the panel**; nothing is hardcoded in the
scripts. Add a row to add a locus, delete a row to drop one.

| target | loci |
|---|---|
| Rickettsia | 16S, gltA, ompA, ompB, sca4, htrA/17kDa, groEL |
| Anaplasma | 16S, gltA, groEL, msp4, msp2/p44 |
| Ehrlichia / Neoehrlichia | 16S, gltA, groEL, p28-omp1 |
| Borrelia / Borreliella | 16S, flaB, glpQ, clpA, ospA, 16S–23S IGS |
| Coxiella | 16S, rpoB, dnaK, groEL, icd, IS1111 |
| Francisella | 16S, sdhA, rpoB, pdpD, groEL |
| Bartonella | 16S, gltA, rpoB, ftsZ, groEL |
| Midichloria / Rickettsiella / Lariskella | 16S, gltA, groEL, ftsZ, rpoB |
| Wolbachia | 16S, ftsZ, coxA, fbpA, gatB, wsp, groEL |
| Spiroplasma | 16S, rpoB, dnaA, ftsZ |
| Arsenophonus / Cardinium | 16S, fbaA, yaeT, gyrB, groEL |
| Babesia / Theileria / Cytauxzoon | 18S, cox1, cytb, beta-tubulin, ITS |
| Hepatozoon / Hemolivia | 18S, cox1 |
| Trypanosoma | 18S, gGAPDH, cox1 |
| Nairoviridae / Phenuiviridae / Peribunyaviridae | RdRp (L), N, GPC |
| Flaviviridae | polyprotein, NS5 (RdRp), NS3 (helicase), VP1, VP2 |
| Chuviridae / Rhabdoviridae | RdRp (L), N, G |
| Orthomyxoviridae | PB1, PB2, NP |
| Reoviridae | RdRp (VP1) |
| Mycoplasma | 16S, rpoB, rnpB |

---

## The steps

| script | does |
|---|---|
| `11_fetch_marker_refs.py` | fetches reference **probes** per locus from NCBI Entrez; cached |
| `12_extract_target_contigs.py` | pulls target contigs out of the 70 GB of assemblies, relabelled |
| `13_find_markers.py` | blastn (rRNA) / tblastn (protein), assembles HSPs into loci, resolves paralogues |
| `14_extract_markers.py` | cuts the loci out, fixes the reading frame, names every sequence |
| `15_align_markers.py` | MAFFT per target × marker; codon-aware for protein-coding loci |
| `16_marker_report.py` | what was recovered, and what was not |

### The reference probes are not your tree references

Step 11 fetches a handful of sequences per locus **only so blastn and tblastn
can find and orient the gene**. Nothing downstream of step 14 uses them, and
none of them appear in an alignment. The reference set for your phylogeny —
type strains, named species, whatever the locus needs — is the thing you add by
hand, and it is deliberately left to you.

Cached under `results/markers/refs/`, with `reference_manifest.tsv` recording
which Entrez query produced each set. A locus whose query returns nothing is
listed in `unresolved_queries.tsv` rather than silently searched with no probe.

### Three things step 13 and 14 do that matter

**Paralogue resolution.** *Rickettsia* `ompA`, `ompB` and `sca4` are all
Sca-family autotransporters and cross-hit each other's probes; so do *Anaplasma*
`msp2`/`p44`, and a viral polyprotein contains its own NS5. Left alone, one
stretch of contig would be exported three times under three gene names and each
would get a tree. Overlapping loci of different markers on one contig are
resolved competitively — best single-probe bitscore takes the region — and the
loser is kept in `marker_loci.tsv` as `paralog_suppressed`, naming what beat it.

**Reading frame.** Each protein-coding locus is trimmed to the forward frame
with the fewest internal stops, so the exported nucleotide sequence is in frame
and a whole number of codons. This is what makes the codon alignment possible.
The frame is recovered from the sequence rather than taken from tblastn's
`sframe`, because a locus assembled from merged HSPs need not begin on a codon
boundary and a minus-strand frame is reported relative to the whole contig.

**Internal stop codons are flagged, not fixed.** A stop that survives frame
selection usually means an assembly frameshift or a pseudogene. Those sequences
go to `*.flagged.fna` and stay out of the alignments — you should see them
before deciding, not have them quietly patched.

### Probe with the relatives that are actually there

The viral rows are probed with the taxa observed in this cohort, not with the
type members of the family, because the first pass showed that matters. The
7,118 contigs on the `Flaviviridae` target are 7,108 contigs of **Cheeloo
Jingmen-like virus** — a *segmented* flavi-like virus that encodes no
polyprotein at all — so classical flavivirus polyprotein probes returned
literally zero HSPs from 6.2 Mb of contigs. Probed with Jingmenvirus segment
proteins instead, VP2 comes back at 83% identity. The `Nairoviridae` and
`Phenuiviridae` rows are likewise probed with Dali Nairo / Huangpi / Taggert and
Okutama / Dabieshan / Kaisodi rather than the CCHF and SFTS groups.

Two results survive that correction and are worth knowing before you read
`marker_gaps.tsv`:

- The **Jingmen-like virus signal is structural segments only.** VP1 and VP2 are
  recovered; NS5 and NS3 give no hit at all even at `-evalue 1` with correct
  Jingmen tick virus probes. The replicase segments were not assembled.
- The 92 **Nairoviridae** contigs yield nothing from any of L, N or GPC, best
  blastx hit *e* = 0.43 against probes from the very viruses their LCA names.
  Whatever those contigs are, they do not encode a detectable nairovirus
  protein — treat the family-level LCA calls on them with suspicion.

---

## Sequence naming

Every sequence, from the contig pull to the final alignment, is named:

```
<sample_id>|<tick_species>|<target>|<marker>|<contig_id>|<len>bp[|copyN]

2331-Apersicus|Argas_persicus|Anaplasma|groEL|NODE_5734_length_1727_cov_2.5|1359bp
```

so a tip label states the sample, the tick it came from and which gene it is
with no lookup. `copyN` appears only when one sample yields the locus on more
than one contig, numbered longest first — `copy1` is always the best sequence
for that sample.

**The tick species is a guess until you say otherwise.** Until
`metadata/sample_metadata.tsv` exists, the name comes from
`tick_species_guess` in the metadata template, derived from the sample-id
abbreviation. `tick_species_source` travels in every output table so any tip
label can be traced back to how its species was assigned.

On the current metadata, **59 of the 141 samples carry no species at all** —
their sample ids hold no abbreviation the guesser recognises — and are labelled
`unknown_species` in every tip. Step 12 prints the list. Fill in
`metadata/sample_metadata.tsv` and re-run steps 12–15 to relabel; nothing else
has to be redone, because the reference probes and the BLAST outputs are cached.

---

## Scope: every target contig, not just the confident ones

Extraction runs over **all** of `target_contig_index.tsv`, including samples
below the presence call in `target_prevalence.tsv`. A recovered marker gene is
itself detection evidence — filtering on read counts first would discard real
detections. The step-10 call travels with every sequence as annotation
(`target_present_in_sample`, `target_evidence`), so you can filter afterwards on
whatever the locus in front of you justifies.

---

## Outputs

```
results/markers/
  refs/          reference probes + reference_manifest.tsv, unresolved_queries.tsv
  contigs/       <target>.contigs.fa, contig_metadata.tsv
  hits/          raw blastn/tblastn HSPs, kept for audit
  marker_loci.tsv        every candidate locus, incl. paralog_suppressed
  sequences/<target>/    <target>__<marker>.fna / .faa
                         <target>__<marker>.flagged.fna  (internal stops)
  alignments/<target>/   <target>__<marker>.aln.fasta      <- nucleotide, codon-aware
                         <target>__<marker>.aa.aln.fasta   <- protein
  marker_sequences.tsv   one row per extracted sequence, fully annotated
  marker_summary.tsv     per target x marker: n, lengths, identity, alignment
  marker_by_sample.tsv   per sample x target x marker
  marker_gaps.tsv        loci that produced nothing, and why
```

Read `marker_gaps.tsv` before concluding a locus is absent from a taxon: it
separates *the target had no contigs in any sample* from *the target had contigs
and this locus was not among them*. Only the second is a biological statement.

---

## Adding your references and building trees

The alignments hold sample sequences only. To add references without
re-aligning everything:

```bash
conda activate tick_markers

mafft --add my_Rickettsia_gltA_refs.fasta --keeplength --auto \
      results/markers/alignments/Rickettsia/Rickettsia__gltA.aln.fasta \
      > Rickettsia_gltA_with_refs.aln.fasta
```

Drop `--keeplength` if your references should be allowed to open new columns —
usually the right choice when they are full-length genes and the sample
sequences are fragments.

For a protein-coding locus, keep working in codon space: align references at the
amino-acid level against `*.aa.aln.fasta`, then back-translate, or add them to
the codon alignment directly if they are already in frame.

`iqtree`, `fasttree` and `trimal` are in the `tick_markers` env if you want
them, but the trees are yours to build — the choice of model, of outgroup and
of which fragmentary sequences deserve to be in the analysis is not something
this pipeline should be making.

### A short sequence is not a species identification

Many of these loci come off contigs of a few hundred bp. A 300 bp `gltA`
fragment places a sequence in a genus and rarely more. `marker_summary.tsv`
carries `median_len_nt` and `median_ref_cov_pct` — how much of the gene was
actually recovered — for exactly this judgement.

---

## Setup

```bash
mamba env create -f env/tick_markers.yml   # mafft, muscle, trimal, iqtree, fasttree
export NCBI_API_KEY=...                    # optional; 3 req/s -> 10 req/s
```

`run_markers.sh` puts the env on `PATH` itself; you do not need to activate it.
blast+, seqkit and samtools come from the base miniforge3 install.

Re-running is cheap: fetched references and BLAST outputs are cached. Use
`--force` on step 13 to redo the searches, `--refresh` on step 11 to refetch.
