# Marker-gene phylogenies for the target taxa

Takes the marker loci recovered by `run_markers.sh`, downloads the named-species
references each locus needs, cuts the locus out of them, aligns everything
together and builds a tree per locus.

```bash
bash run_trees.sh                                     # every locus
bash run_trees.sh --target Rickettsia --target Coxiella
bash run_trees.sh --target Coxiella --marker rpoB
```

Every locus comes out as four files in `results/markers/trees/<target>/`:

```
<target>__<marker>.aln.fna       the alignment the tree was built from
<target>__<marker>.nwk           FastTree output, UNROOTED
<target>__<marker>.rooted.nwk    rooted -- on the outgroup, else midpoint
<target>__<marker>.png           the figure
```

---

## Why this exists

`README_markers.md` ends with "the trees are yours to build -- the choice of
model, of outgroup and of which fragmentary sequences deserve to be in the
analysis is not something this pipeline should be making." That was the right
call for a hand-built final analysis and the wrong one for getting a first look
at 47 loci. Steps 17-19 make those choices explicitly, write them down, and
build the trees.

They are **screening trees**. They answer "is this *Coxiella* contig the
endosymbiont or *C. burnetii*?" well enough to decide which loci deserve a
proper analysis. They are not the figure for the paper. What you are meant to
do with them is read `tree_summary.tsv`, pick the loci that carry real signal,
and rebuild those under a model-tested partitioned analysis -- `iqtree` is in
the same environment and step 25 shows the pattern.

---

## Where the references come from

Two panels, both plain TSV, both meant to be edited.

**`metadata/tree_references.tsv`** -- one row per target x marker. An Entrez
nuccore query for the named species and published tick-derived sequences of
that locus, plus an optional outgroup query. 47 rows, one for every locus that
recovered a sample sequence.

**`metadata/tree_ref_genomes.tsv`** -- one row per target. Complete genomes.

The genome tier is there because several of the loci that actually resolve
these genera have almost no standalone gene records in GenBank -- *Rickettsia*
`sca4` and the 17 kDa antigen, *Coxiella* `icd`/`dnaK`/`IS1111`, *Francisella*
`sdhA`, everything in *Rickettsiella* and *Midichloria*. The named species that
carry them exist only inside complete genomes. Step 18 blasts the locus out of a
genome exactly as step 13 blasts it out of a contig, so **one genome download
per target supplies a reference for every locus of that target**. For the
piroplasms and *Hepatozoon* the "genome" is the mitochondrial genome, which is
where `cox1` and `cytb` live; for the viruses it is the complete genome or
segment record.

### Breadth, not depth

NCBI's relevance sort hands you forty near-identical sequences of the
best-studied species. A tree of forty *Rickettsia felis* `gltA` sequences has no
topology worth reading.

So candidates are screened on **esummary** -- accession, title, length -- before
anything is downloaded. The title carries the organism, so the quota is spent
round-robin across organisms, capped at `trees.curate.max_per_organism` records
each. A query matching 400 records costs one summary call, not 400 downloads,
and what comes back spans species.

### Nothing is trusted until it is blasted

A record matched by title need not carry the locus, and a 1.3 Mb genome is
99.9% not the gene. Step 18 blasts every downloaded record with this locus's
step-11 probes plus the longest few sample sequences, against the same
`min_aln_len` / `min_pident` floors the panel used on the contigs, and keeps
only the region that hit:

- HSPs on one record and one strand are merged into a locus; the best-scoring
  locus per record wins. A genome holds one `groEL` and twenty `IS1111` copies —
  in both cases the best hit is the answer, because a second copy of the same
  locus is not extra phylogenetic information.
- The locus is cut out, reverse-complemented if it was on the minus strand, and
  for a protein-coding marker trimmed to the forward frame with the fewest
  internal stops and to a whole number of codons — the same treatment step 14
  gives a sample sequence. That is what lets references and samples share one
  codon alignment.
- A reference whose locus still carries an internal stop is dropped. It is a
  pseudogene or a mis-annotation and it would distort the tree.

A loose query therefore costs download time, not correctness.
`tree_reference_sequences.tsv` lists every reference that survived, with the
identity and coverage that got it there, and its GenBank defline.

---

## Alignment, tree, root

**Alignment.** rRNA loci: MAFFT `--adjustdirection` on the nucleotides.
Protein-coding loci: MAFFT on the translations, back-translated onto the
in-frame nucleotides — a true codon alignment, gaps on codon boundaries.
Aligning these genes as raw nucleotide would put gaps mid-codon and hand the
substitution model an alignment it does not describe.

Column trimming is **off** by default. Most sample sequences are fragments of
the gene, so occupancy is low everywhere and a gap-threshold trim removes real
signal along with the noise. Set `trees.align.trim` to `gappyout` or `gt:0.05`
to enable trimal; the untrimmed alignment is kept alongside.

**Tree.** FastTree, GTR+CAT, with its SH-like local support at every node. Those
are **not bootstrap values** — they are a local test on each split and they are
optimistic. Treat them as "this split is not obviously an artefact", not as a
posterior.

**Root.** On the outgroup when the panel supplied one and step 18 verified it,
by midpoint otherwise. The figure says on its face which was used, because a
midpoint root is a display convention, not a phylogenetic claim. The outgroup is
the one thing not held to the panel's `min_pident` floor — it is the most
divergent sequence in the tree by construction — so check its realised identity
in `tree_reference_sequences.tsv` before leaning on the root.

---

## Reading the figure

| | |
|---|---|
| **red, bold** | sequence from one of these ticks — `sample_id [tick species]` |
| black, italic | GenBank reference — `organism accession` |
| blue, italic | outgroup |
| grey numbers | FastTree SH-like local support, drawn at or above `trees.figures.min_support_shown` |

The tick species in a red tip is a **guess** until `metadata/sample_metadata.tsv`
exists — see the naming section of `README_markers.md`. `unknown_species` means
the sample id carried no abbreviation the guesser recognised.

### A short sequence is not a species identification

This bears repeating here, because a tree makes a weak result look strong. Many
of these loci come off contigs of a few hundred bp; `tree_summary.tsv` carries
`median_ungapped_nt` and `pct_gaps` for exactly this judgement, and
`marker_summary.tsv` from step 16 carries `median_ref_cov_pct` — how much of the
gene was actually recovered. A 300 bp `gltA` fragment sitting inside a clade
tells you the genus and the clade. It does not tell you the species, however
confidently the branch is drawn.

---

## Outputs

```
results/markers/trees/
  refs/gene/<target>__<marker>.fna      raw downloads, gene-level records
  refs/genome/<target>.fna              raw downloads, whole genomes
  refs/outgroup/<target>__<marker>.fna  raw downloads, outgroup
  refs/sets/<target>__<marker>.refs.fna the locus cut out of them (+ .faa)
  hits/                                 raw blast HSPs against the references
  <target>/<target>__<marker>.aln.fna   alignment  <- sample + reference
            <target>__<marker>.nwk      unrooted FastTree
            <target>__<marker>.rooted.nwk
            <target>__<marker>.png      figure
            <target>__<marker>.aa.aln.faa   protein alignment, CDS loci
  tree_reference_manifest.tsv   what each query returned, per tier
  tree_reference_sequences.tsv  every reference in every tree, annotated
  tree_summary.tsv              per locus: tips, samples, refs, organisms,
                                alignment columns, gap %, whether a tree was built
```

Everything is cached. Re-running `run_trees.sh` re-fetches nothing; use
`--refresh` on step 17 to re-download, `--force` on step 18 to re-blast, and
`--force` on step 19 to realign and rebuild.

---

## Changing what goes in a tree

**Different references:** edit the query in `metadata/tree_references.tsv`, then

```bash
python3 code/analysis/17_fetch_tree_refs.py --target Coxiella --marker rpoB --refresh
python3 code/analysis/18_curate_tree_refs.py --target Coxiella --marker rpoB --force
python3 code/analysis/19_build_trees.py      --target Coxiella --marker rpoB --force
Rscript  code/analysis/19b_tree_figures.R    Coxiella
```

**Your own sequences, no Entrez at all:** append them to
`results/markers/trees/refs/sets/<target>__<marker>.refs.fna` (naming them
`REF|<id>|<organism>|<target>|<marker>|<len>bp`, or `OUTGROUP|...` to root on
them) and re-run steps 19 and 19b with `--force`. For a protein-coding locus
they must be in frame and a whole number of codons, and they must also go into
the matching `.refs.faa`, or step 19 will stop and tell you which sequence has
no translation.

**A locus the panel does not cover:** step 17 falls back to an auto-generated
query and marks it `auto` in the manifest, so a newly added marker still gets a
tree. Write a real row when you care about the result.

**Sequences with internal stop codons** (step 14's `*.flagged.fna`) are left out
of every tree. `--include-flagged` on step 19 puts them in — decide that per
locus, after looking at them.

---

## Setup

```bash
mamba env create -f env/tick_markers.yml   # mafft, fasttree, trimal, blast
mamba env create -f env/tick_micro.yml     # R, ape, phangorn (step 19b)
export NCBI_API_KEY=...                    # 3 req/s -> 10 req/s
```

`run_trees.sh` puts both environments where they are needed itself; you do not
need to activate anything. **Set `NCBI_API_KEY` if you have one** — step 17 makes
hundreds of Entrez requests and the genome tier downloads hundreds of megabytes.

### If step 17 stalls or reports SSL errors

This compute node reaches NCBI through something that intermittently presents
its own certificate — requests fail in bursts with

```
esearch.fcgi failed after 6 tries: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED]
certificate verify failed: self-signed certificate in certificate chain>
```

while the same request from a second process at the same moment succeeds. It is
the network, not the query. The client already backs off up to ~45 s per attempt
to ride an episode out, and a locus that still fails is **not cached**, so the
fix is simply to run step 17 again — everything already downloaded is skipped:

```bash
python3 code/analysis/17_fetch_tree_refs.py      # repeat until the tally stops changing
```

The last line of the run tells you how many loci still have no gene-level
references, and `tree_reference_manifest.tsv` names them.
