# Blood-meal host identification

Which vertebrate(s) each tick fed on, read from the host DNA in its own
assembly, compared with the host recorded in the metadata.

```bash
bash run_bloodmeal.sh
```

Runs independently of `run_all.sh` / `run_mito.sh`. Settings: the `bloodmeal:`
block of `config/analysis_params.yaml`.

---

## The steps

| script | does |
|---|---|
| `30_scan_vertebrate_hits.py` | one pass over each raw nt BLAST (`blast/`, ~480 GB) keeping vertebrate HSPs on mitochondrial / rRNA subjects; cached per sample |
| `31_fetch_vert_mito_refs.py` | RefSeq vertebrate mitogenomes (7,793): FASTA, per-gene coordinates, BLAST db |
| `32_resolve_metadata_hosts.py` | metadata host → NCBI taxon (+ `metadata/host_name_fixes.tsv`); how many nt / RefSeq mito records exist for each host |
| `33_classify_vertebrate_contigs.py` | per contig: mtDNA vs 18S, competitive check, LCA, identity cap, tier; RefSeq gene span |
| `33b_fetch_host_family_refs.py` | all nt mito records of each candidate host family (Entrez); cached |
| `33c_refine_host_contigs.py` | re-names each host mito contig by a like-for-like identity comparison |
| `34_call_bloodmeal.py` | per-sample tiered call, multiple hosts, cross-talk flag, comparison with metadata |
| `35_bloodmeal_figures.R` | concordance, class matrix, per-sample host figure |

---

## 1. Evidence tiers

Per sample, evidence is used hierarchically:

1. **Tier 1: mitochondrial genome / segment.** Contigs whose best vertebrate
   mtDNA HSP is ≥ `tier1_min_aln` (150 bp). These range from 150 bp fragments
   to complete 17 kb host mitogenomes.
2. **Tier 2, only when tier 1 is empty: 18S rRNA and COI.** Vertebrate 18S
   contigs, and COI fragments shorter than the tier-1 floor
   (≥ `tier2_min_aln`, 50 bp). COI is itself mitochondrial, so any COI contig
   ≥ 150 bp is already tier 1; tier 2 can only add the shorter fragments.
   18S separates vertebrate classes, and orders only at ≥ 99.5 % identity, so a
   tier-2 call is a class-level (sometimes order-level) answer by construction.

## 2. Is a contig vertebrate at all?

The taxonomy comes from the samples' **existing nt BLAST**, which is
competitive across all life. The node this was built on has no nt database
mounted and two CPUs, so nothing is re-BLASTed against nt.

- **Strict mtDNA subjects.** nt titles mentioning "mitochondrial" include
  nuclear genes for mitochondrial proteins ("translocase of inner
  mitochondrial membrane 17A", "lon peptidase 1, mitochondrial"), which turn up
  in nearly every sample. `bloodmeal_lib.locus_class` drops RefSeq
  transcripts (XM_/NM_/XR_), mRNAs, PREDICTED and "-like" records, and names
  the mitochondrial gene from the title (COX1, CYTB, 12S, 16S, ND, ATP, CR,
  complete genome).
- **Competitive check.** The best vertebrate HSP must sit inside the contig's
  own LCA window: bitscore ≥ 90 % of the contig's top nt bitscore over all
  taxa (`lca_run.top_percent`). A tick mitogenome contig that also hits a
  vertebrate COX1 at 75 % fails (`status = not_best_hit`).
- **Identity cap.** The contig is named at the LCA of the vertebrate subjects
  in that window, then truncated to the deepest rank the best identity
  supports: species ≥ 98 %, genus ≥ 95 %, family ≥ 90 %, order ≥ 85 %,
  otherwise class. nt holds a mitogenome for few Australian squamates, so a
  91 % hit to the only *Varanus* with one is read as Varanidae, not as that
  species.

## 3. Re-naming host mito contigs (steps 33b/33c)

The nt BLAST kept only 10 subjects per contig (`-max_target_seqs 10`). On a
long host contig, those 10 are complete mitogenomes of related species: a
16 kb alignment at 95 % outscores a 1 kb CYTB record of the right species at
99 %. Example: 8398 (recorded host *Notamacropus rufogriseus*) assembled a
complete 17 kb host mitogenome, but nt has no complete *N. rufogriseus*
mitogenome, only gene records, so all 10 subjects were other *Notamacropus*
species at 95 %.

Step 33b therefore fetches **every nt mitochondrial record** of each candidate
host family. Candidates are the families step 33 named on any usable mtDNA
contig, plus the family of every metadata host. At most 3 records are kept per
species × gene, and 2 complete genomes per species. Adding the recorded host's
family is not circular: it only adds references, and identity decides.

Step 33c re-BLASTs the host contigs against those records + RefSeq and names
each contig like for like:

1. The contig is cut into 300 bp windows. Each HSP (≥ 100 bp) gives its
   taxon an identity in each window it covers.
2. Candidates are the taxa within 0.5 % of the best identity in any window.
3. **Pairwise, on shared windows only**, A eliminates B if A's mean identity
   there is > 0.5 % higher. A congener's complete mitogenome at 95 % is
   eliminated by the true species' 99 % CYTB on the CYTB windows. It cannot
   win by covering 50 other windows where the true species has no sequence.
   (A naive window vote would let it do exactly that.)
4. The surviving taxa name the contig (LCA), capped by identity as above.
   Survivors that never overlap (one known only from 12S, another only from
   CYTB) cannot be separated and are named at their common ancestor.

`tables/refined_contig_calls.tsv` lists the survivors with identity and
number of windows per contig. Step 34 uses the refined name wherever there is
one, and keeps the first-pass nt-LCA name otherwise.

## 4. Per-sample call

- Host taxa = the tier's contigs grouped by name; each must have ≥ 1 contig
  and ≥ 5 mapped reads (`mapping/<s>.coverage.tsv`).
- **Lab DNA** (human, chicken, cow, pig: `call.lab_contaminants`) and **lane
  cross-talk** are reported, but kept out of the primary call unless they are
  the only signal. A family is cross-talk when a sample on the same lane
  (`seq_run`) carries ≥ 100× its reads. Human is a genuine recorded host for
  4 samples, so it is never discarded.
- **Primary call**: from species up to class, the first rank at which one
  taxon holds ≥ 80 % of the read weight resolved to that rank.
- **Hosts**: every leaf taxon with ≥ 10 % of the weight. `multiple_families`
  marks two or more non-lab, non-cross-talk families. That can mean a mixed
  or sequential blood meal, a pooled tube, or contamination.

## 5. Comparison with the metadata

The recorded host is resolved by step 32. Corrections are in
`metadata/host_name_fixes.tsv`: 7 misspellings, 2 genus-only entries, 1
subspecies, the cattle egret, and **one conflict to check**: "Brown boobie
(*Onychoprion fuscatus*)". The common name is the brown booby (*Sula
leucogaster*, Sulidae), but the scientific name is the sooty tern (Laridae).
The scientific name is used.

The call and the recorded host are compared rank by rank from class down:

| outcome | meaning |
|---|---|
| `MATCH_<rank>` | agree at every rank both resolve; `<rank>` is the deepest compared |
| `MISMATCH_<rank>` | first rank at which they differ |
| `no_bloodmeal_detected` / `no_dominant_host` | nothing passed, or no taxon reached 80 % |
| `no_host_recorded` / `questing` | nothing to compare against |

`metadata_host_in_hosts` = the recorded host is consistent with *any*
reported host, not just the primary one. `meta_nt_mito_species` = how many nt
mito records exist for the recorded host species. A species-level mismatch
with 0–5 records is very likely a reference gap, not a wrong host, and the
`note` column says so.

## Outputs (`results/bloodmeal/`)

| file | |
|---|---|
| `tables/bloodmeal_calls.tsv` | **one row per sample**: metadata host, tier, primary call, hosts, comparison, notes |
| `tables/bloodmeal_host_taxa.tsv` | sample × host taxon: reads, contigs, identity, genes, flags |
| `tables/bloodmeal_concordance.tsv` | outcome counts overall and by host class / tick genus / life stage / tier |
| `tables/vertebrate_contigs.tsv` | every vertebrate-hit contig, incl. rejected ones and why |
| `tables/refined_contig_calls.tsv` | step 33c survivors per contig |
| `tables/metadata_hosts.tsv` | metadata host → taxon, reference availability |
| `tables/nuclear_vertebrate_background.tsv` | reads on contigs whose whole-contig nt LCA is vertebrate, by class (context only; never used for a call) |
| `figures/` | `bloodmeal_concordance`, `bloodmeal_class_matrix`, `bloodmeal_per_sample` (png + svg) |

---

## Results on this cohort (run 2026-10-13)

146 samples analysed; 113 have a recorded host.

| | samples |
|---|---|
| host DNA found, tier 1 (mito), ≥ 5 reads | 65 |
| tier 2 only (18S / short COI) | 8 |
| weak (mito contigs, < 5 reads) | 12 |
| no vertebrate mito / 18S contig at all | 61 |

Against the 113 recorded hosts: **42 agree** (19 species, 12 genus, 4
family, 7 order/class), **19 disagree**, 3 show only lab DNA / lane cross-talk,
4 have no dominant host, and 45 have no detectable blood meal.

Disagreements worth reading first (`tables/bloodmeal_calls.tsv`):

- **8336-Argasidae-4L-Boobie**: recorded "Brown boobie (*Onychoprion
  fuscatus*)". The blood meal is *Sula dactylatra* (masked booby), 100 % over 8
  contigs. The common name is right and the scientific name is wrong.
- **8398-Booyal-3-H-B-F**: recorded *Notamacropus rufogriseus*. A complete 17 kb
  host mitogenome is 98.5–99.6 % *N. dorsalis* (black-striped wallaby) and only
  90–92 % to both *N. rufogriseus* complete mitogenomes. This is a genuine
  species disagreement, not a reference gap.
- **8016 / 8016-F1-Iorni** (*Ixodes ornithorhynchi*, the platypus tick):
  recorded bridled tern; host DNA is Monotremata (platypus mtDNA contig
  at 99.3 % + mammalian 18S).
- **8395** (recorded long-nosed bandicoot) → echidna; **B4926** (southern brown
  bandicoot) → *Macropus fuliginosus*; **8220** (red-necked wallaby) →
  Peramelidae; **0895, 5070** (echidna / grey kangaroo) → *Varanus*.
- **8244-A-M1-Alimb** (water python): the python is there (*Liasis*, 99 %), and
  so is a tortoise (Testudinidae, ~92 %; no closer reference in nt), at
  ~45 % each. AZWH samples 8352/8356A/8378 carry the same tortoise signal at low
  reads.

Caveats: no extraction blanks were sequenced, so low-read signals cannot be
separated from contamination by design. The `cross_talk` flag only covers
same-lane bleed. Read the `evidence`, `cross_talk_taxa` and `note` columns
before treating a single-contig call as a blood meal.
