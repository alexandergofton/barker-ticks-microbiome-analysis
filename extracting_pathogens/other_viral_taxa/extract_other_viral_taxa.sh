#!/usr/bin/env bash
# Extract every contig whose LCA is a virus OTHER than Cheeloo Jingmen-like
# virus (3064049) and North shore virus (2485875) -- those two have their own
# directories -- one FASTA per sample, plus every BLAST hit for those contigs
# (all subjects, not only viral ones) so the LCA call can be checked.
#
# viral_taxids.tsv = results/taxonomy/taxon_lineage.tsv rows with
# superkingdom_group == Viruses, minus the two taxids above.
#
# Same layout as ../north_shore_virus/extract_north_shore_virus.sh. Resumable:
# a sample already listed in summary.tsv (written last) is skipped.
set -euo pipefail
ROOT=/data/hb-austick/work/AlexG_data/Barker_ticks/microbe_analysis/old_analysis
OUT=$ROOT/extracting_pathogens/other_viral_taxa
TAXA=$OUT/viral_taxids.tsv
BLAST_HDR=$'qaccver\tsaccver\tstitle\tsseqid\tssciname\tscomname\tstaxid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'

[[ -f "$OUT/lca_assignments.tsv" ]] || printf 'sample\tquery\tlca_taxid\tlca_name\tlca_rank\tn_hits\ttop_bitscore\n' > "$OUT/lca_assignments.tsv"
[[ -f "$OUT/summary.tsv" ]]         || printf 'sample\tn_contigs\tn_fasta_seqs\ttotal_bp\tn_blast_hits\n' > "$OUT/summary.tsv"

for lca in "$ROOT"/lca/*.blast.lca; do
  s=$(basename "$lca" .blast.lca)
  grep -q "^$s	" "$OUT/summary.tsv" && continue
  sed -i "/^$s	/d" "$OUT/lca_assignments.tsv"   # drop a partial earlier attempt
  awk -F'\t' 'NR==FNR{t[$1]; next} FNR>1 && ($2 in t)' "$TAXA" "$lca" > "$OUT/ids/$s.lca.tmp"
  if [[ ! -s "$OUT/ids/$s.lca.tmp" ]]; then rm "$OUT/ids/$s.lca.tmp"; continue; fi
  cut -f1 "$OUT/ids/$s.lca.tmp" > "$OUT/ids/$s.ids"

  # LCA stems can carry a suffix the contigs/BLAST files lack (8334...Turn7):
  # fall back to the longest contigs file name that is a prefix of the sample.
  fa="$ROOT/contigs/$s.contigs.fa"
  [[ -f "$fa" ]] || fa="$ROOT/contigs_with_missing_reads/$s.contigs.fa"
  base=$s
  if [[ ! -f "$fa" ]]; then
    fa=""; best=0
    for c in "$ROOT"/contigs/*.contigs.fa "$ROOT"/contigs_with_missing_reads/*.contigs.fa; do
      b=$(basename "$c" .contigs.fa)
      if [[ "$s" == "$b"* && ${#b} -gt $best ]]; then fa=$c; base=$b; best=${#b}; fi
    done
    if [[ -z "$fa" ]]; then echo "WARNING: no contigs file for $s, skipped" >&2; rm "$OUT/ids/$s.lca.tmp"; continue; fi
    echo "NOTE: $s -> contigs/BLAST files named $base" >&2
  fi
  seqkit grep -f "$OUT/ids/$s.ids" "$fa" > "$OUT/fasta/$s.viral.fasta"

  { echo "$BLAST_HDR"
    mawk -F'\t' 'NR==FNR{k[$1]; next} ($1 in k)' "$OUT/ids/$s.ids" "$ROOT/blast/$base.blast.txt"
  } > "$OUT/blast/$s.viral.blast.tsv"

  awk -v s="$s" '{print s"\t"$0}' "$OUT/ids/$s.lca.tmp" >> "$OUT/lca_assignments.tsv"; rm "$OUT/ids/$s.lca.tmp"
  n=$(wc -l < "$OUT/ids/$s.ids")
  read nf bp < <(seqkit stats -T "$OUT/fasta/$s.viral.fasta" | awk 'NR==2{print $4, $5}')
  nb=$(( $(wc -l < "$OUT/blast/$s.viral.blast.tsv") - 1 ))
  # summary row last: it marks the sample done for a restart
  printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$n" "$nf" "$bp" "$nb" >> "$OUT/summary.tsv"
  echo "$(date +%T) $s contigs=$n fasta=$nf blast_hits=$nb"
done
