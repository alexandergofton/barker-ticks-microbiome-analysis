#!/usr/bin/env bash
# Extract contigs assigned by LCA to North shore virus (taxid 2485875),
# one FASTA per sample, plus every BLAST hit for those contigs (all subjects, not
# only the North shore virus ones) so the LCA call can be checked.
set -euo pipefail
ROOT=/data/hb-austick/work/AlexG_data/Barker_ticks/microbe_analysis/old_analysis
OUT=$ROOT/extracting_pathogens/north_shore_virus
TAXID=2485875
BLAST_HDR=$'qaccver\tsaccver\tstitle\tsseqid\tssciname\tscomname\tstaxid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'

printf 'sample\tquery\tlca_taxid\tlca_name\tlca_rank\tn_hits\ttop_bitscore\n' > "$OUT/lca_assignments.tsv"
printf 'sample\tn_contigs\tn_fasta_seqs\ttotal_bp\tn_blast_hits\n' > "$OUT/summary.tsv"

for lca in "$ROOT"/lca/*.blast.lca; do
  s=$(basename "$lca" .blast.lca)
  awk -F'\t' -v t=$TAXID '$2==t' "$lca" > "$OUT/ids/$s.lca.tmp"
  if [[ ! -s "$OUT/ids/$s.lca.tmp" ]]; then rm "$OUT/ids/$s.lca.tmp"; continue; fi
  awk -v s="$s" '{print s"\t"$0}' "$OUT/ids/$s.lca.tmp" >> "$OUT/lca_assignments.tsv"
  cut -f1 "$OUT/ids/$s.lca.tmp" > "$OUT/ids/$s.ids"; rm "$OUT/ids/$s.lca.tmp"

  # A few LCA files carry a suffix the contigs/BLAST files lack (e.g.
  # 8334-AMLoculosum-BriddledTurn7 vs 8334-AMLoculosum-BriddledTurn), so fall
  # back to the longest contigs file name that is a prefix of the sample name.
  fa="$ROOT/contigs/$s.contigs.fa"
  [[ -f "$fa" ]] || fa="$ROOT/contigs_with_missing_reads/$s.contigs.fa"
  base=$s
  if [[ ! -f "$fa" ]]; then
    fa=""; best=0
    for c in "$ROOT"/contigs/*.contigs.fa "$ROOT"/contigs_with_missing_reads/*.contigs.fa; do
      b=$(basename "$c" .contigs.fa)
      if [[ "$s" == "$b"* && ${#b} -gt $best ]]; then fa=$c; base=$b; best=${#b}; fi
    done
    if [[ -z "$fa" ]]; then echo "WARNING: no contigs file for $s, skipped" >&2; continue; fi
    echo "NOTE: $s -> contigs/BLAST files named $base" >&2
  fi
  seqkit grep -f "$OUT/ids/$s.ids" "$fa" > "$OUT/fasta/$s.north_shore_virus.fasta"

  { echo "$BLAST_HDR"
    awk -F'\t' 'NR==FNR{k[$1]; next} ($1 in k)' "$OUT/ids/$s.ids" "$ROOT/blast/$base.blast.txt"
  } > "$OUT/blast/$s.north_shore_virus.blast.tsv"

  n=$(wc -l < "$OUT/ids/$s.ids")
  read nf bp < <(seqkit stats -T "$OUT/fasta/$s.north_shore_virus.fasta" | awk 'NR==2{print $4, $5}')
  nb=$(( $(wc -l < "$OUT/blast/$s.north_shore_virus.blast.tsv") - 1 ))
  printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$n" "$nf" "$bp" "$nb" >> "$OUT/summary.tsv"
  echo "$(date +%T) $s contigs=$n fasta=$nf blast_hits=$nb"
done
