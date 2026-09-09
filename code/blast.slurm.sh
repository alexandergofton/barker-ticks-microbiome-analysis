#!/bin/bash

#SBATCH --job-name=blastn
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/blast_%A_%a.log
#SBATCH --error=logs/blast_%A_%a.err
#SBATCH --mem=512G
#SBATCH --array=0-7
#SBATCH --account=OD-240875
##SBATCH --qos=express

module load blast+/2.16.0
module load bioref

BASEDIR="/datasets/work/hb-austicks/work/AlexG_data/Barker_ticks/microbe_analysis/old_analysis"

SPADESDIR=${BASEDIR}/contigs
OUTDIR=${BASEDIR}/blast

SAMPLE=( "8392-18-4-24-Amoreliae-1F.contigs.fa" \
  "8394-Btachyg-2M-5Legs-EastEnd.contigs.fa" \
  "8395-Btachyg-2M-3Legs-Damscus.contigs.fa" \
  "8398-Booyal-3-H-B-F.contigs.fa" \
  "8399-F-Benaraby-M-gig-14-4-24.contigs.fa" \
  "8401-Hban-1F-Farnborough.contigs.fa" \
  "8404-Amtrig-M-EtnaCk.contigs.fa" \
  "8405-Bhydrosauri-1F-VIC.contigs.fa" );

OUTNAME=( "8392-18-4-24-Amoreliae-1F.blast.txt" \
  "8394-Btachyg-2M-5Legs-EastEnd.blast.txt" \
  "8395-Btachyg-2M-3Legs-Damscus.blast.txt" \
  "8398-Booyal-3-H-B-F.blast.txt" \
  "8399-F-Benaraby-M-gig-14-4-24.blast.txt" \
  "8401-Hban-1F-Farnborough.blast.txt" \
  "8404-Amtrig-M-EtnaCk.blast.txt" \
  "8405-Bhydrosauri-1F-VIC.blast.txt" );

DBPATH="/data/bioref/blast/ncbi_20250504/nt"

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]; then
  i=${SLURM_ARRAY_TASK_ID}

  #source ~/bin/localntdb.sh DBPATH

  blastn \
  -query ${SPADESDIR}/${SAMPLE[$i]} \
  -db $DBPATH \
  -out ${OUTDIR}/${OUTNAME[$i]}.blastn \
  -outfmt '6 qaccver saccver stitle sseqid ssciname scomname staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -num_threads ${SLURM_CPUS_PER_TASK} \
  -max_target_seqs 10 \
  -evalue 1e-5

else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
