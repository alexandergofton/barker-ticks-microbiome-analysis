#!/bin/bash

#SBATCH --job-name=bowtie2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --output=logs/mapping_%A.log
#SBATCH --error=logs/mapping_%A.err
#SBATCH --mem=200G
##SBATCH --array=0-7
#SBATCH --account=OD-240875
##SBATCH --qos=express

module load bowtie/2.5.4
module load samtools/1.18.0

BASEDIR="/datasets/work/hb-austicks/work/AlexG_data/Barker_ticks/microbe_analysis/old_analysis"

SPADESDIR=${BASEDIR}/contigs
READSDIR=${BASEDIR}/QC
OUTDIR=${BASEDIR}/mapping
CODEDIR=${BASEDIR}/code

for CONTIGFILE in ${SPADESDIR}/*.contigs.fa
do
  ${CODEDIR}/map_reads_to_contigs.sh \
    --sample $(basename $CONTIGFILE .contigs.fa) \
		--contigs $CONTIGFILE \
		--r1 ${READSDIR}/$(basename $CONTIGFILE .contigs.fa).R1.fq.gz \
		--r2 ${READSDIR}/$(basename $CONTIGFILE .contigs.fa).R2.fq.gz \
		--outdir ${OUTDIR} \
		--threads ${SLURM_CPUS_PER_TASK}
done
