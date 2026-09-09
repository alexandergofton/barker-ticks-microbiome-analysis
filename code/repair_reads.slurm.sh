#!/bin/bash

#SBATCH --job-name=repair_fq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=512G
#SBATCH --output=code/logs/repair_fq_%A_%a.log
#SBATCH --error=code/logs/repair_fq_%A_%a.err
#SBATCH --account=OD-240875

module load bbmap/39.43

for x in \
	"B7115_N3_S40.R1.fq.gz" \
    "B7115_N5_S41.R1.fq.gz" \
    "B7116_N3_8_78_S8.R1.fq.gz" \
    "B7119_F1_5_34_S5.R1.fq.gz" \
    "B7203_N1_B3644_S26.R1.fq.gz" \
    "B7217_F1_8_44_S22.R1.fq.gz" \
    "B7229_F1_11_1_S17.R1.fq.gz" \
    "B7380_F1_34_S21.R1.fq.gz"
do
	R1=$(basename $x)
    R2=$(basename $x .R1.fq.gz).R2.fq.gz
    CONTIGS=$(basename $x .R1.fq.gz).contigs.fa
    SAMPLE=$(basename $x .R1.fq.gz)
	SINGLETONS=${SAMPLE}.singletons.fq.gz

    echo ${SAMPLE}
   
	#Repair reads
	repair.sh \
		in=QC/${R1} \
    	in2=QC/${R2} \
    	out=repaired_QC/${R1} \
    	out2=repaired_QC/${R2} \
		outs=repaired_QC/${SINGLETONS}
    	-Xmx512g \
    	-tossbrokenreads

        # ./map_reads_to_contigs.sh \
        #       --sample ${SAMPLE} \
        #       --contigs contigs/${CONTIGS} \
        #       --r1 QC/${R1} \
        #       --r2 QC/${R2} \
        #       --outdir mapping \
        #       --threads 48

   echo ""
done
