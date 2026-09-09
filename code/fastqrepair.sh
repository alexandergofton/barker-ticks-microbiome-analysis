#!/bin/bash

for x in B2490_F1_74_9_S2 B3015_F1_0_49_S28 B3648_M1_14_8_S19 B4092_F1_28_3_S15 B4276_F1_28_7_S20 B5053a_N1_1_14_S10 B6565_F1_1_03_S13 B6658_L1_15_S25 B6852_F1_38_S1 B6894_F1_0_913_S23 B7002_1_F1_5_22_S9 B7017_M2_S22 B7024_M1_21_2_S11 B7106_F1_S21 B7115_N3_S40 B7115_N5_S41 B7116_N3_8_78_S8 B7119_F1_5_34_S5 B7203_N1_B3644_S26 B7217_F1_8_44_S22 B7229_F1_11_1_S17 B7380_F1_34_S21
do

  R1=${x}.R1.fq.gz
  R2=${x}.R2.fq.gz

  echo "------------"; echo ""; echo $x; echo ""

  # repair.sh \
  #   in=repaired_QC/${R1} \
  #   in2=repaired_QC/${R2} \
  #   out=repaired_QC/fastqrepair/${R1} \
  #   out2=repaired_QC/fastqrepair/${R2} \
  #   outs=repaired_QC/fastqrepair/${x}.singletons.fq.gz \
  #   -Xmx16g

  #zcat repaired_QC/fastqrepair/${R1} | wc -l
  #zcat repaired_QC/fastqrepair/${R2} | wc -l

  #zcat repaired_QC/fastqrepair/${R1} | awk 'NR%4==2 || NR%4==0 {if (length($0) == 0) print "Empty data line at line " NR}'
  #zcat repaired_QC/fastqrepair/${R2} | awk 'NR%4==2 || NR%4==0 {if (length($0) == 0) print "Empty data line at line " NR}'

  # reformat.sh \
  #   in1=repaired_QC/fastqrepair/${R1} \
  #   in2=repaired_QC/fastqrepair/${R2} \
  #   out1=repaired_QC/fastqrepair/reformat/${R1} \
  #   out2=repaired_QC/fastqrepair/reformat/${R2} \
  #   minlen=10 \
  #   -Xmx16g

  # code/map_reads_to_contigs.sh \
  # --sample ${x} \
  # --contigs contigs/B1201_M1_1_96_S24.contigs.fa \
  # --r1 repaired_QC/fastqrepair/reformat/${R1} \
  # --r2 repaired_QC/fastqrepair/reformat/${R2} \
  # --outdir mapping \
  # --threads 8

  ls -1 mapping/${x}*
done
