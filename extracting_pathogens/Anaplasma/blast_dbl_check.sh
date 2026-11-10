#!/bin/bash

for x in *.fa
do
	blastn \
		-query $x \
		-db /data/bioref/blast/ncbi_20250504/nt \
		-out blast_dbl_check/${x}.blast \
		-outfmt '6 qaccver saccver stitle sseqid ssciname scomname staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
		-num_threads 24 \
		-max_target_seqs 10 \
		-evalue 1e-5
			
	echo $x
done
