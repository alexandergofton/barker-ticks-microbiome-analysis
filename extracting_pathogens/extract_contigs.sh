#!/bin/bash

# Parse arguments
while getopts "k:" opt; do
  case $opt in
    k)
      keyword="$OPTARG"
      ;;
    *)
      echo "Usage: $0 -k <keyword>"
      exit 1
      ;;
  esac
done

# Check if keyword is provided
if [ -z "$keyword" ]; then
  echo "Error: Keyword not provided."
  echo "Usage: $0 -k <keyword>"
  exit 1
fi

blast_out="${keyword}.blast.txt"
seqids="${keyword}.blast.seqids"
blast_files="${keyword}.blast.files"

# Use grep to search for the keyword in blast.txt
#grep $keyword blast_results/*.blast.txt > $blast_out

# Get seq ids
awk -F "\t" '{print $1}' $blast_out | awk -F ':' '{print $2}' | sort -u > $seqids

# Get files that contain contigs
awk -F "\t" '{print $1}' $blast_out | awk -F ':' '{print $1}' | sort -u > $blast_files

mkdir -p $keyword

# Extract contigs into a new file
for fasta in contigs/*; do
    seqtk subseq $fasta $seqids > ${keyword}/$(basename $fasta)
done

# remove empty files
find $keyword -type f -empty -exec rm {} \;
