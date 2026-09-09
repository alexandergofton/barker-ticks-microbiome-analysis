#!/data/hb-austick/work/miniforge3/bin/python3

"""
Count sequences in FASTA files within a directory.

Finds files ending in ".contigs.fa" in the given directory, counts how
many sequences are in each (by counting header lines starting with ">"),
and writes the results to a CSV file.

count_contigs.py contigs/ -o contig_counts.csv
"""

import argparse
import csv
from pathlib import Path


def count_sequences(fasta_path):
    """Count sequences in a FASTA file by counting header lines."""
    count = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Count sequences in .contigs.fa files in a directory."
    )
    parser.add_argument("directory", help="Directory containing the FASTA files")
    parser.add_argument(
        "-o", "--output",
        default="sequence_counts.csv",
        help="Output CSV filename (default: sequence_counts.csv)"
    )
    parser.add_argument(
        "-r", "--recursive",
        action="store_true",
        help="Also search subdirectories"
    )
    args = parser.parse_args()

    directory = Path(args.directory)
    pattern = "**/*.contigs.fa" if args.recursive else "*.contigs.fa"
    fasta_files = sorted(directory.glob(pattern))

    if not fasta_files:
        print(f"No .contigs.fa files found in {directory}")
        return

    results = []
    for fasta_file in fasta_files:
        seq_count = count_sequences(fasta_file)
        results.append((fasta_file.name, seq_count))
        print(f"{fasta_file.name}: {seq_count} sequences")

    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["filename", "sequence_count"])
        writer.writerows(results)

    print(f"\nResults written to {args.output}")


if __name__ == "__main__":
    main()

