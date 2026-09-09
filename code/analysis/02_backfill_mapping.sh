#!/bin/bash
#SBATCH --job-name=map_backfill
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/map_backfill_%A_%a.log
#SBATCH --error=results/logs/map_backfill_%A_%a.err
#SBATCH --mem=100G
#SBATCH --array=0-23
#SBATCH --account=OD-240875
#
# Step 2 -- backfill read mapping for samples that have an LCA but no coverage.
#
# Thin wrapper around the existing code/map_reads_to_contigs.sh, which already
# does bowtie2 --end-to-end --no-unal | samtools sort | index | idxstats -> TPM.
# None of the mapping or TPM maths is reimplemented here.
#
# Queue: results/audit/read_dir_priority.tsv (from 00_audit_samples.py). It
# resolves each sample to the MOST-REPAIRED read pair available
# (repaired_QC/fastqrepair/reformat > fastqrepair > repaired_QC > QC). That
# matters: the QC/ copies for these samples have pairing defects (see
# code/fastqrepair.sh) and bowtie2 paired-end either fails or mis-pairs on them.
#
# NOTE: the commented-out mapping command in code/fastqrepair.sh hardcodes
# --contigs contigs/B1201_M1_1_96_S24.contigs.fa for EVERY sample, which would
# map every library against one sample's assembly. This script takes the
# per-sample contigs path from the queue instead.
#
#   sbatch code/analysis/02_backfill_mapping.sh                 # cluster (24 tasks)
#   THREADS=8 nohup bash code/analysis/02_backfill_mapping.sh & # serial fallback
#
# Optional env:
#   THREADS        override thread count (default: SLURM_CPUS_PER_TASK, else 4)
#   BOWTIE2_BIN    directory holding bowtie2/samtools if not already on PATH
#
# PORTABILITY: POSIX sh. Slurm at this site runs the batch script with /bin/sh
# regardless of the shebang, so no [[ ]], (( )), mapfile, <<< or arrays.

module load bowtie/2.5.4
module load samtools

set -u

# Locate the project root.
#
# Under Slurm, $0 is the spool copy of this script
# (/cm/local/apps/slurm/var/spool/jobNNNN/slurm_script), NOT its path in the
# repo, so deriving the root from $0 alone breaks on the cluster. Order:
#   1. TICK_ROOT if the user set it
#   2. SLURM_SUBMIT_DIR (Slurm's default --chdir is the submit directory)
#   3. the script's own location (correct for a direct local invocation)
#   4. the current directory
# Whichever wins must contain code/analysis/, or we stop rather than guess.
find_root() {
    for cand in "${TICK_ROOT:-}" "${SLURM_SUBMIT_DIR:-}" \
                "$(cd "$(dirname "$0")/../.." 2>/dev/null && pwd)" "$(pwd)"; do
        [ -n "$cand" ] || continue
        if [ -d "$cand/code/analysis" ]; then printf '%s\n' "$cand"; return 0; fi
    done
    return 1
}

BASEDIR=$(find_root) || {
    echo "ERROR: cannot locate the project root." >&2
    echo "  Submit from the project directory, or set TICK_ROOT, e.g." >&2
    echo "  sbatch --export=ALL,TICK_ROOT=/path/to/old_analysis $0" >&2
    exit 1
}
cd "$BASEDIR" || exit 1
echo "project root: $BASEDIR"

QUEUE="results/audit/read_dir_priority.tsv"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"

mkdir -p results/logs mapping

if [ ! -f "$QUEUE" ]; then
    echo "ERROR: $QUEUE missing. Run 00_audit_samples.py first." >&2
    exit 1
fi

# bowtie2 lives in a conda env on this system rather than on the default PATH.
#if ! command -v bowtie2 >/dev/null 2>&1; then
#    PATH="${BOWTIE2_BIN:-/data/hb-austick/work/miniforge3/envs/bowtie2/bin}:$PATH"
#    export PATH
#fi
#for prog in bowtie2 bowtie2-build samtools; do
#    command -v "$prog" >/dev/null 2>&1 || {
#        echo "ERROR: $prog not on PATH. Set BOWTIE2_BIN or module load it." >&2
#        exit 1
#    }
#done

run_one() {
    row="$1"
    sample=$(printf '%s\n'  "$row" | cut -f1)
    contigs=$(printf '%s\n' "$row" | cut -f2)
    r1=$(printf '%s\n'      "$row" | cut -f3)
    r2=$(printf '%s\n'      "$row" | cut -f4)
    src=$(printf '%s\n'     "$row" | cut -f5)

    [ -n "$sample" ] || { echo "SKIP  empty queue row"; return 0; }

    if [ -s "mapping/${sample}.coverage.tsv" ]; then
        echo "SKIP  $sample - coverage already exists"
        return 0
    fi
    for f in "$contigs" "$r1" "$r2"; do
        if [ ! -f "$f" ]; then
            echo "SKIP  $sample - missing $f" >&2
            return 0
        fi
    done

    echo "==== $sample  reads from ${src}  start $(date -Is)"
    bash code/map_reads_to_contigs.sh \
        --sample  "$sample" \
        --contigs "$contigs" \
        --r1 "$r1" --r2 "$r2" \
        --outdir  mapping \
        --threads "$THREADS"
    if [ $? -ne 0 ]; then
        echo "FAILED mapping: $sample" >&2
        return 1
    fi
    echo "==== $sample  done $(date -Is)"
}

N=$(awk 'NR>1 && NF' "$QUEUE" | wc -l)

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    IDX="$SLURM_ARRAY_TASK_ID"
    if [ "$IDX" -ge "$N" ]; then
        echo "array index $IDX >= queue length $N - nothing to do"
        exit 0
    fi
    ROW=$(awk -v n="$IDX" 'NR>1 && NF { if (i++ == n) { print; exit } }' "$QUEUE")
    run_one "$ROW"
else
    echo "No SLURM_ARRAY_TASK_ID - running $N samples serially (${THREADS} threads)"
    i=0
    while [ "$i" -lt "$N" ]; do
        ROW=$(awk -v n="$i" 'NR>1 && NF { if (j++ == n) { print; exit } }' "$QUEUE")
        run_one "$ROW"
        i=$((i + 1))
    done
fi
