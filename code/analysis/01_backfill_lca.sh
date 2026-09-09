#!/bin/bash
#SBATCH --job-name=lca_backfill
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --output=results/logs/lca_backfill_%A_%a.log
#SBATCH --error=results/logs/lca_backfill_%A_%a.err
#SBATCH --mem=16G
#SBATCH --array=0-11
#SBATCH --account=OD-240875
#
# Step 1 -- backfill LCA for samples that have BLAST output but no .blast.lca.
#
# Queue: results/audit/lca_backfill_queue.tsv (from 00_audit_samples.py),
# largest BLAST file first. Array index N processes queue row N.
#
#   sbatch code/analysis/01_backfill_lca.sh                  # cluster (12 tasks)
#   MAX_GB=2 nohup bash code/analysis/01_backfill_lca.sh &   # serial, small files only
#
# Optional env:
#   MAX_GB   skip any BLAST file larger than this many GB
#   PYTHON   python interpreter to use (default: python3 on PATH)
#
# PORTABILITY: written to POSIX sh. Slurm at this site executes the batch
# script with /bin/sh regardless of the shebang, so bashisms fail at submit
# time -- no [[ ]], no (( )), no mapfile, no <<< here-strings, no arrays, and
# no bc (awk does the float comparison instead).

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

QUEUE="results/audit/lca_backfill_queue.tsv"
NODES="taxdump/nodes.dmp"
NAMES="taxdump/names.dmp"
MERGED="taxdump/merged.filtered.dmp"
MAX_GB=180
PYTHON="${PYTHON:-python3}"

mkdir -p results/logs lca

for f in "$QUEUE" "$NODES" "$NAMES" "$MERGED"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing $f" >&2
        echo "  run: python3 code/analysis/00_audit_samples.py" >&2
        echo "       python3 code/analysis/make_merged_filtered.py" >&2
        exit 1
    fi
done

# gt A B -> exit 0 if A > B. awk, not bc: bc is often absent on compute nodes.
gt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a > b) }'; }

run_one() {
    row="$1"
    sample=$(printf '%s\n' "$row" | cut -f1)
    blast=$(printf '%s\n'  "$row" | cut -f2)
    gb=$(printf '%s\n'     "$row" | cut -f3)

    [ -n "$sample" ] || { echo "SKIP  empty queue row"; return 0; }

    if gt "$gb" "$MAX_GB"; then
        echo "SKIP  $sample - ${gb} GB exceeds MAX_GB=${MAX_GB}"
        return 0
    fi
    if [ -s "lca/${sample}.blast.lca" ]; then
        echo "SKIP  $sample - lca already exists"
        return 0
    fi
    if [ ! -f "$blast" ]; then
        echo "SKIP  $sample - BLAST file not found: $blast" >&2
        return 0
    fi

    echo "==== $sample  ${gb} GB  start $(date -Is)"

    # The duplicate-query guard holds one hash per query and is the only part
    # of lca_stream.py that scales with input size, so it is dropped on the
    # very large files. BLAST emits hits grouped by query, so the guard is a
    # belt-and-braces check rather than a correctness requirement.
    EXTRA=""
    if gt "$gb" 20; then EXTRA="--no-order-check"; fi

    # shellcheck disable=SC2086
    "$PYTHON" code/analysis/lca_stream.py \
        -i "$blast" \
        -o "lca/${sample}.blast.lca" \
        --nodes "$NODES" --names "$NAMES" --merged "$MERGED" \
        --evalue 1e-5 --pident 0 --top-percent 10 \
        $EXTRA
    if [ $? -ne 0 ]; then
        echo "FAILED lca: $sample" >&2
        return 1
    fi

    # Pavian view. NOTE: the kreport counts CONTIGS, not reads, so it is a
    # richness view only. Read-weighted abundance comes from Step 4.
    "$PYTHON" code/lca_to_kreport.py \
        -i "lca/${sample}.blast.lca" \
        -o "lca/${sample}.blast.lca.kreport" \
        --nodes "$NODES" --names "$NAMES" \
      || echo "WARN: kreport failed for $sample" >&2

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
    echo "No SLURM_ARRAY_TASK_ID - running $N samples serially (MAX_GB=$MAX_GB)"
    i=0
    while [ "$i" -lt "$N" ]; do
        ROW=$(awk -v n="$i" 'NR>1 && NF { if (j++ == n) { print; exit } }' "$QUEUE")
        run_one "$ROW"
        i=$((i + 1))
    done
fi
