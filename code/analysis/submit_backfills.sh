#!/bin/bash
# Print (or run with --go) the sbatch commands for the two backfill steps,
# with --array sized to the CURRENT queue length.
#
# The #SBATCH --array in each script is a static default. Queues shrink as
# samples complete, so read the length at submit time instead of trusting it.
# Out-of-range tasks exit cleanly anyway; this just avoids wasting slots.
#
#   bash code/analysis/submit_backfills.sh          # show the commands
#   bash code/analysis/submit_backfills.sh --go     # submit them
set -u
cd "$(dirname "$0")/../.." || exit 1
ROOT=$(pwd)

emit() {
    queue="$1"; script="$2"
    n=$(awk 'NR>1 && NF' "$queue" 2>/dev/null | wc -l)
    if [ "$n" -eq 0 ]; then
        echo "# $(basename "$queue"): empty - nothing to submit"
        return
    fi
    cmd="sbatch --array=0-$((n-1)) --export=ALL,TICK_ROOT=$ROOT $script"
    echo "$cmd"
    if [ "${GO:-0}" = "1" ]; then eval "$cmd"; fi
}

[ "${1:-}" = "--go" ] && GO=1 || GO=0
export GO
echo "# project root: $ROOT"
emit results/audit/lca_backfill_queue.tsv  code/analysis/01_backfill_lca.sh
emit results/audit/read_dir_priority.tsv   code/analysis/02_backfill_mapping.sh
[ "$GO" = "0" ] && echo "# (re-run with --go to submit)"
exit 0
