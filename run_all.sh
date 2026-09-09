#!/bin/bash
# Run the full export/analysis pipeline. Steps 1 and 2 (the compute-heavy
# backfills) are NOT run here -- submit them to Slurm separately, see README.
#
#   bash run_all.sh
#
# Optional env:
#   RSCRIPT   Rscript to use (default: the tick_micro conda env)
set -u
cd "$(dirname "$0")" || exit 1

PY="${PYTHON:-python3}"
RS="${RSCRIPT:-/data/hb-austick/work/miniforge3/envs/tick_micro/bin/Rscript}"

step() { echo; echo "=============== $* ==============="; }
die()  { echo "FAILED: $*" >&2; exit 1; }

step "0  sample reconciliation audit"
$PY code/analysis/00_audit_samples.py            || die "audit"

step "3  taxonomy table"
$PY code/analysis/03_build_taxonomy.py           || die "taxonomy"

step "4  join LCA to mapping, build matrices"
$PY code/analysis/04_build_matrices.py           || die "matrices"

step "4b filter sensitivity"
$PY code/analysis/04b_filter_sensitivity.py 30   || echo "WARN: sensitivity failed (non-fatal)"

step "5  metadata template"
$PY code/analysis/05_make_metadata_template.py   || die "metadata template"

step "verification"
$PY code/analysis/verify.py                      || die "verification"

if [ -f metadata/sample_metadata.tsv ]; then
    step "6  validate metadata"
    $RS code/analysis/06_validate_metadata.R     || die "metadata validation"
else
    echo
    echo "SKIP step 6: metadata/sample_metadata.tsv not found."
    echo "  Fill in metadata/sample_metadata_TEMPLATE.tsv and save it under that"
    echo "  name. Until then the TSE falls back to the id-string guesses and no"
    echo "  biological hypothesis can be tested properly."
fi

step "7  build TreeSummarizedExperiment"
$RS code/analysis/07_build_tse.R                 || die "TSE"

step "8  composition"
$RS code/analysis/08_composition.R               || die "composition"

step "9  diversity"
$RS code/analysis/09_diversity.R                 || die "diversity"

step "10 pathogen / endosymbiont targets"
$PY code/analysis/10_export_targets.py           || die "targets"

echo
echo "=============== done ==============="
echo "outputs under results/"
