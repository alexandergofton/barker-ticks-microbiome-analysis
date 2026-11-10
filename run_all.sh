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
# Maps + report env (sf, leaflet, quarto); see env/tick_report.yml. Run through
# `conda run` because sf (PROJ/GDAL data) and quarto (QUARTO_SHARE_PATH) rely on
# variables set by the env's activation scripts.
CONDA="${CONDA_EXE:-/data/hb-austick/work/miniforge3/bin/conda}"
RS_REPORT="${RSCRIPT_REPORT:-$CONDA run -n tick_report --cwd $PWD Rscript}"

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

if [ -f metadata/barker_microbe_analysis_metadata.csv ]; then
    step "5b import curated metadata"
    $PY code/analysis/05b_import_metadata.py     || die "metadata import"
fi

step "verification"
$PY code/analysis/verify.py                      || die "verification"

if [ -f metadata/sample_metadata.tsv ]; then
    step "6  validate metadata"
    $RS code/analysis/06_validate_metadata.R     || die "metadata validation"
else
    echo
    echo "SKIP step 6: metadata/sample_metadata.tsv not found."
    echo "  It is built by step 5b from metadata/barker_microbe_analysis_metadata.csv."
    echo "  Until then the TSE falls back to the id-string guesses and no"
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

if [ -f metadata/sample_metadata.tsv ]; then
    step "10b tick genus / species comparison"
    $RS code/analysis/10b_tick_group_comparison.R || die "tick group comparison"

    step "10c sampling maps"
    $RS_REPORT code/analysis/10c_sampling_maps.R || die "maps"
fi

if [ "${RENDER_REPORT:-1}" = "1" ]; then
    step "report  (quarto render; RENDER_REPORT=0 to skip)"
    $CONDA run -n tick_report --cwd "$PWD/report" quarto render tick_microbiome_report.qmd \
        || echo "WARN: report render failed (non-fatal)"
fi

echo
echo "=============== done ==============="
echo "outputs under results/"
