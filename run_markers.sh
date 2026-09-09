#!/bin/bash
# Phylogenetic marker gene workflow: find, extract and align the marker genes
# for every target taxon, ready for you to add reference sequences and build
# trees.
#
#   bash run_markers.sh                 # everything
#   bash run_markers.sh --target Rickettsia --target Coxiella
#
# Runs after run_all.sh -- step 10 must have written
# results/targets/target_contig_index.tsv.
#
# The panel of loci is metadata/target_markers.tsv; the settings are the
# `markers:` block of config/analysis_params.yaml. Edit those, not this file.
#
# Optional env:
#   PYTHON        python3 to use
#   MARKER_ENV    prefix of the conda env holding mafft (default: tick_markers)
#   NCBI_API_KEY  raises the Entrez rate limit from 3/s to 10/s
set -u
cd "$(dirname "$0")" || exit 1

PY="${PYTHON:-python3}"
MARKER_ENV="${MARKER_ENV:-/data/hb-austick/work/miniforge3/envs/tick_markers}"

# mafft lives only in the marker env. Prepend it rather than requiring the user
# to activate: every script finds its tools through $PATH and nothing hardcodes
# a path into the env.
if [ -x "$MARKER_ENV/bin/mafft" ]; then
    PATH="$MARKER_ENV/bin:$PATH"; export PATH
else
    echo "WARNING: mafft not found at $MARKER_ENV/bin."
    echo "  mamba env create -f env/tick_markers.yml"
    echo "  Steps 11-14 will still run; step 15 (alignment) will stop."
fi

step() { echo; echo "=============== $* ==============="; }
die()  { echo "FAILED: $*" >&2; exit 1; }

step "11 fetch marker reference probes (cached; Entrez)"
$PY code/analysis/11_fetch_marker_refs.py            || die "marker refs"

step "12 extract target contigs from the assemblies"
$PY code/analysis/12_extract_target_contigs.py "$@"  || die "target contigs"

step "13 locate marker loci (blastn / tblastn)"
$PY code/analysis/13_find_markers.py "$@"            || die "find markers"

step "14 cut out and label the marker sequences"
$PY code/analysis/14_extract_markers.py "$@"         || die "extract markers"

step "15 align each target x marker set"
$PY code/analysis/15_align_markers.py "$@"           || die "align markers"

step "16 recovery report"
$PY code/analysis/16_marker_report.py                || die "marker report"

echo
echo "=============== done ==============="
echo "alignments: results/markers/alignments/<target>/<target>__<marker>.aln.fasta"
echo "next:       add your reference sequences, then build trees."
echo "            see README_markers.md"
