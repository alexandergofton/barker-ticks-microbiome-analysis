#!/bin/bash
# Marker-gene phylogenies: download the reference sequences for every locus the
# samples yielded, cut the locus out of them, align and build a tree.
#
#   bash run_trees.sh                                   # every locus
#   bash run_trees.sh --target Rickettsia --target Coxiella
#   bash run_trees.sh --target Coxiella --marker rpoB
#
# Runs after run_markers.sh -- step 14 must have written
# results/markers/sequences/<target>/<target>__<marker>.fna.
#
# The reference panels are metadata/tree_references.tsv (per target x marker)
# and metadata/tree_ref_genomes.tsv (per target); the settings are the `trees:`
# block of config/analysis_params.yaml. Edit those, not this file.
#
# Every locus ends up as three files in results/markers/trees/<target>/:
#   <target>__<marker>.aln.fna      alignment (codon-aware for protein loci)
#   <target>__<marker>.nwk          FastTree output, unrooted
#   <target>__<marker>.rooted.nwk   rooted -- on the outgroup, else midpoint
#   <target>__<marker>.png          figure
#
# Optional env:
#   PYTHON        python3 to use
#   MARKER_ENV    prefix of the conda env holding mafft/fasttree (default: tick_markers)
#   RSCRIPT       Rscript for step 19b (default: the tick_micro conda env)
#   NCBI_API_KEY  raises the Entrez rate limit from 3/s to 10/s -- worth setting,
#                 step 17 makes hundreds of requests
set -u
cd "$(dirname "$0")" || exit 1

PY="${PYTHON:-python3}"
MARKER_ENV="${MARKER_ENV:-/data/hb-austick/work/miniforge3/envs/tick_markers}"
RS="${RSCRIPT:-/data/hb-austick/work/miniforge3/envs/tick_micro/bin/Rscript}"

# mafft, fasttree and trimal live only in the marker env; blast+ comes from the
# base install. Prepend rather than require an activated env, as run_markers.sh
# does -- every script finds its tools through $PATH.
if [ -x "$MARKER_ENV/bin/mafft" ] && [ -x "$MARKER_ENV/bin/fasttree" ]; then
    PATH="$MARKER_ENV/bin:$PATH"; export PATH
else
    echo "WARNING: mafft/fasttree not found at $MARKER_ENV/bin."
    echo "  mamba env create -f env/tick_markers.yml"
    echo "  Steps 17-18 will still run; step 19 (align + tree) will stop."
fi

step() { echo; echo "=============== $* ==============="; }
die()  { echo "FAILED: $*" >&2; exit 1; }

step "17 download reference sequences for each locus (cached; Entrez)"
$PY code/analysis/17_fetch_tree_refs.py "$@"     || die "fetch tree refs"

step "18 cut the marker locus out of each reference (blastn / tblastn)"
$PY code/analysis/18_curate_tree_refs.py "$@"    || die "curate tree refs"

step "19 align samples + references, build one tree per locus (FastTree)"
$PY code/analysis/19_build_trees.py "$@"         || die "build trees"

step "19b root each tree and draw it (R)"
# The R step takes bare target names, not --target flags, and ignores --marker:
# it draws whatever trees exist for those targets.
targets=""
prev=""
for a in "$@"; do
    [ "$prev" = "--target" ] && targets="$targets $a"
    prev="$a"
done
# shellcheck disable=SC2086
$RS code/analysis/19b_tree_figures.R $targets     || die "tree figures"

echo
echo "=============== done ==============="
echo "trees:   results/markers/trees/<target>/<target>__<marker>.{aln.fna,nwk,rooted.nwk,png}"
echo "tables:  results/markers/trees/tree_summary.tsv"
echo "         results/markers/trees/tree_reference_sequences.tsv"
echo "see README_trees.md"
