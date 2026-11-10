#!/bin/bash
# Tick mitochondrial genome workflow: extract one mitogenome candidate per
# sample, annotate its genes, align each protein-coding gene separately --
# together with ~100 published reference tick mitogenomes
# (metadata/mito_ref_genomes.tsv) -- and build one partitioned ML phylogeny
# with per-gene models in IQ-TREE.
#
#   bash run_mito.sh                      # everything
#   bash run_mito.sh --sample 8404-Amtrig-M-EtnaCk --sample 2331-Apersicus
#
# (the --sample/--gene filters apply to steps 21-23 only; 24-26 always run over
#  whatever is in results/mito/ at the time, since they operate on the whole
#  cohort by construction)
#
# Runs independently of run_all.sh / run_markers.sh -- no shared inputs beyond
# contigs/ and metadata/.
#
# See README_mito.md for the MITOS2 note: the tool is installed (env mitos2)
# but its Zenodo-hosted reference data could not be fetched from this compute
# node, so annotation runs on the homology fallback by default.
#
# Optional env:
#   PYTHON        python3 to use
#   PYTHON_NET    python3 for the NCBI download steps (20, 20b). Default
#                 /usr/bin/python3: the tick_markers env's own python does not
#                 trust this site's TLS-inspecting proxy certificate, and its
#                 bin/ is put first on PATH below, so plain `python3` would
#                 resolve to it.
#   MARKER_ENV    prefix of the conda env holding mafft/iqtree (default: tick_markers)
#   RSCRIPT       Rscript for step 26 (default: the tick_micro conda env)
#   NCBI_API_KEY  raises the Entrez rate limit from 3/s to 10/s
set -u
cd "$(dirname "$0")" || exit 1

PY="${PYTHON:-python3}"
PYNET="${PYTHON_NET:-/usr/bin/python3}"
MARKER_ENV="${MARKER_ENV:-/data/hb-austick/work/miniforge3/envs/tick_markers}"
RS="${RSCRIPT:-/data/hb-austick/work/miniforge3/envs/tick_micro/bin/Rscript}"

if [ -x "$MARKER_ENV/bin/mafft" ]; then
    PATH="$MARKER_ENV/bin:$PATH"; export PATH
else
    echo "WARNING: mafft/iqtree not found at $MARKER_ENV/bin -- steps 23/25 will fail."
    echo "  mamba env create -f env/tick_markers.yml"
fi

step() { echo; echo "=============== $* ==============="; }
die()  { echo "FAILED: $*" >&2; exit 1; }

step "20 fetch mitogenome gene reference probes (cached; Entrez)"
$PYNET code/analysis/20_fetch_mito_refs.py               || die "mito refs"

step "20b download + annotate the published reference mitogenomes (cached; Entrez)"
$PYNET code/analysis/20b_fetch_mito_ref_genomes.py        || die "mito ref genomes"

step "21 find one mitogenome candidate per sample (tblastn/blastn)"
$PY code/analysis/21_find_mito_candidates.py "$@"        || die "mito candidates"

step "22 annotate (MITOS2 if configured, else homology fallback)"
$PY code/analysis/22_annotate_mitogenomes.py "$@"         || die "annotate"

step "23 align each gene separately (codon-aware for PCGs)"
$PY code/analysis/23_align_mito_genes.py                  || die "align genes"

step "24 concatenate into one partitioned supermatrix"
$PY code/analysis/24_concat_partition.py                  || die "concat/partition"

step "25 IQ-TREE: ModelFinder + partitioned ML tree"
$PY code/analysis/25_run_iqtree.py                         || die "iqtree"

step "26 final figure + distance tables (R)"
$RS code/analysis/26_mito_figures.R                        || die "R figures"

step "26b cross-check curated tick identifications against the tree (R)"
$RS code/analysis/26b_mito_vs_metadata.R                   || die "mito vs metadata"

step "26c audit genus-discordant samples: contaminant / mixed / wrong tick (raw BLAST; cached)"
$PY code/analysis/26c_mito_discordance_audit.py            || die "discordance audit"

echo
echo "=============== done ==============="
echo "tree:      results/mito/tree/mito.treefile"
echo "figure:    results/mito/figures/mito_tree.png"
echo "distances: results/mito/tables/mito_distance_{pdistance,patristic}.tsv"
echo "audit:     results/mito/tables/mito_discordance_audit.tsv"
echo "see README_mito.md"
