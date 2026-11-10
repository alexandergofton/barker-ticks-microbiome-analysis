#!/bin/bash
# Blood-meal host identification: which vertebrate(s) each tick fed on, from
# the host DNA in its own assembly, compared with the host recorded in the
# metadata. See README_bloodmeal.md.
#
#   bash run_bloodmeal.sh
#
# Tiered evidence per sample:
#   1  vertebrate mitochondrial genome / segment contigs
#   2  only when tier 1 is empty: vertebrate 18S rRNA contigs + short COI fragments
#
# Runs independently of run_all.sh / run_mito.sh. Needs the main pipeline's
# blast/, lca/, mapping/, contigs/ and results/taxonomy/taxon_lineage.tsv
# (step 03), and metadata/sample_metadata.tsv (step 5b).
#
# Step 30 streams every raw nt BLAST file once (~480 GB, ~2-3 h on this node);
# its per-sample cache makes every later run take minutes.
#
# Optional env:
#   PYTHON        python3 to use
#   PYTHON_NET    python3 for the network steps (31, 32, 33b). Default
#                 /usr/bin/python3 -- see run_mito.sh for why.
#   MARKER_ENV    conda env holding blastn/makeblastdb/seqkit (default: tick_markers)
#   RSCRIPT       Rscript for step 35 (default: the tick_micro conda env)
#   NCBI_API_KEY  raises the Entrez rate limit from 3/s to 10/s
set -u
cd "$(dirname "$0")" || exit 1

PY="${PYTHON:-python3}"
PYNET="${PYTHON_NET:-/usr/bin/python3}"
MARKER_ENV="${MARKER_ENV:-/data/hb-austick/work/miniforge3/envs/tick_markers}"
RS="${RSCRIPT:-/data/hb-austick/work/miniforge3/envs/tick_micro/bin/Rscript}"

if [ -x "$MARKER_ENV/bin/blastn" ]; then
    PATH="$MARKER_ENV/bin:$PATH"; export PATH
else
    echo "WARNING: blastn/seqkit not found at $MARKER_ENV/bin -- steps 31/33/33c will fail."
fi

step() { echo; echo "=============== $* ==============="; }
die()  { echo "FAILED: $*" >&2; exit 1; }

step "30 scan raw nt BLAST for vertebrate HSPs (cached per sample)"
$PY code/analysis/30_scan_vertebrate_hits.py              || die "scan"

step "31 RefSeq vertebrate mitogenomes + gene table (cached; FTP)"
$PYNET code/analysis/31_fetch_vert_mito_refs.py           || die "refseq mito"

step "32 resolve metadata hosts to NCBI taxa (Entrez counts cached)"
$PYNET code/analysis/32_resolve_metadata_hosts.py         || die "metadata hosts"

step "33 classify vertebrate contigs (tiers, LCA, identity cap, RefSeq genes)"
$PY code/analysis/33_classify_vertebrate_contigs.py       || die "classify"

step "33b fetch nt mito records of candidate host families (cached; Entrez)"
$PYNET code/analysis/33b_fetch_host_family_refs.py        || die "family refs"

step "33c re-name host mito contigs by identity vote"
$PY code/analysis/33c_refine_host_contigs.py              || die "refine"

step "34 per-sample call + comparison with metadata"
$PY code/analysis/34_call_bloodmeal.py                    || die "call"

step "35 figures"
$RS code/analysis/35_bloodmeal_figures.R                  || die "figures"

echo
echo "=============== done ==============="
echo "outputs under results/bloodmeal/"
