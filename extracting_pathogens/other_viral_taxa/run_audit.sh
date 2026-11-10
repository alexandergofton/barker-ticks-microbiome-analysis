#!/usr/bin/env bash
# Master script for the "are the other viral contigs real?" audit.
#
#   cd extracting_pathogens/other_viral_taxa
#   nohup ./run_audit.sh > audit.log 2>&1 &
#   tail -f audit.log                      # progress; safe to disconnect
#
# Safe to re-run after an interruption: every expensive step is cached or
# resumable, so a re-run picks up where the last one stopped.
#   extraction  (~3 h, 460 GB scan)  resumes per sample (summary.tsv)
#   02 NCBI fetch                    fetches only accessions not yet cached
#   03 protein DB                    skipped if protdb.dmnd is newer than its inputs
#   04 DIAMOND  (~1 h)               cached in are_these_real/diamond.tsv
#   08 remote NCBI BLAST             cached in are_these_real/remote_blast.json
# The other steps take minutes and always re-run. To force a step, delete its
# cache file. Finally the report is rendered (are_these_real.html).
#
# Memory: this node has ~3 GB. Steps run one at a time; DIAMOND peaks at
# ~0.5 GB with the settings in scripts/04. Do not run run_all.sh at the same
# time (step 03/R steps of the main pipeline need ~1-1.5 GB).
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
export PATH=/data/hb-austick/work/miniforge3/bin:/usr/bin:$PATH
CONDA=/data/hb-austick/work/miniforge3/bin/conda
W=are_these_real

# one run at a time
exec 9> .run_audit.lock
flock -n 9 || { echo "another run_audit.sh is already running"; exit 1; }

step() { echo; echo "$(date '+%F %T') ==== $*"; }

step "extraction (resumable)"
done_n=$(( $(wc -l < summary.tsv 2>/dev/null || echo 1) - 1 ))
echo "samples already extracted: $done_n"
./extract_other_viral_taxa.sh >> extract.log 2>&1
grep -E "NOTE|WARN" extract.log | sort -u || true

step "01 contig evidence"
python3 scripts/01_contig_evidence.py

step "02 fetch hit records (cached)"
python3 scripts/02_fetch_hit_records.py 2>&1 | grep -v '^    [0-9]' || true

step "03 protein database"
if [[ -s $W/protdb.dmnd && $W/protdb.dmnd -nt $W/ncbi/hitrec_prot.fa \
      && $W/protdb.dmnd -nt scripts/03_build_protein_db.py ]]; then
  echo "up to date, skipped"
else
  python3 scripts/03_build_protein_db.py
  rm -f $W/diamond.tsv        # a new database invalidates the DIAMOND cache
fi

step "04 DIAMOND blastx + protein evidence (cached)"
[[ -s $W/diamond.tsv ]] && echo "diamond.tsv cached" || echo "running DIAMOND (~1 h)"
rm -f $W/diamond.tsv.part
python3 scripts/04_protein_evidence.py

step "05 reference hotspots";      python3 scripts/05_reference_hotspots.py
step "06 classify (pass 1)";       python3 scripts/06_classify.py
step "07 elements";                python3 scripts/07_elements.py
step "06 classify (pass 2)";       python3 scripts/06_classify.py
step "08 remote BLAST (cached)";   python3 scripts/08_remote_blast_elements.py
step "09 reading frames";          python3 scripts/09_reading_frames.py

step "render are_these_real.qmd"
# own Quarto cache: the shared ~/.cache/quarto/sass/sass.kv failed to open (2026-10-13)
XDG_CACHE_HOME=$PWD/.quarto_cache $CONDA run -n tick_report --cwd "$PWD" quarto render are_these_real.qmd

step "done"
cut -f12 $W/contig_verdicts.tsv | tail -n +2 | sort | uniq -c | sort -rn
