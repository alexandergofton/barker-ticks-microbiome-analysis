#!/usr/bin/env bash
# map_reads_to_contigs.sh
#
# Map Illumina reads back to SPAdes de novo contigs with Bowtie2,
# then produce a per-contig read-count table.
#
# OUTPUTS
#   <outdir>/<sample>.bam            sorted, indexed BAM
#   <outdir>/<sample>.coverage.tsv   per-contig table:
#                                      contig_id | contig_len | mapped_reads |
#                                      pct_of_mapped | tpm | pct_length_adj
#
# USAGE — single sample
#   bash map_reads_to_contigs.sh \
#       --sample  SAMPLEID \
#       --contigs path/to/SAMPLE.contigs.fa[.gz] \
#       --r1      path/to/SAMPLE_R1.fastq.gz \
#       --r2      path/to/SAMPLE_R2.fastq.gz \   # omit for single-end
#       --outdir  mapping/
#
# USAGE — batch (all samples, reads and contigs share a naming prefix)
#   bash map_reads_to_contigs.sh \
#       --batch \
#       --contigs-dir  contigs/ \
#       --reads-dir    reads/ \
#       --r1-suffix    _R1.fastq.gz \   # default: _R1.fastq.gz
#       --r2-suffix    _R2.fastq.gz \   # default: _R2.fastq.gz  (omit for single-end)
#       --outdir       mapping/
#
# DEPENDENCIES  conda install -c bioconda bowtie2 samtools
#
# NOTES
#   • pct_of_mapped (col 4) sums to 100% = share of all mapped READS. This is
#     length-BIASED: a long contig collects more reads at equal abundance. Use it
#     for read-fraction statements ("host DNA = 99% of reads").
#   • tpm (col 5) = Transcripts-Per-Million-style length normalisation:
#         RPK_i    = mapped_i / (len_i / 1000)
#         TPM_i    = RPK_i / (sum(RPK) / 1e6)
#     TPM removes contig-length bias, so it sums to 1,000,000 and reflects relative
#     abundance of the underlying molecules/genomes — the better basis for comparing
#     taxonomic diversity across contigs/samples.
#   • pct_length_adj (col 6) = tpm/10000, i.e. TPM expressed as a percentage that
#     sums to 100 — a length-adjusted analogue of pct_of_mapped.
#   • Contigs with 0 mapped reads are included (all metrics = 0).
#   • Index files are written to a temp directory and removed afterwards.
#   • Compressed contigs (.fa.gz) are handled automatically.

set -euo pipefail

# ─── defaults ─────────────────────────────────────────────────────────────────
THREADS=8
R1_SUFFIX="_R1.fastq.gz"
R2_SUFFIX="_R2.fastq.gz"
OUTDIR="mapping"
BATCH=false
SAMPLE=""
CONTIGS=""
READS_R1=""
READS_R2=""
CONTIGS_DIR=""
READS_DIR=""

# ─── parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)       SAMPLE="$2";      shift 2 ;;
        --contigs)      CONTIGS="$2";     shift 2 ;;
        --r1)           READS_R1="$2";    shift 2 ;;
        --r2)           READS_R2="$2";    shift 2 ;;
        --outdir)       OUTDIR="$2";      shift 2 ;;
        --threads)      THREADS="$2";     shift 2 ;;
        --batch)        BATCH=true;       shift   ;;
        --contigs-dir)  CONTIGS_DIR="$2"; shift 2 ;;
        --reads-dir)    READS_DIR="$2";   shift 2 ;;
        --r1-suffix)    R1_SUFFIX="$2";   shift 2 ;;
        --r2-suffix)    R2_SUFFIX="$2";   shift 2 ;;
        -h|--help)
            sed -n '2,30p' "$0" | grep "^#" | sed 's/^# \?//'
            exit 0 ;;
        *)  echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ─── dependency check ─────────────────────────────────────────────────────────
for prog in bowtie2 bowtie2-build samtools; do
    if ! command -v "$prog" &>/dev/null; then
        echo "ERROR: '$prog' not found on PATH." >&2
        echo "Install with: conda install -c bioconda bowtie2 samtools" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

# ─── core function: map one sample ────────────────────────────────────────────
map_sample() {
    local sample="$1"
    local contigs="$2"
    local r1="$3"
    local r2="${4:-}"   # empty string = single-end

    local bam="${OUTDIR}/${sample}.bam"
    local table="${OUTDIR}/${sample}.coverage.tsv"

    echo "" >&2
    echo "════════════════════════════════════════" >&2
    echo "  Sample : ${sample}" >&2
    echo "  Contigs: ${contigs}" >&2
    echo "  R1     : ${r1}" >&2
    [[ -n "$r2" ]] && echo "  R2     : ${r2}" >&2 || echo "  Mode   : single-end" >&2
    echo "════════════════════════════════════════" >&2

    # Validate inputs
    for f in "$contigs" "$r1"; do
        [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; return 1; }
    done
    [[ -n "$r2" && ! -f "$r2" ]] && { echo "ERROR: R2 not found: $r2" >&2; return 1; }

    # ── 1. Build Bowtie2 index ──────────────────────────────────────────────
    local tmpdir
    tmpdir=$(mktemp -d)
    local idx="${tmpdir}/idx"

    echo "[1/4] Building Bowtie2 index..." >&2

    # Handle gzip-compressed contigs (bowtie2-build supports .gz natively in ≥2.3.5)
    if [[ "$contigs" == *.gz ]]; then
        bowtie2-build --threads "$THREADS" "$contigs" "$idx" 2>/dev/null \
            || { echo "  NOTE: bowtie2-build failed on .gz — decompressing first..." >&2
                 gunzip -c "$contigs" > "${tmpdir}/contigs.fa"
                 bowtie2-build --threads "$THREADS" "${tmpdir}/contigs.fa" "$idx" 2>/dev/null; }
    else
        bowtie2-build --threads "$THREADS" "$contigs" "$idx" 2>/dev/null
    fi

    # ── 2. Map reads ────────────────────────────────────────────────────────
    echo "[2/4] Mapping reads with Bowtie2 (${THREADS} threads)..." >&2

    if [[ -n "$r2" ]]; then
        # Paired-end
        bowtie2 \
            -x "$idx" \
            -1 "$r1" \
            -2 "$r2" \
            -p "$THREADS" \
            --no-unal \
            --end-to-end \
            2>"${OUTDIR}/${sample}.bowtie2.log" \
        | samtools sort -@ "$THREADS" -o "$bam" -
    else
        # Single-end
        bowtie2 \
            -x "$idx" \
            -U "$r1" \
            -p "$THREADS" \
            --no-unal \
            --end-to-end \
            2>"${OUTDIR}/${sample}.bowtie2.log" \
        | samtools sort -@ "$THREADS" -o "$bam" -
    fi

    # ── 3. Index BAM ────────────────────────────────────────────────────────
    echo "[3/4] Indexing BAM..." >&2
    samtools index -@ "$THREADS" "$bam"

    # Report mapping rate from Bowtie2 log
    if [[ -f "${OUTDIR}/${sample}.bowtie2.log" ]]; then
        grep -E "overall alignment rate|reads; of these" "${OUTDIR}/${sample}.bowtie2.log" \
            | sed 's/^/    /' >&2 || true
    fi

    # ── 4. Per-contig read counts → coverage table ──────────────────────────
    echo "[4/4] Generating per-contig coverage table..." >&2

    # samtools idxstats columns: contig_name, length, mapped_reads, unmapped_reads
    # Last row is "*" for unplaced reads — skip it.
    # We compute total mapped across all contigs first, then calculate proportions.

    samtools idxstats "$bam" \
    | awk '
        BEGIN { OFS="\t" }

        # First pass: accumulate mapped reads and reads-per-kilobase (skip * row).
        # rpk = reads normalised by contig length, the basis for TPM.
        $1 != "*" {
            contigs[NR]=$1; lengths[NR]=$2; mapped[NR]=$3; total += $3
            rpk[NR]     = ($2 > 0) ? $3 / ($2 / 1000) : 0   # guard zero-length contigs
            rpk_total  += rpk[NR]
        }

        END {
            # TPM per-million scaling factor: divides each RPK so TPM sums to 1e6.
            scaling = (rpk_total > 0) ? rpk_total / 1000000 : 0

            # Print header
            #   pct_of_mapped  = share of TOTAL READS (length-biased; = "% of reads")
            #   tpm            = length-normalised abundance, sums to 1,000,000
            #   pct_length_adj = tpm as a percentage (sums to 100; length-adjusted "%")
            print "contig_id", "contig_len", "mapped_reads", "pct_of_mapped", "tpm", "pct_length_adj"

            # Print one row per contig (including those with 0 reads)
            for (i=1; i<=length(contigs); i++) {
                pct     = (total   > 0) ? (mapped[i] / total) * 100 : 0
                tpm     = (scaling > 0) ? rpk[i] / scaling          : 0
                pct_len = tpm / 10000   # tpm/1e6*100 → length-adjusted percentage
                printf "%s\t%d\t%d\t%.6f\t%.4f\t%.6f\n", \
                       contigs[i], lengths[i], mapped[i], pct, tpm, pct_len
            }
        }
    ' > "$table"

    # Sanity check: proportions should sum to 100, TPM to 1,000,000
    local n_contigs mapped_total pct_sum tpm_sum pctlen_sum
    n_contigs=$(tail -n +2 "$table" | wc -l | tr -d ' ')
    mapped_total=$(tail -n +2 "$table" | awk '{sum+=$3} END{print sum}')
    pct_sum=$(tail -n +2 "$table" | awk '{sum+=$4} END{printf "%.2f", sum}')
    tpm_sum=$(tail -n +2 "$table" | awk '{sum+=$5} END{printf "%.0f", sum}')
    pctlen_sum=$(tail -n +2 "$table" | awk '{sum+=$6} END{printf "%.2f", sum}')

    echo "" >&2
    echo "  Output BAM    : ${bam}" >&2
    echo "  Coverage table: ${table}" >&2
    echo "  Contigs       : ${n_contigs}" >&2
    echo "  Total mapped  : ${mapped_total} reads" >&2
    echo "  Pct sum check : ${pct_sum}%  (should be ~100.00)" >&2
    echo "  TPM sum check : ${tpm_sum}  (should be ~1000000)" >&2
    echo "  TPM-pct check : ${pctlen_sum}%  (should be ~100.00)" >&2

    # Clean up index
    rm -rf "$tmpdir"
}


# ─── single-sample mode ───────────────────────────────────────────────────────
if [[ "$BATCH" == "false" ]]; then
    [[ -z "$SAMPLE"  ]] && { echo "ERROR: --sample required in single-sample mode" >&2; exit 1; }
    [[ -z "$CONTIGS" ]] && { echo "ERROR: --contigs required" >&2; exit 1; }
    [[ -z "$READS_R1" ]] && { echo "ERROR: --r1 required" >&2; exit 1; }

    map_sample "$SAMPLE" "$CONTIGS" "$READS_R1" "$READS_R2"

    echo "" >&2
    echo "Done." >&2
    exit 0
fi


# ─── batch mode ───────────────────────────────────────────────────────────────
[[ -z "$CONTIGS_DIR" ]] && { echo "ERROR: --contigs-dir required in batch mode" >&2; exit 1; }
[[ -z "$READS_DIR"   ]] && { echo "ERROR: --reads-dir required in batch mode" >&2; exit 1; }

# Discover samples from contig files (support .fa and .fa.gz)
shopt -s nullglob
contig_files=( "${CONTIGS_DIR}"/*.contigs.fa "${CONTIGS_DIR}"/*.contigs.fa.gz )
shopt -u nullglob

if [[ ${#contig_files[@]} -eq 0 ]]; then
    echo "ERROR: no *.contigs.fa or *.contigs.fa.gz found in ${CONTIGS_DIR}" >&2
    exit 1
fi

n_total=${#contig_files[@]}
n_done=0
n_fail=0

echo "Batch mode: ${n_total} samples found in ${CONTIGS_DIR}" >&2

for contig_file in "${contig_files[@]}"; do
    # Derive sample ID by stripping directory and .contigs.fa[.gz] suffix
    fname=$(basename "$contig_file")
    sample="${fname%.contigs.fa.gz}"
    sample="${sample%.contigs.fa}"

    r1="${READS_DIR}/${sample}${R1_SUFFIX}"
    r2="${READS_DIR}/${sample}${R2_SUFFIX}"

    # Check if R1 exists
    if [[ ! -f "$r1" ]]; then
        echo "WARNING: R1 not found for ${sample} (expected: ${r1}) — skipping" >&2
        (( n_fail++ ))
        continue
    fi

    # R2 is optional
    if [[ ! -f "$r2" ]]; then
        echo "  NOTE: no R2 found for ${sample} — running single-end" >&2
        r2=""
    fi

    (( n_done++ ))
    echo "" >&2
    echo "── Sample ${n_done}/${n_total}: ${sample}" >&2

    map_sample "$sample" "$contig_file" "$r1" "$r2" || (( n_fail++ ))
done

echo "" >&2
echo "════════════════════════════════════════" >&2
echo "Batch complete." >&2
echo "  Processed : $((n_done - n_fail)) / ${n_total} samples" >&2
[[ $n_fail -gt 0 ]] && echo "  Failed    : ${n_fail}" >&2
echo "  Output dir: ${OUTDIR}/" >&2
echo "" >&2
echo "Coverage tables: ${OUTDIR}/*.coverage.tsv" >&2
echo "BAM files      : ${OUTDIR}/*.bam" >&2

