#!/bin/bash
# SLURM job: where do real indels sit, and are they VARIANTS or SLIPPAGE ERRORS?
#
# WHY THIS EXISTS. The realism panel measures cand_per_mb 57.5 on real HCC1395 normal reads
# and exactly 0 on simulated ones (job 21639409), and the 318 real candidate sites sit in
# sequence with a mean longest homopolymer of 9.71 bp against 4.30 at controls 5 kb away.
# Clipping happens where the aligner cannot confidently place a gap, which is inside a
# homopolymer or STR, so the boundaries agree across reads and a caller clusters them.
#
# THE FORK THIS SETTLES. Indels land in homopolymers by two different mechanisms, and they
# point at different code:
#
#   * replication slippage in the GENOME -> a germline or somatic indel VARIANT, present in
#     roughly half the reads at its locus. Fixing that is issue #378, indel placement
#     conditioned on sequence context.
#   * slippage during PCR or sequencing -> an ERROR, present in a handful of reads. Fixing
#     that is the sequencing error model, which is also context-neutral, and #378 would be
#     the wrong thing to build.
#
# Whichever of the `high` / `low` columns carries the enrichment is the answer. Do not
# implement #378 before reading it.
#
# Prereqs: a real BAM, the reference it was aligned to, and a regions BED -- reuse the one a
# realism_panel run wrote ($OUTDIR/regions.bed) so the two measurements cover the same loci.
#
# Usage:
#   REAL_BAM=$SCRATCH/neat_data/hcc1395_chr20_21_22/normal.bam \
#   REFERENCE=$SCRATCH/neat_data/GRCh38_chr20_21_22.fa \
#   REGIONS=$SCRATCH/realism_21639409/regions.bed \
#     sbatch scripts/delta/indel_context.sbatch
#
#SBATCH --job-name=eidolon-indelctx
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
exec 2>&1
echo "[selftest] stderr reaches this log (if this line is missing, diagnostics are being lost)"

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

# samtools comes from a MODULE on Delta, not from the bioinf conda env -- see
# sv_pipeline.sbatch and stage_hcc1395.sh, which both load it explicitly. Activating bioinf
# alone left `samtools: command not found` at the first read (job 21656030), a minute into
# a job that had already printed a clean header.
module load samtools/1.22-cce19.0.0 2>/dev/null || true
module load htslib/1.22-gcc13.3.1  2>/dev/null || true
setup_conda
conda_activate bioinf

# Preflight, so a missing tool fails in a second with a clear message rather than partway
# through the read loop. A job that inherits its tools from whatever shell submitted it is
# a job that works until someone submits it from a clean one.
missing=""
for t in samtools awk; do command -v "$t" >/dev/null 2>&1 || missing="$missing $t"; done
if [[ -n "$missing" ]]; then
    echo "ERROR: required tool(s) not found:$missing" >&2
    echo "       samtools is a module on Delta:  module load samtools/1.22-cce19.0.0" >&2
    echo "       This job loads it itself; if it is still missing the module name has" >&2
    echo "       changed or this is not Delta." >&2
    exit 1
fi

REAL_BAM="${REAL_BAM:?REAL_BAM is required}"
REFERENCE="${REFERENCE:?REFERENCE is required}"
REGIONS="${REGIONS:?REGIONS is required -- reuse the regions.bed a realism_panel run wrote}"
OUTDIR="${OUTDIR:-$SCRATCH/indelctx_${SLURM_JOB_ID:-local}}"
# Optional: a realism_panel --dump-candidates TSV. Given one, the job also reports whether
# the clip boundaries have a VARIANT near them or nothing at all -- see
# indel_context_join.awk for why that distinction picks the fix.
CANDIDATES="${CANDIDATES:-}"
JOIN_WIN="${JOIN_WIN:-25}"        # how far from a clip boundary an indel still counts
MAX_RUN="${MAX_RUN:-10}"          # run lengths above this are pooled into one bin
HIGH_FRAC="${HIGH_FRAC:-0.25}"    # >= this share of local depth reads as a variant
LOW_FRAC="${LOW_FRAC:-0.10}"      # <  this share reads as a slippage error

for f in "$REAL_BAM" "$REFERENCE" "$REGIONS"; do
    [[ -s "$f" ]] || { echo "ERROR: not found or empty: $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

echo "=== eidolon indel context: ${SLURM_JOB_ID:-local} ==="
echo "BAM:        $REAL_BAM"
echo "Reference:  $REFERENCE"
echo "Regions:    $REGIONS  ($(wc -l < "$REGIONS") loci)"
echo "Bins:       run length 1..$MAX_RUN (pooled above)"
echo "Support:    high >= $HIGH_FRAC of depth, low < $LOW_FRAC"

# ── 1. ONE pass over the reads ──────────────────────────────────────────────
# Cached as (contig, pos, cigar) rather than re-queried. The previous version called
# `samtools view` for the indels and `samtools depth` for the denominator, reading every
# read twice; on 10 x 400 kb at 46x that did not finish on a login node.
CACHE="$OUTDIR/reads.tsv"
echo; echo "[1/4] reading the BAM once..."
: > "$CACHE"
while read -r c s e _; do
    [[ -n "${c:-}" ]] || continue
    samtools view -F 0x904 -q 11 "$REAL_BAM" "$c:$((s+1))-$e" \
      | awk -v OFS='\t' '{ print $3, $4, $6 }' >> "$CACHE"
done < "$REGIONS"
NREADS=$(wc -l < "$CACHE")
[[ "$NREADS" -gt 0 ]] || { echo "ERROR: no reads passed the filter over these regions" >&2; exit 1; }
echo "      $NREADS reads cached"

# ── 2. indel positions, from the cached CIGARs ──────────────────────────────
echo "[2/4] extracting indel operations..."
awk -f "$REPO_ROOT/scripts/delta/indel_context_extract.awk" "$CACHE" \
  | sort -k1,1 -k2,2n > "$OUTDIR/indels.tsv"
NIND=$(wc -l < "$OUTDIR/indels.tsv")
[[ "$NIND" -gt 0 ]] || { echo "ERROR: no indel operations found -- nothing to measure" >&2; exit 1; }
echo "      $NIND distinct indel positions"

# ── 3. depth AT THOSE POSITIONS ONLY, from the same cache ───────────────────
# A depth array over every reference base would be ~4M awk entries; only the indel
# positions need a denominator.
echo "[3/4] depth at indel positions..."
awk -v OFS='\t' '
    FNR == NR { want[$1 SUBSEP $2] = 1; next }
    {
        pos = $2; cig = $3; n = ""
        span = 0
        for (i = 1; i <= length(cig); i++) {
            ch = substr(cig, i, 1)
            if (ch ~ /[0-9]/) { n = n ch; continue }
            len = n + 0; n = ""
            if (ch ~ /[MDN=X]/) span += len
        }
        for (p = pos; p < pos + span; p++) if (($1 SUBSEP p) in want) d[$1 SUBSEP p]++
    }
    END { for (k in d) { split(k, a, SUBSEP); print a[1], a[2], d[k] } }' \
    "$OUTDIR/indels.tsv" "$CACHE" > "$OUTDIR/depth.tsv"
echo "      $(wc -l < "$OUTDIR/depth.tsv") positions with depth"

# ── 4. reference context, one faidx per region, exact background ────────────
echo "[4/4] reference homopolymer context..."
: > "$OUTDIR/ctx.tsv"; : > "$OUTDIR/bg.tsv"
while read -r c s e _; do
    [[ -n "${c:-}" ]] || continue
    samtools faidx "$REFERENCE" "$c:$((s+1))-$e" | tail -n +2 | tr -d '\n' | tr a-z A-Z \
      | awk -v c="$c" -v off="$((s+1))" -v OFS='\t' \
            -v ctx="$OUTDIR/ctx.tsv" -v bg="$OUTDIR/bg.tsv" '
        FNR == NR { if ($1 == c) want[$2 + 0] = 1; next }
        {
            L = length($0); i = 1
            while (i <= L) {
                b = substr($0, i, 1); j = i + 1
                while (j <= L && substr($0, j, 1) == b) j++
                R = j - i
                # N runs are unsequenced reference, not homopolymers, and belong in
                # neither the observed nor the background side.
                if (b != "N") {
                    bgh[R] += R
                    for (p = i; p < j; p++) { gp = off + p - 1; if (gp in want) print c, gp, R >> ctx }
                }
                i = j
            }
            for (R in bgh) print R, bgh[R] >> bg
        }' "$OUTDIR/indels.tsv" -
done < "$REGIONS"

awk -v mx="$MAX_RUN" -v hf="$HIGH_FRAC" -v lf="$LOW_FRAC" \
    -f "$REPO_ROOT/scripts/delta/indel_context_summarise.awk" \
    "$OUTDIR/indels.tsv" "$OUTDIR/depth.tsv" "$OUTDIR/ctx.tsv" "$OUTDIR/bg.tsv" \
  | tee "$OUTDIR/summary.txt"

if [[ -n "$CANDIDATES" ]]; then
    if [[ -s "$CANDIDATES" ]]; then
        cp -f "$CANDIDATES" "$OUTDIR/candidates.tsv"
        awk -v win="$JOIN_WIN" -v hf="$HIGH_FRAC" -v lf="$LOW_FRAC" \
            -f "$REPO_ROOT/scripts/delta/indel_context_join.awk" \
            "$OUTDIR/indels.tsv" "$OUTDIR/depth.tsv" "$OUTDIR/candidates.tsv" \
          | tee -a "$OUTDIR/summary.txt"
    else
        # Refusing rather than skipping: a run asked to make this comparison and silently
        # not making it looks identical to one that was never asked.
        echo "ERROR: CANDIDATES is set to '$CANDIDATES', which is missing or empty." >&2
        exit 1
    fi
fi

rm -f "$CACHE"   # the read cache is large and reproducible; the summary is not
archive_run "indelctx" "$OUTDIR" "summary.txt" "indels.tsv" || true
