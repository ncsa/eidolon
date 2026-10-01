#!/usr/bin/env bash
# Fit a substitution matrix from one library's mate-overlap disagreements (#779), so the
# shipped default can be an average over several libraries rather than one sample's fit.
#
# Only the matrix is wanted. gen-seq-error-model also fits a quality model, which needs a FASTQ,
# so a small one is pulled from the BAM itself; that half of the output is not used.
#
# INPUTS, all required but VCF:
#   BAM   aligned, coordinate-sorted, with MD tags; duplicates marked if possible (they are
#         skipped). Mates must overlap: fragments shorter than the two reads combined.
#   TAG   a short name for this library, used for the output directory
#   VCF   the sample's own variants [optional]. Without it, errors that land exactly on a
#         variant are not masked; that was 0.24% of HG002's disagreements. Germline variants
#         themselves cancel either way, since both mates carry them.
#
# REBUILD FIRST: step 0 refuses a binary without #783's per-cycle reweighting.
#   cd $REPO && git pull && CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc \
#     CARGO_TARGET_DIR=$SCRATCH/cargo-target/eidolon cargo build --release
#
# USAGE
#   BAM=/path/normal.bam TAG=hcc1395_normal sbatch scripts/delta/fit_overlap_matrix.sbatch
# then check the exit state, not just the log:
#   sacct -j <id> --format=JobID,State,ExitCode,Elapsed

#SBATCH --job-name=eidolon-overlapfit
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=overlapfit-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"

BAM="${BAM:?set BAM to the aligned library}"
TAG="${TAG:?set TAG to a short library name}"
VCF="${VCF:-}"
MIN_MAPQ="${MIN_MAPQ:-20}"
N_READS="${N_READS:-200000}"
WORK="${WORK:-$SCRATCH/overlap_fits/$TAG}"
EIDOLON="${EIDOLON:-$SCRATCH/cargo-target/eidolon/release/eidolon}"
# #783: the per-cycle reweighting. The fit depends on it.
REQUIRED_COMMIT="${REQUIRED_COMMIT:-6f3c6cd}"

module load samtools/1.22-cce19.0.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH after module load" >&2; exit 1; }
[[ -s "$BAM" ]] || { echo "FATAL: missing BAM: $BAM" >&2; exit 1; }
[[ -z "$VCF" || -s "$VCF" ]] || { echo "FATAL: missing VCF: $VCF" >&2; exit 1; }
mkdir -p "$WORK"

ver="$("$EIDOLON" --version 2>/dev/null | tr -d '\r\n')" || { echo "FATAL: no eidolon at $EIDOLON" >&2; exit 1; }
bin_sha="${ver##*+}"
bin_sha="${bin_sha%-dirty}"
echo "=== fit_overlap_matrix: $TAG ==="
echo "repo:    $REPO ($(git -C "$REPO" rev-parse --short HEAD))"
echo "binary:  $ver"
echo "bam:     $BAM"
echo "vcf:     ${VCF:-none}"
if ! git -C "$REPO" merge-base --is-ancestor "$REQUIRED_COMMIT" "$bin_sha" 2>/dev/null; then
    echo "FATAL: binary $bin_sha does not contain $REQUIRED_COMMIT (#783). Rebuild (see header)." >&2
    exit 1
fi

# The quality half's input: the first N_READS primary reads. `head` closes the pipe early, so
# samtools takes SIGPIPE; pipefail comes off for exactly this pipeline.
FQ="$WORK/reads.fastq.gz"
if [[ ! -s "$FQ" ]]; then
    set +o pipefail
    samtools fastq -F 0x900 "$BAM" 2>/dev/null | head -n $((N_READS * 4)) | gzip -c > "$FQ"
    set -o pipefail
fi
n=$(zcat -f -- "$FQ" | awk 'END { print NR / 4 }')
[[ "${n%.*}" -ge 1000 ]] || { echo "FATAL: only $n reads extracted from $BAM" >&2; exit 1; }

MODEL="$WORK/seqerr_$TAG.json.gz"
CFG="$WORK/seqerr.yml"
printf 'fastq_file: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\nbam_file: %s\nbam_method: overlap\nbam_min_mapq: %s\n' \
    "$FQ" "$MODEL" "$BAM" "$MIN_MAPQ" > "$CFG"
[[ -n "$VCF" ]] && printf 'known_variants_vcf: %s\n' "$VCF" >> "$CFG"

rc=0
"$EIDOLON" gen-seq-error-model -c "$CFG" > "$WORK/seqerr.log" 2>&1 || rc=$?
sed 's/\x1b\[[0-9;]*m//g' "$WORK/seqerr.log" | grep -E "Overlap fit|Masking|WARN|ERROR" || true
if [[ "$rc" -ne 0 || ! -s "$MODEL" ]]; then
    echo "FATAL: fit failed (rc $rc). Last 20 lines:" >&2
    tail -20 "$WORK/seqerr.log" >&2
    exit 1
fi

archive_run "overlap_fits/$TAG" "$WORK" \
    "seqerr_$TAG.json.gz" "seqerr_$TAG.json.gz.overlap_bins.tsv" seqerr.yml seqerr.log
echo "DONE."
