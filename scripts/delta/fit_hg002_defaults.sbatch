#!/usr/bin/env bash
# Fit the germline defaults that can be fitted from GIAB HG002, so the shipped set describes
# ONE library rather than the best available source per component (#752).
#
#   step  tool                    from                   replaces
#   1     gen-mut-model           GIAB v4.2.1 truth VCF  NEAT2 NA12878 mutation/indel/trinuc
#   2     gen-frag-length-model   HG002 BAM              HCC1395 NovaSeq fragment model
#   3     gen-gc-bias-model       HG002 BAM              nothing: shipped as an OPTIONAL model;
#                                                        the uniform default stays
#   4     gen-seq-error-model     HG002 BAM + the reads  NEAT2 substitution matrix
#                                  job 22233888 fitted
#
# Step 4 refits the shipped quality model from the same subsample as job 22233888, adding only
# `bam_file:`. Same reads and settings, so its error_rate and degraded fractions must come back
# exactly as shipped: 0.003774, R1 11.02%, R2 25.98%. Anything else means the inputs differ,
# and the model is not the one it claims to extend.
#
# What this job does NOT do is judge the fitted files. They get vetted on the workstation when
# they land in model_data/, as tests CI keeps running.
#
# REBUILD FIRST. The binary is a build artifact, not the checkout; step 0 refuses one that
# predates #765, whose gen-mut-model fixes (empty trinucleotide rows, all-skipped SNPs) this
# fit depends on:
#   cd $REPO && git pull && CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc \
#     CARGO_TARGET_DIR=$SCRATCH/cargo-target/eidolon cargo build --release
#
# INPUTS, all staged already: the BAM on /work/nvme (the 21x alignment of the same
# NIST_Illumina_2x250bps library), the pair subsample by fit_hg002_pair.sbatch (job 22233888),
# and the truth VCF and its BED by fetch_validation_data.sh.
#
# USAGE
#   sbatch scripts/delta/fit_hg002_defaults.sbatch
# then check the exit state, not just the log:
#   sacct -j <id> --format=JobID,State,ExitCode,Elapsed
#
# Each step skips when its output exists, so a resubmit resumes. Delete a model to refit it.

#SBATCH --job-name=eidolon-hg002defaults
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=hg002defaults-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"

DATA_DIR="${DATA_DIR:-$SCRATCH/neat_data/hg002}"
# The aligned library lives on nvme, not beside the VCF (#752).
BAM="${BAM:-/work/nvme/bhrd/jallen17/hg002/hg002.bam}"
VCF="${VCF:-$DATA_DIR/HG002_GRCh38_v4.2.1_benchmark.vcf.gz}"
BED="${BED:-$DATA_DIR/HG002_GRCh38_v4.2.1_highconf.bed}"
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
PAIR_WORK="${PAIR_WORK:-$SCRATCH/hg002_pair}"
WORK="${WORK:-$SCRATCH/hg002_defaults}"
EIDOLON="${EIDOLON:-$SCRATCH/cargo-target/eidolon/release/eidolon}"
# The first #765 commit. The fit needs the binary to contain it.
REQUIRED_COMMIT="${REQUIRED_COMMIT:-bcfc623}"

# samtools is for the fragment validator, which runs AFTER the fragment fit, so check it here
# rather than find out hours in. Same module the other Delta scripts load.
module load samtools/1.22-cce19.0.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH after module load" >&2; exit 1; }

for f in "$BAM" "$VCF" "$BED" "$REFERENCE" "$PAIR_WORK/sub_r1.fastq.gz" "$PAIR_WORK/sub_r2.fastq.gz"; do
    [[ -s "$f" ]] || { echo "FATAL: missing input: $f" >&2; exit 1; }
done
mkdir -p "$WORK"

# ── 0. the binary contains the fixes this fit depends on ────────────────────────
ver="$("$EIDOLON" --version 2>/dev/null | tr -d '\r\n')" || { echo "FATAL: no eidolon at $EIDOLON" >&2; exit 1; }
bin_sha="${ver##*+}"
bin_sha="${bin_sha%-dirty}"
echo "=== fit_hg002_defaults ==="
echo "repo:    $REPO ($(git -C "$REPO" rev-parse --short HEAD))"
echo "binary:  $ver"
if ! git -C "$REPO" merge-base --is-ancestor "$REQUIRED_COMMIT" "$bin_sha" 2>/dev/null; then
    echo "FATAL: binary $bin_sha does not contain $REQUIRED_COMMIT (#765). Rebuild (see header)." >&2
    exit 1
fi
echo

# Runs one builder, keeping its log, and fails the job loudly if the builder does.
fit() {  # <label> <subcommand> <config> <output>
    local label="$1" sub="$2" cfg="$3" out="$4"
    if [[ -s "$out" ]]; then
        echo "[$label] $out exists; skipping. Delete it to refit."
        return 0
    fi
    echo "[$label] $sub ..."
    local rc=0
    "$EIDOLON" "$sub" -c "$cfg" > "$WORK/$label.log" 2>&1 || rc=$?
    if [[ "$rc" -ne 0 || ! -s "$out" ]]; then
        echo "FATAL: $label failed (rc $rc). Last 20 lines:" >&2
        tail -20 "$WORK/$label.log" >&2
        exit 1
    fi
    sed 's/\x1b\[[0-9;]*m//g' "$WORK/$label.log" | grep -E "WARN|ERROR|INFO.*(rate|SNP|rows|reads|pairs|mismatch)" | tail -15 || true
    echo
}

# ── 1. mutation model from the truth VCF ────────────────────────────────────────
# The BED is required, not optional: a benchmark VCF only calls variants inside its confident
# regions, so dividing by the whole genome would understate the rate.
printf 'reference: %s\nvcf_file: %s\nbed_file: %s\noutput_file: %s\noverwrite_output: true\n' \
    "$REFERENCE" "$VCF" "$BED" "$WORK/mut_hg002.json.gz" > "$WORK/mut.yml"
fit mut gen-mut-model "$WORK/mut.yml" "$WORK/mut_hg002.json.gz"

# ── 2. fragment length, same method as the model it replaces ────────────────────
printf 'input_file: %s\noutput_file: %s\noverwrite_output: true\ndistribution: discrete\nmin_reads: 100\n' \
    "$BAM" "$WORK/frag_hg002.json.gz" > "$WORK/frag.yml"
fit frag gen-frag-length-model "$WORK/frag.yml" "$WORK/frag_hg002.json.gz"
# Known-answer check against the BAM itself, measured with samtools rather than eidolon.
echo "[frag] validating against the BAM ..."
bash "$REPO/scripts/delta/validate_frag_model.sh" "$BAM" "$WORK/frag_hg002.json.gz" 2>&1 | tee "$WORK/frag_validate.log"
echo

# ── 3. GC bias, an optional model; the shipped default stays uniform ────────────
printf 'reference: %s\nbam_file: %s\nmin_mapq: 0\nbed_file: .\noutput_file: %s\noverwrite_output: true\nwindow_size: 100\nwindow_stride: 100\nmin_windows_per_bin: 10\n' \
    "$REFERENCE" "$BAM" "$WORK/gc_hg002.json.gz" > "$WORK/gc.yml"
fit gc gen-gc-bias-model "$WORK/gc.yml" "$WORK/gc_hg002.json.gz"

# ── 4. the shipped quality model plus a BAM-fitted substitution matrix ──────────
# Settings are job 22233888's fit_pair.yml, with bam_file added.
printf 'fastq_file: %s\nfastq_file_r2: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\nfit_quality_degradation: true\ndegradation_tail_cut: 25\ndegradation_tail_window: 50\nbam_file: %s\n' \
    "$PAIR_WORK/sub_r1.fastq.gz" "$PAIR_WORK/sub_r2.fastq.gz" "$WORK/seqerr_hg002.json.gz" "$BAM" > "$WORK/seqerr.yml"
fit seqerr gen-seq-error-model "$WORK/seqerr.yml" "$WORK/seqerr_hg002.json.gz"
echo "[seqerr] must match the shipped model: error_rate 0.003774, degraded R1 11.02%, R2 25.98%"
sed 's/\x1b\[[0-9;]*m//g' "$WORK/seqerr.log" | grep -E "error_rate|[Dd]egraded population|mismatches" || true
echo

archive_run hg002_defaults "$WORK" \
    mut_hg002.json.gz frag_hg002.json.gz gc_hg002.json.gz seqerr_hg002.json.gz \
    mut.yml frag.yml gc.yml seqerr.yml mut.log frag.log frag_validate.log gc.log seqerr.log
echo "DONE."
