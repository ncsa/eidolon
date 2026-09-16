#!/usr/bin/env bash
# Fit the two-population quality model (#694) on a real library and measure whether the reads
# it generates carry the degraded population the library has.
#
# WHY THIS EXISTS. The two-population model is validated only on a synthetic fixture whose
# degraded population was constructed: 9.64% generated against 9.84% planted, where a
# single-population fit reads 1.80%. That shows the mechanism works end to end. It does not
# show it reproduces a real library, and `fit_quality_degradation` stays off by default until
# it does. This job is what would justify flipping that default.
#
# WHAT IT DOES. Subsamples the library, fits it twice (with and without the degraded
# population), generates reads from each model, and measures the collapsed-tail rate of all
# three: the real subsample, and the two simulated arms. The control arm matters as much as the
# test arm -- without it, a number from the test arm has nothing to sit against.
#
# REBUILD THE BINARY FIRST. $EIDOLON is a build artifact in $SCRATCH, not the repo
# checkout: pulling the branch does not rebuild it, and gen-seq-error-model reads its
# config by key lookup, so an older binary ignores fit_quality_degradation rather than
# rejecting it -- both arms then silently become the control.
#
# WHY IT SUBSAMPLES HERE, WHICH IS INTERIM. `gen-seq-error-model`'s own `max_reads` takes the
# FIRST N records, and Illumina FASTQs are written in flowcell order: every one of the first
# 200,000 records in HG002 R1 is lane 1, tile 1101. Capping the fit without spanning the file
# fits one corner of one lane. That is a defect in the fitter, tracked in #721, not a property
# this script should be routing around -- the stride below goes away once the fitter can sample
# for itself, and this whole step with it.
#
# USAGE
#   sbatch scripts/delta/validate_quality_degradation.sh
#   FASTQ=/path/hg002_R2.fastq.gz NAME=r2 sbatch scripts/delta/validate_quality_degradation.sh
#
# Defaults to HG002 R1. Run R2 as well: it is 2.3x worse on collapsed tails and is the fit #695
# has been waiting on.

#SBATCH --job-name=eidolon-qualdegfit
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=qualdegfit-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
FASTQ="${FASTQ:-${DATA_DIR:-/work/nvme/bhrd/jallen17/hg002}/hg002_R1.fastq.gz}"
# Delta does not export $SCRATCH to jobs and this runs under `set -u`: lib_report.sh
# resolves it (and RESULTS_DIR) the way every other job here does. Source before any use.
source "$REPO/scripts/delta/lib_report.sh"

NAME="${NAME:-r1}"
STRIDE="${STRIDE:-65}"          # ~2M reads from a 130M-record library
READ_LEN="${READ_LEN:-250}"
COVERAGE="${COVERAGE:-4}"
TAIL_CUT="${TAIL_CUT:-25}"
WORK="${WORK:-$SCRATCH/qualdeg_$NAME}"
EIDOLON="${EIDOLON:-$SCRATCH/cargo-target/eidolon/release/eidolon}"
REFERENCE="${REFERENCE:-$REPO/eidolon/test_data/references/ecoli.fa}"

[[ -x "$EIDOLON" ]] || { echo "FATAL: no eidolon binary at $EIDOLON. Build it first:" >&2
    echo "  cd $REPO && bash scripts/delta/setup.sh && cargo build --release" >&2; exit 1; }
[[ -f "$FASTQ" ]]   || { echo "FATAL: no FASTQ at $FASTQ" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "FATAL: no reference at $REFERENCE" >&2; exit 1; }

mkdir -p "$WORK"
echo "=== validate_quality_degradation ($NAME) ==="
echo "repo:     $REPO ($(cd "$REPO" && git rev-parse --short HEAD))"
echo "binary:   $EIDOLON"
echo "fastq:    $FASTQ   stride $STRIDE"
echo "work:     $WORK"
echo

# ── 1. strided subsample ────────────────────────────────────────────────────────
SUB="$WORK/subsample.fastq.gz"
if [[ ! -s "$SUB" ]]; then
    echo "[1/5] subsampling: keeping 1 record in $STRIDE ..."
    # Rename on success only: a killed job leaves a truncated but non-empty .gz, which the
    # `-s` cache above would accept. pipefail off around the pipe: see CLAUDE.md on SIGPIPE.
    set +o pipefail
    rc=0
    zcat -f -- "$FASTQ" \
      | awk -v s="$STRIDE" 'NR%4==1 { rec++ } ((rec - 1) % s) == 0' \
      | gzip -c > "$SUB.part" || rc=$?
    set -o pipefail
    [[ "$rc" -eq 0 ]] || { echo "FATAL: subsample failed (rc $rc)" >&2; rm -f "$SUB.part"; exit 1; }
    mv -f "$SUB.part" "$SUB"
fi
n_sub=$(zcat -f -- "$SUB" | awk 'END { print NR/4 }')
[[ "${n_sub%.*}" -gt 1000 ]] || { echo "FATAL: subsample holds only $n_sub reads" >&2; exit 1; }
echo "      $n_sub reads"

# ── 2. fit both models ──────────────────────────────────────────────────────────
for mode in plain degraded; do
    flag=false; [[ "$mode" == degraded ]] && flag=true
    cfg="$WORK/fit_$mode.yml"
    printf 'fastq_file: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\nfit_quality_degradation: %s\ndegradation_tail_cut: %s\n' \
        "$SUB" "$WORK/model_$mode.json.gz" "$flag" "$TAIL_CUT" > "$cfg"
    echo "[2/5] fitting $mode ..."
    "$EIDOLON" gen-seq-error-model -c "$cfg" 2>&1 | grep -E "Degraded population|Separability|error_rate|read length" || true
done

# ── 3. generate from each ───────────────────────────────────────────────────────
for mode in plain degraded; do
    cfg="$WORK/gen_$mode.yml"
    printf 'reference: %s\noutput_dir: %s\noutput_filename: gen_%s\noverwrite_output: true\nread_len: %s\ncoverage: %s\nproduce_fastq: true\nproduce_bam: false\nproduce_vcf: false\npaired_ended: false\nrng_seed: qualdeg\nsequence_error_model: %s\n' \
        "$REFERENCE" "$WORK" "$mode" "$READ_LEN" "$COVERAGE" "$WORK/model_$mode.json.gz" > "$cfg"
    echo "[3/5] generating from $mode ..."
    "$EIDOLON" gen-reads -c "$cfg" >/dev/null 2>&1
done

# ── 4. measure all three arms ───────────────────────────────────────────────────
echo "[4/5] measuring ..."
measure() {  # <label> <fastq>
    FASTQ="$2" MAX_READS=0 OUT="$WORK/meas_$1.txt" \
        bash "$REPO/scripts/delta/measure_quality_degradation.sh" >/dev/null
    awk -v l="$1" '/tail Q<25/{a=$5} /tail Q<20/{b=$5} /mean quality  *collapsed/{c=$4}
        END { printf "  %-22s Q<25 %-8s Q<20 %-8s head-collapsed %s\n", l, a, b, c }' "$WORK/meas_$1.txt"
}

# ── 5. the comparison ───────────────────────────────────────────────────────────
echo
echo "[5/5] === RESULT ($NAME) ==="
measure "real (subsample)"      "$SUB"
measure "sim, one population"   "$WORK/gen_plain_r1.fastq.gz"
measure "sim, two populations"  "$WORK/gen_degraded_r1.fastq.gz"
echo
echo "The two-population arm should approach the real row. The one-population arm is the"
echo "control and is expected to read far below it; without that row the test arm means little."
echo
echo "Outputs in $WORK. Paste the RESULT block back."
