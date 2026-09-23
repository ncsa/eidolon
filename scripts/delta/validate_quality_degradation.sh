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
# One value, two consumers. TAIL_CUT is `degradation_tail_cut` for the fitter AND the
# collapsed-tail threshold the measurement below reports at; TAIL_WINDOW likewise is
# `degradation_tail_window` and the measurement's window. Letting them drift is how a run at
# TAIL_CUT=30 would present a Q<25 rate as though it said something about a model fitted at
# Q<30. The measurement helper defaults to 25/50 on its own, which is why this was invisible
# at the defaults and only wrong when overridden.
TAIL_CUT="${TAIL_CUT:-25}"
TAIL_WINDOW="${TAIL_WINDOW:-50}"
# The second, stricter rate the measurement reports. Tracks TAIL_CUT five Q down, matching that
# script's own default, so the defaults stay the 25/20 pair #694 quotes.
TAIL_CUT_DEEP="${TAIL_CUT_DEEP:-$((TAIL_CUT - 5))}"
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
echo "cut:      Q<$TAIL_CUT over the last $TAIL_WINDOW bases (fit and measurement both)"
echo "binary:   $EIDOLON"
echo "fastq:    $FASTQ   stride $STRIDE"
echo "work:     $WORK"
echo

# ── 1. strided subsample ────────────────────────────────────────────────────────
SUB="$WORK/subsample.fastq.gz"
if [[ ! -s "$SUB" ]]; then
    echo "[1/5] subsampling: keeping 1 record in $STRIDE ..."
    # Rename on success only: a killed job leaves a truncated but non-empty .gz, which the
    # `-s` cache above would accept. pipefail STAYS ON: awk and gzip both read to EOF, so
    # there is no early consumer and no SIGPIPE here. With it off, a truncated input makes
    # zcat fail while awk and gzip succeed, the pipeline reports 0, and the script caches and
    # fits the partial sample. Measured on a half-truncated fixture: rc 0, 2448.25 reads
    # cached -- a fractional count, so even the final record was cut in half.
    rc=0
    zcat -f -- "$FASTQ" \
      | awk -v s="$STRIDE" 'NR%4==1 { rec++ } ((rec - 1) % s) == 0' \
      | gzip -c > "$SUB.part" || rc=$?
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
    printf 'fastq_file: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\nfit_quality_degradation: %s\ndegradation_tail_cut: %s\ndegradation_tail_window: %s\n' \
        "$SUB" "$WORK/model_$mode.json.gz" "$flag" "$TAIL_CUT" "$TAIL_WINDOW" > "$cfg"
    echo "[2/5] fitting $mode ..."
    # Status first, filter second. Ending the pipeline with `| grep ... || true` hid a failed
    # fit, and because WORK is reused the previous run's model_$mode.json.gz survives it, so
    # step 3 would generate from a stale model and step 5 would print a plausible, invalid
    # comparison. The `|| true` on the grep is still needed: no matching lines is exit 1.
    if ! "$EIDOLON" gen-seq-error-model -c "$cfg" > "$WORK/fit_$mode.log" 2>&1; then
        echo "FATAL: fitting $mode failed. Last 20 lines of $WORK/fit_$mode.log:" >&2
        tail -20 "$WORK/fit_$mode.log" >&2
        exit 1
    fi
    grep -E "Degraded population|Separability|error_rate|read length" "$WORK/fit_$mode.log" || true
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
    # The SAME cut and window the fit used. Reporting a rate at a threshold the model was not
    # fitted at is not a result about that model.
    FASTQ="$2" MAX_READS=0 OUT="$WORK/meas_$1.txt" \
        TAIL_CUT="$TAIL_CUT" TAIL_CUT_DEEP="$TAIL_CUT_DEEP" TAIL_WINDOW="$TAIL_WINDOW" \
        bash "$REPO/scripts/delta/measure_quality_degradation.sh" >/dev/null
    # Match on the cut VALUE rather than a literal Q<25. This is also the check that the
    # override reached the helper: if it were dropped, the helper would report Q<25/Q<20 rows,
    # no row would match, and this aborts. Without the END guard an unmatched row prints as a
    # blank column, which reads as "measured, and small" rather than "not measured".
    awk -v l="$1" -v cut="$TAIL_CUT" -v deep="$TAIL_CUT_DEEP" '
        $1 == "tail" && $2 == "Q<" cut  { a = $5 }
        $1 == "tail" && $2 == "Q<" deep { b = $5 }
        /mean quality  *collapsed/      { h = $4 }
        END {
            if (a == "" || b == "") {
                printf "FATAL: no Q<%s / Q<%s row in the measurement of %s\n", cut, deep, l > "/dev/stderr"
                exit 1
            }
            printf "  %-22s Q<%-3s %-8s Q<%-3s %-8s head-collapsed %s\n", l, cut, a, deep, b, h
        }' "$WORK/meas_$1.txt"
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
