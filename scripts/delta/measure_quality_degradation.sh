#!/usr/bin/env bash
# Measure the per-read quality-degradation statistics that #694's model change needs.
#
# WHY THIS EXISTS
#   #694 says eidolon's quality model cannot represent degraded reads: 12.22% of real reads
#   have a collapsed tail (mean Q < 25 over their last 50 bases) against ~0.36% simulated.
#   The proposed fix is a two-component model -- a healthy population and a degraded one with
#   a per-read onset -- and fitting that needs numbers this repo does not have:
#
#     * how often a healthy read dips below the threshold and RECOVERS (transient excursions),
#     * how long those excursions last,
#     * how deep they go,
#     * where a collapse that does NOT recover begins.
#
#   The local fixture reproducing #694 currently invents all four. This script replaces them
#   with measurements. It is read-only and touches nothing but the input FASTQ.
#
# WHAT IT REPORTS
#   Section 1: the #694 headline rates, recomputed here so this script's view of the file is
#              known to agree with the issue before anything else is believed.
#   Section 2: low-run statistics -- every maximal run of bases below LOW_Q, classified as
#              TRANSIENT (recovers before the read ends) or TERMINAL (runs to the read end),
#              with length and depth distributions for each.
#   Section 3: onset positions of terminal runs -- the per-read collapse point the two-
#              component model needs a distribution for.
#   Section 4: mean quality by position, to compare against #694's profile table.
#
# USAGE
#   bash scripts/delta/measure_quality_degradation.sh
#   FASTQ=/path/to/reads.fastq.gz MAX_READS=2000000 bash scripts/delta/measure_quality_degradation.sh
#
# Defaults to HG002 R1, the library #694 and #695 were measured on. Set FASTQ to do R2, which
# #695 explicitly asks for and which is systematically worse in paired Illumina data.
#
# Runtime is one pass, no alignment, no eidolon build required -- minutes on a login node for
# a few million reads. It does NOT need to be a SLURM job.

set -euo pipefail

FASTQ="${FASTQ:-${DATA_DIR:-$SCRATCH/neat_data/hg002}/hg002_R1.fastq.gz}"
MAX_READS="${MAX_READS:-2000000}"   # 0 = whole file
QUAL_OFFSET="${QUAL_OFFSET:-33}"
LOW_Q="${LOW_Q:-25}"                # a base at or below this counts as "low"
TAIL_WINDOW="${TAIL_WINDOW:-50}"    # #694's window
OUT="${OUT:-quality_degradation.txt}"

[[ -f "$FASTQ" ]] || { echo "FASTQ not found: $FASTQ" >&2; exit 1; }

echo "=== measure_quality_degradation ==="
echo "fastq:       $FASTQ"
echo "max_reads:   $MAX_READS (0 = all)"
echo "qual_offset: $QUAL_OFFSET   low_q: $LOW_Q   tail_window: $TAIL_WINDOW"
echo "writing:     $OUT"
echo

# `zcat | awk` with `set -o pipefail` makes zcat take SIGPIPE if awk exits early, and `set -e`
# then aborts a step that succeeded. awk here always drains its input, but the guard costs
# nothing and this footgun has bitten this repo before (see CLAUDE.md).
set +o pipefail
zcat -f -- "$FASTQ" | awk -v offset="$QUAL_OFFSET" \
                          -v max_reads="$MAX_READS" \
                          -v low_q="$LOW_Q" \
                          -v win="$TAIL_WINDOW" '
    NR % 4 != 0 { next }                       # quality line only
    max_reads > 0 && reads >= max_reads { exit }
    {
        reads++
        n = length($0)
        if (n < win) { too_short++; next }
        if (n > maxlen) maxlen = n

        # ---- section 1: the #694 headline rates ----
        tail = 0
        for (i = n - win + 1; i <= n; i++) {
            q = index_of(substr($0, i, 1))
            tail += q
        }
        tail /= win
        if (tail < 25) lt25++
        if (tail < 20) lt20++

        # ---- sections 2-4: walk the read once ----
        run_len = 0; run_sum = 0; run_start = 0
        for (i = 1; i <= n; i++) {
            q = index_of(substr($0, i, 1))
            pos_sum[i] += q; pos_n[i]++
            if (q <= low_q) {
                if (run_len == 0) run_start = i
                run_len++; run_sum += q
            } else if (run_len > 0) {
                record_run(run_len, run_sum, run_start, 0, n)
                run_len = 0; run_sum = 0
            }
        }
        if (run_len > 0) record_run(run_len, run_sum, run_start, 1, n)
        total_bases += n
    }

    function index_of(c) { return index(chars, c) - 1 + 0 }

    # terminal = the run reaches the end of the read and never recovers
    function record_run(len, sum, start, terminal, readlen,   bucket) {
        if (terminal) {
            term_n++; term_len_sum += len; term_depth_sum += sum / len
            bucket = int((start - 1) * 10 / readlen); if (bucket > 9) bucket = 9
            onset[bucket]++
            if (len > term_len_max) term_len_max = len
        } else {
            trans_n++; trans_len_sum += len; trans_depth_sum += sum / len
            if (len <= 2) tl[1]++; else if (len <= 5) tl[2]++
            else if (len <= 10) tl[3]++; else if (len <= 25) tl[4]++; else tl[5]++
            if (len > trans_len_max) trans_len_max = len
        }
        total_low_bases += len
    }

    BEGIN { chars = "" ; for (i = 0; i < 94; i++) chars = chars sprintf("%c", i + offset) }

    END {
        if (reads == 0) { print "FATAL: no reads read"; exit 1 }
        printf "reads scanned:            %d\n", reads
        printf "reads shorter than window: %d\n", too_short + 0
        printf "read length (max):        %d\n\n", maxlen

        print "--- 1. #694 headline rates (expect ~12.22% and ~5.61% on HG002 R1) ---"
        printf "tail Q<25 (last %d):       %.2f%%  (%d reads)\n", win, 100 * lt25 / reads, lt25
        printf "tail Q<20 (last %d):       %.2f%%  (%d reads)\n\n", win, 100 * lt20 / reads, lt20

        print "--- 2. low runs (base <= Q" low_q "), per read ---"
        printf "TRANSIENT runs (recover):  %d   (%.3f per read)\n", trans_n, trans_n / reads
        if (trans_n > 0) {
            printf "  mean length:             %.2f bases   max: %d\n", trans_len_sum / trans_n, trans_len_max
            printf "  mean depth:              Q%.1f\n", trans_depth_sum / trans_n
            printf "  length 1-2 / 3-5 / 6-10 / 11-25 / 26+:  %d / %d / %d / %d / %d\n",
                   tl[1]+0, tl[2]+0, tl[3]+0, tl[4]+0, tl[5]+0
        }
        printf "TERMINAL runs (to read end): %d   (%.2f%% of reads)\n", term_n, 100 * term_n / reads
        if (term_n > 0) {
            printf "  mean length:             %.2f bases   max: %d\n", term_len_sum / term_n, term_len_max
            printf "  mean depth:              Q%.1f\n", term_depth_sum / term_n
        }
        printf "low bases overall:         %.2f%% of %d bases scanned\n\n", 100 * total_low_bases / total_bases, total_bases

        print "--- 3. onset of TERMINAL runs, by decile of read length ---"
        if (term_n > 0) for (b = 0; b <= 9; b++)
            printf "  %3d-%3d%%: %7d  (%.2f%%)\n", b*10, b*10+10, onset[b]+0, 100 * (onset[b]+0) / term_n
        print ""

        print "--- 4. mean quality by position (every 20th) ---"
        for (i = 1; i <= maxlen; i += 20)
            if (pos_n[i] > 0) printf "  pos %3d: Q%.1f\n", i, pos_sum[i] / pos_n[i]
        if (pos_n[maxlen] > 0) printf "  pos %3d: Q%.1f\n", maxlen, pos_sum[maxlen] / pos_n[maxlen]
    }
' | tee "$OUT"
set -o pipefail

echo
echo "=== done. Paste $OUT back. ==="
