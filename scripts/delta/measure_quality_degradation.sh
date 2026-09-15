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
#   sbatch scripts/delta/measure_quality_degradation.sh                        # a real library
#   FASTQ=/path/reads.fastq.gz STRIDE=65 sbatch scripts/delta/measure_quality_degradation.sh
#   MAX_READS=20000 bash scripts/delta/measure_quality_degradation.sh          # smoke, inline
#
# SUBMIT IT. Login nodes cap at 30 minutes and this pass does not fit: the full library is
# ~65 billion base iterations at stride 1. Two runs were killed partway learning that, and the
# second reported "done" over an empty file, which is why the guard below exists.
#
# `bash` still works and is right for a smoke run over a small MAX_READS. For a real library
# use sbatch, with MAX_READS=0 and a STRIDE large enough to cover the file -- a head-sample is
# one flowcell tile, see STRIDE below.
#
# Defaults to HG002 R1, the library #694 and #695 were measured on. Set FASTQ to do R2, which
# #695 explicitly asks for and which is systematically worse in paired Illumina data.
#
# One pass, no alignment, no eidolon build required. Read-only; it touches nothing but the
# input FASTQ.

#SBATCH --job-name=eidolon-qualdeg
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=qualdeg-%j.log

set -euo pipefail

FASTQ="${FASTQ:-${DATA_DIR:-$SCRATCH/neat_data/hg002}/hg002_R1.fastq.gz}"
MAX_READS="${MAX_READS:-2000000}"   # 0 = whole file
# Take every STRIDE-th record rather than the first MAX_READS.
#
# WHY THIS IS NOT OPTIONAL FOR A REAL LIBRARY: Illumina FASTQs are written in flowcell order.
# Measured on HG002 R1, every one of the first 200,000 records is lane 1, tile 1101 --
#   zcat R1 | awk 'NR%4==1' | head -200000 | cut -d: -f4,5 | uniq -c
#   200000 1:1101
# so a head-sample characterizes one corner of one lane, not the library. With MAX_READS=50000
# and STRIDE=1 this script reported a 11.77% collapsed-tail rate that is a tile-1101 number.
#
# STRIDE=1 (the default) preserves the old behavior for small fixtures and tests, where the
# input is synthetic and order carries no meaning.
STRIDE="${STRIDE:-1}"
QUAL_OFFSET="${QUAL_OFFSET:-33}"
LOW_Q="${LOW_Q:-25}"                # a base at or below this counts as "low"
TAIL_WINDOW="${TAIL_WINDOW:-50}"    # #694's window
# Head of the read, used to ask whether reads that END badly were ALREADY worse at the start.
# 100 sits well before the decline becomes visible in the mean profile, so a difference there
# cannot be the tail effect leaking backwards.
HEAD_WINDOW="${HEAD_WINDOW:-100}"
OUT="${OUT:-quality_degradation.txt}"

[[ -f "$FASTQ" ]] || { echo "FASTQ not found: $FASTQ" >&2; exit 1; }

echo "=== measure_quality_degradation ==="
echo "fastq:       $FASTQ"
echo "max_reads:   $MAX_READS (0 = all)   stride: $STRIDE"
echo "qual_offset: $QUAL_OFFSET   low_q: $LOW_Q   tail_window: $TAIL_WINDOW"
echo "writing:     $OUT"
echo

# `zcat | awk` with `set -o pipefail` makes zcat take SIGPIPE if awk exits early, and `set -e`
# then aborts a step that succeeded. awk here always drains its input, but the guard costs
# nothing and this footgun has bitten this repo before (see CLAUDE.md).
set +o pipefail
# `awk_status=$?` on the line AFTER the pipeline does not work: `set -e` aborts on the failing
# pipeline first, so the guard below never runs and the script dies silently. CLAUDE.md
# documents this exact footgun and prescribes `cmd || rc=$?`, which is what this is.
awk_status=0
zcat -f -- "$FASTQ" | awk -v offset="$QUAL_OFFSET" \
                          -v max_reads="$MAX_READS" \
                          -v stride="$STRIDE" \
                          -v low_q="$LOW_Q" \
                          -v win="$TAIL_WINDOW" \
                          -v head_win="$HEAD_WINDOW" '
    NR % 4 != 0 { next }                       # quality line only
    {
        rec++
        if (stride > 1 && (rec % stride) != 1) next
        if (max_reads > 0 && reads >= max_reads) exit
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

        collapsed = (tail < 25)
        if (collapsed) coll_reads++; else heal_reads++

        # ---- walk the read once for everything else ----
        run_len = 0; run_sum = 0; run_start = 0
        longest = 0; longest_start = 0
        head_sum = 0; head_low = 0; head_runs = 0; head_n = 0
        for (i = 1; i <= n; i++) {
            q = index_of(substr($0, i, 1))
            pos_sum[i] += q; pos_n[i]++

            # section 5: the head of the read, well before any tail effect
            if (i <= head_win) {
                head_n++; head_sum += q
                if (q <= low_q) {
                    head_low++
                    if (run_len == 0) head_runs++
                }
            }
            if (i > n - win && q <= low_q) tail_low++   # low bases inside the tail window

            if (q <= low_q) {
                if (run_len == 0) run_start = i
                run_len++; run_sum += q
            } else if (run_len > 0) {
                record_run(run_len, run_sum, run_start)
                if (run_len > longest) { longest = run_len; longest_start = run_start }
                run_len = 0; run_sum = 0
            }
        }
        if (run_len > 0) {
            record_run(run_len, run_sum, run_start)
            if (run_len > longest) { longest = run_len; longest_start = run_start }
        }

        # section 2: the longest low run this read carries, and where it starts. Replaces the
        # old "terminal run" metric, which counted any run touching the final base -- with mean
        # quality at position 250 around Q24, that caught 43% of reads and conflated "ends on a
        # low base" with "collapsed".
        if (longest > 0) {
            lr_n++; lr_sum += longest
            if (longest > lr_max) lr_max = longest
            if (longest <= 2) lr[1]++; else if (longest <= 5) lr[2]++
            else if (longest <= 10) lr[3]++; else if (longest <= 25) lr[4]++
            else if (longest <= 50) lr[5]++; else lr[6]++
            b = int((longest_start - 1) * 10 / n); if (b > 9) b = 9
            onset[b]++
            if (collapsed) { c_lr_n++; c_lr_sum += longest } else { h_lr_n++; h_lr_sum += longest }
        }

        # section 5: head-of-read statistics, split by whether the read ends up collapsed.
        if (head_n > 0) {
            if (collapsed) {
                c_head_q += head_sum / head_n; c_head_low += head_low; c_head_runs += head_runs
                c_head_bases += head_n
            } else {
                h_head_q += head_sum / head_n; h_head_low += head_low; h_head_runs += head_runs
                h_head_bases += head_n
            }
        }
        total_bases += n
    }

    function index_of(c) { return index(chars, c) - 1 + 0 }

    function record_run(len, sum, start) {
        run_n++; run_len_sum += len; run_depth_sum += sum / len
        if (len <= 2) tl[1]++; else if (len <= 5) tl[2]++
        else if (len <= 10) tl[3]++; else if (len <= 25) tl[4]++; else tl[5]++
        if (len > run_len_max) run_len_max = len
        total_low_bases += len
    }

    BEGIN { chars = "" ; for (i = 0; i < 94; i++) chars = chars sprintf("%c", i + offset) }

    END {
        if (reads == 0) { print "FATAL: no reads read"; exit 1 }
        printf "reads scanned:            %d of %d records (stride %d)\n", reads, rec, stride
        printf "reads shorter than window: %d\n", too_short + 0
        printf "read length (max):        %d\n\n", maxlen

        print "--- 1. #694 headline rates (expect ~12.22% and ~5.61% on HG002 R1) ---"
        printf "tail Q<25 (last %d):       %.2f%%  (%d reads)\n", win, 100 * lt25 / reads, lt25
        printf "tail Q<20 (last %d):       %.2f%%  (%d reads)\n\n", win, 100 * lt20 / reads, lt20

        print "--- 2. low runs (base <= Q" low_q ") ---"
        printf "all runs:                  %d   (%.3f per read)\n", run_n, run_n / reads
        if (run_n > 0) {
            printf "  mean length:             %.2f bases   max: %d\n", run_len_sum / run_n, run_len_max
            printf "  mean depth:              Q%.1f\n", run_depth_sum / run_n
            printf "  length 1-2 / 3-5 / 6-10 / 11-25 / 26+:  %d / %d / %d / %d / %d\n",
                   tl[1]+0, tl[2]+0, tl[3]+0, tl[4]+0, tl[5]+0
        }
        printf "reads with any low run:    %d  (%.2f%%)\n", lr_n, 100 * lr_n / reads
        if (lr_n > 0) {
            printf "  LONGEST run per read, mean %.2f bases   max: %d\n", lr_sum / lr_n, lr_max
            printf "  longest 1-2 / 3-5 / 6-10 / 11-25 / 26-50 / 51+:  %d / %d / %d / %d / %d / %d\n",
                   lr[1]+0, lr[2]+0, lr[3]+0, lr[4]+0, lr[5]+0, lr[6]+0
            if (c_lr_n > 0) printf "  longest run, COLLAPSED reads %.2f bases\n", c_lr_sum / c_lr_n
            if (h_lr_n > 0) printf "  longest run, healthy reads   %.2f bases\n", h_lr_sum / h_lr_n
        }
        printf "low bases overall:         %.2f%% of %d bases scanned\n", 100 * total_low_bases / total_bases, total_bases
        printf "low bases in the last %d:   %.2f%%\n\n", win, 100 * tail_low / (reads * win)

        print "--- 3. start of each read LONGEST low run, by decile of read length ---"
        if (lr_n > 0) for (b = 0; b <= 9; b++)
            printf "  %3d-%3d%%: %7d  (%.2f%%)\n", b*10, b*10+10, onset[b]+0, 100 * (onset[b]+0) / lr_n
        print ""

        print "--- 5. DOES A COLLAPSED READ START BADLY? first " head_win " positions ---"
        print "  Onset model: a read is normal until it collapses -> these should MATCH."
        print "  Propensity:  collapsed reads are noisier throughout -> they should DIFFER."
        printf "  collapsed reads %d   healthy %d\n", coll_reads+0, heal_reads+0
        if (coll_reads > 0 && heal_reads > 0) {
            printf "  mean quality   collapsed Q%.2f   healthy Q%.2f   difference %+.2f\n",
                   c_head_q / coll_reads, h_head_q / heal_reads,
                   (c_head_q / coll_reads) - (h_head_q / heal_reads)
            printf "  low bases      collapsed %.3f%%   healthy %.3f%%   ratio %.2fx\n",
                   100 * c_head_low / c_head_bases, 100 * h_head_low / h_head_bases,
                   (h_head_low > 0) ? (c_head_low / c_head_bases) / (h_head_low / h_head_bases) : 0
            printf "  low runs/read  collapsed %.3f    healthy %.3f    ratio %.2fx\n",
                   c_head_runs / coll_reads, h_head_runs / heal_reads,
                   (h_head_runs > 0) ? (c_head_runs / coll_reads) / (h_head_runs / heal_reads) : 0
        }
        print ""

        print "--- 4. mean quality by position (every 20th) ---"
        for (i = 1; i <= maxlen; i += 20)
            if (pos_n[i] > 0) printf "  pos %3d: Q%.1f\n", i, pos_sum[i] / pos_n[i]
        if (pos_n[maxlen] > 0) printf "  pos %3d: Q%.1f\n", maxlen, pos_sum[maxlen] / pos_n[maxlen]
    }
' > "$OUT" || awk_status=$?
set -o pipefail

# A killed awk writes nothing and its END block -- where the "no reads" guard lives -- never
# runs. Without this check the script printed "done" over an empty file, which is a measurement
# tool reporting success having measured nothing.
#
# This is how it happened: 260,990,178 reads at stride 1 is ~65 billion base iterations, and a
# login node killed it partway. Use sbatch, or a stride, or both.
# One check, not two: a non-zero status and a truncated file are the same failure seen from
# two sides, and as separate guards each masked the other under mutation -- neither was pinned.
if [[ "$awk_status" -ne 0 ]] || ! grep -q "headline rates" "$OUT" 2>/dev/null; then
    echo "FATAL: the measurement pass did not reach its report (awk exit $awk_status)." >&2
    echo "  Nothing in $OUT can be trusted -- do not read numbers off a partial file." >&2
    echo "  Most likely killed for running too long: a login node will not carry a" >&2
    echo "  full-file pass over a large library. 260,990,178 reads at stride 1 is about" >&2
    echo "  65 billion base iterations. Submit with sbatch, and/or raise STRIDE." >&2
    exit 1
fi

cat "$OUT"
echo
echo "=== done. Paste $OUT back. ==="
