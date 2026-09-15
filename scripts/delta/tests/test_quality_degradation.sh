#!/usr/bin/env bash
# Tests for measure_quality_degradation.sh.
#
# It exists to replace three invented fixture parameters in #694 with measured ones, so a
# wrong number here would be worse than no number: it would be believed. Every expectation
# below is computed by hand from a fixture small enough to check by eye.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT:-$HERE/../measure_quality_degradation.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$SCRIPT" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        # An unapplied mutation and a surviving one produce identical output. Assert the edit
        # landed before believing anything about the verdict.
        if cmp -s "$SCRIPT" "$WORK/mutant.sh"; then
            printf '  ERROR    %-50s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if SCRIPT="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-50s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTS'
tail threshold 25 -> 35@if (tail < 25) lt25++@if (tail < 35) lt25++
collapsed classified by the wrong threshold@collapsed = (tail < 25)@collapsed = (tail < 35)
longest-run start bucketed by global max@int((longest_start - 1) * 10 / n)@int((longest_start - 1) * 10 / maxlen)
low-base denominator back to reads*maxlen@total_low_bases / total_bases@total_low_bases / (reads * maxlen)
run continues across a recovery@} else if (run_len > 0) {@} else if (0) {
stride ignored, head-sampling restored@if (stride > 1 && (rec % stride) != 1) next@if (0) next
longest run never updated@if (run_len > longest) { longest = run_len; longest_start = run_start }@if (0) { longest = run_len; longest_start = run_start }
head window ignored, whole read counted as head@if (i <= head_win) {@if (1) {
collapsed and healthy populations swapped@if (collapsed) coll_reads++; else heal_reads++@if (!collapsed) coll_reads++; else heal_reads++
tail-window low bases counted everywhere@if (i > n - win && q <= low_q) tail_low++@if (q <= low_q) tail_low++
a killed or truncated pass reports success@if [[ "$awk_status" -ne 0 ]] || ! grep -q "headline rates" "$OUT" 2>/dev/null; then@if false; then
status swallowed by set -e instead of captured@' > "$OUT" || awk_status=$?@' > "$OUT"; awk_status=$?
MUTS
    echo
    [[ "$survived" -eq 0 ]] && { echo "all mutations caught"; exit 0; } || { echo "$survived survived"; exit 1; }
fi

PASS=0; FAIL=0
# Floor on how many assertions must execute. Raise it when adding tests; if it ever reads low,
# an assertion stopped running rather than started failing.
MIN_ASSERTIONS=36
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
hasw(){ local h n; h="$(printf '%s' "$2" | tr -d ' ')"; n="$(printf '%s' "$3" | tr -d ' ')"
        case "$h" in *"$n"*) ok "$1";; *) bad "$1" "contains (spaces removed): $3" "$2";; esac; }

# Build a FASTQ from a list of "<count> <qual-string>" pairs on stdin.
make_fastq() {
    local out="$1" i=0 count qual
    : > "$WORK/plain.fastq"
    while read -r count qual; do
        [[ -n "$count" ]] || continue
        for ((c = 0; c < count; c++)); do
            printf '@r%d\n%s\n+\n%s\n' "$i" "$(printf 'A%.0s' $(seq 1 ${#qual}))" "$qual" >> "$WORK/plain.fastq"
            i=$((i + 1))
        done
    done
    gzip -cf "$WORK/plain.fastq" > "$out"
}

# Q40 = 'I', Q10 = '+', Q30 = '?'  (Phred+33)
rep() { printf "%${2}s" | tr ' ' "$1"; }

echo "=== known-answer fixture: 4 healthy, 3 with one transient dip, 3 that collapse ==="
# 100-base reads.
#   4 x flat Q40                                 -> no low runs
#   3 x Q40, a 5-base Q10 dip at 21-25, then Q40 -> 3 TRANSIENT runs, length 5, depth Q10
#   3 x Q40 to 60 then Q10 to the end            -> 3 TERMINAL runs, length 40, onset 61
# Last-50 means: healthy 40, dipped 40 (dip is outside the window), collapsed 10.
#   tail Q<25 = 3/10 = 30.00%   tail Q<20 = 3/10 = 30.00%
#   low bases = (3*5 + 3*40) / 1000 = 13.50%
make_fastq "$WORK/known.fastq.gz" <<EOF
4 $(rep I 100)
3 $(rep I 20)$(rep + 5)$(rep I 75)
3 $(rep I 60)$(rep + 40)
EOF
# HEAD_WINDOW=50 so the head statistics cover positions 1-50 only; at the default 100 the
# "head" would be the whole 100-base read and could not distinguish a start from an end.
OUT_KNOWN="$(FASTQ="$WORK/known.fastq.gz" MAX_READS=0 HEAD_WINDOW=50 OUT="$WORK/known.txt" bash "$SCRIPT" 2>&1)"

hasw "scans all 10 reads"                  "$OUT_KNOWN" "reads scanned:            10"
hasw "read length reported as 100"         "$OUT_KNOWN" "read length (max):        100"
hasw "tail Q<25 is 30.00%"                 "$OUT_KNOWN" "tail Q<25 (last 50):       30.00%  (3 reads)"
hasw "tail Q<20 is 30.00%"                 "$OUT_KNOWN" "tail Q<20 (last 50):       30.00%  (3 reads)"
# 3 dips of 5 + 3 collapses of 40 = 6 runs; mean (15+120)/6 = 22.5; buckets 3-5 and 26+.
hasw "6 low runs in total"                 "$OUT_KNOWN" "all runs:                  6   (0.600 per read)"
hasw "mean run length 22.50"               "$OUT_KNOWN" "mean length:             22.50 bases   max: 40"
hasw "run depth Q10.0"                     "$OUT_KNOWN" "mean depth:              Q10.0"
hasw "run length buckets"                  "$OUT_KNOWN" "1-2 / 3-5 / 6-10 / 11-25 / 26+:  0 / 3 / 0 / 0 / 3"
hasw "6 of 10 reads carry a low run"       "$OUT_KNOWN" "reads with any low run:    6  (60.00%)"
hasw "longest-run mean is 22.50"           "$OUT_KNOWN" "LONGEST run per read, mean 22.50 bases   max: 40"
hasw "longest-run buckets"                 "$OUT_KNOWN" "longest 1-2 / 3-5 / 6-10 / 11-25 / 26-50 / 51+:  0 / 3 / 0 / 0 / 3 / 0"
# The split is the point: collapsed reads carry the 40-base run, healthy ones the 5-base dip.
hasw "collapsed reads longest run 40.00"   "$OUT_KNOWN" "longest run, COLLAPSED reads 40.00 bases"
hasw "healthy reads longest run 5.00"      "$OUT_KNOWN" "longest run, healthy reads   5.00 bases"
hasw "low bases 13.50% of 1000"            "$OUT_KNOWN" "low bases overall:         13.50% of 1000 bases scanned"
# Only the collapses fall inside the last 50: 3 x 40 of 10 x 50.
hasw "low bases in the tail window 24.00%" "$OUT_KNOWN" "low bases in the last 50:   24.00%"
# Longest run starts: dips at 21 -> decile 20-30%, collapses at 61 -> decile 60-70%.
hasw "dips start in the 20-30% decile"     "$OUT_KNOWN" "20- 30%:       3  (50.00%)"
hasw "collapses start in the 60-70% decile" "$OUT_KNOWN" "60- 70%:       3  (50.00%)"
# Position 21 is Q10 in 3 reads and Q40 in 7: (7*40 + 3*10)/10 = 31.0
hasw "mean at position 21 is Q31.0"        "$OUT_KNOWN" "pos  21: Q31.0"
hasw "mean at position 1 is Q40.0"         "$OUT_KNOWN" "pos   1: Q40.0"

echo "=== must not fire: a file with no low bases reports no runs ==="
make_fastq "$WORK/healthy.fastq.gz" <<EOF
10 $(rep I 100)
EOF
OUT_HEALTHY="$(FASTQ="$WORK/healthy.fastq.gz" MAX_READS=0 OUT="$WORK/healthy.txt" bash "$SCRIPT" 2>&1)"
hasw "healthy file: no low runs at all" "$OUT_HEALTHY" "all runs:                  0"
hasw "healthy file: no read carries one" "$OUT_HEALTHY" "reads with any low run:    0"
hasw "healthy file: tail Q<25 is 0.00%" "$OUT_HEALTHY" "tail Q<25 (last 50):       0.00%"

echo "=== the tail threshold is Q25, not merely 'some threshold' ==="
# A read whose tail sits BETWEEN the thresholds is what discriminates them. Flat Q30 must not
# count at Q<25; it would at Q<35. Without this the threshold is unpinned -- found by mutation.
make_fastq "$WORK/between.fastq.gz" <<EOF
1 $(rep I 100)
1 $(rep ? 100)
EOF
OUT_BETWEEN="$(FASTQ="$WORK/between.fastq.gz" MAX_READS=0 OUT="$WORK/between.txt" bash "$SCRIPT" 2>&1)"
hasw "a flat-Q30 read does not count as a collapsed tail" \
     "$OUT_BETWEEN" "tail Q<25 (last 50):       0.00%"
hasw "and it has no low runs either, being above Q25" \
     "$OUT_BETWEEN" "all runs:                  0"
# Q30 sits between the thresholds, so this is what pins the CLASSIFIER rather than just the
# reported rate: at Q<25 neither read is collapsed; at Q<35 the flat-Q30 one would be.
hasw "neither read is classified as collapsed at the Q25 threshold" \
     "$OUT_BETWEEN" "collapsed reads 0   healthy 2"

echo "=== onset decile uses THIS read's length, not the file maximum ==="
# One 200-base read collapsing at 101 (decile 50-60%) and one 100-base read collapsing at 51
# (also 50-60%). Bucketing by the file maximum would put the short read's onset at 25-30%.
make_fastq "$WORK/mixed_len.fastq.gz" <<EOF
1 $(rep I 100)$(rep + 100)
1 $(rep I 50)$(rep + 50)
EOF
OUT_MIXED="$(FASTQ="$WORK/mixed_len.fastq.gz" MAX_READS=0 OUT="$WORK/mixed.txt" bash "$SCRIPT" 2>&1)"
hasw "both onsets land in the same decile despite different read lengths" \
     "$OUT_MIXED" "50- 60%:       2  (100.00%)"
# 150 low bases of 300 scanned. Dividing by reads*maxlen (2*200) would report 37.50%.
hasw "low-base rate divides by bases scanned, not reads x longest read" \
     "$OUT_MIXED" "low bases overall:         50.00% of 300 bases scanned"

echo "=== section 5 separates an onset population from a propensity one ==="
# ONSET shape: collapsed reads are pristine until they collapse at 61, so their first 50
# positions must be indistinguishable from a healthy read.
make_fastq "$WORK/onset.fastq.gz" <<EOF
6 $(rep I 100)
4 $(rep I 60)$(rep + 40)
EOF
OUT_ONSET="$(FASTQ="$WORK/onset.fastq.gz" MAX_READS=0 HEAD_WINDOW=50 OUT="$WORK/onset.txt" bash "$SCRIPT" 2>&1)"
hasw "onset fixture: head quality is identical, difference +0.00" \
     "$OUT_ONSET" "mean quality   collapsed Q40.00   healthy Q40.00   difference +0.00"

# PROPENSITY shape: collapsed reads are noisier from the very start. Same collapse, but their
# head carries dips the healthy reads do not.
make_fastq "$WORK/propensity.fastq.gz" <<EOF
6 $(rep I 100)
4 $(rep I 10)$(rep + 10)$(rep I 10)$(rep + 10)$(rep I 20)$(rep + 40)
EOF
OUT_PROP="$(FASTQ="$WORK/propensity.fastq.gz" MAX_READS=0 HEAD_WINDOW=50 OUT="$WORK/prop.txt" bash "$SCRIPT" 2>&1)"
# Head of a collapsed read: 30 x Q40 + 20 x Q10 = Q28.0 against Q40.0 healthy.
hasw "propensity fixture: collapsed reads are already worse by 12 Q" \
     "$OUT_PROP" "mean quality   collapsed Q28.00   healthy Q40.00   difference -12.00"
hasw "propensity fixture: and carry low runs the healthy ones do not" \
     "$OUT_PROP" "low runs/read  collapsed 2.000    healthy 0.000"

echo "=== STRIDE samples across the file, not the head ==="
# Alternating healthy / collapsed, so the two halves are separable by position. A head-sample
# of a real FASTQ is one flowcell tile (measured: 200,000 of 200,000 head reads on HG002 R1 are
# tile 1101), which is what STRIDE exists to avoid. Records are 1-indexed: stride 2 takes
# records 1,3,5,7,9 -- every healthy one, and no collapsed one.
make_fastq "$WORK/alternating.fastq.gz" <<EOF
1 $(rep I 100)
1 $(rep I 50)$(rep + 50)
1 $(rep I 100)
1 $(rep I 50)$(rep + 50)
1 $(rep I 100)
1 $(rep I 50)$(rep + 50)
1 $(rep I 100)
1 $(rep I 50)$(rep + 50)
1 $(rep I 100)
1 $(rep I 50)$(rep + 50)
EOF
OUT_ALL="$(FASTQ="$WORK/alternating.fastq.gz" MAX_READS=0 OUT="$WORK/all.txt" bash "$SCRIPT" 2>&1)"
hasw "stride 1 sees both populations: half the reads collapse" \
     "$OUT_ALL" "tail Q<25 (last 50):       50.00%"
OUT_STRIDE="$(FASTQ="$WORK/alternating.fastq.gz" MAX_READS=0 STRIDE=2 OUT="$WORK/stride.txt" bash "$SCRIPT" 2>&1)"
hasw "stride 2 selects only the healthy records" \
     "$OUT_STRIDE" "tail Q<25 (last 50):       0.00%"
hasw "and reports what it sampled and out of how many" \
     "$OUT_STRIDE" "reads scanned:            5 of 10 records (stride 2)"

echo "=== a pass that does not reach its report is a hard failure ==="
# A killed awk writes nothing and its END guard never runs. This drives the same path: an
# input with no records means the pass exits non-zero, and the script must refuse rather than
# print "done" over an empty file. Measured on the real thing -- 261M reads at stride 1 on a
# login node was killed partway and the script reported success having measured nothing.
: | gzip -c > "$WORK/empty.fastq.gz"
if OUT_EMPTY="$(FASTQ="$WORK/empty.fastq.gz" MAX_READS=0 OUT="$WORK/empty.txt" bash "$SCRIPT" 2>&1)"; then
    bad "an empty pass must exit non-zero" "non-zero exit" "exit 0"
else
    ok "an empty pass exits non-zero"
fi
case "$OUT_EMPTY" in
    *FATAL*) ok "and says FATAL rather than done" ;;
    *) bad "and says FATAL rather than done" "FATAL in output" "$OUT_EMPTY" ;;
esac
case "$OUT_EMPTY" in
    *"=== done"*) bad "and does NOT claim done" "no done banner" "$OUT_EMPTY" ;;
    *) ok "and does NOT claim done" ;;
esac

echo
printf 'assertions: %d passed, %d failed\n' "$PASS" "$FAIL"
if [[ $((PASS + FAIL)) -lt $MIN_ASSERTIONS ]]; then
    echo "FAIL: only $((PASS + FAIL)) assertions ran, expected at least $MIN_ASSERTIONS" >&2
    exit 1
fi
[[ "$FAIL" -eq 0 ]] || exit 1
