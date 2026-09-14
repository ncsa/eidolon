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
terminal runs counted as transient@if (terminal) {@if (0) {
tail threshold 25 -> 35@if (tail < 25) lt25++@if (tail < 35) lt25++
onset bucketed by global max, not read length@int((start - 1) * 10 / readlen)@int((start - 1) * 10 / maxlen)
low-base denominator back to reads*maxlen@total_low_bases / total_bases@total_low_bases / (reads * maxlen)
run continues across a recovery@} else if (run_len > 0) {@} else if (0) {
MUTS
    echo
    [[ "$survived" -eq 0 ]] && { echo "all mutations caught"; exit 0; } || { echo "$survived survived"; exit 1; }
fi

PASS=0; FAIL=0
# Floor on how many assertions must execute. Raise it when adding tests; if it ever reads low,
# an assertion stopped running rather than started failing.
MIN_ASSERTIONS=21
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
OUT_KNOWN="$(FASTQ="$WORK/known.fastq.gz" MAX_READS=0 OUT="$WORK/known.txt" bash "$SCRIPT" 2>&1)"

hasw "scans all 10 reads"                  "$OUT_KNOWN" "reads scanned:            10"
hasw "read length reported as 100"         "$OUT_KNOWN" "read length (max):        100"
hasw "tail Q<25 is 30.00%"                 "$OUT_KNOWN" "tail Q<25 (last 50):       30.00%  (3 reads)"
hasw "tail Q<20 is 30.00%"                 "$OUT_KNOWN" "tail Q<20 (last 50):       30.00%  (3 reads)"
hasw "3 transient runs"                    "$OUT_KNOWN" "TRANSIENT runs (recover):  3"
hasw "transient mean length 5.00"          "$OUT_KNOWN" "mean length:             5.00 bases   max: 5"
hasw "transient depth Q10.0"               "$OUT_KNOWN" "mean depth:              Q10.0"
hasw "transient length bucket 3-5 holds 3" "$OUT_KNOWN" "1-2 / 3-5 / 6-10 / 11-25 / 26+:  0 / 3 / 0 / 0 / 0"
hasw "3 terminal runs, 30% of reads"       "$OUT_KNOWN" "TERMINAL runs (to read end): 3   (30.00% of reads)"
hasw "terminal mean length 40.00"          "$OUT_KNOWN" "mean length:             40.00 bases   max: 40"
hasw "low bases 13.50% of 1000"            "$OUT_KNOWN" "low bases overall:         13.50% of 1000 bases scanned"
hasw "onset lands in the 60-70% decile"    "$OUT_KNOWN" "60- 70%:       3  (100.00%)"
# Position 21 is Q10 in 3 reads and Q40 in 7: (7*40 + 3*10)/10 = 31.0
hasw "mean at position 21 is Q31.0"        "$OUT_KNOWN" "pos  21: Q31.0"
hasw "mean at position 1 is Q40.0"         "$OUT_KNOWN" "pos   1: Q40.0"

echo "=== must not fire: a file with no low bases reports no runs ==="
make_fastq "$WORK/healthy.fastq.gz" <<EOF
10 $(rep I 100)
EOF
OUT_HEALTHY="$(FASTQ="$WORK/healthy.fastq.gz" MAX_READS=0 OUT="$WORK/healthy.txt" bash "$SCRIPT" 2>&1)"
hasw "healthy file: no transient runs" "$OUT_HEALTHY" "TRANSIENT runs (recover):  0"
hasw "healthy file: no terminal runs"  "$OUT_HEALTHY" "TERMINAL runs (to read end): 0"
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
     "$OUT_BETWEEN" "TRANSIENT runs (recover):  0"

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

echo
printf 'assertions: %d passed, %d failed\n' "$PASS" "$FAIL"
if [[ $((PASS + FAIL)) -lt $MIN_ASSERTIONS ]]; then
    echo "FAIL: only $((PASS + FAIL)) assertions ran, expected at least $MIN_ASSERTIONS" >&2
    exit 1
fi
[[ "$FAIL" -eq 0 ]] || exit 1
