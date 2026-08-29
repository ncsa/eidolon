#!/usr/bin/env bash
# Regression tests for the CIGAR-op (#589) and DEL-depth (#590/#221) gates in
# sv_pipeline.sbatch. Functions are extracted verbatim; stubs remove bcftools/samtools.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PIPELINE:-$HERE/../sv_pipeline.sbatch}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$PIPELINE" "$WORK/mutant.sbatch"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sbatch"
        if cmp -s "$PIPELINE" "$WORK/mutant.sbatch"; then
            printf '  ERROR   %-54s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if PIPELINE="$WORK/mutant.sbatch" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-53s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
max op reads the wrong CIGAR field@n = split($6, parts, /[0-9]+/)@n = split($10, parts, /[0-9]+/)
max op ignores the requested op kind@if (parts[i] == want && lens[i - 1] + 0 > max)@if (lens[i - 1] + 0 > max)
a missing clip is not counted@n_missing=$(( n_missing + 1 ))@:
clip threshold ignored, any clip counts@if (ops[i] == "S" && lens[i - 1] + 0 >= min)@if (ops[i] == "S")
hard clips counted as soft@if (ops[i] == "S" && lens[i - 1] + 0 >= min)@if ((ops[i] == "S" || ops[i] == "H") && lens[i - 1] + 0 >= min)
control window ignored, bare threshold used@if [[ "${at_locus:-0}" -le "${background:-0}" ]]; then@if [[ "${at_locus:-0}" -le 0 ]]; then
background takes the larger control@[[ "${ctl_r:-0}" -lt "${background:-0}" ]] && background="$ctl_r"@[[ "${ctl_r:-0}" -gt "${background:-0}" ]] && background="$ctl_r"
dead flank is scored as a good deletion@if awk -v f="$flank" 'BEGIN { exit !(f < 1) }'; then@if false; then
undeleted span passes the ratio test@-v t="${DEL_DEPTH_MAX_RATIO:-0.95}" 'BEGIN { exit !(x > t) }'@-v t="${DEL_DEPTH_MAX_RATIO:-0.95}" 'BEGIN { exit !(x > 99) }'
depth measurement failure is swallowed@d="$(samtools depth -a -r "$2:$3-$4" "$1" 2>/dev/null@d="$(samtools depth -a -r "$2:$3-$4" "$1" 2>/dev/null || true
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

extract() { awk "/^$1\(\)/,/^}\$/" "$PIPELINE"; }
for fn in max_cigar_op count_clipped_reads region_mean_depth verify_sv_cigar_ops verify_del_depth; do
    src="$(extract "$fn")"
    [[ -n "$src" ]] || { echo "FATAL: could not extract $fn"; exit 2; }
    eval "$src"
done

PASS=0; FAIL=0
ok() { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
is() { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$2" "$3"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

echo "=== max_cigar_op (known answers) ==="
sam() { printf 'r\t0\tchr1\t1\t60\t%s\t*\t0\t0\tACGT\t*\n' "$1"; }
is "single I op"            300 "$(sam '50M300I50M' | max_cigar_op I)"
is "largest of several I"   300 "$(printf '%s' "$(sam '10I50M'; sam '50M300I')" | max_cigar_op I)"
is "D ops are not I ops"      0 "$(sam '50M300D50M' | max_cigar_op I)"
is "D op is found"          300 "$(sam '50M300D50M' | max_cigar_op D)"
is "soft clip is not an I"    0 "$(sam '50S100M' | max_cigar_op I)"
is "pure M yields zero"       0 "$(sam '150M' | max_cigar_op I)"
is "no records yields zero"   0 "$(printf '' | max_cigar_op I)"

echo "=== count_clipped_reads (known answers) ==="
is "soft clip at the end counts"        1 "$(sam '100M51S' | count_clipped_reads 20)"
is "soft clip at the start counts"      1 "$(sam '51S100M' | count_clipped_reads 20)"
is "clip below the threshold does not"  0 "$(sam '140M11S' | count_clipped_reads 20)"
is "hard clip is not a soft clip"       0 "$(sam '100M51H' | count_clipped_reads 20)"
is "an I op is not a clip"              0 "$(sam '50M51I50M' | count_clipped_reads 20)"
is "pure M is not a clip"               0 "$(sam '151M' | count_clipped_reads 20)"
is "each read counted once, not twice"  1 "$(sam '30S91M30S' | count_clipped_reads 20)"
is "no records yields zero"             0 "$(printf '' | count_clipped_reads 20)"

echo "=== verify_sv_cigar_ops: clip support against a flanking control ==="
mkdir -p "$WORK/d"; : > "$WORK/d/truth_sv_INS.vcf.gz"; : > "$WORK/bam"
READ_LEN=151
bcftools() { [[ "${1:-}" == query ]] || return 2; printf 'chr1\t1000\tA\tA%s\t.\n' "$(printf 'C%.0s' {1..200})"; }
# The locus window is centred on POS=1000; the controls sit ~5 kb away. Keying on the region
# string lets the stub answer each independently.
samtools() {
    [[ "${1:-}" == view ]] || return 2
    [[ "${CIG_MODE:-clipped}" == fail ]] && return 2
    case "$3" in
        chr1:700-1300)
            case "${CIG_MODE:-clipped}" in
                clipped)   for i in 1 2 3 4 5; do sam '100M51S'; done ;;
                pureM)     for i in 1 2 3 4 5; do sam '151M'; done ;;
                iop_only)  for i in 1 2 3 4 5; do sam '50M2I99M'; done ;;
                noisy)     for i in 1 2 3; do sam '100M51S'; done ;;
                asym)      for i in 1 2 3 4 5; do sam '100M51S'; done ;;
            esac ;;
        chr1:1-301)        # LEFT flanking control
            [[ "${CIG_MODE:-clipped}" == noisy ]] && for i in 1 2 3 4; do sam '100M51S'; done
            ;;
        chr1:5900-6500)    # RIGHT flanking control
            case "${CIG_MODE:-clipped}" in
                noisy) for i in 1 2 3 4; do sam '100M51S'; done ;;
                # A neighbouring SV lands in this window only. The clean side is the honest
                # background, so the minimum is what must be used.
                asym)  for i in 1 2 3 4 5 6 7 8 9; do sam '100M51S'; done ;;
            esac
            ;;
        *) ;;
    esac
    return 0
}

CIGAR_MISSING=0
verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "clipped reads above control passes" 0 "$?"
is "and are not counted missing"        0 "$CIGAR_MISSING"
has "reports locus and control counts"  "$(<"$WORK/o")" "5 clipped (control 0)"

# THE #624 REGRESSION: an insertion the aligner represents as pure M must FAIL. Under the old
# check this passed on any stray 1-2 bp I op from the small-indel model.
CIGAR_MISSING=0
CIG_MODE=pureM verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "pure-M insertion is counted missing (#589)" 1 "$CIGAR_MISSING"
has "says clip support is absent" "$(<"$WORK/o")" "NO CLIP SUPPORT"

# MUST NOT FIRE ON NOISE: a small I op with no excess clipping is background from the
# mutation model, not evidence of a 200 bp insertion. This is exactly what job 21575385's
# "87 of 87" was reading.
CIGAR_MISSING=0
CIG_MODE=iop_only verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "a stray I op alone does NOT count as support (#624)" 1 "$CIGAR_MISSING"

# MUST NOT FIRE: clipping that is no higher than the surrounding reference is not evidence.
CIGAR_MISSING=0
CIG_MODE=noisy verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "clipping at or below the control is not support" 1 "$CIGAR_MISSING"
has "reports the control it lost to" "$(<"$WORK/o")" "vs 4 in flanking control"

# An SV in ONE control window must not fail a perfectly good locus. Background is the
# minimum of the two flanks for exactly this reason: 5 clipped here, 0 on the clean side,
# 9 on the contaminated side. Taking the maximum would report this locus as unsupported.
CIGAR_MISSING=0
CIG_MODE=asym verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "a contaminated control flank does not fail a good locus" 0 "$CIGAR_MISSING"

CIG_MODE=fail verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "samtools failure is non-zero" 1 "$?"
has "failure is not read as missing support" "$(<"$WORK/o")" "was not evaluated"

echo "=== verify_del_depth ==="
bcftools() { [[ "${1:-}" == query ]] || return 2; printf 'chr1\t1000\tA\t<DEL>\t1500\n'; }
samtools() {
    [[ "${1:-}" == depth ]] || return 2
    [[ "${DEPTH_MODE:-}" == fail ]] && return 2
    local reg="$4" v
    case "$reg" in
        chr1:1001-1500) v="${INSIDE:-0}" ;;
        *)              v="${FLANK:-60}" ;;
    esac
    local i; for i in 1 2 3; do printf 'chr1\t%d\t%s\n' "$i" "$v"; done
}
DEL_DEPTH_BAD=0; DEL_DEPTH_DEAD=0
verify_del_depth x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "removed coverage passes"        0 "$DEL_DEPTH_BAD"
is "live flank is not called dead"  0 "$DEL_DEPTH_DEAD"

DEL_DEPTH_BAD=0; DEL_DEPTH_DEAD=0
INSIDE=60 verify_del_depth x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "undeleted span is caught"       1 "$DEL_DEPTH_BAD"
has "reports the ratio"  "$(<"$WORK/o")" "coverage NOT removed"

# The first #590 attempt: inside depth 0.00 looked like success while the contig was empty.
DEL_DEPTH_BAD=0; DEL_DEPTH_DEAD=0
INSIDE=0 FLANK=0 verify_del_depth x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "dead flank is not scored as a deletion" 0 "$DEL_DEPTH_BAD"
is "dead flank is counted separately"       1 "$DEL_DEPTH_DEAD"
has "dead flank says uninterpretable" "$(<"$WORK/o")" "no local baseline"

DEPTH_MODE=fail verify_del_depth x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "depth failure is non-zero" 1 "$?"
has "depth failure is not read as a deletion" "$(<"$WORK/o")" "was not evaluated"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
