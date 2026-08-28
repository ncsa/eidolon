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
a missing CIGAR op is not counted@n_missing=$(( n_missing + 1 ))@:
dead flank is scored as a good deletion@if awk -v f="$flank" 'BEGIN { exit !(f < 1) }'; then@if false; then
undeleted span passes the ratio test@-v t="${DEL_DEPTH_MAX_RATIO:-0.95}" 'BEGIN { exit !(x > t) }'@-v t="${DEL_DEPTH_MAX_RATIO:-0.95}" 'BEGIN { exit !(x > 99) }'
depth measurement failure is swallowed@d="$(samtools depth -a -r "$2:$3-$4" "$1" 2>/dev/null@d="$(samtools depth -a -r "$2:$3-$4" "$1" 2>/dev/null || true
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

extract() { awk "/^$1\(\)/,/^}\$/" "$PIPELINE"; }
for fn in max_cigar_op region_mean_depth verify_sv_cigar_ops verify_del_depth; do
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

echo "=== verify_sv_cigar_ops ==="
mkdir -p "$WORK/d"; : > "$WORK/d/truth_sv_INS.vcf.gz"; : > "$WORK/bam"
bcftools() { [[ "${1:-}" == query ]] || return 2; printf 'chr1\t1000\tA\tA%s\t.\n' "$(printf 'C%.0s' {1..200})"; }
samtools() {
    [[ "${1:-}" == view ]] || return 2
    case "${CIG_MODE:-ins}" in
        fail) return 2 ;;
        ins)  printf 'r\t0\tchr1\t900\t60\t50M200I50M\t*\t0\t0\tACGT\t*\n' ;;
        pureM) printf 'r\t0\tchr1\t900\t60\t150M\t*\t0\t0\tACGT\t*\n' ;;
    esac
}
CIGAR_MISSING=0
verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "insertion with an I op passes" 0 "$?"; is "not counted missing" 0 "$CIGAR_MISSING"
has "reports the denominator" "$(<"$WORK/o")" "1 of 1 planted"

CIGAR_MISSING=0
CIG_MODE=pureM verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "pure-M insertion is counted missing (#589)" 1 "$CIGAR_MISSING"
has "names the signature" "$(<"$WORK/o")" "NO I OP IN ANY READ"

CIG_MODE=fail verify_sv_cigar_ops x "$WORK/bam" "$WORK/d" > "$WORK/o" 2>&1
is "samtools failure is non-zero" 1 "$?"
has "failure is not read as missing" "$(<"$WORK/o")" "was not evaluated"

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
