#!/usr/bin/env bash
# Regression tests for the INS read-evidence gate in sv_pipeline.sbatch (#540).
#
# Functions are extracted verbatim from the production pipeline. Command stubs make this
# independent of bcftools, samtools, and a BAM; the cases under test are whether every probe was
# evaluated and whether failures are distinguishable from zero read support.
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
            printf '  ERROR   %-50s mutation did not apply\n' "$label"
            survived=$((survived + 1))
            continue
        fi
        if PIPELINE="$WORK/mutant.sbatch" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-49s <- nothing caught this\n' "$label"
            survived=$((survived + 1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
probe matcher drops zero-support rows@END { for (i = 1; i <= n; i++) print chrom[i], pos[i], len[i], hits[i] }@END { for (i = 1; i <= n; i++) if (hits[i] > 0) print chrom[i], pos[i], len[i], hits[i] }
samtools failures are swallowed@if ! samtools view@if false; then
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]
    exit $?
fi

extract() { awk "/^$1\(\)/,/^}\$/" "$PIPELINE"; }
for fn in count_probe_hits verify_planted_ins; do
    src="$(extract "$fn")"
    [[ -n "$src" ]] || { echo "FATAL: could not extract $fn from $PIPELINE"; exit 2; }
    eval "$src"
    [[ "$fn" == "count_probe_hits" ]] && count_probe_hits_src="$src"
done

PASS=0
FAIL=0
ok() { PASS=$((PASS + 1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL + 1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
is() { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$2" "$3"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

# The stub emits one literal insertion with a 60 bp novel sequence. Every successful samtools
# read contains the derived middle probe in SAM field 10.
PROBE="$(printf 'C%.0s' {1..30})"
bcftools() {
    [[ "${1:-}" == "query" ]] || return 2
    printf 'chr1\t100\tA\tA%s\n' "$(printf 'C%.0s' {1..60})"
}
samtools() {
    case "${SAMTOOLS_MODE:-support}" in
        fail) return 2 ;;
        empty) return 0 ;;
        support)
            printf 'read1\t0\tchr1\t1\t60\t30M\t*\t0\t0\t%s\t*\n' "$PROBE"
            return 0
            ;;
    esac
    return 2
}

truth="$WORK/truth.vcf.gz"
bam="$WORK/tumor.bam"
: > "$truth"
: > "$bam"

echo "=== supported probe ==="
INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
rc=$?
out="$(<"$WORK/out")"
is "supported probe evaluates successfully" 0 "$rc"
has "supported probe reports read evidence" "$out" "2 read(s)"
is "supported probe is not marked unsupported" 0 "$INS_UNSUPPORTED"

echo "=== evaluated probe with no matching reads ==="
SAMTOOLS_MODE=empty INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
rc=$?
out="$(<"$WORK/out")"
is "no-support probe still evaluates successfully" 0 "$rc"
has "no-support probe is reported explicitly" "$out" "NO READ SUPPORT"
is "no-support probe is counted" 1 "$INS_UNSUPPORTED"

echo "=== probe matcher failure ==="
count_probe_hits() { return 2; }
SAMTOOLS_MODE=support INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
rc=$?
out="$(<"$WORK/out")"
is "probe matcher failure is non-zero" 1 "$rc"
has "probe matcher failure is not reported as no support" "$out" "probe matching failed"
eval "$count_probe_hits_src"

echo "=== tool failure ==="
SAMTOOLS_MODE=fail INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
rc=$?
out="$(<"$WORK/out")"
is "samtools failure is non-zero" 1 "$rc"
has "samtools failure is not reported as no support" "$out" "verification was not evaluated"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
