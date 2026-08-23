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
only one probe per insertion@for (k = 1; k <= 5; k++) {@for (k = 3; k <= 3; k++) {
support floor never trips@if [[ "$pct" -lt "$min_pct" ]]; then@if false; then
probes taken from the head@s = int(L * k / 6) - 14@s = 1 + 0 * k
unsupported insertions are not counted@unsupported=$((unsupported + 1))@unsupported=$((unsupported + 0))
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

# ── multi-probe behaviour (five interior probes per insertion) ────────────────────────
# A 300 bp insertion built from distinguishable blocks, so a stub can expose part of it and
# leave the rest absent. Probes land at int(300*k/6)-14 for k=1..5 -> 36, 86, 136, 186, 236.
INS300="$(awk 'BEGIN{ b="ACGTTGCAAGGCTTACCGGATCAGTCAGGT"; for(i=0;i<10;i++){
    printf "%s%s", substr("ACGT", (i%4)+1, 1) substr("TGCA", (i%4)+1, 1), substr(b, 3) } }')"
INS300="$(printf '%s' "$INS300" | cut -c1-300)"
bcftools() { [[ "${1:-}" == "query" ]] || return 2; printf 'chr1\t100\tA\tA%s\n' "$INS300"; }

# EXPOSE is the slice(s) of the insertion the stub puts into reads.
samtools() {
    case "${SAMTOOLS_MODE:-support}" in
        fail) return 2 ;;
        empty) return 0 ;;
        *) local e; for e in $EXPOSE; do
               printf 'r\t0\tchr1\t1\t60\t*\t*\t0\t0\t%s\t*\n' "$(printf '%s' "$INS300" | cut -c"$e")"
           done; return 0 ;;
    esac
}

echo "=== interior probes rescue an insertion whose exact midpoint has no read ==="
# Head and tail present, a hole straight through the middle probe at offset 136.
EXPOSE="1-120 190-300" SAMTOOLS_MODE=partial INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "midpoint hole still evaluates successfully" 0 "$rc"
is "midpoint hole is NOT counted unsupported" 0 "$INS_UNSUPPORTED"
has "midpoint hole reports partial probe support" "$out" "/5 probes"
case "$out" in *"NO READ SUPPORT"*) bad "midpoint hole is not called unsupported" "no NO READ SUPPORT" "$out";; *) ok "midpoint hole is not called unsupported";; esac

echo "=== an insertion with no interior anywhere is still fatal ==="
EXPOSE="1-40" SAMTOOLS_MODE=partial INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
out="$(<"$WORK/out")"
is "absent interior is counted unsupported" 1 "$INS_UNSUPPORTED"
has "absent interior is reported" "$out" "NO READ SUPPORT"
has "absent interior trips the support floor" "$out" "below the ${INS_SUPPORT_MIN_PCT:-95}% floor"

echo "=== a thin minority stays under the floor and is reported, not fatal ==="
# 20 insertions, one of them with no interior support: 95% >= the 95% floor.
bcftools() {
    [[ "${1:-}" == "query" ]] || return 2
    local i; for i in $(seq 1 20); do printf 'chr1\t%s\tA\tA%s\n' "$((100 * i))" "$INS300"; done
}
# Every locus fetch returns the full insertion, so all 20 are supported.
EXPOSE="1-300" SAMTOOLS_MODE=partial INS_UNSUPPORTED=0
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
out="$(<"$WORK/out")"
has "20 supported insertions are all counted" "$out" "20 of 20 planted insertion(s)"
is "20 supported insertions are not unsupported" 0 "$INS_UNSUPPORTED"

echo "=== the floor is configurable and enforced ==="
EXPOSE="1-40" SAMTOOLS_MODE=partial INS_UNSUPPORTED=0 INS_SUPPORT_MIN_PCT=100
verify_planted_ins "$truth" "$bam" "$WORK" > "$WORK/out" 2>&1
out="$(<"$WORK/out")"
has "a 100% floor rejects any unsupported insertion" "$out" "below the 100% floor"
unset INS_SUPPORT_MIN_PCT

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
