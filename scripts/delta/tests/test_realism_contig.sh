#!/usr/bin/env bash
# Tests for select_contig — the realism panel's BAM/reference agreement pre-flight.
#
# These build REAL BAMs with samtools rather than faking idxstats output, because the
# failure being guarded against is a real BAM disagreeing with a real reference, and a
# hand-written fixture can agree with a wrong understanding of the format.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB="${LIB:-$HERE/../lib_realism_contig.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if ! command -v samtools >/dev/null 2>&1; then
    source "$HERE/../lib_report.sh" 2>/dev/null || true
    setup_conda 2>/dev/null || true
    conda_activate bioinf 2>/dev/null || true
fi
command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools unavailable"; exit 0; }

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$LIB" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        if cmp -s "$LIB" "$WORK/mutant.sh"; then
            printf '  ERROR   %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if LIB="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-51s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
length disagreement is tolerated@    if [[ "$bam_len" != "$ref_len" ]]; then@    if false; then
a contig absent from the reference is accepted@    if [[ -z "$ref_len" ]]; then@    if false; then
the reference is ignored when auto-picking@$3>0 && ($1 in ref) {print $1; exit}@$3>0 {print $1; exit}
a contig with no reads is accepted@$1==c && $3>0 {found=1}@$1==c {found=1}
an empty pick is treated as success@        if [[ -z "$contig" ]]; then@        if false; then
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

source "$LIB"

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  FAIL  %s\n' "$1"; [[ -n "${2:-}" ]] && printf '        %s\n' "$2"; }
eq()   { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "want [$3] got [$2]"; }

# Build a BAM whose header declares the given contigs, with one mapped read on each named
# in $with_reads. Real samtools throughout: header lengths and idxstats come from the file.
make_bam() {  # <out.bam> <"name:len name:len..."> <"names with reads">
    local out="$1" contigs="$2" reads="$3" c n
    # Separate statement: a command substitution in the SAME `local` runs in a subshell
    # that cannot see the locals being declared beside it, and under `set -u` that is an
    # unbound-variable error which does not stop the function.
    local sam="$WORK/$(basename "$out").sam"
    : > "$sam"
    for c in $contigs; do printf '@SQ\tSN:%s\tLN:%s\n' "${c%%:*}" "${c##*:}" >> "$sam"; done
    for n in $reads; do
        printf 'r_%s\t0\t%s\t1000\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n' "$n" "$n" >> "$sam"
    done
    samtools view -b -o "$out" "$sam" 2>/dev/null
    samtools index "$out" 2>/dev/null
    # Without this, a fixture that failed to build makes every "is rejected" assertion
    # below pass for the wrong reason — select_contig would be failing on a missing file
    # rather than on the condition under test.
    [[ -s "$out" ]] || { echo "FIXTURE ERROR: $out was not built" >&2; exit 1; }
    [[ "$(samtools view -c "$out" 2>/dev/null)" -eq "$(set -- $reads; echo $#)" ]] \
        || { echo "FIXTURE ERROR: $out has the wrong read count" >&2; exit 1; }
}
idxstats() { samtools idxstats "$1" > "$2" 2>/dev/null; }
make_fai() { : > "$2"; local c; for c in $1; do printf '%s\t%s\t0\t60\t61\n' "${c%%:*}" "${c##*:}" >> "$2"; done; }

echo "=== it picks a contig present on both sides ==="
# The load-bearing case: a whole-genome BAM against a single-contig reference. Picking from
# the BAM alone returns chr1 here, which is the bug.
make_bam "$WORK/wgs.bam" "chr1:248956422 chr22:50818468" "chr1 chr22"
idxstats "$WORK/wgs.bam" "$WORK/wgs.idx"
make_fai "chr22:50818468" "$WORK/chr22.fai"
got="$(select_contig "$WORK/wgs.bam" "$WORK/chr22.fai" "$WORK/wgs.idx" 2>"$WORK/err")"
eq "a whole-genome BAM against a chr22 reference picks chr22, not chr1" \
   "$got" "$(printf 'chr22\t50818468')"

echo "=== it refuses what it cannot compare ==="
# Length mismatch = different builds. This is the one that would otherwise print a clean
# table built from coordinates that mean different things on each side.
make_fai "chr22:50818469" "$WORK/wrong.fai"
select_contig "$WORK/wgs.bam" "$WORK/wrong.fai" "$WORK/wgs.idx" >/dev/null 2>"$WORK/err" \
  && bad "a one-base length disagreement is rejected" "it was accepted" \
  || ok "a one-base length disagreement is rejected"
grep -q "different builds" "$WORK/err" \
  && ok "the length error says what a length disagreement means" \
  || bad "the length error says what a length disagreement means" "got: $(cat "$WORK/err")"

# chr22 vs 22 — the naming-convention split, which produces zero common contigs.
make_fai "22:50818468" "$WORK/nochr.fai"
select_contig "$WORK/wgs.bam" "$WORK/nochr.fai" "$WORK/wgs.idx" >/dev/null 2>"$WORK/err" \
  && bad "a naming-convention mismatch is rejected" "it was accepted" \
  || ok "a naming-convention mismatch is rejected"
grep -q "chr22 vs 22" "$WORK/err" \
  && ok "the error names the convention split as a likely cause" \
  || bad "the error names the convention split as a likely cause" "got: $(cat "$WORK/err")"
grep -q "Reference contigs" "$WORK/err" \
  && ok "the error lists both sides so the mismatch is visible" \
  || bad "the error lists both sides so the mismatch is visible" "got: $(cat "$WORK/err")"

echo "=== a declared-but-unsequenced contig is not a match ==="
# The distinction rule 4 exists for: present in the header is not the same as measurable.
# A contig with zero reads yields an empty region, which reads as "no artifacts here".
make_bam "$WORK/empty22.bam" "chr1:248956422 chr22:50818468" "chr1"
idxstats "$WORK/empty22.bam" "$WORK/empty22.idx"
select_contig "$WORK/empty22.bam" "$WORK/chr22.fai" "$WORK/empty22.idx" >/dev/null 2>"$WORK/err" \
  && bad "a contig with no reads is rejected, not silently measured" "it was accepted" \
  || ok "a contig with no reads is rejected, not silently measured"

echo "=== an explicit CONTIG is honoured, and still checked ==="
make_fai "chr1:248956422 chr22:50818468" "$WORK/both.fai"
got="$(select_contig "$WORK/wgs.bam" "$WORK/both.fai" "$WORK/wgs.idx" chr1 2>"$WORK/err")"
eq "an explicit CONTIG overrides the auto-pick" "$got" "$(printf 'chr1\t248956422')"
select_contig "$WORK/empty22.bam" "$WORK/both.fai" "$WORK/empty22.idx" chr22 >/dev/null 2>"$WORK/err" \
  && bad "an explicit CONTIG with no reads is still rejected" "it was accepted" \
  || ok "an explicit CONTIG with no reads is still rejected"

# An explicit contig that HAS reads but is absent from the reference. The length check
# would also reject this — with "different builds", which is the wrong story and sends
# whoever reads it looking for a build mismatch that isn't there. Assert the diagnostic,
# not just the exit status: a mutation removing this branch survives an exit-status-only
# test, because the next check happens to fail too.
select_contig "$WORK/wgs.bam" "$WORK/chr22.fai" "$WORK/wgs.idx" chr1 >/dev/null 2>"$WORK/err" \
  && bad "an explicit CONTIG absent from the reference is rejected" "it was accepted" \
  || ok "an explicit CONTIG absent from the reference is rejected"
grep -q "not present in the reference" "$WORK/err" \
  && ok "that rejection says 'not in the reference', not 'different builds'" \
  || bad "that rejection says 'not in the reference', not 'different builds'" "got: $(cat "$WORK/err")"

echo "=== it does not fire on data that is fine ==="
# The must-not-fire case. A check that rejects everything also passes every test above.
got="$(select_contig "$WORK/wgs.bam" "$WORK/both.fai" "$WORK/wgs.idx" 2>"$WORK/err")"
eq "a fully matching BAM and reference are accepted" "$got" "$(printf 'chr1\t248956422')"
[[ ! -s "$WORK/err" ]] && ok "a clean pairing emits no diagnostic" \
                       || bad "a clean pairing emits no diagnostic" "got: $(cat "$WORK/err")"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
