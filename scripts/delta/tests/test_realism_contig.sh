#!/usr/bin/env bash
# Tests for select_contigs — the realism panel's BAM/reference agreement pre-flight.
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
the reference is ignored when auto-picking@($1 in seq) { print $1 }@{ print $1 }
a contig with no reads is accepted@$1==c && $3>0 {found=1}@$1==c {found=1}
an empty pick is treated as success@        if [[ -z "$contigs" ]]; then@        if false; then
the N mask is ignored when placing@if (s < me[i] && e > ms[i]) { bad = 1; break }@if (0) { bad = 1; break }
an all-N line is not recognised as N@alln = ($0 ~ /^[Nn]+$/)@alln = ($0 !~ /^[Nn]+$/)
a shortfall is placed anyway rather than reported@                n = (NR < k) ? NR : k@                n = k
the shortfall is not named@        [[ "$nkept" -lt "$want" ]] && \@        [[ 1 -eq 0 ]] && \
only the first shared contig is used@($1 in seq) { print $1 }@($1 in seq) { print $1; exit }
a short contig is silently dropped@            echo "  NOTE: $c ($l bp) is too short@            true "  NOTE: $c ($l bp) is too short
zero eligible contigs is not fatal@    if [[ "$n_elig" -eq 0 ]]; then@    if false; then
the remainder is never distributed@        if [[ "$idx" -lt "$extra" ]]; then want=$(( want + 1 )); fi@        if false; then want=$(( want + 1 )); fi
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
    # below pass for the wrong reason — select_contigs would be failing on a missing file
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
got="$(select_contigs "$WORK/wgs.bam" "$WORK/chr22.fai" "$WORK/wgs.idx" 2>"$WORK/err")"
eq "a whole-genome BAM against a chr22 reference picks chr22, not chr1" \
   "$got" "$(printf 'chr22\t50818468')"

echo "=== it refuses what it cannot compare ==="
# Length mismatch = different builds. This is the one that would otherwise print a clean
# table built from coordinates that mean different things on each side.
make_fai "chr22:50818469" "$WORK/wrong.fai"
select_contigs "$WORK/wgs.bam" "$WORK/wrong.fai" "$WORK/wgs.idx" >/dev/null 2>"$WORK/err" \
  && bad "a one-base length disagreement is rejected" "it was accepted" \
  || ok "a one-base length disagreement is rejected"
grep -q "different builds" "$WORK/err" \
  && ok "the length error says what a length disagreement means" \
  || bad "the length error says what a length disagreement means" "got: $(cat "$WORK/err")"

# chr22 vs 22 — the naming-convention split, which produces zero common contigs.
make_fai "22:50818468" "$WORK/nochr.fai"
select_contigs "$WORK/wgs.bam" "$WORK/nochr.fai" "$WORK/wgs.idx" >/dev/null 2>"$WORK/err" \
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
select_contigs "$WORK/empty22.bam" "$WORK/chr22.fai" "$WORK/empty22.idx" >/dev/null 2>"$WORK/err" \
  && bad "a contig with no reads is rejected, not silently measured" "it was accepted" \
  || ok "a contig with no reads is rejected, not silently measured"

echo "=== an explicit CONTIG is honoured, and still checked ==="
make_fai "chr1:248956422 chr22:50818468" "$WORK/both.fai"
got="$(select_contigs "$WORK/wgs.bam" "$WORK/both.fai" "$WORK/wgs.idx" chr1 2>"$WORK/err")"
eq "an explicit CONTIG overrides the auto-pick" "$got" "$(printf 'chr1\t248956422')"
select_contigs "$WORK/empty22.bam" "$WORK/both.fai" "$WORK/empty22.idx" chr22 >/dev/null 2>"$WORK/err" \
  && bad "an explicit CONTIG with no reads is still rejected" "it was accepted" \
  || ok "an explicit CONTIG with no reads is still rejected"

# An explicit contig that HAS reads but is absent from the reference. The length check
# would also reject this — with "different builds", which is the wrong story and sends
# whoever reads it looking for a build mismatch that isn't there. Assert the diagnostic,
# not just the exit status: a mutation removing this branch survives an exit-status-only
# test, because the next check happens to fail too.
select_contigs "$WORK/wgs.bam" "$WORK/chr22.fai" "$WORK/wgs.idx" chr1 >/dev/null 2>"$WORK/err" \
  && bad "an explicit CONTIG absent from the reference is rejected" "it was accepted" \
  || ok "an explicit CONTIG absent from the reference is rejected"
grep -q "not present in the reference" "$WORK/err" \
  && ok "that rejection says 'not in the reference', not 'different builds'" \
  || bad "that rejection says 'not in the reference', not 'different builds'" "got: $(cat "$WORK/err")"

echo "=== it does not fire on data that is fine ==="
# The must-not-fire case. A check that rejects everything also passes every test above.
got="$(select_contigs "$WORK/wgs.bam" "$WORK/both.fai" "$WORK/wgs.idx" 2>"$WORK/err")"
eq "a fully matching BAM and reference yield BOTH contigs, not just the first" \
   "$got" "$(printf 'chr1\t248956422\nchr22\t50818468')"
[[ ! -s "$WORK/err" ]] && ok "a clean pairing emits no diagnostic" \
                       || bad "a clean pairing emits no diagnostic" "got: $(cat "$WORK/err")"

echo "=== regions spread across contigs instead of piling onto the first ==="
# The simulation covers the whole reference whatever gets measured, so measuring one contig
# of a three-chromosome reference pays for three and reads one.
printf 'chr20\t64444167\nchr21\t46709983\nchr22\t50818468\n' > "$WORK/three.tsv"
res="$(place_regions "$WORK/three.tsv" 9 200000 2000000 "$WORK/r.bed" 2>"$WORK/err")"
eq "9 regions over 3 contigs places 9 across 3" "$res" "$(printf '9\t3')"
eq "each contig gets an equal share" \
   "$(cut -f1 "$WORK/r.bed" | sort | uniq -c | awk '{printf "%s:%s ", $2, $1}')" \
   "chr20:3 chr21:3 chr22:3 "

echo "=== an uneven split still places exactly what was asked for ==="
# This is the case that exposed `[[ ... ]] && want=$((want+1))`: that list returns 1 when
# the test is false, and under `set -e` a bare failing list exits the job — so every run
# with more contigs than remainder would have died here.
res="$(place_regions "$WORK/three.tsv" 10 200000 2000000 "$WORK/r.bed" 2>"$WORK/err")"
eq "10 over 3 places 10, not 9" "$res" "$(printf '10\t3')"
eq "the remainder goes to the earliest contig" \
   "$(cut -f1 "$WORK/r.bed" | sort | uniq -c | awk '{printf "%s:%s ", $2, $1}')" \
   "chr20:4 chr21:3 chr22:3 "
eq "the BED has exactly the rows it reported" "$(wc -l < "$WORK/r.bed")" "10"

echo "=== regions stay inside the contig, clear of the ends ==="
# An interval running past the contig end measures nothing there and reports it as zero.
bad_rows="$(awk -F'\t' -v m=2000000 '
    $1=="chr20" {L=64444167} $1=="chr21" {L=46709983} $1=="chr22" {L=50818468}
    $2 < m || $3 > L - m {n++} END {print n+0}' "$WORK/r.bed")"
eq "no region crosses a margin or a contig end" "$bad_rows" "0"
eq "no region is empty or inverted" \
   "$(awk -F'\t' '$3 <= $2 {n++} END {print n+0}' "$WORK/r.bed")" "0"

echo "=== short contigs are reported, not silently dropped ==="
printf 'chr20\t64444167\ntiny\t100000\n' > "$WORK/mixed.tsv"
res="$(place_regions "$WORK/mixed.tsv" 4 200000 2000000 "$WORK/r.bed" 2>"$WORK/err")"
eq "the eligible count excludes the short contig" "$res" "$(printf '4\t1')"
grep -q "tiny (100000 bp) is too short" "$WORK/err" \
  && ok "the short contig is named in the output" \
  || bad "the short contig is named in the output" "got: $(cat "$WORK/err")"
eq "all 4 regions land on the contig that can hold them" \
   "$(cut -f1 "$WORK/r.bed" | sort -u | tr '\n' ' ')" "chr20 "

echo "=== placing nothing is fatal, not an empty table ==="
# "No artifacts here" and "nothing was measured here" are indistinguishable in the table
# this feeds, which is the whole point of the job (rule 4).
printf 'tiny\t100000\n' > "$WORK/short.tsv"
place_regions "$WORK/short.tsv" 4 200000 2000000 "$WORK/r.bed" >/dev/null 2>"$WORK/err" \
  && bad "a reference with no usable contig is fatal" "it succeeded" \
  || ok "a reference with no usable contig is fatal"
grep -q "no contig is long enough" "$WORK/err" \
  && ok "the fatal error says why nothing could be placed" \
  || bad "the fatal error says why nothing could be placed" "got: $(cat "$WORK/err")"

echo "=== a single contig still works ==="
# The must-not-fire case for the distribution logic: generalising to N contigs must not
# break the 1-contig reference every other Delta job uses.
printf 'chr22\t50818468\n' > "$WORK/one.tsv"
res="$(place_regions "$WORK/one.tsv" 5 200000 2000000 "$WORK/r.bed" 2>"$WORK/err")"
eq "5 regions on 1 contig places 5" "$res" "$(printf '5\t1')"
eq "all 5 are distinct start positions" "$(cut -f2 "$WORK/r.bed" | sort -u | wc -l)" "5"
[[ ! -s "$WORK/err" ]] && ok "a single usable contig emits no diagnostic" \
                       || bad "a single usable contig emits no diagnostic" "got: $(cat "$WORK/err")"

echo "=== N blocks: the reference is read, not assumed ==="
# A REAL fasta through samtools, not a hand-written mask. The defect being guarded against
# (job 21622644) was placement disagreeing with what the reference actually contains, and a
# hand-made mask can agree with a wrong understanding just as easily as the code can.
#
# Shaped like an acrocentric chromosome: a leading N block, then sequence. 30000 is exactly
# 500 lines of 60, so the boundary is checkable to the base despite the line-width scan.
mkfa() {   # <path> <contig> <n_len> <total_len>
    local p="$1" c="$2" n="$3" t="$4"
    { echo ">$c"
      awk -v n="$n" -v t="$t" 'BEGIN{
          s=""; for(i=0;i<n;i++) s=s "N"
          b="ACGT"; for(i=n;i<t;i++) s=s substr(b, (i%4)+1, 1)
          for(i=1;i<=length(s);i+=60) print substr(s,i,60)
      }'
    } > "$p"
    samtools faidx "$p"
}
mkfa "$WORK/acro.fa" acro 30000 100000
printf 'acro\t100000\n' > "$WORK/acro.tsv"

got="$(n_mask "$WORK/acro.fa" "$WORK/acro.tsv" 500)"
eq "n_mask finds the leading N block at its exact bounds" "$got" "$(printf 'acro\t0\t30000')"

# Known answer, computed without the code under test: 30000 N's are in the file.
eq "the fixture really does hold 30000 N (known answer, not the code's opinion)" \
   "$(grep -v '^>' "$WORK/acro.fa" | tr -d '\n' | tr -cd 'N' | wc -c)" "30000"

eq "a run shorter than min_run is not masked" \
   "$(n_mask "$WORK/acro.fa" "$WORK/acro.tsv" 40000 | wc -l)" "0"

echo "=== placement avoids the hole ==="
# Without the mask this is the exact failure from job 21622644: the first windows land
# inside the N block and the panel refuses the whole run.
n_mask "$WORK/acro.fa" "$WORK/acro.tsv" 500 > "$WORK/acro.mask"
res="$(place_regions "$WORK/acro.tsv" 4 1000 5000 "$WORK/nomask.bed" 2>/dev/null)"
in_n="$(awk -F'\t' '$2 < 30000 {n++} END{print n+0}' "$WORK/nomask.bed")"
[[ "$in_n" -gt 0 ]] && ok "WITHOUT a mask, windows do land in the N block (the bug reproduces)" \
                    || bad "WITHOUT a mask, windows do land in the N block" "none did — fixture is wrong"

res="$(place_regions "$WORK/acro.tsv" 4 1000 5000 "$WORK/masked.bed" "$WORK/acro.mask" 2>"$WORK/err")"
eq "with a mask, all 4 regions are still placed" "$res" "$(printf '4\t1')"
eq "not one region overlaps the N block" \
   "$(awk -F'\t' '$2 < 30000 {n++} END{print n+0}' "$WORK/masked.bed")" "0"
eq "regions still respect the margin and the contig end" \
   "$(awk -F'\t' '$2 < 5000 || $3 > 100000 - 5000 {n++} END{print n+0}' "$WORK/masked.bed")" "0"
eq "regions are distinct" "$(cut -f2 "$WORK/masked.bed" | sort -u | wc -l)" "4"

echo "=== the must-not-fire case: a reference with no N must not be perturbed ==="
# A mask that rejects everything would pass every assertion above.
mkfa "$WORK/clean.fa" clean 0 100000
printf 'clean\t100000\n' > "$WORK/clean.tsv"
eq "n_mask emits nothing for a reference with no N" \
   "$(n_mask "$WORK/clean.fa" "$WORK/clean.tsv" 500 | wc -l)" "0"
: > "$WORK/empty.mask"
res="$(place_regions "$WORK/clean.tsv" 4 1000 5000 "$WORK/clean.bed" "$WORK/empty.mask" 2>"$WORK/err")"
eq "an empty mask still places the full count" "$res" "$(printf '4\t1')"
eq "an empty mask discards nothing" "$(grep -c 'discarded' "$WORK/err" || true)" "0"

echo "=== under-delivery is reported, not hidden ==="
# Rule 4: the median of 2 loci prints just as confidently as the median of 4.
mkfa "$WORK/mostly.fa" mostly 85000 100000
printf 'mostly\t100000\n' > "$WORK/mostly.tsv"
n_mask "$WORK/mostly.fa" "$WORK/mostly.tsv" 500 > "$WORK/mostly.mask"
res="$(place_regions "$WORK/mostly.tsv" 4 1000 5000 "$WORK/mostly.bed" "$WORK/mostly.mask" 2>"$WORK/err")"
eq "it places what it can rather than what it was asked for" "${res%%$'\t'*}" \
   "$(wc -l < "$WORK/mostly.bed" | tr -d ' ')"
[[ "${res%%$'\t'*}" -lt 4 ]] && ok "a mostly-N contig under-delivers instead of faking the count" \
                             || bad "a mostly-N contig under-delivers" "got $res"
grep -q "clear of N" "$WORK/err" \
  && ok "the shortfall is named on stderr" \
  || bad "the shortfall is named on stderr" "got: $(cat "$WORK/err")"
eq "nothing it did place is in the N block" \
   "$(awk -F'\t' '$2 < 85000 {n++} END{print n+0}' "$WORK/mostly.bed")" "0"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
