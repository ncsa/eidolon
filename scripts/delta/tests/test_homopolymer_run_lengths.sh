#!/usr/bin/env bash
# Tests for homopolymer_run_lengths.sh and its three awk programs.
#
# KNOWN ANSWER, computed independently of the code under test. The fixture is a synthetic
# contig with a 25 bp poly-A run at a known position, and reads built by slicing that
# reference — so how many reads carry 23, 25 or 27 A's is fixed by construction, not by what
# the measurement reports.
#
# THE CIGARs ARE DELIBERATELY WRONG. Every run-length read is given a plain full-length `M`
# CIGAR even when its sequence carries a 2 bp indel. That is the whole point of the tool: it
# reads run length out of the bases and must not consult the alignment. A measurement that
# quietly fell back to the CIGAR would report every read as matching the reference, and this
# fixture is what fails instead.
#
# THE CLIP CASES COVER ALL THREE BANDS, which is what makes the tool able to answer #692's
# question rather than just produce a number:
#   +50  a deletion the aligner declined to gap
#     0  a clip that would have mapped with no indel at all
#   -10  an insertion
#
#   scripts/delta/tests/test_homopolymer_run_lengths.sh            # run
#   scripts/delta/tests/test_homopolymer_run_lengths.sh --mutate   # prove non-vacuity
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRIVER="${DRIVER:-$HERE/../homopolymer_run_lengths.sh}"
LOCI_AWK="${LOCI_AWK:-$HERE/../homopolymer_loci.awk}"
MEASURE_AWK="${MEASURE_AWK:-$HERE/../homopolymer_measure.awk}"
SUMMARISE_AWK="${SUMMARISE_AWK:-$HERE/../homopolymer_summarise.awk}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# ── mutation mode ────────────────────────────────────────────────────────────
#
# Each mutant is verified to have APPLIED before its verdict is believed: a pattern that does
# not match produces a "survivor" that is really an unmodified program (rule 2). Files are
# copied, never round-tripped through a shell variable — `$(cat f)` strips trailing newlines.
if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r target label from to; do
        [[ -n "$label" ]] || continue
        case "$target" in
            measure)   src="$MEASURE_AWK";   var=MEASURE_AWK;;
            loci)      src="$LOCI_AWK";      var=LOCI_AWK;;
            summarise) src="$SUMMARISE_AWK"; var=SUMMARISE_AWK;;
            *) printf '  ERROR    %-52s unknown target %s\n' "$label" "$target"; survived=$((survived+1)); continue;;
        esac
        cp "$src" "$WORK/mutant.awk"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.awk"
        if cmp -s "$src" "$WORK/mutant.awk"; then
            printf '  ERROR    %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if env "$var=$WORK/mutant.awk" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-52s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
measure@a run reaching the read end is reported as a measurement@skip("run reaches the read end (left anchor)")@printf "%s\tRUN\t%d\n", locus, n - lref
measure@the run length is not compared to the reference@printf "%s\tRUN\t%d\n", locus, n - lref@printf "%s\tRUN\t%d\n", locus, n
measure@the clip offset is off by one@(endpos + tc - tail + 1)@(endpos + tc - tail)
measure@the 5' offset sign is inverted@(pos - lc) - (winstart + index(win, t) - 1)@(winstart + index(win, t) - 1) - (pos - lc)
measure@a clip tail absent from the window is measured anyway@if (occ == 1) {@if (occ >= 0) {
measure@deletions stop consuming reference@if (c == "M" || c == "D" || c == "N" || c == "=" || c == "X") span += num + 0@if (c == "M" || c == "N" || c == "=" || c == "X") span += num + 0
loci@the locus position is off by one@run_abs = win_start + rs - 1@run_abs = win_start + rs
loci@a run shorter than min_run is kept@if (BEST_LEN < min_run) { skipped["no run >= min_run"]++; return }@if (BEST_LEN < 1) { return }
summarise@the nonzero-delta count is wrong@if (d != 0) n_run_nonzero++@if (d != 0) n_run_nonzero += 2
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]
    exit $?
fi

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

for f in "$DRIVER" "$LOCI_AWK" "$MEASURE_AWK" "$SUMMARISE_AWK"; do
    [[ -f "$f" ]] || { echo "FATAL: missing $f" >&2; exit 1; }
done
if ! command -v samtools >/dev/null 2>&1; then
    echo "SKIP: samtools not on PATH (conda activate bioinf)"
    exit 0
fi

# ── fixture ──────────────────────────────────────────────────────────────────
#
# Park-Miller rather than awk's rand(): srand() seeding differs between gawk and mawk, and a
# fixture that is not reproducible across implementations is not one anyone can re-run.
# 16807 * (2^31-2) stays inside awk's exact-integer range; a glibc-style multiplier does not.
RUN_START=801       # 1-based
RUN_LEN=25
awk -v run_start="$RUN_START" -v run_len="$RUN_LEN" -v total=2000 '
function urand() { lcg = (16807 * lcg) % 2147483647; return lcg / 2147483647 }
BEGIN {
    lcg = 12345
    for (i = 0; i < 8; i++) urand()
    split("A|C|G|T", b, "|")
    s = ""
    for (i = 1; i <= total; i++) {
        if (i >= run_start && i < run_start + run_len) s = s "A"
        else s = s b[int(urand() * 4) + 1]
    }
    print ">ctg1"
    for (i = 1; i <= length(s); i += 60) print substr(s, i, 60)
    print s > "/dev/stderr"
}' > "$WORK/ref.fa" 2>"$WORK/refseq.txt"
samtools faidx "$WORK/ref.fa"
REFSEQ="$(cat "$WORK/refseq.txt")"
eq "the fixture reference is the expected length" "${#REFSEQ}" "2000"
eq "the planted run is exactly $RUN_LEN A's" \
   "$(printf '%s' "${REFSEQ:$((RUN_START-1)):$RUN_LEN}")" "$(printf 'A%.0s' $(seq 1 $RUN_LEN))"

sub() { printf '%s' "${REFSEQ:$(( $1 - 1 )):$2}"; }   # 1-based start, length
qual() { printf 'I%.0s' $(seq 1 "$1"); }
rec()  { printf '%s\t0\tctg1\t%s\t60\t%s\t*\t0\t0\t%s\t%s\n' "$1" "$2" "$3" "$4" "$(qual ${#4})"; }

{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ctg1\tLN:2000\n'
  # ── run-length reads: ref[780..800] + A x N + ref[826..885] ──────────────
  # CIGAR is a plain full-length M for ALL of them, including the ones whose sequence
  # carries an indel. The measurement must ignore it.
  pre="$(sub 780 21)"; post="$(sub 826 60)"
  for i in 1 2 3 4 5; do   # 5 reads at the reference length
      s="${pre}$(printf 'A%.0s' $(seq 1 25))${post}"; rec "same_$i" 780 "${#s}M" "$s"
  done
  for i in 1 2 3; do       # 3 reads two bases SHORT
      s="${pre}$(printf 'A%.0s' $(seq 1 23))${post}"; rec "short_$i" 780 "${#s}M" "$s"
  done
  for i in 1 2; do         # 2 reads two bases LONG
      s="${pre}$(printf 'A%.0s' $(seq 1 27))${post}"; rec "long_$i" 780 "${#s}M" "$s"
  done
  # A read that ENDS inside the run: length is a lower bound, not a measurement.
  s="${pre}$(printf 'A%.0s' $(seq 1 10))"; rec "truncated_1" 780 "${#s}M" "$s"

  # ── clip reads: aligned ref[780..879] (spans the run), then a 50 bp clip ──
  aln="$(sub 780 100)"
  s="${aln}$(sub 930 50)"; rec "clipdel_1"  780 "100M50S" "$s"   # 50 bp deletion  -> +50
  s="${aln}$(sub 880 50)"; rec "clipzero_1" 780 "100M50S" "$s"   # no indel        ->   0
  s="${aln}GGGGGGGGGG$(sub 880 40)"; rec "clipins_1" 780 "100M50S" "$s"  # 10 bp ins -> -10
  # A clip that is not downstream reference at all: must be refused, not assigned an offset.
  s="${aln}$(printf 'G%.0s' $(seq 1 50))"; rec "clipnovel_1" 780 "100M50S" "$s"
  # A 5' clip: ref[680..729] then an alignment at 780 skips ref[730..779], a 50 bp deletion.
  s="$(sub 680 50)${aln}"; rec "clip5del_1" 780 "50S100M" "$s"
  # A CIGAR carrying a D, so the offset depends on deletions consuming reference. Aligned
  # ref[780..829] + ref[850..899] (20 bp deleted), clip from ref[930..979] -> +30, and only
  # if the D is counted; otherwise it reads +50.
  s="$(sub 780 50)$(sub 850 50)$(sub 930 50)"; rec "dread_1" 780 "50M20D50M50S" "$s"
} > "$WORK/reads.sam"

samtools view -b "$WORK/reads.sam" > "$WORK/reads.bam" 2>/dev/null
samtools index "$WORK/reads.bam"
printf 'ctg1\t%d\t%d\n' "$((RUN_START-1))" "$((RUN_START+RUN_LEN-1))" > "$WORK/regions.bed"
# A second interval, covered by the same reads, whose longest run is well under min_run. It
# must not become a locus -- min_run is what keeps short runs out of the distribution.
printf 'ctg1\t700\t780\n' >> "$WORK/regions.bed"

echo "=== the driver runs and finds exactly one locus ==="
PAD=400 FLANK=12 MIN_RUN=12 TAIL=15 MIN_CLIP=20 \
LOCI_AWK="$LOCI_AWK" MEASURE_AWK="$MEASURE_AWK" SUMMARISE_AWK="$SUMMARISE_AWK" \
  bash "$DRIVER" --bam "$WORK/reads.bam" --reference "$WORK/ref.fa" \
    --regions "$WORK/regions.bed" --outdir "$WORK/out" > "$WORK/stdout.txt" 2>"$WORK/stderr.txt"
rc=$?
eq "the driver exits 0" "$rc" "0"
eq "exactly one locus was found" "$(grep -c . "$WORK/out/loci.tsv" 2>/dev/null || echo 0)" "1"
eq "the locus reports the planted run length" "$(cut -f3 "$WORK/out/loci.tsv")" "$RUN_LEN"
eq "the locus reports the planted base"       "$(cut -f4 "$WORK/out/loci.tsv")" "A"
eq "the locus reports the planted position"   "$(cut -f2 "$WORK/out/loci.tsv")" "$RUN_START"

M="$WORK/out/measurements.tsv"
echo "=== run lengths are read from the bases, not the CIGAR ==="
# 5 same + 3 clip reads (whose aligned part carries the reference run) = 8 at delta 0.
eq "5 reference-length reads and 6 clip reads report delta 0" \
   "$(awk -F'\t' '$2=="RUN" && $3==0' "$M" | wc -l | tr -d ' ')" "11"
eq "3 reads two bases short report delta -2" \
   "$(awk -F'\t' '$2=="RUN" && $3==-2' "$M" | wc -l | tr -d ' ')" "3"
eq "2 reads two bases long report delta +2" \
   "$(awk -F'\t' '$2=="RUN" && $3==2' "$M" | wc -l | tr -d ' ')" "2"
eq "no other run-length value appears" \
   "$(awk -F'\t' '$2=="RUN" && $3!=0 && $3!=-2 && $3!=2' "$M" | wc -l | tr -d ' ')" "0"
# Must not fire: a read ending inside the run is a lower bound and must be excluded.
eq "the read ending inside the run yields no measurement" \
   "$(awk -F'\t' '$2=="RUN"' "$M" | wc -l | tr -d ' ')" "16"
has "and it is reported as a skip, not dropped silently" \
    "$(awk -F'\t' '$2=="SKIP"{print $3}' "$M")" "run reaches the read end"

echo "=== implied indels are recovered from the clipped bases ==="
o3() { awk -F'\t' -v v="$1" '$2=="OFF" && $3=="3prime" && $4==v' "$M" | wc -l | tr -d ' '; }
eq "a 3' clip over a 50 bp deletion reads +50"   "$(o3 50)"  "1"
eq "a 3' clip needing no indel reads 0"          "$(o3 0)"   "1"
eq "a 3' clip over a 10 bp insertion reads -10"  "$(o3 -10)" "1"
# Depends on the D in 50M20D50M50S consuming reference; drop that and this reads +50.
eq "a clip after a CIGAR deletion reads +30"     "$(o3 30)"  "1"
# The 5' sign convention is mirrored, so it needs its own case: an inverted sign would
# otherwise read as a perfectly plausible number.
eq "a 5' clip over a 50 bp deletion also reads +50" \
   "$(awk -F'\t' '$2=="OFF" && $3=="5prime" && $4==50' "$M" | wc -l | tr -d ' ')" "1"
eq "exactly five clips were measured" \
   "$(awk -F'\t' '$2=="OFF"' "$M" | wc -l | tr -d ' ')" "5"
# Must not fire: a clip that is not downstream reference gets no offset at all. Assigning it
# one would invent an indel out of a read that simply does not belong there.
has "a clip absent from the reference window is refused" \
    "$(awk -F'\t' '$2=="SKIP"{print $3}' "$M")" "3' clip tail is not in the reference window"

echo "=== the report states its denominators ==="
R="$WORK/out/report.txt"
has "the report counts run-length measurements" "$(cat "$R")" "run-length measurements"
has "the report counts clip offsets"            "$(cat "$R")" "clip-offset measurements"
has "the report counts skips"                   "$(cat "$R")" "skipped read/measurement pairs"
has "the slippage band is named"                "$(cat "$R")" "1-3 (slippage)"
has "the no-indel band is named"                "$(cat "$R")" "0 (no indel implied)"
eq "5 of 16 reads are reported as differing from the reference" \
   "$(grep -oE '^  [0-9]+ of [0-9]+ reads' "$R" | head -1 | awk '{print $1" "$3}')" "5 16"

echo "=== failures are refused rather than reported as zero ==="
# Rule 4: an empty result and an unmeasurable one must not look alike.
printf 'nosuchcontig\t100\t200\n' > "$WORK/bad.bed"
if bash "$DRIVER" --bam "$WORK/reads.bam" --reference "$WORK/ref.fa" \
     --regions "$WORK/bad.bed" --outdir "$WORK/bad" >/dev/null 2>"$WORK/bad.err"; then
    bad "a BED naming an absent contig exits non-zero" "failure" "exit 0"
else
    ok "a BED naming an absent contig exits non-zero"
fi
# A region with no homopolymer must fail, not report "no slippage".
printf 'ctg1\t100\t200\n' > "$WORK/norun.bed"
if MIN_RUN=12 bash "$DRIVER" --bam "$WORK/reads.bam" --reference "$WORK/ref.fa" \
     --regions "$WORK/norun.bed" --outdir "$WORK/norun" >/dev/null 2>"$WORK/norun.err"; then
    bad "a region with no homopolymer exits non-zero" "failure" "exit 0"
else
    ok "a region with no homopolymer exits non-zero"
fi
has "and says zero loci is a failure, not a result" "$(cat "$WORK/norun.err")" "no loci survived"
if bash "$DRIVER" --bam "$WORK/reads.bam" --reference "$WORK/ref.fa" >/dev/null 2>&1; then
    bad "missing arguments are refused" "failure" "exit 0"
else
    ok "missing arguments are refused"
fi

# Floor on how many assertions must execute. If this reads low, an assertion stopped running
# rather than started failing. Raise it when adding tests.
MIN_ASSERTIONS=30
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
