#!/usr/bin/env bash
# Tests for indel_context.sbatch's two decision-making programs.
#
# Both are separate .awk files so this runs THE SAME code the job runs. No samtools: the
# runner installs only bcftools and tabix, and what needs testing is the CIGAR cursor
# arithmetic and the binning -- not whether samtools can read a BAM.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT="${EXTRACT:-$HERE/../indel_context_extract.awk}"
SUMMARISE="${SUMMARISE:-$HERE/../indel_context_summarise.awk}"
JOIN="${JOIN:-$HERE/../indel_context_join.awk}"
PIPELINE="${PIPELINE:-$HERE/../indel_context.sbatch}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    run_muts() {  # <file-under-test> <var-name>
        local target="$1" var="$2"
        while IFS='@' read -r label from to; do
            [[ -n "$label" ]] || continue
            cp "$target" "$WORK/mutant.awk"
            FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.awk"
            if cmp -s "$target" "$WORK/mutant.awk"; then
                printf '  ERROR   %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
            fi
            if env "$var=$WORK/mutant.awk" bash "$0" >/dev/null 2>&1; then
                printf '  SURVIVED %-51s <- nothing caught this\n' "$label"; survived=$((survived+1))
            else
                printf '  caught   %s\n' "$label"
            fi
        done
    }
    run_muts "$EXTRACT" EXTRACT <<'M1'
a deletion does not advance the reference cursor@            if (ch == "D") pos += len@            if (0) pos += len
soft clips advance the cursor@        } else if (ch ~ /[MN=X]/) {@        } else if (ch ~ /[MN=XS]/) {
insertions are not recorded@        if (ch == "I" || ch == "D") {@        if (ch == "D") {
support is not counted per position@            c[$1 SUBSEP pos]++@            c[$1 SUBSEP pos] = 1
M1
    run_muts "$JOIN" JOIN <<'M3'
a candidate with no nearby indel is counted as having one@    if (best < 0)          { none++ }@    if (0)                 { none++ }
the join runs without being told which file is which@    if (f_ind == "" || f_dep == "" || f_cand == "") {@    if (0) {
the search window is ignored@    for (d = 0; d <= win; d++) {@    for (d = 0; d <= 100000; d++) {
the header line is counted as a candidate@    if ($2 !~ /^[0-9]+$/) next@    if (0) next
M3
    run_muts "$SUMMARISE" SUMMARISE <<'M2'
the background is not binned like the observed side@FILENAME == f_bg  { b = $1 + 0; if (b > mx) b = mx; bgn[b] += $2; bgtot += $2; next }@FILENAME == f_bg  { bgn[$1 + 0] += $2; bgtot += $2; next }
files are matched by substring instead of exact path@FILENAME == f_ctx { hp[$1 SUBSEP $2]  = $3; next }@FILENAME ~ /ctx/ { hp[$1 SUBSEP $2]  = $3; next }
an empty background is not fatal@    if (bgtot == 0) {@    if (0) {
support class thresholds are inverted@        if (f >= hf)     { hi[h]++; thi++ }@        if (f < hf)      { hi[h]++; thi++ }
enrichment ignores the background@        e_all  = (bs > 0 && tot > 0) ? (n[h]  + 0) / tot / bs : 0@        e_all  = (bs > 0 && tot > 0) ? (n[h]  + 0) / tot : 0
the low-support curve is normalized by the pooled total@        e_low  = (bs > 0 && tlo > 0) ? (lo[h] + 0) / tlo / bs : 0@        e_low  = (bs > 0 && tlo > 0) ? (lo[h] + 0) / tot / bs : 0
the low-support curve counts every indel, not just slippage@        e_low  = (bs > 0 && tlo > 0) ? (lo[h] + 0) / tlo / bs : 0@        e_low  = (bs > 0 && tlo > 0) ? (n[h] + 0) / tlo / bs : 0
the variant curve is fed the slippage counts@        e_high = (bs > 0 && thi > 0) ? (hi[h] + 0) / thi / bs : 0@        e_high = (bs > 0 && thi > 0) ? (lo[h] + 0) / thi / bs : 0
M2
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

# Floor on how many assertions must execute. Raise it when adding tests; if it ever reads
# low, an assertion stopped running rather than started failing.
MIN_ASSERTIONS=67
PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
# Padding-insensitive. Asserting on a %5.1f field means the assertion changes when a value
# crosses a width boundary -- "2 ( 66.7%)" has a leading space and "3 (100.0%)" does not, so
# the same assertion passed for one and failed for the other. Spaces are DELETED, not
# squeezed: squeezing leaves the single space inside "( 66.7%)" and the needle would still
# have to guess the padding.
hasw() { local h n; h="$(printf '%s' "$2" | tr -d ' ')"; n="$(printf '%s' "$3" | tr -d ' ')"
         case "$h" in *"$n"*) ok "$1";; *) bad "$1" "contains (spaces removed): $3" "$2";; esac; }


extract() { awk -f "$EXTRACT" "$@"; }

echo "=== the reference cursor: only M/D/N/=/X consume reference ==="
# Hand-computed. A read at POS=100 with 10S20M1D30M:
#   10S -> cursor stays 100
#   20M -> cursor 100..119, ends at 120
#   1D  -> deletion recorded AT 120, then cursor 121
#   30M -> 121..150
printf 'chr1\t100\t10S20M1D30M\n' > "$WORK/one.tsv"
eq "a deletion is recorded at the cursor, after the soft clip and match" \
   "$(extract "$WORK/one.tsv")" "$(printf 'chr1\t120\t1')"

# An insertion does NOT advance the cursor, so a second op after it stays in frame:
#   POS=200, 5M2I5M1D  -> I at 205, cursor still 205, 5M -> 210, D at 210
printf 'chr1\t200\t5M2I5M1D\n' > "$WORK/two.tsv"
eq "an insertion is recorded without consuming reference" \
   "$(extract "$WORK/two.tsv" | sort -k2,2n)" "$(printf 'chr1\t205\t1\nchr1\t210\t1')"

# TWO deletions, because the cursor advance is only observable in what comes AFTER it.
# Every single-deletion fixture above passes whether or not `D` advances -- there is
# nothing downstream to shift -- and mutating the advance away survived until this existed.
#   POS=100, 10M5D10M3D10M:
#     10M -> 110 ; 5D recorded AT 110, cursor 115 ; 10M -> 125 ; 3D recorded AT 125
#   Without the advance the second deletion would land at 120.
printf 'chr1\t100\t10M5D10M3D10M\n' > "$WORK/twodel.tsv"
eq "a second deletion sits where the first one advanced the cursor to" \
   "$(extract "$WORK/twodel.tsv" | sort -k2,2n)" "$(printf 'chr1\t110\t1\nchr1\t125\t1')"

echo "=== hard clips are query-only too ==="
printf 'chr1\t300\t8H10M1D10M\n' > "$WORK/three.tsv"
eq "a hard clip does not shift the deletion position" \
   "$(extract "$WORK/three.tsv")" "$(printf 'chr1\t310\t1')"

echo "=== support is the number of reads agreeing at one position ==="
{ for _ in 1 2 3 4 5; do printf 'chr1\t100\t20M1D30M\n'; done
  printf 'chr1\t100\t20M2D30M\n'; } > "$WORK/sup.tsv"
# All six put a deletion at 120 -- differing LENGTH is still the same junction position.
eq "six reads at one junction report support 6" \
   "$(extract "$WORK/sup.tsv")" "$(printf 'chr1\t120\t6')"

echo "=== a read with no indel contributes nothing ==="
printf 'chr1\t100\t150M\nchr1\t400\t100M50S\n' > "$WORK/none.tsv"
eq "clip-only and match-only reads yield no positions" "$(extract "$WORK/none.tsv" | wc -l)" "0"

echo "=== the summariser: binning, support classes, enrichment ==="
# Known answer, all four inputs hand-written.
#   two indels in a run-12 context: one high support (20/40), one low (2/40)
#   one indel in a run-1 context, mid support (6/40 = 0.15)
printf 'chr1\t1000\t20\nchr1\t2000\t2\nchr1\t3000\t6\n' > "$WORK/indels.tsv"
printf 'chr1\t1000\t40\nchr1\t2000\t40\nchr1\t3000\t40\n' > "$WORK/depth.tsv"
printf 'chr1\t1000\t12\nchr1\t2000\t12\nchr1\t3000\t1\n'  > "$WORK/ctx.tsv"
# background: 9000 bases in run-1, 1000 in run-12
printf '1\t9000\n12\t1000\n' > "$WORK/bg.tsv"
summarise() { awk -v mx=10 -v hf=0.25 -v lf=0.10 \
      -v f_ind="$WORK/indels.tsv" -v f_dep="$WORK/depth.tsv" \
      -v f_ctx="$WORK/ctx.tsv"    -v f_bg="$WORK/bg.tsv" \
      -f "$SUMMARISE" "$WORK/indels.tsv" "$WORK/depth.tsv" "$WORK/ctx.tsv" "$WORK/bg.tsv"; }
out="$(summarise)"
has "run 12 is pooled into the >=10 bin"          "$out" ">=10"
has "totals split the three support classes"       "$out" "1 high, 1 mid, 1 low"
has "background reports its own denominator"       "$out" "10000 reference bases"
# 2 of 3 observed in a bin holding 1000/10000 = 10% of the reference -> 6.67x
has "enrichment is observed share over background share" "$out" "6.67x"
# and the run-1 bin: 1 of 3 over 9000/10000 = 90% -> 0.37x
has "an under-represented bin reads below 1x"      "$out" "0.37x"

echo "=== the three enrichment columns are computed from different populations ==="
# The pooled column feeds NO model: enr_low is the sequencing-error curve (#661/#662) and
# enr_high is the variant curve (#378). A fixture where all three agree cannot show they
# were computed separately, so this one INVERTS them: slippage sits in the long runs,
# variants in the short ones.
#
#   3 low-support indels  (2/40 = 0.05)  in a run-12 context
#   3 high-support indels (20/40 = 0.50) in a run-1 context
#   background: 9000 bases at run 1, 1000 at run 12  ->  bg_share 0.9 / 0.1
#
# Hand-computed, independent of the awk:
#   run >=10:  enr_all = (3/6)/0.1 = 5.00x   enr_low = (3/3)/0.1 = 10.00x   enr_high = 0.00x
#   run 1:     enr_all = (3/6)/0.9 = 0.56x   enr_low = 0.00x                enr_high = (3/3)/0.9 = 1.11x
{ printf 'chr1\t1000\t2\nchr1\t1100\t2\nchr1\t1200\t2\n'
  printf 'chr1\t2000\t20\nchr1\t2100\t20\nchr1\t2200\t20\n'; } > "$WORK/s2_indels.tsv"
{ for p in 1000 1100 1200 2000 2100 2200; do printf 'chr1\t%s\t40\n' "$p"; done; } > "$WORK/s2_depth.tsv"
{ printf 'chr1\t1000\t12\nchr1\t1100\t12\nchr1\t1200\t12\n'
  printf 'chr1\t2000\t1\nchr1\t2100\t1\nchr1\t2200\t1\n'; } > "$WORK/s2_ctx.tsv"
printf '1\t9000\n12\t1000\n' > "$WORK/s2_bg.tsv"
summarise2() { awk -v mx=10 -v hf=0.25 -v lf=0.10 \
      -v f_ind="$WORK/s2_indels.tsv" -v f_dep="$WORK/s2_depth.tsv" \
      -v f_ctx="$WORK/s2_ctx.tsv"    -v f_bg="$WORK/s2_bg.tsv" \
      -f "$SUMMARISE" "$WORK/s2_indels.tsv" "$WORK/s2_depth.tsv" "$WORK/s2_ctx.tsv" "$WORK/s2_bg.tsv"; }
out2="$(summarise2)"
row10="$(printf '%s' "$out2" | awk '$1 == ">=10"')"
row1="$(printf '%s' "$out2" | awk '$1 == "1"')"
hasw "the header names all three enrichment columns" "$out2" "enr_all"
hasw "the header names the low-support curve"        "$out2" "enr_low"
hasw "the header names the variant curve"            "$out2" "enr_high"
# The whole point: at run >=10 the three columns must read 5.00 / 10.00 / 0.00.
eq "long runs: pooled enrichment is 5.00x"   "$(printf '%s' "$row10" | awk '{print $7}')" "5.00x"
eq "long runs: the slippage curve is 10.00x" "$(printf '%s' "$row10" | awk '{print $8}')" "10.00x"
eq "long runs: the variant curve is 0.00x"   "$(printf '%s' "$row10" | awk '{print $9}')" "0.00x"
# Inverted at run 1, so neither column can be a copy of the other or of the pooled one.
eq "short runs: pooled enrichment is 0.56x"  "$(printf '%s' "$row1" | awk '{print $7}')" "0.56x"
eq "short runs: the slippage curve is 0.00x" "$(printf '%s' "$row1" | awk '{print $8}')" "0.00x"
eq "short runs: the variant curve is 1.11x"  "$(printf '%s' "$row1" | awk '{print $9}')" "1.11x"
# Non-vacuity: if any two columns were the same expression these would coincide.
[[ "$(printf '%s' "$row10" | awk '{print $7}')" != "$(printf '%s' "$row10" | awk '{print $8}')" ]] \
  && ok "pooled and slippage columns are not the same expression" \
  || bad "pooled and slippage columns are not the same expression" "different values" "identical"
hasw "the legend says which column #662 must fit" "$out2" "enr_low  IS the sequencing-error context curve"
hasw "the legend warns the pooled column is neither" "$out2" "enr_all  is neither"

echo "=== the archived job 21674484 reproduces all three published curves ==="
# Known answer against REAL data, not a synthetic fixture. The counts below are job
# 21674484's summary.txt verbatim -- the run the shipped curve was derived from. Rebuilding
# its four inputs from those counts and re-summarising must land on all three published
# numbers at once:
#
#   enr_all  -> the enrichment column that job actually printed  (52.06x at >=10)
#   enr_low  -> the curve shipped in sequencing_error_model.rs   (39.20x at >=10)
#   enr_high -> the variant figure #378 cites                    (60.44x at >=10)
#
# One fixture, three independently-sourced targets. Before this file printed enr_low, the
# only visible column was 52.06x while the model shipped 39.20x, and nothing connected them.
# Support is encoded by depth 40 with 20 / 2 / 6 reads -> 0.50 high, 0.05 low, 0.15 mid.
: > "$WORK/j_ind.tsv"; : > "$WORK/j_dep.tsv"; : > "$WORK/j_ctx.tsv"; : > "$WORK/j_bg.tsv"
jpos=1000
# run  n    high  low  mid  bg_share
while read -r run n hi lo mid share; do
    [[ -n "$run" ]] || continue
    [[ $((hi + lo + mid)) -eq "$n" ]] || { bad "job fixture run $run splits to n" "$n" "$((hi+lo+mid))"; }
    for pair in "20 $hi" "2 $lo" "6 $mid"; do
        set -- $pair
        for _ in $(seq 1 "$2"); do
            printf 'chr1\t%s\t%s\n' "$jpos" "$1" >> "$WORK/j_ind.tsv"
            printf 'chr1\t%s\t40\n'  "$jpos"      >> "$WORK/j_dep.tsv"
            printf 'chr1\t%s\t%s\n' "$jpos" "$run" >> "$WORK/j_ctx.tsv"
            jpos=$((jpos + 10))
        done
    done
    awk -v r="$run" -v s="$share" 'BEGIN{printf "%s\t%d\n", r, int(s * 3999990 + 0.5)}' >> "$WORK/j_bg.tsv"
done <<'JOB'
1 920 253 562 105 0.506731
2 507 113 348 46 0.264043
3 287 81 184 22 0.130416
4 179 54 108 17 0.056547
5 98 22 65 11 0.023781
6 51 19 25 7 0.007890
7 67 32 33 2 0.003392
8 56 22 28 6 0.001334
9 82 39 39 4 0.000932
10 777 270 334 173 0.004936
JOB
jout="$(awk -v mx=10 -v hf=0.25 -v lf=0.10 \
    -v f_ind="$WORK/j_ind.tsv" -v f_dep="$WORK/j_dep.tsv" \
    -v f_ctx="$WORK/j_ctx.tsv" -v f_bg="$WORK/j_bg.tsv" \
    -f "$SUMMARISE" "$WORK/j_ind.tsv" "$WORK/j_dep.tsv" "$WORK/j_ctx.tsv" "$WORK/j_bg.tsv")"
jrow() { printf '%s' "$jout" | awk -v r="$1" '$1 == r'; }
hasw "the archived job's support split is reproduced" "$jout" "3024 positions -- 905 high, 393 mid, 1726 low"
# The shipped curve, every bucket. These are the ten numbers in
# DEFAULT_INDEL_CONTEXT_CURVE; if this row ever disagrees, one of the two moved.
i=0
for expect in 0.64 0.76 0.82 1.11 1.58 1.84 5.64 12.16 24.24 39.20; do
    i=$((i+1)); label=$i; [[ "$i" -eq 10 ]] && label=">=10"
    eq "enr_low at run $label is the shipped curve's ${expect}x" \
       "$(jrow "$label" | awk '{print $8}')" "${expect}x"
done
# The other two columns, at the bucket where they diverge most -- so a single-column bug
# cannot hide behind the curve being right.
eq "enr_all at >=10 is the pooled 52.06x the job printed" "$(jrow ">=10" | awk '{print $7}')" "52.06x"
eq "enr_high at >=10 is the 60.44x variant figure"        "$(jrow ">=10" | awk '{print $9}')" "60.44x"

echo "=== an OUTDIR whose name contains a file marker does not break the matching ==="
# Job 21671697 ran with OUTDIR=/scratch/.../indelctx_21671697. bg.tsv's full path contained
# "ctx", the old `FILENAME ~ /ctx/` rule matched it before the /bg/ rule, and the background
# was read as context data -- bgtot 0, every enrichment 0.00x, and a table that looked
# complete. The local fixture used `mktemp -d`, whose path has no "ctx" in it, so the test
# could not see the bug.
CTXDIR="$WORK/indelctx_21671697"; mkdir -p "$CTXDIR"
cp "$WORK/indels.tsv" "$WORK/depth.tsv" "$WORK/ctx.tsv" "$WORK/bg.tsv" "$CTXDIR/"
out="$(awk -v mx=10 -v hf=0.25 -v lf=0.10 \
      -v f_ind="$CTXDIR/indels.tsv" -v f_dep="$CTXDIR/depth.tsv" \
      -v f_ctx="$CTXDIR/ctx.tsv"    -v f_bg="$CTXDIR/bg.tsv" \
      -f "$SUMMARISE" "$CTXDIR/indels.tsv" "$CTXDIR/depth.tsv" "$CTXDIR/ctx.tsv" "$CTXDIR/bg.tsv")"
hasw "the background still accumulates" "$out" "10000 reference bases"
hasw "and the enrichment is a number"   "$out" "6.67x"

echo "=== an empty background is fatal, not a table of 0.00x ==="
: > "$WORK/empty_bg.tsv"
if awk -v mx=10 -v hf=0.25 -v lf=0.10 \
     -v f_ind="$WORK/indels.tsv" -v f_dep="$WORK/depth.tsv" \
     -v f_ctx="$WORK/ctx.tsv"    -v f_bg="$WORK/empty_bg.tsv" \
     -f "$SUMMARISE" "$WORK/indels.tsv" "$WORK/depth.tsv" "$WORK/ctx.tsv" "$WORK/empty_bg.tsv" \
     >"$WORK/eb.out" 2>&1; then
  bad "an empty background exits non-zero" "non-zero" "it exited 0"
else ok "an empty background exits non-zero"; fi
hasw "and says the enrichments mean nothing" "$(cat "$WORK/eb.out")" "0.00x by"

echo "=== the summariser refuses to guess which file is which ==="
if awk -v mx=10 -f "$SUMMARISE" "$WORK/indels.tsv" >/dev/null 2>&1; then
  bad "missing -v paths is an error" "non-zero" "it exited 0"
else ok "missing -v paths is an error"; fi

echo "=== a zero-depth position is reported, not silently dropped ==="
printf 'chr1\t5000\t3\n' >> "$WORK/indels.tsv"
printf 'chr1\t5000\t1\n' >> "$WORK/ctx.tsv"
out="$(summarise)"
has "a position with no depth is counted as such" "$out" "1 without depth"
has "and still counts toward the total"           "$out" "totals: 4 positions"

echo "=== the job states the fork it exists to settle ==="
has "it names variant placement as one branch" "$(sed -n '/HOW TO READ IT/,+5p' "$SUMMARISE")" "#378"
has "it names the error model as the other"    "$(sed -n '/HOW TO READ IT/,+5p' "$SUMMARISE")" "error model"
guard="$(grep -c 'MEASURES' "$SUMMARISE")"
eq "and refuses to call the result a verdict" "$guard" "1"

echo "=== the join: is there a VARIANT at the clip boundary, or nothing? ==="
# Known answer, hand-written. Three candidates:
#   1000 -> an indel at 1005 with 20/40 support   -> HIGH  (a variant)
#   2000 -> an indel at 2003 with 2/40  support   -> LOW   (slippage)
#   3000 -> nothing within the window             -> NONE  (alignment difficulty)
printf 'chr1\t1005\t20\nchr1\t2003\t2\nchr1\t9999\t9\n' > "$WORK/j_indels.tsv"
printf 'chr1\t1005\t40\nchr1\t2003\t40\nchr1\t9999\t40\n' > "$WORK/j_depth.tsv"
{ printf 'contig\tpos\tside\tsupport\tdepth\tmapq0\tsupport_frac\tmapq0_frac\n'
  printf 'chr1\t1000\tL\t6\t40\t0\t0.15\t0.00\n'
  printf 'chr1\t2000\tR\t5\t40\t0\t0.12\t0.00\n'
  printf 'chr1\t3000\tL\t4\t40\t0\t0.10\t0.00\n'; } > "$WORK/j_candidates.tsv"
join_run() { awk -v win="$1" -v hf=0.25 -v lf=0.10 \
      -v f_ind="$JD/j_indels.tsv" -v f_dep="$JD/j_depth.tsv" -v f_cand="$JD/j_candidates.tsv" \
      -f "$JOIN" "$JD/j_indels.tsv" "$JD/j_depth.tsv" "$JD/j_candidates.tsv"; }
# Run it out of a directory whose name contains a marker, so substring matching would break.
JD="$WORK/indelctx_21671697"; mkdir -p "$JD"
cp "$WORK"/j_*.tsv "$JD/"
out="$(join_run 25)"
hasw "the header row is not counted as a candidate" "$out" "3 candidate sites"
hasw "two of three have an indel in range"         "$out" "2 (66.7%) have an indel"
hasw "and one has none"                            "$out" "1 (33.3%) have NONE"
hasw "the support classes are split"               "$out" "1 high-support, 0 mid, 1 low"

echo "=== the join refuses to guess which file is which ==="
# The summariser's inputs were once matched by substring, and an OUTDIR named "indelctx"
# made bg.tsv match the /ctx/ rule -- an empty background and a table of 0.00x. Requiring
# the paths explicitly is what removes the whole class.
if awk -v win=25 -f "$JOIN" "$JD/j_indels.tsv" >/dev/null 2>&1; then
  bad "missing -v paths is an error" "non-zero" "it exited 0"
else ok "missing -v paths is an error"; fi

echo "=== an indel just outside the window does not count ==="
# 3000 -> 9999 is far; narrowing the window to 2 also drops 1005 (distance 5) and
# 2003 (distance 3), leaving nothing.
out="$(join_run 2)"
hasw "a tight window finds no indels at all" "$out" "3 (100.0%) have NONE"

echo "=== the job names all three outcomes and what each implies ==="
has "high -> placement"   "$(grep -A4 'HOW TO READ IT' "$JOIN")" "#378"
has "low -> error model"  "$(grep -A4 'HOW TO READ IT' "$JOIN")" "error model"
has "none -> neither"     "$(grep -A6 'HOW TO READ IT' "$JOIN")" "Neither"

echo "=== a CANDIDATES path that does not exist is refused, not skipped ==="
guard="$(sed -n '/CANDIDATES is set to/,+2p' "$PIPELINE")"
has "the refusal says why silence would be worse" "$(sed -n '/Refusing rather than skipping/,+2p' "$PIPELINE")" "never asked"

echo "=== the sbatch passes the file paths both awks now require ==="
# The awks refuse to guess which input is which, so a caller that omits the -v flags gets
# exit 2 -- and the unit tests would not notice, because they invoke the awks directly.
# A rebase silently dropped these flags from the summariser call while the suite stayed
# green, which is exactly the gap this closes.
# Each call is isolated to its own command line -- a check against the whole file would
# pass on the OTHER awk's flags and assert nothing.
sum_call="$(awk '/indel_context_summarise\.awk/{for(i=NR-5;i<NR;i++) print a[i]; print} {a[NR]=$0}' "$PIPELINE")"
for v in f_ind f_dep f_ctx f_bg; do
    has "summariser -v $v is on its command line" "$sum_call" "$v="
done
join_call="$(awk '/indel_context_join\.awk/{for(i=NR-5;i<NR;i++) print a[i]; print} {a[NR]=$0}' "$PIPELINE")"
for v in f_ind f_dep f_cand; do
    has "join -v $v is on its command line" "$join_call" "$v="
done

echo "=== the job loads its own tools and preflights them ==="
# Job 21656030 printed a clean header and then died on `samtools: command not found` a
# minute in. samtools is a MODULE on Delta, not part of the bioinf conda env.
has "it loads the samtools module itself"  "$(grep -v '^[[:space:]]*#' "$PIPELINE")" "module load samtools"
has "and fails fast when a tool is absent" "$(sed -n '/required tool(s) not found/,+3p' "$PIPELINE")" "module load samtools"
guard="$(grep -vE '^[[:space:]]*#' "$PIPELINE" | grep -c 'command -v "\$t"')"
eq "the preflight actually checks each tool" "$guard" "1"

echo "=== the sbatch reads the BAM once, not twice ==="
# The first version called `samtools view` for the indels and `samtools depth` for the
# denominator, reading every read twice; it did not finish on a login node.
# Code lines only. The first version of this assertion counted the COMMENTS that mention
# those commands and reported 1 and 2 -- a fixture measuring the prose beside the code.
code_lines() { grep -vE '^[[:space:]]*#' "$PIPELINE" | grep -c "$1" || true; }
eq "no samtools depth pass" "$(code_lines 'samtools depth')" "0"
eq "one samtools view loop"  "$(code_lines 'samtools view')" "1"

# Guard against an assertion that never ran. A helper typo, a missing function, a `case`
# that silently matched nothing -- all of them leave PASS+FAIL short while the suite prints
# a clean summary.
RAN=$((PASS + FAIL))
if [[ "$RAN" -lt "$MIN_ASSERTIONS" ]]; then
    printf '  FAIL  only %d assertions ran, expected at least %d -- some did not execute\n' \
           "$RAN" "$MIN_ASSERTIONS"
    FAIL=$((FAIL + 1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
