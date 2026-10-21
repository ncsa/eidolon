#!/usr/bin/env bash
# Tests for the indel_probability measurement (#746): indel_context_subs.awk, which finds
# mismatches by walking each read's CIGAR against the reference, and indel_probability.awk,
# which classifies indel and substitution positions by the same support rule and turns the
# low-support (error) counts into indel_probability.
#
# Both are separate .awk files, so this runs THE SAME code the job runs. The fixtures are
# small enough that every expected number is counted by hand in the comments.
#
# Usage:
#   scripts/delta/tests/test_indel_probability.sh            # run the suite
#   scripts/delta/tests/test_indel_probability.sh --mutate   # prove the suite is non-vacuous
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBS="${SUBS:-$HERE/../indel_context_subs.awk}"
PROB="${PROB:-$HERE/../indel_probability.awk}"
EXTRACT="$HERE/../indel_context_extract.awk"
PIPELINE="$HERE/../indel_context.sbatch"
RUST_MODEL="$HERE/../../../eidolon-core/src/models/sequencing_error_model.rs"
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
    run_muts "$SUBS" SUBS <<'M1'
soft clips do not consume query@        } else if (ch == "S") {@        } else if (ch == "X_NEVER") {
insertions do not consume query@        } else if (ch == "I") {@        } else if (ch == "I_NEVER") {
deletions do not advance the reference cursor@        } else if (ch == "D" || ch == "N") {@        } else if (ch == "N") {
an N in the read counts as a mismatch@    if (rb == "N" || refb == "N" || refb == "") return@    if (refb == "") return
support is not counted per read@k = c SUBSEP p; mm[k]++@k = c SUBSEP p; mm[k] = 1
aligned bases are not counted@    bases++@    bases += 0
M1
    run_muts "$PROB" PROB <<'M2'
mid-support substitutions are counted as errors@else { smid++; smev += sup_s[k] }@else { smid++; smev += sup_s[k]; sev += sup_s[k] }
the ambiguous fit leaves out ambiguous substitutions@share_amb = (iev + imev) / (iev + imev + sev + smev)@share_amb = (iev + imev) / (iev + imev + sev)
indel errors are counted by position, not by read@else if (f < lf) { ilo++; iev += sup_i[k] }@else if (f < lf) { ilo++; iev += 1 }
substitution errors are counted by position, not by read@else if (f < lf) { slo++; sev += sup_s[k] }@else if (f < lf) { slo++; sev += 1 }
the curve mean is not weighted by the background@cm += bgn[h] * cw[h]@cm += cw[h]
the fit ignores the context curve@fit = share / cmean@fit = share
a zero substitution count is not fatal@if (sev == 0 || iev == 0 || bases == 0) {@if (0) {
the pooled substitution count drops the mid class@spool += sup_s[k]@if (f >= hf || f < lf) spool += sup_s[k]
M2
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
row() { awk -F'\t' -v k="$1" '$1 == k { print $2 }' "$2"; }   # key/value lookup

# ── fixture: 60 bp of ACGT repeats, so every base is a homopolymer run of length 1 ──
REF=""; for i in $(seq 1 15); do REF="${REF}ACGT"; done
printf 'chr1\t0\t%s\n' "$REF" > "$WORK/refseq.tsv"
refsub() { echo "${REF:$(($1 - 1)):$2}"; }                     # 1-based, length
alt() { case "$1" in A) echo C;; C) echo G;; G) echo T;; T) echo A;; esac; }
mutate_at() {  # <seq> <0-based offset> -> seq with that base changed
    local s="$1" o="$2"; echo "${s:0:$o}$(alt "${s:$o:1}")${s:$((o + 1))}"
}

# 20 reads at POS 21 covering 21-30. Hand counts:
#   pos 23: 1 read mismatched  -> 1/20 = 0.05 -> low  (substitution ERROR, 1 event)
#   pos 26: reads 1-10 mismatched -> 10/20 = 0.50 -> high (substitution VARIANT)
#   pos 28: reads 11-13 mismatched -> 3/20 = 0.15 -> mid (neither)
#   read 14: 4M1D5M -> 1 bp deletion at 25 -> 1/20 -> low (indel ERROR, 1 event)
#   reads 15-20: 5M1I5M -> insertion at 26 -> 6/20 = 0.30 -> high (indel VARIANT)
BASE10="$(refsub 21 10)"
: > "$WORK/cache.tsv"
for r in $(seq 1 20); do
    seq="$BASE10"; cig="10M"
    if   [[ $r -le 10 ]]; then seq="$(mutate_at "$seq" 5)"                 # pos 26
    elif [[ $r -le 13 ]]; then seq="$(mutate_at "$seq" 7)"                 # pos 28
    elif [[ $r -eq 14 ]]; then cig="4M1D5M"; seq="$(refsub 21 4)$(refsub 26 5)"
    else cig="5M1I5M"; seq="$(refsub 21 5)T$(refsub 26 5)"
    fi
    [[ $r -eq 1 ]] && seq="$(mutate_at "$seq" 2)"                          # pos 23
    printf 'chr1\t21\t%s\t%s\n' "$cig" "$seq" >> "$WORK/cache.tsv"
done

echo "=== mismatches: found by walking the CIGAR against the reference ==="
awk -v f_ref="$WORK/refseq.tsv" -v f_bases="$WORK/bases.tsv" -f "$SUBS" \
    "$WORK/refseq.tsv" "$WORK/cache.tsv" | sort -k2,2n > "$WORK/subs.tsv"
eq "three mismatch positions"          "$(wc -l < "$WORK/subs.tsv")" "3"
eq "pos 23 has 1 supporting read"      "$(awk '$2==23{print $3}' "$WORK/subs.tsv")" "1"
eq "pos 26 has 10 supporting reads"    "$(awk '$2==26{print $3}' "$WORK/subs.tsv")" "10"
eq "pos 28 has 3 supporting reads"     "$(awk '$2==28{print $3}' "$WORK/subs.tsv")" "3"
# 13 reads x 10M + read 14's 9 aligned bases + 6 reads x 10 aligned (the I is not aligned).
eq "aligned bases counted by hand"     "$(row aligned_bases "$WORK/bases.tsv")" "199"

echo "=== the cursor: soft clips and insertions consume query, not reference ==="
# 2S8M at 21: query 3 is ref 21. Mismatch at query 5 -> ref 23.
s="$(refsub 19 10)"; s="$(mutate_at "$s" 4)"
printf 'chr1\t21\t2S8M\t%s\n' "$s" > "$WORK/c_clip.tsv"
eq "a soft clip shifts query, not reference" \
   "$(awk -v f_ref="$WORK/refseq.tsv" -v f_bases="$WORK/b.tsv" -f "$SUBS" "$WORK/refseq.tsv" "$WORK/c_clip.tsv" | cut -f2)" "23"
# 3M2I5M at 21: query 1-3 = ref 21-23, query 4-5 inserted, query 6-10 = ref 24-28.
# Mismatch at query 7 -> ref 25.
s="$(refsub 21 3)GG$(refsub 24 5)"; s="$(mutate_at "$s" 6)"
printf 'chr1\t21\t3M2I5M\t%s\n' "$s" > "$WORK/c_ins.tsv"
eq "an insertion shifts query, not reference" \
   "$(awk -v f_ref="$WORK/refseq.tsv" -v f_bases="$WORK/b.tsv" -f "$SUBS" "$WORK/refseq.tsv" "$WORK/c_ins.tsv" | cut -f2)" "25"
# 3M2D5M at 21: query 1-3 = ref 21-23, 24-25 deleted, query 4-8 = ref 26-30. Mismatch at
# query 5 -> ref 27.
s="$(refsub 21 3)$(refsub 26 5)"; s="$(mutate_at "$s" 4)"
printf 'chr1\t21\t3M2D5M\t%s\n' "$s" > "$WORK/c_del.tsv"
eq "a deletion shifts reference, not query" \
   "$(awk -v f_ref="$WORK/refseq.tsv" -v f_bases="$WORK/b.tsv" -f "$SUBS" "$WORK/refseq.tsv" "$WORK/c_del.tsv" | cut -f2)" "27"

echo "=== an N is a no-call, not an error (must-not-fire) ==="
s="$(refsub 21 10)"; s="${s:0:4}N${s:5}"
printf 'chr1\t21\t10M\t%s\n' "$s" > "$WORK/c_n.tsv"
eq "a read N is not a mismatch" \
   "$(awk -v f_ref="$WORK/refseq.tsv" -v f_bases="$WORK/b.tsv" -f "$SUBS" "$WORK/refseq.tsv" "$WORK/c_n.tsv" | wc -l)" "0"

echo "=== indel_probability: errors are the low-support class, counted by read ==="
awk -f "$EXTRACT" "$WORK/cache.tsv" | sort -k2,2n > "$WORK/indels.tsv"
# Depth at every measured position: all 20 reads span 21-30 (a D consumes reference).
for p in 23 25 26 28; do printf 'chr1\t%d\t20\n' "$p"; done > "$WORK/depth.tsv"
printf '1\t60\n' > "$WORK/bg.tsv"              # 60 reference bases, all in runs of 1
SHIPPED="0.64,0.76,0.82,1.11,1.58,1.84,5.64,12.16,24.24,39.20"
awk -v hf=0.25 -v lf=0.10 -v mx=10 -v curve="$SHIPPED" \
    -v f_ind="$WORK/indels.tsv" -v f_sub="$WORK/subs.tsv" -v f_dep="$WORK/depth.tsv" \
    -v f_bg="$WORK/bg.tsv" -v f_bases="$WORK/bases.tsv" -v f_out="$WORK/prob.tsv" \
    -f "$PROB" "$WORK/indels.tsv" "$WORK/subs.tsv" "$WORK/depth.tsv" "$WORK/bg.tsv" \
    "$WORK/bases.tsv" > "$WORK/prob.txt"
eq "exit 0 with a usable measurement"  "$?" "0"
P="$WORK/prob.tsv"
eq "indel error events"                "$(row indel_error_events "$P")" "1"
eq "substitution error events"         "$(row sub_error_events "$P")" "1"
eq "substitution events, every class"  "$(row sub_events_pooled "$P")" "14"
eq "substitution variant events"       "$(row sub_variant_events "$P")" "10"
eq "indel variant events"              "$(row indel_variant_events "$P")" "6"
eq "observed indel share of errors"    "$(row indel_share "$P")" "0.500000"
# Every base is a run of 1, so the background mean of the shipped curve is its first entry.
eq "background mean of the curve"      "$(row curve_bg_mean "$P")" "0.640000"
eq "indel_probability = share / mean"  "$(row indel_probability "$P")" "0.781250"
# Ambiguous: pos 28's 3 substitution reads, no indels. 1 / (1 + 0 + 1 + 3) = 0.2, / 0.64.
eq "ambiguous substitution events"     "$(row sub_ambiguous_events "$P")" "3"
eq "ambiguous indel events"            "$(row indel_ambiguous_events "$P")" "0"
eq "fit counting ambiguous as errors"  "$(row indel_probability_with_ambiguous "$P")" "0.312500"
eq "aligned bases carried through"     "$(row aligned_bases "$P")" "199"
eq "indel errors per aligned base"     "$(row indel_error_rate "$P")" "0.00502513"
has "the report names the denominator" "$(cat "$WORK/prob.txt")" "199 aligned bases"

echo "=== events are reads, not positions: one error position carried by several reads ==="
# An indel error at 3 of 50 reads (0.06) and a substitution error at 2 of 50 (0.04): both
# low, one position each. By reads the share is 3/(3+2) = 0.6; counting positions instead
# would give 1/2, and either kind alone by position 1/3 or 3/4.
printf 'chr1\t100\t3\t-1\n' > "$WORK/i_multi.tsv"
printf 'chr1\t200\t2\n'     > "$WORK/s_multi.tsv"
printf 'chr1\t100\t50\nchr1\t200\t50\n' > "$WORK/d_multi.tsv"
printf 'aligned_bases\t1000\n' > "$WORK/b_multi.tsv"
awk -v hf=0.25 -v lf=0.10 -v mx=10 -v curve="$SHIPPED" \
    -v f_ind="$WORK/i_multi.tsv" -v f_sub="$WORK/s_multi.tsv" -v f_dep="$WORK/d_multi.tsv" \
    -v f_bg="$WORK/bg.tsv" -v f_bases="$WORK/b_multi.tsv" -v f_out="$WORK/p_multi.tsv" \
    -f "$PROB" "$WORK/i_multi.tsv" "$WORK/s_multi.tsv" "$WORK/d_multi.tsv" "$WORK/bg.tsv" \
    "$WORK/b_multi.tsv" > /dev/null
eq "indel events are the 3 reads"      "$(row indel_error_events "$WORK/p_multi.tsv")" "3"
eq "substitution events are the 2 reads" "$(row sub_error_events "$WORK/p_multi.tsv")" "2"
eq "share is 3/(3+2)"                  "$(row indel_share "$WORK/p_multi.tsv")" "0.600000"

echo "=== nothing to measure is fatal, not a ratio of zeros ==="
: > "$WORK/subs_empty.tsv"
awk -v hf=0.25 -v lf=0.10 -v mx=10 -v curve="$SHIPPED" \
    -v f_ind="$WORK/indels.tsv" -v f_sub="$WORK/subs_empty.tsv" -v f_dep="$WORK/depth.tsv" \
    -v f_bg="$WORK/bg.tsv" -v f_bases="$WORK/bases.tsv" -v f_out="$WORK/p2.tsv" \
    -f "$PROB" "$WORK/indels.tsv" "$WORK/subs_empty.tsv" "$WORK/depth.tsv" "$WORK/bg.tsv" \
    "$WORK/bases.tsv" > "$WORK/p2.txt" 2>&1
rc=$?
eq "no substitution errors exits non-zero" "$([[ $rc -ne 0 ]] && echo yes || echo no)" "yes"
has "and says why" "$(cat "$WORK/p2.txt")" "FATAL"

echo "=== the harness's copy of the shipped curve matches the Rust default ==="
rust_curve="$(sed -n '/DEFAULT_INDEL_CONTEXT_CURVE: \[f64; 10\] = \[/,/^\];/p' "$RUST_MODEL" \
              | grep -oE '[0-9]+\.[0-9]+' | paste -sd, -)"
harness_curve="$(grep -oE 'INDEL_CURVE:-[0-9.,]+' "$PIPELINE" | cut -d- -f2)"
eq "Rust curve parsed"                  "$(echo "$rust_curve" | tr ',' '\n' | wc -l)" "10"
eq "harness curve parsed"               "$(echo "$harness_curve" | tr ',' '\n' | grep -c .)" "10"
eq "harness default equals the Rust default" "$harness_curve" "$rust_curve"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 && "$PASS" -ge 32 ]]
