#!/usr/bin/env bash
# Automated tests for the BND geometry functions in sv_pipeline.sbatch.
#
# WHY THIS FILE EXISTS: scripts/delta/*.sh determine numbers the ACCESS report cites and
# were exercised by nothing automatic (#466). Every quiet failure found in this project
# so far has been a harness reporting a metric it never checked, so a harness with no
# tests of its own is the least defensible thing in the repo.
#
# The functions are EXTRACTED VERBATIM from sv_pipeline.sbatch rather than copied here.
# A copy would drift, and then this file would be testing something that is not what runs
# on Delta — which is the same defect it exists to catch.
#
# Every fixture below has an answer computable BY HAND, independently of the code under
# test, and every assertion is accompanied by a mutation that must break it (--mutate).
#
# Usage:
#   scripts/delta/tests/test_bnd_geometry.sh            # run the suite
#   scripts/delta/tests/test_bnd_geometry.sh --mutate   # prove the suite is non-vacuous
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PIPELINE:-$HERE/../sv_pipeline.sbatch}"   # overridden by --mutate
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# ── Mutation driver ────────────────────────────────────────────────────────────
# A suite that has never been seen to fail is decoration. Each entry below breaks the
# production code in a plausible way and names the assertion that must catch it; a
# mutation that SURVIVES means that assertion is not testing what it claims.
#
# The mutation is applied to a COPY of sv_pipeline.sbatch and verified to have actually
# changed the file before the suite runs — a silently no-op edit would otherwise be
# reported as "caught" when nothing was ever mutated.
if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    # label @ literal-to-find @ replacement
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$PIPELINE" "$WORK/mutant.sbatch"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sbatch"
        if cmp -s "$PIPELINE" "$WORK/mutant.sbatch"; then
            printf '  ERROR   %-42s mutation did not apply — anchor text not found\n' "$label"
            survived=$((survived+1)); continue
        fi
        if PIPELINE="$WORK/mutant.sbatch" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-41s <- nothing caught this\n' "$label"
            survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
mates_ok accepts any pairing@function mates_ok(a,b) {@function mates_ok(a,b) { return 1;
unpaired mates not counted@if (!(t in g)) { u@if (!(t in g)) { u0
unparsable ALTs not counted@{ bad++; badex@{ badex
single-geometry NOTE always fires@if (nd==1 && parsed>=24)@if (nd>=1 && parsed>=1)
geometry check always exits 0@exit (bad+unpaired+mismatch > 0) ? 1 : 0@exit 0
scoring loop goes back to globbing@for svt in "${SV_TYPES_SCORED_IN_LOOP[@]}"; do@for typed in "$OUTDIR"/truth_sv_*.vcf.gz; do
a derived artifact stops being scored separately@BND|CNV) continue ;;@BND) continue ;;
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]
    exit $?
fi

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
is()   { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$2" "$3"; }
has()  { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
hasnt(){ case "$2" in *"$3"*) bad "$1" "must NOT contain: $3" "$2";; *) ok "$1";; esac; }

# Skipping is a convenience for a workstation without bcftools, never for CI. If the
# install step ever breaks, exiting 0 here would turn "tested nothing" into a green
# check — the false-pass shape this whole suite exists to prevent.
if ! command -v bcftools >/dev/null; then
    if [[ -n "${CI:-}" ]]; then
        echo "FATAL: bcftools is not on PATH and CI is set. Refusing to report success" >&2
        echo "  for a suite that would run zero assertions." >&2
        exit 2
    fi
    echo "SKIP: bcftools not on PATH (set CI=1 to make this fatal)"
    exit 0
fi

# ── Load the functions under test, verbatim ────────────────────────────────────
extract() { awk "/^$1\(\)/,/^}\$/" "$PIPELINE"; }
for fn in check_bnd_geometry; do
    src="$(extract "$fn")"
    [[ -n "$src" ]] || { echo "FATAL: could not extract $fn from $PIPELINE"; exit 2; }
    eval "$src"
done
OUTDIR="$WORK"   # check_bnd_geometry writes its scratch files here

# ── Fixture builder ────────────────────────────────────────────────────────────
vcf_header() {
    printf '##fileformat=VCFv4.2\n##contig=<ID=c1,length=1000000>\n'
    printf '##ALT=<ID=BND,Description="Breakend">\n'
    printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
}
rec() { printf 'c1\t%s\t.\t%s\t%s\t60\tPASS\tSVTYPE=BND\n' "$1" "$2" "$3"; }

# KNOWN-ANSWER FIXTURE. Hand-counted, independent of the implementation:
#   2 direct junctions      -> 2x t[p[ and 2x ]p]t
#   3 head-to-head junctions-> 6x t]p]
#   1 tail-to-tail junction -> 2x [p[t
#   ------------------------------------------------
#   12 records, 6 junctions; inversion-oriented = 6+2 = 8, direct = 2+2 = 4
build_mixed() {
    { vcf_header
      rec 1000 T "T[c1:2000["   ; rec 2000 A "]c1:1000]A"    # direct #1
      rec 3000 G "G[c1:4000["   ; rec 4000 C "]c1:3000]C"    # direct #2
      rec 5000 T "T]c1:6000]"   ; rec 6000 A "A]c1:5000]"    # head-to-head #1
      rec 7000 G "G]c1:8000]"   ; rec 8000 C "C]c1:7000]"    # head-to-head #2
      rec 9000 T "T]c1:10000]"  ; rec 10000 A "A]c1:9000]"   # head-to-head #3
      rec 11000 G "[c1:12000[G" ; rec 12000 C "[c1:11000[C"  # tail-to-tail #1
    } > "$WORK/mixed.vcf"
}

echo "=== check_bnd_geometry: known-answer fixture (12 records, 6 junctions) ==="
build_mixed
out="$(check_bnd_geometry "$WORK/mixed.vcf" 2>&1)"; rc=$?
is "accepts a well-formed truth (rc)"          "0"  "$rc"
is "counts t[p[ (hand-count: 2)"   "2" "$(sed -n 's/.*t\[p\[ *direct\/deletion-like *\([0-9]*\).*/\1/p' <<<"$out")"
is "counts t]p] (hand-count: 6)"   "6" "$(sed -n 's/.*t\]p\] *head-to-head *\([0-9]*\).*/\1/p' <<<"$out")"
is "counts [p[t (hand-count: 2)"   "2" "$(sed -n 's/.*\[p\[t *tail-to-tail *\([0-9]*\).*/\1/p' <<<"$out")"
is "counts ]p]t (hand-count: 2)"   "2" "$(sed -n 's/.*\]p\]t *direct\/duplication-like *\([0-9]*\).*/\1/p' <<<"$out")"
has "reports full reciprocity"     "$out" "0 unpaired, 0 mispaired"
has "reports parse coverage"       "$out" "12 record(s), 12 parsed into a form, 0 unparsable"
hasnt "does NOT flag a mixed truth as single-geometry" "$out" "NOTE:"

echo "=== check_bnd_geometry: each way the truth can be wrong ==="
# Direct forms pair two DIFFERENT shapes, so a mate given the same shape is invalid.
# A naive "both sides match" check would pass this, which is why it is tested.
sed 's/\]c1:1000\]A/A]c1:1000]/' "$WORK/mixed.vcf" > "$WORK/mispaired.vcf"
out="$(check_bnd_geometry "$WORK/mispaired.vcf" 2>&1)"; rc=$?
is  "rejects a mate with the wrong form (rc)" "1" "$rc"
has "names the mispaired junction" "$out" "first mispaired:"

grep -v $'\t6000\t' "$WORK/mixed.vcf" > "$WORK/unpaired.vcf"
out="$(check_bnd_geometry "$WORK/unpaired.vcf" 2>&1)"; rc=$?
is  "rejects a junction whose mate is absent (rc)" "1" "$rc"
has "names the unpaired junction" "$out" "(no record there)"

sed 's/\[c1:12000\[G/<BND>/' "$WORK/mixed.vcf" > "$WORK/unparsable.vcf"
out="$(check_bnd_geometry "$WORK/unparsable.vcf" 2>&1)"; rc=$?
is  "rejects an ALT that is not breakend notation (rc)" "1" "$rc"
has "quotes the unparsable ALT" "$out" "first unparsable ALT: <BND>"

echo "=== check_bnd_geometry: the pre-v3.1.0 shape is consistent but must be flagged ==="
# 14 head-to-head junctions = 28 records: internally reciprocal, so every check above
# passes. This is exactly what eidolon emitted before geometry sampling, and it must
# not be able to regress silently.
{ vcf_header
  for i in $(seq 1 14); do
      a=$((i*1000)); b=$((i*1000+500))
      rec "$a" A "A]c1:$b]"; rec "$b" C "C]c1:$a]"
  done | sort -k2,2n
} > "$WORK/all_h2h.vcf"
out="$(check_bnd_geometry "$WORK/all_h2h.vcf" 2>&1)"; rc=$?
is  "a single-geometry truth is still internally consistent (rc)" "0" "$rc"
has "flags the single-geometry regression" "$out" "share ONE geometry"

echo "=== cross-component: a derived truth artifact must not become an SVTYPE ==="
# REGRESSION GUARD, job 20745149. The scoring loop globbed truth_sv_*.vcf.gz and relied
# on a hand-maintained exclusion list, so writing truth_sv_BNDinv.vcf.gz silently
# promoted it to an SVTYPE. The harness then printed, for the same label:
#     manta_BNDinv recall=0.000  (caller emitted no BNDinv records)
#     manta_BNDinv       recall=0.738 ... TP=354  FN=126
# Neither split_bnd_by_geometry nor the scoring loop is wrong on its own — they simply
# disagreed about what a truth_sv_* file means. That is an invariant BETWEEN two
# components, and it needs its own assertion or it is nobody's test.
# Anchored to non-comment lines: the fix's own comment quotes the old glob verbatim,
# and a check that cannot tell code from a comment would fail on the correct code.
is "scoring loop does not glob the truth directory" "0" \
   "$(grep -cE '^[^#]*for +typed +in .*truth_sv_\*' "$PIPELINE")"
is "scoring loop iterates the canonical type list" "1" \
   "$(grep -cF 'for svt in "${SV_TYPES_SCORED_IN_LOOP[@]}"' "$PIPELINE")"

# The loop must cover exactly the types not handled separately. Hardcoded here on
# purpose: if someone edits SV_TYPES, this must be reviewed rather than silently follow.
SV_TYPES=(); SV_TYPES_SCORED_IN_LOOP=()
eval "$(awk '/^SV_TYPES=\(/,/^done$/' "$PIPELINE")"
is "canonical SVTYPEs"            "DEL DUP INV INS BND CNV" "${SV_TYPES[*]}"
is "types scored by the generic loop" "DEL DUP INV INS"     "${SV_TYPES_SCORED_IN_LOOP[*]}"

# Every truth_sv_<X> file the script writes must be a canonical SVTYPE or a KNOWN
# derived artifact. A new artifact that is neither shows up here, not on Delta.
# BNDinv/BNDdirect retired with the split (#507); if either reappears as a truth_sv_*
# artifact it is a genuine unclassified artifact and this must fail.
KNOWN_DERIVED=" BNDspan scoreable "
unclassified=""
# Comment lines are excluded, for the same reason the glob assertion above is anchored to
# non-comment lines: the retirement of BNDinv (#507) leaves a historical note quoting
# `truth_sv_BNDinv.vcf.gz` verbatim, and a check that cannot tell code from a comment fails
# on correct code. Only artifacts the script actually WRITES are in scope.
for nm in $(grep -vE '^[[:space:]]*#' "$PIPELINE" \
            | grep -oE 'truth_sv_[A-Za-z]+' | sed 's/^truth_sv_//' | sort -u); do
    case " ${SV_TYPES[*]} " in *" $nm "*) continue ;; esac
    case "$KNOWN_DERIVED"     in *" $nm "*) continue ;; esac
    unclassified="$unclassified $nm"
done
is "no unclassified truth_sv_* artifact" "" "$unclassified"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
