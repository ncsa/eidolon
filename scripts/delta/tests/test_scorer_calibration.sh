#!/usr/bin/env bash
# Automated tests for the scorer CALIBRATION controls in sv_pipeline.sbatch:
# truvari_metric, selftest_truvari (positive control) and decoy_truvari (negative control).
#
# WHY THIS FILE EXISTS: the decoy is the only thing in the pipeline that distinguishes
# "the caller found the variant" from "the matching window is loose enough to accept
# anything". It was also the one control that could not fail a run — it always returned 0,
# its call sites ignored the status, and a MISSING summary.json was indistinguishable from
# a legitimate recall of 0, so a truvari that crashed printed
#   decoy    scoreable: PASS (shifted truth recall=0, expected 0)
# A negative control that passes because it never ran is worse than no control, because it
# is affirmatively reported as evidence. Job 20853511 is the case history: it archived an
# empty truvari_decoy_scoreable/ and emitted no decoy verdict at all.
#
# The functions are EXTRACTED VERBATIM from sv_pipeline.sbatch rather than copied here,
# for the same reason test_bnd_geometry.sh does it: a copy drifts, and then this file tests
# something that is not what runs on Delta.
#
# truvari is NOT installed on a workstation and is not needed: it is STUBBED, which is the
# only way to exercise the paths that matter — a crash that writes no summary, and a decoy
# that matches. Those cannot be provoked reliably from a real truvari. bcftools IS used, so
# the decoy VCF is really built and really sorted.
#
# Usage:
#   scripts/delta/tests/test_scorer_calibration.sh            # run the suite
#   scripts/delta/tests/test_scorer_calibration.sh --mutate   # prove the suite is non-vacuous
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PIPELINE:-$HERE/../sv_pipeline.sbatch}"   # overridden by --mutate
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# ── Mutation driver ────────────────────────────────────────────────────────────
# Each entry breaks the production code in a plausible way and names the assertion that
# must catch it. A mutation that SURVIVES means that assertion is not testing what it
# claims. `\n` in an anchor is expanded to a real newline, so an anchor can span lines —
# needed because `if ! r=$(truvari_metric ...)` appears in BOTH controls and only the
# surrounding lines tell them apart.
if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$PIPELINE" "$WORK/mutant.sbatch"
        FROM="$(printf '%b' "$from")" TO="$(printf '%b' "$to")" \
            perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sbatch"
        if cmp -s "$PIPELINE" "$WORK/mutant.sbatch"; then
            printf '  ERROR   %-46s mutation did not apply — anchor not found\n' "$label"
            survived=$((survived+1)); continue
        fi
        if PIPELINE="$WORK/mutant.sbatch" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-45s <- nothing caught this\n' "$label"
            survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
absent summary silently reads as 0@[[ -f "$1" ]] || return 1@[[ -f "$1" ]] || { echo 0; return 0; }
decoy always passes@if awk -v r="$r" 'BEGIN{exit !(r < 0.001)}'; then@if true; then
decoy tolerates a control that did not run@    if ! r=$(truvari_metric "$out/summary.json" recall); then\n        echo "ERROR: decoy control for $label DID NOT RUN@    if ! r=$(truvari_metric "$out/summary.json" recall); then\n        r=0; if false; then echo "ERROR: decoy control for $label DID NOT RUN
selftest tolerates a control that did not run@    if ! r=$(truvari_metric "$out/summary.json" recall); then\n        echo "ERROR: scorer selftest for $label DID NOT RUN@    if ! r=$(truvari_metric "$out/summary.json" recall); then\n        r=1.0; if false; then echo "ERROR: scorer selftest for $label DID NOT RUN
decoy match is only a warning@echo "ERROR: decoy control for $label matched@return 0; echo "ERROR: decoy control for $label matched
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

# Skipping is a convenience for a workstation without bcftools, never for CI — exiting 0
# here would turn "tested nothing" into a green check, the false-pass shape this suite
# exists to prevent.
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
for fn in truvari_metric selftest_truvari decoy_truvari check_denominator; do
    src="$(extract "$fn")"
    [[ -n "$src" ]] || { echo "FATAL: could not extract $fn from $PIPELINE"; exit 2; }
    eval "$src"
done

OUTDIR="$WORK"
REFERENCE="$WORK/ref.fa"       # only ever handed to the stub
TRUVARI_SIZEMAX=100000000
BND_DIST=500
: > "$REFERENCE"

# ── The truvari stub ───────────────────────────────────────────────────────────
# Behaviour is driven by $STUB_MODE so a single stub covers every path:
#   crash   -> exits non-zero and writes NO summary.json  (the regression under test)
#   empty   -> creates the output dir but no summary.json  (what job 20853511 archived)
#   json    -> writes summary.json from $STUB_RECALL / $STUB_TPBASE / $STUB_FN
truvari() {
    local out="" prev=""
    for a in "$@"; do [[ "$prev" == "-o" ]] && out="$a"; prev="$a"; done
    case "$STUB_MODE" in
        crash) echo "truvari: simulated failure" >&2; return 2 ;;
        empty) mkdir -p "$out"; echo "truvari: simulated failure" >&2; return 2 ;;
        json)
            mkdir -p "$out"
            printf '{"recall": %s, "precision": 1.0, "f1": 1.0, "TP-base": %s, "FN": %s, "FP": 0}\n' \
                "$STUB_RECALL" "${STUB_TPBASE:-0}" "${STUB_FN:-0}" > "$out/summary.json"
            return 0 ;;
    esac
}

# ── Fixtures ──────────────────────────────────────────────────────────────────
# A 4-record truth whose answers are hand-countable: 4 records, so a selftest that scores
# all of them has TP-base=4, FN=0, and check_denominator must agree.
#
# Three records are SMALLER than the 2 Mb displacement floor and one is much larger — that
# split is deliberate, because the shift has two branches and #508 only added the second.
# The 6,111,196 bp DUP is the real one from job 20853511 (chr1:51454985-57566181), the size
# class that defeated the old flat 2 Mb shift by still overlapping itself ~67%.
build_truth() {
    { printf '##fileformat=VCFv4.2\n##contig=<ID=c1,length=100000000>\n'
      printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
      printf '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="l">\n'
      printf '##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
      printf 'c1\t1000\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-500;END=1500\n'
      printf 'c1\t20000\t.\tG\t<DUP>\t60\tPASS\tSVTYPE=DUP;SVLEN=3000;END=23000\n'
      printf 'c1\t50000\t.\tG\t<INV>\t60\tPASS\tSVTYPE=INV;SVLEN=800;END=50800\n'
      printf 'c1\t51454985\t.\tG\t<DUP>\t60\tPASS\tSVTYPE=DUP;SVLEN=6111196;END=57566181\n'
    } > "$WORK/truth.vcf"
    bgzip -f -c "$WORK/truth.vcf" > "$WORK/truth.vcf.gz"
    bcftools index -f -t "$WORK/truth.vcf.gz"
}
build_truth
TRUTH="$WORK/truth.vcf.gz"

echo "── truvari_metric: an absent measurement is not a zero ──"
rm -f "$WORK/none.json"
out="$(truvari_metric "$WORK/none.json" recall)"; rc=$?
is "missing summary.json returns non-zero" 1 "$rc"
is "missing summary.json prints nothing"    "" "$out"

printf '{"recall": 0.25}\n' > "$WORK/s.json"
out="$(truvari_metric "$WORK/s.json" recall)"; rc=$?
is "present summary returns 0"        0 "$rc"
is "present summary prints the value" 0.25 "$out"

# truvari writes null (not 0) for a metric with no comparisons; that IS a real measured
# zero and must be readable, which is exactly why absence needs a different channel.
printf '{"recall": null}\n' > "$WORK/null.json"
out="$(truvari_metric "$WORK/null.json" recall)"; rc=$?
is "null metric returns 0"       0 "$rc"
is "null metric reads as 0.0"    0.0 "$out"

printf 'not json at all' > "$WORK/bad.json"
truvari_metric "$WORK/bad.json" recall >/dev/null 2>&1; rc=$?
is "unparseable summary returns non-zero" 1 "$rc"

echo "── decoy_truvari: THE REGRESSION — a control that did not run must FAIL ──"
for mode in crash empty; do
    STUB_MODE="$mode"; STUB_RECALL=0
    outp="$(decoy_truvari "$TRUTH" "dc_$mode" 2>&1)"; rc=$?
    is  "decoy($mode) returns non-zero"      1 "$rc"
    hasnt "decoy($mode) does not claim PASS"  "$outp" "PASS"
    has "decoy($mode) says it did not run"   "$outp" "DID NOT RUN"
    has "decoy($mode) surfaces truvari stderr" "$outp" "simulated failure"
done

echo "── decoy_truvari: a real measured zero still passes ──"
STUB_MODE=json; STUB_RECALL=0.0
outp="$(decoy_truvari "$TRUTH" dc_zero 2>&1)"; rc=$?
is  "decoy(recall=0) returns 0"    0 "$rc"
has "decoy(recall=0) reports PASS" "$outp" "PASS"

echo "── decoy_truvari: a matching decoy is a failure, not a warning ──"
# 0.202 is the real figure from job 20824321's flat-shift decoy.
STUB_MODE=json; STUB_RECALL=0.202
outp="$(decoy_truvari "$TRUTH" dc_match 2>&1)"; rc=$?
is   "decoy(recall=0.202) returns non-zero" 1 "$rc"
hasnt "decoy(recall=0.202) does not say PASS" "$outp" "PASS"
has  "decoy(recall=0.202) reports the value"  "$outp" "0.202"

echo "── decoy displacement exceeds the event (#508), and END travels with POS ──"
STUB_MODE=json; STUB_RECALL=0.0
decoy_truvari "$TRUTH" geom >/dev/null 2>&1
DECOY="$WORK/.decoy_geom.vcf.gz"
q() { bcftools query -f '%POS\t%INFO/END\n' "$DECOY" 2>/dev/null | awk -F'\t' -v p="$1" '$1==p'; }

# Branch 1, the 2 Mb FLOOR. The 3000 bp DUP at POS=20000 is far smaller than the floor, so
# it moves by 2000000, not by 2*3000: POS 2020000, END 23000+2000000=2023000.
is "small DUP moves by the 2Mb floor, END carried" "2020000	2023000" "$(q 2020000)"

# Branch 2, the SVLEN SCALING that #508 added, and the branch the old flat shift got wrong.
# 6111196 bp DUP: 2*6111196 = 12222392 exceeds the floor, so POS 51454985+12222392=63677377
# and END 57566181+12222392=69788573. Under the old flat 2 Mb shift this record would have
# moved to 53454985 while keeping END=57566181 — still 67% self-overlapping, hence matched.
is "large DUP moves by 2x SVLEN, END carried" "63677377	69788573" "$(q 63677377)"

neg="$(bcftools query -f '%POS\t%INFO/END\n' "$DECOY" 2>/dev/null | awk -F'\t' '$2 < $1' | wc -l)"
is "no decoy record has END before POS" 0 "$neg"
n_decoy="$(bcftools view -H "$DECOY" 2>/dev/null | wc -l)"
is "every truth record is represented in the decoy" 4 "$n_decoy"

echo "── selftest_truvari: did-not-run is distinguished from scored-zero ──"
STUB_MODE=crash; STUB_RECALL=0
outp="$(selftest_truvari "$TRUTH" st_crash 2>&1)"; rc=$?
is  "selftest(crash) returns non-zero"       1 "$rc"
has "selftest(crash) says it did not run"    "$outp" "DID NOT RUN"
hasnt "selftest(crash) does not claim PASS"  "$outp" "PASS"

STUB_MODE=json; STUB_RECALL=1.0; STUB_TPBASE=4; STUB_FN=0
outp="$(selftest_truvari "$TRUTH" st_ok 2>&1)"; rc=$?
is  "selftest(perfect, full denominator) returns 0" 0 "$rc"
has "selftest(perfect) reports PASS"                "$outp" "PASS"

# recall=1.000 over a SUBSET must still fail: a filter drops a record from base and comp
# alike, so the subset matches itself perfectly. Only check_denominator sees it.
STUB_MODE=json; STUB_RECALL=1.0; STUB_TPBASE=3; STUB_FN=0
outp="$(selftest_truvari "$TRUTH" st_subset 2>&1)"; rc=$?
is  "selftest(perfect over 3 of 4) returns non-zero" 1 "$rc"
has "selftest(subset) names the wrong denominator"   "$outp" "WRONG DENOMINATOR"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
