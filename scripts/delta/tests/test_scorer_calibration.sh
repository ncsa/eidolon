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
empty stratum silently dropped (pre-fix)@TYPES_UNMEASURED+=("$svt")@:
coverage-hole summary never printed@echo "  SV types with an EMPTY denominator@echo "  suppressed
empty stratum reported for every type@if [[ "$n" -eq 0 ]]; then@if [[ "$n" -ge 0 ]]; then
empty type leaves a stray truth file@rm -f "$typed" "${typed}.tbi"@:
commit check accepts any commit@elif [[ "$ver" != *"+$want"* ]]; then@elif false; then
unstamped binary accepted@elif [[ "$ver" == *"+unknown"* ]]; then@elif false; then
dirty tree not detected@[[ -n "$(git -C "$root" status --porcelain --untracked-files=no 2>/dev/null)" ]] && dirty="-dirty"@:
provenance failure is advisory only@    [[ "$rc" -eq 0 ]] && return 0@    return 0; [[ "$rc" -eq 0 ]] && return 0
prune deletes BAMs of an unscored run@    if [[ ! -f "$dir/truvari_manta_overall/summary.json" ]]; then@    if false; then
prune ignores PRUNE_BAM=0@[[ "${PRUNE_BAM:-1}" == "1" ]] || { echo "[prune] keeping BAMs (PRUNE_BAM=0)"; return 0; }@:
prune leaves the .bai files behind@rm -f "$dir"/normal.bam "$dir"/normal.bam.bai "$dir"/tumor.bam "$dir"/tumor.bam.bai@rm -f "$dir"/normal.bam "$dir"/tumor.bam
quota parser ignores the over-quota marker@gsub(/\\*/, "", used); gsub(/\\*/, "", hard)@;
quota parser accepts non-numeric values@if (used ~ /^[0-9]+$/ && hard ~ /^[0-9]+$/) { print used, hard; found = 1 }@{ print used, hard; found = 1 }
quota parser reads the soft limit as the cap@$1 == fs {@$1 == fs { $4 = $3 }
query coverage never escalates@    if [[ "$scored" -eq 0 ]]; then@    if false; then
query coverage assumes --passonly@    [[ "$nonpass" -gt 0 ]] && why="--passonly ($nonpass non-PASS)"@    why="--passonly ($nonpass non-PASS)"
query coverage chatters on a healthy stratum@    [[ "$scored" -eq "$n_query" ]] && return 0    # fully accounted for@    :
query coverage computes a negative drop@    if [[ "$scored" -gt "$n_query" ]]; then@    if false; then
peak is a constant again@        gb = 214.0 * (bp / 3100000000.0) * (cov / 30.0)@        gb = 214.0
peak ignores coverage@ * (cov / 30.0)@ * 1.0
unknown reference size reads as small@        if (bp <= 0 || cov <= 0) { print 214; exit }@        if (bp <= 0 || cov <= 0) { print 5; exit }
reference_bp ignores the .fai@    if [[ -s "$ref.fai" ]]; then@    if false; then
probe matcher ignores the reverse strand@{ for (i = 1; i <= n; i++) if (index($0, fwd[i]) || index($0, rev[i])) hits[i]++ }@{ for (i = 1; i <= n; i++) if (index($0, fwd[i])) hits[i]++ }
probe matcher drops zero-support probes@END { for (i = 1; i <= n; i++) print chrom[i], pos[i], len[i], hits[i] }@END { for (i = 1; i <= n; i++) if (hits[i] > 0) print chrom[i], pos[i], len[i], hits[i] }
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
for fn in truvari_metric selftest_truvari decoy_truvari check_denominator split_truth_by_type check_binary_provenance prune_bams parse_lfs_quota_kb count_probe_hits report_query_coverage scaled_replicate_peak_gb reference_bp; do
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

echo "── split_truth_by_type: an empty stratum is named, not dropped ──"
# The known answer is hand-counted from build_truth: DEL 1, DUP 2, INV 1, and INS/BND/CNV
# absent entirely. So exactly INS BND CNV must be reported unmeasured, in SV_TYPES order.
SV_TYPES=(DEL DUP INV INS BND CNV)
TYPES_UNMEASURED=()
mkdir -p "$WORK/split"
# NOT in $( ): the function sets a global, and a command substitution runs it in a subshell
# where that assignment is discarded — which made this assertion pass vacuously at first.
split_truth_by_type "$TRUTH" "$WORK/split" >"$WORK/split.log" 2>&1
outp="$(cat "$WORK/split.log")"
is  "planted types are counted"            "1" "$(grep -cF 'truth DEL: 1' <<<"$outp")"
is  "a type with 2 records counts 2"       "1" "$(grep -cF 'truth DUP: 2' <<<"$outp")"
is  "empty strata are exactly INS BND CNV" "INS BND CNV" "${TYPES_UNMEASURED[*]}"
has "an empty stratum is called out"       "$outp" "NOT MEASURED: 0 INS planted"
has "it is not framed as a pass"           "$outp" "not a passing stratum"
has "the count of holes is reported"       "$outp" "EMPTY denominator: INS BND CNV (3 of 6)"
# A zero-count type must leave NO file behind: the scoring loop skips on `[[ -f ]]`, and a
# stray empty file would make it run truvari against an empty truth instead.
is "no file is left for an empty type" "0" \
   "$(ls "$WORK/split"/truth_sv_INS.vcf.gz 2>/dev/null | wc -l)"
is "files remain for planted types"    "1" \
   "$(ls "$WORK/split"/truth_sv_DEL.vcf.gz 2>/dev/null | wc -l)"

echo "── split_truth_by_type: must NOT fire when every type is planted ──"
# The negative case. Most defects in this repo were things firing when they should not.
{ bcftools view -h "$TRUTH"
  printf 'c1\t70000\t.\tG\t<INS>\t60\tPASS\tSVTYPE=INS;SVLEN=300;END=70000\n'
  printf 'c1\t80000\t.\tG\t<CNV>\t60\tPASS\tSVTYPE=CNV;SVLEN=5000;END=85000\n'
  printf 'c1\t90000\t.\tN\tN[c1:95000[\t60\tPASS\tSVTYPE=BND\n'
  bcftools view -H "$TRUTH"; } > "$WORK/full.vcf"
bcftools sort -O z -o "$WORK/full.vcf.gz" "$WORK/full.vcf" 2>/dev/null
bcftools index -f -t "$WORK/full.vcf.gz"
TYPES_UNMEASURED=(); mkdir -p "$WORK/split_full"
split_truth_by_type "$WORK/full.vcf.gz" "$WORK/split_full" >"$WORK/split_full.log" 2>&1
outp="$(cat "$WORK/split_full.log")"
is    "no type reported unmeasured"        "" "${TYPES_UNMEASURED[*]}"
hasnt "no NOT MEASURED line is emitted"    "$outp" "NOT MEASURED"
hasnt "no coverage-hole summary is emitted" "$outp" "EMPTY denominator"

echo "── check_binary_provenance: a stale binary within one version must be caught ──"
# Build a throwaway repo so the function has a real git HEAD to compare against. The
# defect (#513) is that the semver half passes whenever Cargo.toml is unchanged, which is
# every commit in this project — so every case below holds the version CONSTANT and varies
# only the commit.
PROV_REPO="$WORK/provrepo"
mkdir -p "$PROV_REPO"
printf "[workspace]\nversion = '3.1.0'\n" > "$PROV_REPO/Cargo.toml"
git -C "$PROV_REPO" init -q
git -C "$PROV_REPO" add Cargo.toml
git -C "$PROV_REPO" -c user.email=t@example.com -c user.name=t commit -q -m init
PROV_SHA="$(git -C "$PROV_REPO" rev-parse --short=7 HEAD)"

prov() { ALLOW_STALE_BIN="${2:-0}" check_binary_provenance "$1" "$PROV_REPO" >/dev/null 2>&1; echo $?; }

is "matching version and commit passes"        0 "$(prov "eidolon 3.1.0+$PROV_SHA")"
is "same version, DIFFERENT commit fails"      1 "$(prov "eidolon 3.1.0+deadbee")"
is "unstamped binary (+unknown) fails"         1 "$(prov "eidolon 3.1.0+unknown")"
is "version mismatch still fails"              1 "$(prov "eidolon 3.0.0+$PROV_SHA")"
is "no stamp at all fails"                     1 "$(prov "eidolon 3.1.0")"
is "ALLOW_STALE_BIN=1 overrides a bad commit"  0 "$(prov "eidolon 3.1.0+deadbee" 1)"

# The message must name the actual defect, not just say "stale" — the semver wording sent
# a previous investigation looking for a version problem that did not exist.
outp="$(ALLOW_STALE_BIN=0 check_binary_provenance "eidolon 3.1.0+deadbee" "$PROV_REPO" 2>&1)"
has "commit mismatch names the commit"  "$outp" "DIFFERENT COMMIT"
has "commit mismatch shows both sides"  "$outp" "checkout: $PROV_SHA"

# An UNSTAMPED binary and a WRONG-COMMIT binary both fail, so asserting the exit status
# alone does not distinguish them — a mutation disabling the +unknown branch survived on
# that basis, because control fell through to the commit check and failed anyway. The two
# point at different fixes ("rebuild where git is available" vs "rebuild this checkout"),
# so the diagnosis is part of the contract.
outp="$(ALLOW_STALE_BIN=0 check_binary_provenance "eidolon 3.1.0+unknown" "$PROV_REPO" 2>&1)"
has   "unstamped binary is diagnosed as unstamped" "$outp" "no git stamp"
hasnt "unstamped is not misreported as wrong commit" "$outp" "DIFFERENT COMMIT"

# A dirty tree means the stamp does not describe what was compiled, so the binary built
# FROM that dirty tree must match and a clean-stamped binary must not.
printf "[workspace]\nversion = '3.1.0'\n# edited\n" > "$PROV_REPO/Cargo.toml"
is "dirty tree accepts the -dirty stamp"        0 "$(prov "eidolon 3.1.0+${PROV_SHA}-dirty")"
is "dirty tree rejects a clean stamp"           1 "$(prov "eidolon 3.1.0+$PROV_SHA")"
git -C "$PROV_REPO" checkout -q -- Cargo.toml

# MUST NOT FIRE: outside a git repo there is nothing to compare, so warn and pass rather
# than failing every conda/tarball install that never claimed a commit.
NOGIT="$WORK/nogit"; mkdir -p "$NOGIT"
printf "[workspace]\nversion = '3.1.0'\n" > "$NOGIT/Cargo.toml"
is "non-git checkout warns but passes" 0 \
   "$(ALLOW_STALE_BIN=0 check_binary_provenance "eidolon 3.1.0" "$NOGIT" >/dev/null 2>&1; echo $?)"

echo "── prune_bams: reclaim a scored replicate, never an unscored one ──"
# Job 20904141 died at rep 5 with the project quota full because the external cleanup jobs
# never ran at all. This moved inline, so the destructive half needs its own assertions.
mk_rep() {  # <dir> <scored: yes|no>
    rm -rf "$1"; mkdir -p "$1"
    head -c 1048576 /dev/zero > "$1/normal.bam"; : > "$1/normal.bam.bai"
    head -c 2097152 /dev/zero > "$1/tumor.bam";  : > "$1/tumor.bam.bai"
    if [[ "$2" == "yes" ]]; then
        mkdir -p "$1/truvari_manta_overall"
        printf '{"recall": 0.9}\n' > "$1/truvari_manta_overall/summary.json"
    fi
}
n_bams() { ls "$1"/normal.bam "$1"/tumor.bam 2>/dev/null | wc -l; }

mk_rep "$WORK/rep_scored" yes
outp="$(PRUNE_BAM=1 prune_bams "$WORK/rep_scored" 2>&1)"
is  "scored replicate loses its BAMs"     0 "$(n_bams "$WORK/rep_scored")"
is  "the .bai files go too"               0 "$(ls "$WORK/rep_scored"/*.bai 2>/dev/null | wc -l)"
has "it reports what it reclaimed"        "$outp" "reclaimed"
is  "the summary it keyed on survives"    1 \
    "$(ls "$WORK/rep_scored/truvari_manta_overall/summary.json" 2>/dev/null | wc -l)"

# THE DESTRUCTIVE MUST-NOT-FIRE. A run that died before stage 5 has no summary; its FASTQ is
# already pruned, so deleting the BAMs makes it permanently un-re-callable.
mk_rep "$WORK/rep_unscored" no
outp="$(PRUNE_BAM=1 prune_bams "$WORK/rep_unscored" 2>&1)"
is  "UNSCORED replicate keeps its BAMs"   2 "$(n_bams "$WORK/rep_unscored")"
has "and says why"                        "$outp" "KEEPING BAMs"

mk_rep "$WORK/rep_optout" yes
outp="$(PRUNE_BAM=0 prune_bams "$WORK/rep_optout" 2>&1)"
is  "PRUNE_BAM=0 keeps them"              2 "$(n_bams "$WORK/rep_optout")"
has "and says so"                         "$outp" "PRUNE_BAM=0"

echo "── parse_lfs_quota_kb: the real output format, including over-quota ──"
# VERBATIM from Delta (job 20904141's aftermath). The over-quota row is the one that
# matters and the one that is not guessable: values carry a trailing '*' and a grace
# column appears, so naive field-splitting mangles it. -h is NOT used in production, so
# the numbers are plain KB; the -h form below only exists to prove it is REJECTED rather
# than silently misread as kilobytes.
OVER=$'Disk quotas for prj 21649 (pid 21649):\n      Filesystem    used   bquota  blimit  bgrace   files   iquota  ilimit  igrace\n        /scratch 576458752*  524288000 576716800 6d6h28m43s   21516*  850000  935000       -'
UNDER=$'Disk quotas for prj 21649 (pid 21649):\n      Filesystem    used   bquota  blimit  bgrace   files   iquota  ilimit  igrace\n        /scratch 255852544  524288000 576716800       -   18131  850000  935000       -'
HUMAN=$'      Filesystem    used   bquota  blimit  bgrace   files\n        /scratch  549.9G*    500G    550G 6d6h28m43s   21516*'

is "over-quota row parses used and hard limit" "576458752 576716800" \
   "$(parse_lfs_quota_kb "$OVER" /scratch)"
is "under-quota row parses"                    "255852544 576716800" \
   "$(parse_lfs_quota_kb "$UNDER" /scratch)"
# MUST NOT FIRE: a human-readable row must be REJECTED, not read as kilobytes. "549.9G"
# silently taken as 549 KB would report ~550 TB free and wave through a full filesystem —
# a check that confidently returns the wrong answer, which is worse than none.
is "human-readable (-h) output is rejected"    1 \
   "$(parse_lfs_quota_kb "$HUMAN" /scratch >/dev/null 2>&1; echo $?)"
is "a filesystem not in the output is rejected" 1 \
   "$(parse_lfs_quota_kb "$OVER" /work/nvme >/dev/null 2>&1; echo $?)"
is "empty input is rejected"                   1 \
   "$(parse_lfs_quota_kb "" /scratch >/dev/null 2>&1; echo $?)"

# The arithmetic the gate depends on: 576458752 - ... is NEGATIVE here (over hard limit),
# and free must clamp to 0 rather than wrapping to a huge positive number.
kb="$(parse_lfs_quota_kb "$OVER" /scratch)"
is "over-hard-limit free space clamps to 0" 0 \
   "$(( ( ${kb% *} > ${kb#* } ? 0 : ${kb#* } - ${kb% *} ) / 1048576 ))"
kb="$(parse_lfs_quota_kb "$UNDER" /scratch)"
is "under-quota free space is 306 GB" 306 \
   "$(( ( ${kb% *} > ${kb#* } ? 0 : ${kb#* } - ${kb% *} ) / 1048576 ))"

echo "── count_probe_hits: read-level support for planted insertions ──"
# KNOWN ANSWER, hand-built. Probe P1 appears forward in one sequence, P2 only reverse
# complemented, P3 not at all. #516 is exactly the P3 case, so a probe reading zero must
# report zero rather than being lost.
P1=AAACCCGGGTTTAAACCCGGGTTTAAACCC          # 30bp
P2=CGCGATATCGCGATATCGCGATATCGCGAT          # 30bp
P3=TTTTGGGGCCCCAAAATTTTGGGGCCCCAA          # 30bp, absent
rc() { printf '%s' "$1" | rev | tr ACGT TGCA; }
printf 'chr1\t100\t500\t%s\nchr2\t200\t900\t%s\nchr3\t300\t1200\t%s\n' "$P1" "$P2" "$P3" \
  > "$WORK/probes.tsv"
{ printf 'GGGG%sGGGG\n' "$P1"
  printf 'TTTT%sTTTT\n' "$P1"
  printf 'AAAA%sAAAA\n' "$(rc "$P2")"
  printf 'ACGTACGTACGTACGTACGTACGTACGTACGT\n'; } > "$WORK/seqs.txt"

out="$(count_probe_hits "$WORK/probes.tsv" < "$WORK/seqs.txt")"
is "forward-orientation probe counts both reads" "chr1	100	500	2" \
   "$(grep -P '^chr1\t' <<<"$out")"
# MUST NOT MISS: half of all reads are reverse complemented, so a matcher that only checks
# the forward strand would silently halve every support count and could report 0 for a
# genuinely present insertion.
is "reverse-complemented probe is found"        "chr2	200	900	1" \
   "$(grep -P '^chr2\t' <<<"$out")"
# MUST NOT FIRE: an absent probe must be reported with 0, not dropped — a dropped row would
# shrink the denominator and hide exactly the #516 case this exists to catch.
is "absent probe is reported as zero, not dropped" "chr3	300	1200	0" \
   "$(grep -P '^chr3\t' <<<"$out")"
is "every probe appears in the output"          3 "$(wc -l <<<"$out")"
is "no sequences at all still reports all probes" 3 \
   "$(count_probe_hits "$WORK/probes.tsv" < /dev/null | wc -l)"

echo "── report_query_coverage: precision must not rest on an unverified denominator ──"
# THE MOTIVATING CASE (#511), campaign 20925151. manta_INS reported TP=0 FP=0, which reads as
# "the caller emitted nothing". Manta had emitted THREE INS calls, one a true detection
# (chr7:3198935 SVLEN 61 vs truth chr7:3198936 SVLEN 61), all non-PASS and so removed by
# --passonly before comparison.
outp="$(report_query_coverage manta_INS 3 0 3 2>&1)"; rc=$?
is   "wholly-filtered stratum returns 1"        1 "$rc"
has  "says NONE were scored"                    "$outp" "scored NONE of its 3 query record(s)"
has  "names --passonly with the count"          "$outp" "--passonly (3 non-PASS)"
has  "says the number describes the filter"     "$outp" "FILTER, not the caller"
has  "names it as the pipeline's own choice"    "$outp" "THIS PIPELINE'S choice"

# Wholly filtered, but NOT by --passonly: the attribution must change, not be assumed.
outp="$(report_query_coverage manta_DEL 5 0 0 2>&1)"; rc=$?
is   "wholly filtered without non-PASS still returns 1" 1 "$rc"
has  "attributes to size/other, not --passonly"  "$outp" "size/other filters"
hasnt "does not blame --passonly"                "$outp" "non-PASS"

# PARTIAL loss — the original 60-vs-58 case from job 20853511. Reported, but not escalated:
# precision is still meaningful over 58 records, it just was not 60.
outp="$(report_query_coverage manta_BND 60 58 2 2>&1)"; rc=$?
is   "partial loss returns 0"                    0 "$rc"
has  "reports both numbers"                      "$outp" "58 of 60 record(s) scored, 2 dropped"

# MUST NOT FIRE: a fully-accounted query must be SILENT. Chattering on every healthy
# comparison is how a real warning gets skimmed past.
outp="$(report_query_coverage manta_DUP 37 37 0 2>&1)"; rc=$?
is "fully-scored query returns 0"  0 "$rc"
is "fully-scored query is silent"  "" "$outp"

# MUST NOT FIRE: an empty query is the harness's own doing (it skips truvari), not a filter
# problem, and the caller-emitted-nothing message is printed elsewhere.
outp="$(report_query_coverage manta_INS 0 0 0 2>&1)"; rc=$?
is "empty query returns 0"  0 "$rc"
is "empty query is silent"  "" "$outp"

# Impossible arithmetic must be reported, never turned into a negative "dropped" count.
outp="$(report_query_coverage manta_BND 10 12 0 2>&1)"; rc=$?
is  "more scored than supplied returns 0"  0 "$rc"
has "flags it as not understood"           "$outp" "MORE than were"
hasnt "does not print a negative drop"     "$outp" "-2"

echo "── scaled_replicate_peak_gb: the gate must not over-demand on a small reference ──"
# ANCHOR, measured: 214 GB for GRCh38 (3.1 Gbp) at 30x. Everything else scales from it, because
# both FASTQ and BAM grow with genome size x coverage.
is "GRCh38 at 30x reproduces the anchor" 214 "$(scaled_replicate_peak_gb 3100000000 30)"
is "double the coverage doubles it"      428 "$(scaled_replicate_peak_gb 3100000000 60)"
is "half the genome halves it"           107 "$(scaled_replicate_peak_gb 1550000000 30)"

# THE REGRESSION THIS EXISTS FOR. A constant 214 demanded 61x what a chr22 replicate needs and
# would have refused every smoke run once /scratch passed 336 GB used — breaking the
# smoke-first workflow precisely when disk is tight.
is "chr22 at 30x is bounded by the floor, not 214" 5 "$(scaled_replicate_peak_gb 50818468 30)"
is "a tiny reference still gets the floor"          5 "$(scaled_replicate_peak_gb 100000 30)"

# MUST NOT FIRE: an unmeasurable reference must assume the WORST, not the smallest. "Unknown"
# reading as "small" would wave through a run that cannot fit — the same rule as the
# calibration controls.
is "unknown size falls back to the full anchor" 214 "$(scaled_replicate_peak_gb 0 30)"
is "unknown coverage falls back too"           214 "$(scaled_replicate_peak_gb 3100000000 0)"

echo "── reference_bp: exact from .fai, proxy from the FASTA ──"
printf 'c1\t1000\t0\t60\t61\nc2\t2500\t0\t60\t61\n' > "$WORK/r.fa.fai"
: > "$WORK/r.fa"
is "sums contig lengths from the .fai" 3500 "$(reference_bp "$WORK/r.fa")"
# No .fai (the index is not built until stage 2), so fall back to the byte count.
rm -f "$WORK/r.fa.fai"
head -c 10000 /dev/zero > "$WORK/r.fa"
is "falls back to ~0.98x the FASTA bytes" 9800 "$(reference_bp "$WORK/r.fa")"
is "a missing reference reports 0 (-> worst case)" 0 "$(reference_bp "$WORK/nope.fa")"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
