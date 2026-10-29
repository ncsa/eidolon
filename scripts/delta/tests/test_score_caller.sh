#!/usr/bin/env bash
# Tests for score_caller in sv_pipeline.sbatch (#825).
#
# WHY THIS FILE EXISTS: score_caller assembles every SV recall figure quoted in ACCESS
# §3.5–3.7, and until now nothing tested it. It makes seven decisions, each of which has
# been wrong at least once in a real run:
#   - the aggregate query drops BND and CNV calls, because the aggregate truth has none
#     (job 20675480: manta_overall FP=24 from a 24-record truth, mostly Manta's BNDs);
#   - it drops them with `-e`, so a call with no SVTYPE is kept, not silently lost;
#   - a failed query build is counted as a coverage failure, never scored unfiltered;
#   - the per-type loop follows the canonical type list, never the files that happen to
#     exist (job 20745149: a derived truth_sv_BNDinv.vcf.gz was scored as a fake type);
#   - a type the caller never emitted reads recall=0.000 instead of being scored against
#     an empty query;
#   - the BND query comes from the RAW calls, with --bnddist (the converted file has had
#     its breakends consumed);
#   - BNDspan is scored only when its calibration passed, and always with -t.
#
# score_caller is EXTRACTED VERBATIM from sv_pipeline.sbatch, as in the other suites, so
# this tests what runs on Delta. run_truvari_both is replaced by a recorder: each call is
# logged with its label, its truth file, the IDs of the records in its query and its extra
# arguments. That is exactly the decision under test -- WHAT gets scored, under which
# label, from which records -- without needing truvari. bcftools is real: the subsets are
# really built.
#
# Usage:
#   scripts/delta/tests/test_score_caller.sh            # run the suite
#   scripts/delta/tests/test_score_caller.sh --mutate   # prove the suite is non-vacuous
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PIPELINE:-$HERE/../sv_pipeline.sbatch}"   # overridden by --mutate
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Note on the "no SVTYPE" mutant: `-i 'INFO/SVTYPE!="BND" && INFO/SVTYPE!="CNV"'` is NOT a
# useful mutant. bcftools 1.19 treats a missing SVTYPE as unequal to "BND", so it keeps x1
# exactly as `-e` does (checked 2026-10-28): an equivalent mutant, not a gap. The regression
# the code's comment guards against is a positive type list, which does drop x1.
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
aggregate keeps BND calls@-e 'INFO/SVTYPE="BND" || INFO/SVTYPE="CNV"'@-e 'INFO/SVTYPE="CNV"'
aggregate keeps CNV calls@-e 'INFO/SVTYPE="BND" || INFO/SVTYPE="CNV"'@-e 'INFO/SVTYPE="BND"'
aggregate drops calls with no SVTYPE@bcftools view -e 'INFO/SVTYPE="BND" || INFO/SVTYPE="CNV"'@bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV" || INFO/SVTYPE="INS"'
a failed query build is not counted@            SCORE_COVERAGE_FAIL=$((SCORE_COVERAGE_FAIL + 1))@            :
type loop follows the files on disk@    for svt in "${SV_TYPES_SCORED_IN_LOOP[@]}"; do@    for svt in $(ls "$OUTDIR" | sed -n 's/^truth_sv_\\(.*\\)\\.vcf\\.gz$/\\1/p'); do
an empty type is scored anyway@        if [[ -f "$q" ]] && [[ "$(bcftools view -H "$q" 2>/dev/null | wc -l)" -gt 0 ]]; then@        if [[ -f "$q" ]]; then
BND query built from converted calls@bcftools view -i 'INFO/SVTYPE="BND"' -O z -o "$qb" "$raw_comp"@bcftools view -i 'INFO/SVTYPE="BND"' -O z -o "$qb" "$comp"
BND scored without --bnddist@                --bnddist "$BND_DIST"@
BNDspan scored after a failed calibration@"$OUTDIR/truth_sv_BNDspan.vcf.gz" && "$BNDSPAN_VALID" -eq 1 ]]@"$OUTDIR/truth_sv_BNDspan.vcf.gz" ]]
BNDspan scored without --typeignore@"${caller}_BNDspan" "$raw_comp" -t@"${caller}_BNDspan" "$raw_comp"
missing calls are scored anyway@    [[ -f "$comp" ]] || { echo "  ($caller: no calls to score)"; return 0; }@    :
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

for tool in bcftools bgzip; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        if [[ -n "${CI:-}" ]]; then
            echo "FATAL: $tool is not on PATH and CI is set; refusing to report success" >&2
            exit 2
        fi
        echo "SKIP: $tool not on PATH (set CI=1 to make this fatal)"
        exit 0
    fi
done

# ── Load the code under test, verbatim ─────────────────────────────────────────
extract() { awk "/^$1\(\)/,/^}\$/" "$PIPELINE"; }
for fn in score_caller comp_sample; do
    src="$(extract "$fn")"
    [[ -n "$src" ]] || { echo "FATAL: could not extract $fn from $PIPELINE"; exit 2; }
    eval "$src"
done
# The canonical type list and its loop subset, from the pipeline itself, so a new type in
# SV_TYPES reaches this suite exactly as it reaches the scoring loop.
types_src="$(awk '/^SV_TYPES=\(/{f=1} f{print} f && /^done$/{exit}' "$PIPELINE")"
[[ -n "$types_src" ]] || { echo "FATAL: could not extract SV_TYPES from $PIPELINE"; exit 2; }
eval "$types_src"

BND_DIST=500
# The recorder: one line per scoring call -- label | truth file | query IDs | extra args.
run_truvari_both() {
    local truth="$1" label="$2" q="$3"; shift 3
    local ids
    # One ID per line, joined here: bcftools ends each record with a newline even when the
    # format string has none, which would break this one-line-per-call log.
    ids="$(bcftools query -f '%ID\n' "$q" 2>/dev/null | tr '\n' ',')"
    printf '%s|%s|%s|%s\n' "$label" "$(basename "$truth")" "$ids" "$*" >> "$CALLS"
}
call_for() { grep "^$1|" "$CALLS" 2>/dev/null; }

# A VCF from "CHROM POS ID ALT INFO" lines; INFO may lack SVTYPE on purpose.
mkvcf() {
    local out="$1"; shift
    {
        printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n'
        printf '##contig=<ID=chr2,length=1000000>\n'
        printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n'
        local rec
        for rec in "$@"; do
            set -- $rec
            printf '%s\t%s\t%s\tN\t%s\t.\tPASS\t%s\tGT\t0/1\n' "$1" "$2" "$3" "$4" "$5"
        done
    } | bcftools sort -O z -o "$out" 2>/dev/null
    bcftools index -f -t "$out"
}
# One record of a type, for truth files whose content is never read by the recorder.
truth1() { mkvcf "$1" "chr1 50 t1 <$2> SVTYPE=$2;END=60"; }

# ── Fixture ─────────────────────────────────────────────────────────────────────
# Converted calls: one of each loop type, a BND (b3) and a CNV that the aggregate must
# drop, and x1 with no SVTYPE, which it must keep. No INS call at all.
setup() {
    OUTDIR="$WORK/run$1"; mkdir -p "$OUTDIR"
    CALLS="$OUTDIR/calls.log"; : > "$CALLS"
    SCORE_COVERAGE_FAIL=0; BNDSPAN_VALID=1
    TRUTH_SCOREABLE="$OUTDIR/truth_sv_scoreable.vcf.gz"
    COMP="$OUTDIR/manta.vcf.gz"; RAW="$OUTDIR/manta_raw.vcf.gz"
    mkvcf "$COMP" \
        "chr1 100 d1 <DEL> SVTYPE=DEL;END=500" \
        "chr1 1000 u1 <DUP> SVTYPE=DUP;END=1500" \
        "chr1 2000 i1 <INV> SVTYPE=INV;END=2500" \
        "chr1 3000 x1 <INS> END=3000" \
        "chr1 4000 b3 N[chr2:900[ SVTYPE=BND" \
        "chr1 5000 c1 <CNV> SVTYPE=CNV;END=6000"
    # Raw calls: the breakends the converter consumed (b1, b2), and NOT b3, so a BND
    # query built from the wrong file is told apart by its IDs.
    mkvcf "$RAW" \
        "chr1 100 d1 <DEL> SVTYPE=DEL;END=500" \
        "chr1 7000 b1 N[chr2:100[ SVTYPE=BND" \
        "chr1 8000 b2 ]chr2:200]N SVTYPE=BND"
    truth1 "$TRUTH_SCOREABLE" DEL
    local t
    for t in DEL DUP INV INS BND BNDspan; do truth1 "$OUTDIR/truth_sv_$t.vcf.gz" DEL; done
    # The artifact that became a fake type in job 20745149. It must never be scored.
    truth1 "$OUTDIR/truth_sv_BNDinv.vcf.gz" BND
}

echo "=== a full scoring pass: what is scored, under which label, from which records ==="
setup 1
out="$(score_caller manta "$COMP" "$RAW" 2>&1)"
overall="$(call_for manta_overall)"
has   "the aggregate is scored"                      "$overall" "manta_overall|truth_sv_scoreable.vcf.gz|"
for id in d1 u1 i1; do has "the aggregate keeps $id" "$overall" "$id,"; done
has   "the aggregate keeps a call with no SVTYPE (x1)" "$overall" "x1,"
hasnt "the aggregate drops the BND call (b3)"        "$overall" "b3,"
hasnt "the aggregate drops the CNV call (c1)"        "$overall" "c1,"
has   "it reports what it excluded"                  "$out" "aggregate query: 4 of 6 manta call(s); 2 BND/CNV call(s) excluded"
is    "DEL is scored from DEL calls only"            "$(call_for manta_DEL | cut -d'|' -f3)" "d1,"
is    "DUP is scored from DUP calls only"            "$(call_for manta_DUP | cut -d'|' -f3)" "u1,"
is    "INV is scored from INV calls only"            "$(call_for manta_INV | cut -d'|' -f3)" "i1,"
is    "an INS the caller never emitted is not scored" "$(call_for manta_INS)" ""
has   "it reads recall=0 for INS instead"            "$out" "manta_INS recall=0.000  (caller emitted no INS records)"
[[ ! -e "$OUTDIR/manta_calls_INS.vcf.gz" ]] && ok "and leaves no empty INS query behind" \
    || bad "and leaves no empty INS query behind" "no file" "$OUTDIR/manta_calls_INS.vcf.gz exists"
is    "a stray truth_sv_BNDinv file is never scored" "$(call_for manta_BNDinv)" ""
bnd="$(call_for manta_BND)"
is    "BND is scored from the RAW breakends"         "$(cut -d'|' -f3 <<< "$bnd")" "b1,b2,"
has   "with --bnddist"                               "$(cut -d'|' -f4 <<< "$bnd")" "--bnddist 500"
span="$(call_for manta_BNDspan)"
is    "BNDspan is scored from the unfiltered raw calls" "$(cut -d'|' -f3 <<< "$span")" "d1,b1,b2,"
is    "with --typeignore"                            "$(cut -d'|' -f4 <<< "$span")" "-t"
is    "no coverage failure on a clean pass"          "$SCORE_COVERAGE_FAIL" "0"

echo "=== BNDspan is not scored when its calibration failed ==="
setup 2; BNDSPAN_VALID=0
out="$(score_caller manta "$COMP" "$RAW" 2>&1)"
is    "no BNDspan scoring call"                      "$(call_for manta_BNDspan)" ""
has   "and it says why"                              "$out" "manta_BNDspan: not scored (representation bridge failed calibration)"

echo "=== no scoreable truth: the aggregate is skipped, the per-type scores still run ==="
setup 3; rm -f "$TRUTH_SCOREABLE" "$TRUTH_SCOREABLE.tbi"
out="$(score_caller manta "$COMP" "$RAW" 2>&1)"
is    "no aggregate scoring call"                    "$(call_for manta_overall)" ""
has   "and it says so"                               "$out" "(manta: no truvari-scoreable truth types — aggregate skipped)"
is    "DEL is still scored"                          "$(call_for manta_DEL | cut -d'|' -f3)" "d1,"

echo "=== a failed aggregate query build is counted, never scored unfiltered ==="
setup 4
mkdir -p "$OUTDIR/manta_calls_scoreable.vcf.gz"   # a directory: bcftools cannot write it
score_caller manta "$COMP" "$RAW" > "$OUTDIR/out.txt" 2>&1
out="$(cat "$OUTDIR/out.txt")"
is    "no aggregate scoring call"                    "$(call_for manta_overall)" ""
is    "the failure is counted"                       "$SCORE_COVERAGE_FAIL" "1"
has   "and named"                                    "$out" "ERROR: manta: could not build the aggregate query subset"

echo "=== a caller with no call set is reported, not scored ==="
setup 5
out="$(score_caller delly "$OUTDIR/delly_missing.vcf.gz" 2>&1)"; rc=$?
is    "it returns 0"                                 "$rc" "0"
has   "and says there was nothing to score"          "$out" "(delly: no calls to score)"
is    "no scoring call at all"                       "$(cat "$CALLS")" ""

# Floor on assertions executed: a subshell or an early exit that skips a block would
# otherwise pass with a smaller count. Raise it when adding tests.
MIN_ASSERTIONS=31
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi
printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
