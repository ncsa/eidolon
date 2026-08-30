#!/usr/bin/env bash
# Regression tests for the DEL read-evidence gate in sv_pipeline.sbatch (#590, #607).
#
# Functions are extracted verbatim from the production pipeline. Command stubs make this
# independent of bcftools, samtools, a BAM and a reference. The cases under test are the
# junction coordinate arithmetic (a known answer), the three-way realized/unrealized/uncovered
# split, and whether a tool failure stays distinguishable from "no read support".
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
            printf '  ERROR   %-52s mutation did not apply\n' "$label"
            survived=$((survived + 1))
            continue
        fi
        if PIPELINE="$WORK/mutant.sbatch" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-51s <- nothing caught this\n' "$label"
            survived=$((survived + 1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
searches only the anchor breakpoint@for lo in "$p" "$(( p + dlen ))"; do@for lo in "$p"; do
junction probe joins the wrong flank@printf '%s\t%s\t%s\t%s\n' "$c" "$p" "$dlen" "$left$right"   >> "$jout"@printf '%s\t%s\t%s\t%s\n' "$c" "$p" "$dlen" "$left$delhead" >> "$jout"
left flank is off by one@"$c:$(( p - 14 ))-$p"@"$c:$(( p - 15 ))-$p"
right flank starts at END not END+1@"$c:$(( end + 1 ))-$(( end + 15 ))"@"$c:$end-$(( end + 14 ))"
uncovered locus is scored as a failure@elif [[ "${chits:-0}" -eq 0 ]]; then@elif false; then
unprobeable records are not counted@DEL_SKIPPED=$(( DEL_SKIPPED + 1 ))@:
samtools failures are swallowed@if ! samtools view "$bam" "$c:$lo-$end"@if samtools view "$bam" "$c:$lo-$end"
row misalignment is ignored@if [[ "$misaligned" -ne 0 ]]; then@if false; then
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]
    exit $?
fi

extract() { awk "/^$1\(\)/,/^}\$/" "$PIPELINE"; }
for fn in count_probe_hits build_del_probes verify_planted_del; do
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
hasnt() { case "$2" in *"$3"*) bad "$1" "does NOT contain: $3" "$2";; *) ok "$1";; esac; }

# Hand-computable fixture. One symbolic <DEL>: anchor POS=1000, END=2000, so the deleted span
# is 1001..2000 and the three flanks the prober must fetch are fully determined:
#     left    chr1:986-1000    (15 bp ending AT the anchor)
#     right   chr1:2001-2015   (15 bp starting after END)
#     delhead chr1:1001-1015   (first 15 deleted bases)
LEFT="ACGTACGTACGTACG"
RIGHT="TTTTGGGGCCCCAAA"
DELHEAD="GATTACAGATTACAG"
JUNCTION="$LEFT$RIGHT"
CONTROL="$LEFT$DELHEAD"

bcftools() {
    [[ "${1:-}" == "query" ]] || return 2
    printf '%s\n' "${BCF_ROWS-$'chr1\t1000\tA\t<DEL>\t2000'}"
}
samtools() {
    if [[ "${1:-}" == "faidx" ]]; then
        [[ "${FAIDX_MODE:-ok}" == "fail" ]] && return 2
        printf '%s\n' "$3 $4 $5" >> "$WORK/faidx_regions"
        printf '>%s\n%s\n>%s\n%s\n>%s\n%s\n' "$3" "$LEFT" "$4" "$RIGHT" "$5" "$DELHEAD"
        return 0
    fi
    if [[ "${1:-}" == "view" ]]; then
        # $3 is the region. The anchor window covers POS=1000, the far window covers END=2000.
        local win="near"
        case "${3:-}" in *:1*-3000|*:1000-*) win="far" ;; esac
        [[ "${3:-}" == *"-3000" ]] && win="far"
        case "${VIEW_MODE:-junction}" in
            fail)      return 2 ;;
            uncovered) return 0 ;;
            junction)  printf 'r1\t0\tchr1\t1\t60\t30M\t*\t0\t0\tCCC%sCCC\t*\n' "$JUNCTION" ;;
            through)   printf 'r1\t0\tchr1\t1\t60\t30M\t*\t0\t0\tCCC%sCCC\t*\n' "$CONTROL" ;;
            # Job 21532276: for a large deletion the junction-bearing read anchors at END, and
            # nothing is at POS. The probe searched POS only and called three good deletions
            # unrealized. Must pass.
            far_only)
                if [[ "$win" == "far" ]]; then
                    printf 'r1\t0\tchr1\t1\t60\t30M\t*\t0\t0\tCCC%sCCC\t*\n' "$JUNCTION"
                else
                    printf 'r2\t0\tchr1\t1\t60\t30M\t*\t0\t0\tCCC%sCCC\t*\n' "$CONTROL"
                fi
                ;;
        esac
        return 0
    fi
    return 2
}

truth="$WORK/truth.vcf.gz"; : > "$truth"
bam="$WORK/tumor.bam";      : > "$bam"
ref="$WORK/ref.fa";         : > "$ref"

echo "=== probe construction (known answer) ==="
: > "$WORK/faidx_regions"
build_del_probes "$truth" "$ref" "$WORK/j.tsv" "$WORK/c.tsv"
is "junction probe is left flank + post-END flank" "$(printf 'chr1\t1000\t1000\t%s' "$JUNCTION")" "$(cat "$WORK/j.tsv")"
is "control probe is left flank + first deleted bases" "$(printf 'chr1\t1000\t1000\t%s' "$CONTROL")" "$(cat "$WORK/c.tsv")"
is "flank coordinates" "chr1:986-1000 chr1:2001-2015 chr1:1001-1015" "$(cat "$WORK/faidx_regions")"

echo "=== literal DEL derives END from REF length ==="
: > "$WORK/faidx_regions"
BCF_ROWS="$(printf 'chr1\t1000\t%s\tA\t.' "$(printf 'A%.0s' {1..1001})")" \
    build_del_probes "$truth" "$ref" "$WORK/j2.tsv" "$WORK/c2.tsv"
is "literal DEL yields the same span as the symbolic one" "chr1:986-1000 chr1:2001-2015 chr1:1001-1015" "$(cat "$WORK/faidx_regions")"

echo "=== realized deletion ==="
DEL_UNREALIZED=0
VIEW_MODE=junction verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "realized deletion evaluates successfully" 0 "$rc"
has "realized deletion reports junction evidence" "$out" "1 realized"
is "realized deletion is not counted unrealized" 0 "$DEL_UNREALIZED"

echo "=== reads read THROUGH the deletion (must fire) ==="
DEL_UNREALIZED=0
VIEW_MODE=through verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
has "unrealized deletion is reported explicitly" "$out" "read THROUGH the deletion"
is "unrealized deletion is counted" 1 "$DEL_UNREALIZED"
has "unrealized deletion names the issue" "$out" "#590"

echo "=== junction only at END, nothing at POS (must NOT fire) ==="
# Regression for job 21532276. Past ~a read length there is no read spanning the deletion;
# bwa splits, and the record carrying the junction can anchor on either side. Searching only
# the anchor reported three genuinely-realized deletions as unrealized.
DEL_UNREALIZED=0
VIEW_MODE=far_only verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "junction found at the far breakpoint" 0 "$rc"
is "far-anchored junction is NOT called unrealized" 0 "$DEL_UNREALIZED"
has "far-anchored junction counts as realized" "$out" "1 realized"

echo "=== locus with no spanning reads (must NOT fire) ==="
DEL_UNREALIZED=0
VIEW_MODE=uncovered verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
out="$(<"$WORK/out")"
has "uncovered locus is reported as uninterpretable" "$out" "NO SPANNING READS"
is "uncovered locus is NOT counted as unrealized" 0 "$DEL_UNREALIZED"
hasnt "uncovered locus does not claim reads read through" "$out" "read THROUGH"

echo "=== unprobeable records are counted, not dropped silently ==="
DEL_SKIPPED=0
BCF_ROWS="$(printf 'chr1\t1000\tA\t<DEL>\t1005\nchr1\t8\tA\t<DEL>\t900')" \
    build_del_probes "$truth" "$ref" "$WORK/j3.tsv" "$WORK/c3.tsv"
is "too-short DEL and contig-start DEL are both skipped" 2 "$DEL_SKIPPED"
is "neither produced a probe" 0 "$(wc -l < "$WORK/j3.tsv")"

echo "=== tool failures stay distinguishable from no support ==="
DEL_UNREALIZED=0
VIEW_MODE=fail verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "samtools view failure is non-zero" 1 "$rc"
has "samtools view failure is not reported as no support" "$out" "verification was not evaluated"

DEL_UNREALIZED=0
FAIDX_MODE=fail verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "faidx failure is non-zero" 1 "$rc"

DEL_UNREALIZED=0
BCF_ROWS="" verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "no probeable DEL is not an error" 0 "$rc"
has "no probeable DEL says so" "$out" "nothing to verify"

echo "=== junction/control rows describing different loci (must not be compared) ==="
# The two hit lists come from two independent matcher runs. Nothing but this guard makes them
# describe the same loci, and a silent misalignment would compare unrelated positions and
# report confident nonsense. Force it: return the control rows in reverse order.
count_probe_hits() {
    if [[ "$1" == *control* ]]; then
        tac "$1" | awk -F'\t' 'BEGIN { OFS = "\t" } { print $1, $2, $3, 1 }'
    else
        awk -F'\t' 'BEGIN { OFS = "\t" } { print $1, $2, $3, 1 }' "$1"
    fi
}
DEL_UNREALIZED=0
BCF_ROWS="$(printf 'chr1\t1000\tA\t<DEL>\t2000\nchr1\t3000\tA\t<DEL>\t4000')" \
    verify_planted_del "$truth" "$bam" "$ref" "$WORK" > "$WORK/out" 2>&1
rc=$?; out="$(<"$WORK/out")"
is "misaligned probe rows are a hard failure" 1 "$rc"
has "misalignment is named, not reported as no support" "$out" "describe different loci"
eval "$count_probe_hits_src"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
