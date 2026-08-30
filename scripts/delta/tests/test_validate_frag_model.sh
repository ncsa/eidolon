#!/usr/bin/env bash
# Tests for validate_frag_model.sh — the check that a built fragment model reproduces the
# BAM it came from.
#
# The tool lives in scripts/delta/, not here: it needs a real BAM and real model files, so
# it cannot run argument-free the way this directory's suites must. This is its test.
#
# samtools is STUBBED, the way test_scorer_calibration.sh stubs truvari. The runner installs
# only bcftools and tabix, and what needs testing is the tool's own logic — the RNEXT filter,
# the histogram, the restriction to a model's support, the verdict — not whether samtools
# can read a BAM. A stub also lets the fixture state its own known answer.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="${TOOL:-$HERE/../validate_frag_model.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

command -v jq >/dev/null 2>&1 || { echo "SKIP: jq unavailable"; exit 0; }

# ── fixture: a right-skewed library plus a handful of megabase-scale discordant pairs,
#    which is the shape that made the first real run read 99.92% off ───────────────
mk_sam() {  # writes SAM records to stdout
    awk 'BEGIN {
        # library: dense mode 300-420, thin tail to 700
        for (l = 300; l <= 420; l++) for (i = 0; i < 30; i++) print l
        for (l = 421; l <= 700; l++) { c = 30 - int((l - 421) / 10); if (c < 1) c = 1
                                       for (i = 0; i < c; i++) print l }
        # discordant: mates megabases apart, same contig. Real BAMs carry these; the
        # builder trims them, so the check must too.
        for (i = 0; i < 12; i++) print 3000000 + i * 250000
    }' | awk '{ printf "r%d\t67\tchr1\t%d\t60\t100M\t=\t%d\t%d\t%s\t%s\n",
                NR, NR * 7 + 1, NR * 7 + $1, $1, "A", "I" }'
    # INTER-CHROMOSOMAL pairs, with a TLEN that lands INSIDE the library range. The builder
    # skips these (`ref_id != mate_ref_id`) and so must the check. Putting them at 650 bp
    # rather than at some absurd length is deliberate: an out-of-range value would be
    # discarded by the support restriction anyway, so the RNEXT filter could be deleted and
    # nothing would notice -- which is exactly what happened the first time this fixture
    # was written.
    awk 'BEGIN { for (i = 0; i < 2000; i++)
        printf "x%d\t67\tchr1\t%d\t60\t100M\tchr2\t%d\t650\tA\tI\n", i, i * 11 + 1, i * 11 + 650 }'
}
mk_sam > "$WORK/reads.sam"
LIB_N=$(awk '$9 < 1000000' "$WORK/reads.sam" | wc -l)
DISC_N=$(awk '$9 >= 1000000' "$WORK/reads.sam" | wc -l)

# stub samtools: the tool only ever calls `samtools view <flags> <bam>`
mkdir -p "$WORK/bin"
cat > "$WORK/bin/samtools" <<STUB
#!/usr/bin/env bash
cat "$WORK/reads.sam"
STUB
chmod +x "$WORK/bin/samtools"
export PATH="$WORK/bin:$PATH"
: > "$WORK/fake.bam"

# ── model fixtures: cumulative weights, the way the real files store them ────────
mk_discrete() {  # <out.json.gz> <lo> <hi> [shift]
    local out="$1" lo="$2" hi="$3" shift_by="${4:-0}"
    awk -v lo="$lo" -v hi="$hi" -v sh="$shift_by" '
        BEGIN {
            n = 0
            for (l = lo; l <= hi; l++) {
                c = (l <= 420) ? 30 : 30 - int((l - 421) / 10); if (c < 1) c = 1
                v[++n] = l + sh; w[n] = c; tot += c
            }
            printf "{\"Discrete\":{\"distribution\":{\"values\":["
            for (i = 1; i <= n; i++) printf "%s%d", (i > 1 ? "," : ""), v[i]
            printf "],\"weights\":["
            acc = 0
            for (i = 1; i <= n; i++) { acc += w[i] / tot; printf "%s%.10f", (i > 1 ? "," : ""), acc }
            printf "]}}}"
        }' | gzip > "$out"
}
mk_discrete "$WORK/good.json.gz" 300 700
mk_discrete "$WORK/wrong.json.gz" 300 700 400      # same shape, shifted 400 bp
printf '{"Normal":{"mean":404.0,"st_dev":69.0}}' | gzip > "$WORK/normal.json.gz"

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$TOOL" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        if cmp -s "$TOOL" "$WORK/mutant.sh"; then
            printf '  ERROR   %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if TOOL="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-51s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
a mismatched model is accepted@    awk -v v="$dm" 'BEGIN{exit !(v>2.0)}'  && { echo "      FAIL mean off by ${dm}% (>2%)";  FAIL=1; }@    if false; then echo; fi
the truth is never restricted to the model support@    awk -v lo="$RLO" -v hi="$RHI" '$1>=lo && $1<=hi' "$WORK/truth.tsv" > "$WORK/truth_r.tsv"@    cp "$WORK/truth.tsv" "$WORK/truth_r.tsv"
discordant pairs are counted as library fragments@  | awk '$7=="=" && $9>0 { c[$9]++ } END { for (l in c) print l, c[l] }' \@  | awk '$9>0 { c[$9]++ } END { for (l in c) print l, c[l] }' \
a zero-pair measurement is not fatal@[[ "${NPAIRS:-0}" -gt 0 ]] || { echo "FATAL: zero pairs passed the filter -- nothing to compare against" >&2; exit 1; }@true
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     %s\n' "$1" "${2:-}"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "expected to contain: $3";; esac; }

echo "=== a model built from this data is accepted ==="
out="$(bash "$TOOL" "$WORK/fake.bam" "$WORK/good.json.gz" 2>&1)"; rc=$?
[[ "$rc" -eq 0 ]] && ok "the matching model passes" || bad "the matching model passes" "exit $rc: $out"
has "it reproduces the distribution" "$out" "reproduces the distribution it was built from"

echo "=== the discordant pairs are excluded, and SAID to be ==="
# Known answer: the fixture plants exactly $DISC_N discordant pairs out of $((LIB_N+DISC_N)).
has "the raw row is flagged as outlier-dominated" "$out" "discordant pairs dominate"
has "the excluded mass is reported, not assumed" "$out" "excludes"
has "the restricted row is labelled"              "$out" "REAL BAM in Discrete range"

echo "=== a model from DIFFERENT data is rejected (must-fire) ==="
out="$(bash "$TOOL" "$WORK/fake.bam" "$WORK/wrong.json.gz" 2>&1)"; rc=$?
[[ "$rc" -ne 0 ]] && ok "a shifted model fails" || bad "a shifted model fails" "it exited 0"
has "and says the model does not match"  "$out" "does not match its input"
has "and names the mean as the offender" "$out" "FAIL mean"

echo "=== a Normal is reported but NOT asserted against ==="
# It is EXPECTED to miss the skew. Asserting on it would only re-measure what is known,
# and would make the tool fail on a correct run.
out="$(bash "$TOOL" "$WORK/fake.bam" "$WORK/normal.json.gz" 2>&1)"; rc=$?
[[ "$rc" -eq 0 ]] && ok "a Normal model does not fail the check" || bad "a Normal model does not fail the check" "exit $rc"
has "its skew error is still reported" "$out" "skew off by"

echo "=== a BAM with no passing pairs is fatal, not an empty table ==="
cat > "$WORK/bin/samtools" <<'STUB'
#!/usr/bin/env bash
exit 0
STUB
chmod +x "$WORK/bin/samtools"
out="$(bash "$TOOL" "$WORK/fake.bam" "$WORK/good.json.gz" 2>&1)"; rc=$?
[[ "$rc" -ne 0 ]] && ok "zero pairs is fatal" || bad "zero pairs is fatal" "it exited 0"
has "and says why" "$out" "nothing to compare against"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
