#!/usr/bin/env bash
# Tests for benchmark.sbatch's configuration knobs and the regression suite's use of them.
#
# WHY THIS EXISTS. `REPS=3` was a plain assignment that ignored the environment, while the
# script's own usage header advertised `REPS=2 sbatch ...` as a supported override. Nothing
# failed; the override was simply discarded. That is the worst kind of harness bug — a knob
# that looks configurable sends whoever needs it hunting for the cause somewhere else, and
# the regression suite had no way to trim a perf arm that was timing out (job 21766280).
#
# The config block is EVALUATED OUT OF THE REAL FILE rather than copied here, so this tests
# the same lines the job runs. A copy would keep passing after the original drifted.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="${BENCH:-$HERE/../benchmark.sbatch}"
SUITE="${SUITE:-$HERE/../regression_suite.sh}"

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

# Pull the GENOMES array assignment and the REPS line straight out of the sbatch, and run
# them with nothing else. `sed -n '/^GENOMES=(/,/^})/p'` captures the multi-line array.
genomes_block() { sed -n '/^GENOMES=(/,/^})/p' "$BENCH"; }
reps_line()     { grep -m1 '^REPS=' "$BENCH"; }

probe_genomes() { # env passthrough; prints "count|elem|elem|..."
    bash -c "DATA_DIR=/data
$(genomes_block)
printf '%s' \"\${#GENOMES[@]}\"; for g in \"\${GENOMES[@]}\"; do printf '|%s' \"\$g\"; done"
}
# The printf goes on its OWN LINE: the extracted REPS= line ends in a trailing `#` comment,
# which on one line would comment out everything after it — including the probe. That cost a
# round of "the override still does not work" before the test itself turned out to be wrong.
probe_reps() { bash -c "$(reps_line)
printf '%s' \"\$REPS\""; }

echo "=== the config block is extractable, so the rest of this file is not vacuous ==="
[[ -n "$(genomes_block)" ]] && ok "GENOMES array block found in benchmark.sbatch" \
  || bad "GENOMES array block found in benchmark.sbatch" "a multi-line GENOMES=( ... )" "nothing matched"
[[ -n "$(reps_line)" ]] && ok "REPS assignment found in benchmark.sbatch" \
  || bad "REPS assignment found in benchmark.sbatch" "a REPS= line" "nothing matched"

echo "=== REPS honours the environment (it did not, and was documented as if it did) ==="
eq "REPS defaults to 3"            "$(probe_reps)"          "3"
eq "REPS=1 from the environment wins"  "$(REPS=1 probe_reps)"   "1"
eq "REPS=2, the value the usage header advertises, wins" "$(REPS=2 probe_reps)" "2"

echo "=== GENOMES defaults to the full head-to-head set ==="
out="$(probe_genomes)"
eq "four genomes by default"       "${out%%|*}"             "4"
# Quote handling: the defaults are written with double quotes INSIDE a ${VAR:-...}
# expansion, which is exactly where a literal quote can survive into the value and produce
# a label like `"ecoli`. Assert the elements are clean.
has "the default ecoli entry is unquoted and path-expanded"  "$out" "|ecoli:/data/ecoli.fa"
has "the default chr22 entry is unquoted and path-expanded"  "$out" "|chr22:/data/chr22.fa"
has "the default human entry is unquoted and path-expanded"  "$out" "|human:/data/GRCh38.fa"
case "$out" in *'"'*) bad "no literal double quote survives into a genome entry" "no quotes" "$out";;
             *) ok "no literal double quote survives into a genome entry";; esac

echo "=== GENOMES honours the environment ==="
out2="$(GENOMES="ecoli:/data/ecoli.fa chr22:/data/chr22.fa" probe_genomes)"
eq "an override of two genomes yields exactly two"  "${out2%%|*}"  "2"
has "the overridden set keeps ecoli"                "$out2" "|ecoli:/data/ecoli.fa"
has "the overridden set keeps chr22"                "$out2" "|chr22:/data/chr22.fa"
# Must-not-fire: the expensive arms are GONE, not merely reordered. This is the whole point
# of the override — GRCh38 is ~3 h of the timed-out job and reaches no gate.
case "$out2" in *human*) bad "an override DROPS the genomes it omits" "no human/GRCh38" "$out2";;
              *) ok "an override DROPS the genomes it omits";; esac
case "$out2" in *yeast*) bad "an override drops yeast too" "no yeast" "$out2";;
              *) ok "an override drops yeast too";; esac
eq "a single-genome override yields one" \
   "$(GENOMES="ecoli:/data/ecoli.fa" probe_genomes | cut -d'|' -f1)" "1"

echo "=== the regression suite trims the perf arm to what its own gate reads back ==="
# do_collect extracts only ecoli/chr22, tool==eidolon, threads==1. Anything else the job
# computes is discarded, so the suite must not ask for it.
perf_line="$(sed -n '/sub perf/,/^$/p' "$SUITE")"
has "the perf arm restricts GENOMES"            "$perf_line" "GENOMES="
has "the perf arm asks for ecoli"               "$perf_line" "ecoli"
has "the perf arm asks for chr22"               "$perf_line" "chr22"
has "the perf arm drops the NEAT 4 comparison"  "$perf_line" "NEAT_MAX_GENOME_MB=0"
has "the perf arm collapses the thread sweep"   "$perf_line" "SCALING_THREAD_MODES=1"
case "$perf_line" in *GRCh38*|*human*) bad "the perf arm does not request GRCh38" "no human arm" "$perf_line";;
                   *) ok "the perf arm does not request GRCh38";; esac
# The gate reads exactly these two genomes; if do_collect ever changes, this pairing breaks.
collect_line="$(grep -m1 'for gg in' "$SUITE")"
eq "the gate still collects exactly ecoli and chr22" \
   "$(printf '%s' "$collect_line" | tr -s ' ' | sed 's/^ *//')" "for gg in ecoli chr22; do"

echo "=== both scripts still parse ==="
bash -n "$BENCH"  && ok "benchmark.sbatch parses"   || bad "benchmark.sbatch parses" "exit 0" "syntax error"
bash -n "$SUITE"  && ok "regression_suite.sh parses" || bad "regression_suite.sh parses" "exit 0" "syntax error"

# Floor on how many assertions must execute. Raise it when adding tests; if it ever reads
# low, an assertion stopped running rather than started failing.
MIN_ASSERTIONS=25
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
