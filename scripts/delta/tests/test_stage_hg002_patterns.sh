#!/usr/bin/env bash
# Tests for stage_hg002.sh's read-URL construction.
#
# WHY THIS EXISTS. The default patterns were written as `${R1_PATTERN:-D1_S1_{lane}_R1_{id}...}`.
# A `${VAR:-default}` ends at the FIRST unescaped `}`, so the default was cut off at the brace
# in `{lane}`: the variable came out as `D1_S1_{lane` with the literal tail
# `_R1_{id}.fastq.gz}` appended. `{lane}` then never existed to be substituted, `{id}` survived
# and WAS substituted, and a stray `}` landed on the end. The job fetched
# `D1_S1_{lane_R1_001.fastq.gz}` and got a 404 -- which reads exactly like a dataset that moved,
# and cost a submit round trip to tell apart.
#
# Nothing here reaches the network. `PRINT_URLS=1` makes the script print what it would fetch
# and exit, which is the seam these tests use.
#
# KNOWN ANSWER. The expected filenames are not this script's opinion: the GIAB directory listing
# carries D1_S1_L001_R{1,2}_001.fastq.gz, and both URLs were confirmed HTTP 200 (5.6 GB and
# 6.0 GB) on 2026-09-19.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT:-$HERE/../stage_hg002.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --mutate: break the script under test and confirm this suite notices. CI runs this for
# every suite, because a suite that cannot fail is decoration.
if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$SCRIPT" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        # An unapplied mutation and a surviving one produce identical output. Compare BYTES,
        # not a pattern: a grep whose own escaping is wrong reports "unchanged" for both.
        if cmp -s "$SCRIPT" "$WORK/mutant.sh"; then
            printf '  ERROR    %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if SCRIPT="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-52s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTS'
the original bug: brace-truncated ${VAR:-default}@[[ -n "${R1_PATTERN:-}" ]] || R1_PATTERN='D1_S1_{lane}_R1_{id}.fastq.gz'@R1_PATTERN="${R1_PATTERN:-D1_S1_{lane}_R1_{id}.fastq.gz}"
{lane} is never substituted@p="${p//\{lane\}/$2}"@p="$p"
{id} is never substituted@p="${p//\{id\}/$3}"@p="$p"
the unsubstituted-placeholder guard is disabled@*[{}]*)@*ZZNOMATCHZZ)
MUTS
    echo
    [[ "$survived" -eq 0 ]] && echo "all mutations caught" || echo "$survived mutation(s) survived"
    exit $(( survived > 0 ? 1 : 0 ))
fi

fails=0
check() {  # <label> <expected> <actual>
    if [[ "$2" == "$3" ]]; then
        printf '  ok       %s\n' "$1"
    else
        printf '  FAIL     %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"
        fails=$((fails + 1))
    fi
}

urls() { PRINT_URLS=1 SCRATCH="$WORK" bash "$SCRIPT" 2>&1; }

BASE_DEFAULT="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads"

# 1. The default library. This is the case that was broken.
out="$(urls)"
check "default R1 URL" \
    "$BASE_DEFAULT/D1_S1_L001_R1_001.fastq.gz" "$(echo "$out" | sed -n 1p)"
check "default R2 URL" \
    "$BASE_DEFAULT/D1_S1_L001_R2_001.fastq.gz" "$(echo "$out" | sed -n 2p)"
check "default emits exactly two URLs" "2" "$(echo "$out" | grep -c .)"

# 2. No brace survives into any URL. The specific symptom, asserted directly, because a
#    pattern could be mangled in some other way and still produce two plausible lines.
check "no unsubstituted brace in the default URLs" "" "$(echo "$out" | grep -o '[{}]' | tr -d '\n')"

# 3. Several chunks, and a lane that is not L001 -- the substitution that never ran.
out="$(CHUNKS="L001:001 L002:003" urls)"
check "two chunks give four URLs" "4" "$(echo "$out" | grep -c .)"
check "second chunk substitutes both lane and id" \
    "$BASE_DEFAULT/D1_S1_L002_R1_003.fastq.gz" "$(echo "$out" | sed -n 3p)"

# 4. Another library, which is what the patterns are parameters FOR (#730).
out="$(BASE=https://example.org/reads \
       R1_PATTERN='PRE_{lane}_1_{id}.fq.gz' R2_PATTERN='PRE_{lane}_2_{id}.fq.gz' \
       CHUNKS="L007:042" urls)"
check "custom pattern R1" "https://example.org/reads/PRE_L007_1_042.fq.gz" "$(echo "$out" | sed -n 1p)"
check "custom pattern R2" "https://example.org/reads/PRE_L007_2_042.fq.gz" "$(echo "$out" | sed -n 2p)"

# 5. MUST FIRE: a placeholder that is not {lane} or {id} is refused, rather than becoming a
#    URL the server will 404 on. This is the guard that would have named the original bug.
out="$(R1_PATTERN='X_{lane}_{nope}.gz' urls)"; rc=$?
check "an unknown placeholder is refused (exit)" "1" "$rc"
case "$out" in
    *"unsubstituted placeholder"*) printf '  ok       refusal names the problem\n' ;;
    *) printf '  FAIL     refusal message missing; got: %s\n' "$out"; fails=$((fails + 1)) ;;
esac

# 6. MUST NOT FIRE: a pattern with no placeholders at all is legal -- a single-file dataset.
out="$(BASE=https://example.org R1_PATTERN='only_1.fq.gz' R2_PATTERN='only_2.fq.gz' urls)"
check "a pattern without placeholders is accepted" \
    "https://example.org/only_1.fq.gz" "$(echo "$out" | sed -n 1p)"

echo
if [[ "$fails" -eq 0 ]]; then
    echo "test_stage_hg002_patterns: all checks passed"
else
    echo "test_stage_hg002_patterns: $fails check(s) FAILED"
fi
exit $(( fails > 0 ? 1 : 0 ))
