#!/usr/bin/env bash
# The harness must be able to say why it failed (#583).
#
# Job 21532276 exited 1 after reporting "3 deletion(s) not realized" and NOTHING explaining
# which three ever reached a file. `exec` with redirections and no command applies them to the
# current shell permanently, so `exec 8<>lock 2>/dev/null` in index_reference_locked pointed
# fd 2 at /dev/null for the rest of the job — from the bwa step onward. The three loci were
# recoverable only because the probe files happened to survive on disk.
#
# This suite runs the real function and asserts stderr still works afterwards.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB="${LIB:-$HERE/../lib_report.sh}"
PIPELINE="${PIPELINE:-$HERE/../sv_pipeline.sbatch}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label file from to; do
        [[ -n "$label" ]] || continue
        src="$HERE/../$file"
        cp "$src" "$WORK/orig.$file"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$src"
        if cmp -s "$src" "$WORK/orig.$file"; then
            printf '  ERROR   %-50s mutation did not apply\n' "$label"; survived=$((survived+1))
            cp "$WORK/orig.$file" "$src"; continue
        fi
        if bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-49s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
        cp "$WORK/orig.$file" "$src"
    done <<'MUTATIONS'
exec redirect leaks to the whole shell@lib_report.sh@if { exec 8<>"${ref}.index.lock"; } 2>/dev/null; then@if exec 8<>"${ref}.index.lock" 2>/dev/null; then
pipeline drops its stderr merge@sv_pipeline.sbatch@exec 2>&1@:
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

PASS=0; FAIL=0
ok() { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
is() { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$2" "$3"; }

echo "=== stderr survives index_reference_locked ==="
# Run the REAL function in a subshell whose stderr we capture, with bwa-mem2 stubbed out and
# an index already present so nothing is actually built.
(
    extract() { awk "/^$1\(\)/,/^}\$/" "$LIB"; }
    eval "$(extract index_reference_locked)"
    bwa-mem2() { :; }
    ref="$WORK/ref.fa"; : > "$ref"; printf 'x' > "$ref.bwt.2bit.64"
    echo "MARKER-BEFORE" >&2
    index_reference_locked "$ref" >/dev/null
    echo "MARKER-AFTER" >&2
) 2>"$WORK/err" >/dev/null

# Presence, not exact equality: the function shells out to samtools, whose own stderr
# legitimately lands here too. MARKER-AFTER is the load-bearing one — it is the line the
# permanent redirect swallowed.
grep -q '^MARKER-BEFORE$' "$WORK/err" \
    && ok "stderr works before the locked index" \
    || bad "stderr works before the locked index" "MARKER-BEFORE present" "$(cat "$WORK/err")"
grep -q '^MARKER-AFTER$' "$WORK/err" \
    && ok "stderr SURVIVES the locked index" \
    || bad "stderr SURVIVES the locked index" "MARKER-AFTER present" "$(cat "$WORK/err")"

echo "=== the pipeline merges stderr into stdout ==="
grep -q '^exec 2>&1$' "$PIPELINE" && ok "pipeline sets exec 2>&1" || bad "pipeline sets exec 2>&1" "present" "absent"
grep -q '\[selftest\] stderr reaches this log' "$PIPELINE" \
    && ok "pipeline emits a stderr selftest marker" \
    || bad "pipeline emits a stderr selftest marker" "present" "absent"

echo "=== no bare 'exec <fd> ... <redirect>' anywhere in scripts/delta ==="
# The footgun itself: `exec` with no command makes EVERY redirection on the line permanent.
offenders="$(grep -rnE "^[[:space:]]*(if[[:space:]]+)?!?[[:space:]]*exec[[:space:]]+[0-9]+[<>][^|;]*[[:space:]]+[0-9]*[<>]" \
    "$HERE/.." --include=*.sh --include=*.sbatch 2>/dev/null || true)"
is "no unguarded exec redirect" "" "$offenders"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
