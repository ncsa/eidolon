#!/usr/bin/env bash
# Tests that regression_suite.sh's collector can actually read its harnesses' output.
#
# WHY THIS EXISTS. `d0cd060` renamed the reproducibility harness's verdict labels
# (`multithread_vs_1thread` -> `thread_invariant`, `contig_order_independent` ->
# `contig_name_invariant`) and did not update the collector that greps for them. Neither
# file was wrong on its own. Nothing asserted they had to agree, so:
#
#   - `determinism_thread_invariant` kept its local `ti=FAIL` default and reported FAIL on
#     every run without ever reading a measurement. Its baseline is PASS, so every gate
#     verdict after that commit showed a false regression.
#   - `contig_order_independent` collected an EMPTY value, which an `exact` gate
#     string-compares to a FAIL, making unmeasured indistinguishable from measured-bad.
#
# That is the two-component invariant CLAUDE.md requires be pinned by a test rather than
# left to two files happening to hold the same literal. This test is that pin.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE="${SUITE:-$HERE/../regression_suite.sh}"
ORDER="${ORDER:-$HERE/../run_order_independence.sbatch}"
GATE="${GATE:-$HERE/../regression_gate.sh}"
BASELINE="${BASELINE:-$HERE/../baseline_metrics.tsv}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

# Labels the harness prints, from its own printf statements.
harness_labels() { grep -oE '"[a-z_0-9]+:"' "$ORDER" | tr -d '":' | sort -u; }
# Labels the collector greps for, from its `get <label>` calls in the order) arm.
collector_labels() { sed -n '/^          order)/,/;;/p' "$SUITE" | grep -oE '\$\(get [a-z_0-9]+\)' \
                       | sed 's/\$(get //; s/)//' | sort -u; }

echo "=== the collector's labels must all exist in the harness output ==="
hl="$(harness_labels)"; cl="$(collector_labels)"
[[ -n "$hl" ]] && ok "harness labels were extracted" \
  || bad "harness labels were extracted" "a label list" "nothing matched"
[[ -n "$cl" ]] && ok "collector labels were extracted" \
  || bad "collector labels were extracted" "a label list" "nothing matched"
# THE assertion. Every label the collector asks for must be one the harness prints.
orphans="$(comm -23 <(printf '%s\n' "$cl") <(printf '%s\n' "$hl"))"
if [[ -z "$orphans" ]]; then
    ok "every collector label is printed by the harness"
else
    bad "every collector label is printed by the harness" "no orphans" "orphaned: $(echo $orphans)"
fi
# Non-vacuity: the comparison must be capable of finding an orphan at all.
fake="$(comm -23 <(printf 'a_label_that_cannot_exist\n') <(printf '%s\n' "$hl"))"
eq "the orphan check can detect an orphan" "$fake" "a_label_that_cannot_exist"

echo "=== the specific labels d0cd060 broke ==="
has "the collector reads thread_invariant, the current harness label"    "$cl" "thread_invariant"
has "the collector reads contig_name_invariant, the current label"       "$cl" "contig_name_invariant"
# Outside a COMMENT: the collector deliberately names the retired label when explaining
# why this broke, and a test that forbade the word would push out the explanation.
eq  "the retired multithread_vs_1thread label is gone from executable code" \
    "$(grep -n 'multithread_vs_1thread' "$SUITE" | grep -vE ':[[:space:]]*#' | wc -l | tr -d ' ')" "0"
# Non-vacuity: that filter must still be able to see the word somewhere, or it proves
# nothing about where the word is.
[[ "$(grep -c 'multithread_vs_1thread' "$SUITE")" -gt 0 ]] \
  && ok "the retired label is still explained in a comment" \
  || bad "the retired label is still explained in a comment" "a comment mentioning it" "absent"
eq  "and the harness never printed it either" \
    "$(grep -c 'multithread_vs_1thread' "$ORDER")" "0"

echo "=== every metric the collector emits must have a baseline row, and vice versa ==="
emitted() { sed -n '/^          order)/,/;;/p' "$SUITE" | grep -oE "printf '%s.t%s.n' [a-z_0-9]+" \
              | awk '{print $NF}' | sort -u; }
for m in $(emitted); do
    if grep -qE "^${m}\b" "$BASELINE"; then ok "emitted metric $m has a baseline row"
    else bad "emitted metric $m has a baseline row" "a row in baseline_metrics.tsv" "absent"; fi
done

echo "=== an absent label yields MISSING, never an empty string ==="
# Extract the real get() and run it against output that lacks the label. An empty value is
# what made the drift invisible: `exact` gates compare strings, so "" reads as FAIL.
get_fn() { sed -n '/^            get(){/,/^            }/p' "$SUITE"; }
[[ -n "$(get_fn)" ]] && ok "the collector's get() helper was extracted" \
  || bad "the collector's get() helper was extracted" "a function body" "nothing matched"
printf '  thread_invariant:            PASS  (1thr abc vs 8thr abc)\n' > "$WORK/out.txt"
probe_get() { bash -c "f='$WORK/out.txt'
$(get_fn)
get '$1'" 2>/dev/null; }
eq "a present label returns its verdict"        "$(probe_get thread_invariant)" "PASS"
eq "an absent label returns MISSING"            "$(probe_get contig_name_invariant)" "MISSING"
[[ -n "$(probe_get contig_name_invariant)" ]] && ok "an absent label is never the empty string" \
  || bad "an absent label is never the empty string" "MISSING" "(empty)"
# And it must say so on stderr rather than failing quietly.
warn="$(bash -c "f='$WORK/out.txt'
$(get_fn)
get nonexistent_label" 2>&1 >/dev/null)"
has "an absent label warns on stderr" "$warn" "drifted"

echo "=== the gate reports unmeasured separately from measured-bad ==="
# A FAIL that was never measured must not read like a regression. Baseline row is an
# `exact` gate, which is the case that had no missing-data guard at all.
printf 'determinism_thread_invariant\tMISSING\n' > "$WORK/cand.tsv"
gout="$(bash "$GATE" "$BASELINE" "$WORK/cand.tsv" 2>&1 || true)"
has "a MISSING candidate is called out as NOT MEASURED" "$gout" "NOT MEASURED"
has "and names label drift as the likely cause"          "$gout" "label drift"
# Must-not-fire: a genuinely measured FAIL must still read as an ordinary failure.
printf 'determinism_thread_invariant\tFAIL\n' > "$WORK/cand2.tsv"
gout2="$(bash "$GATE" "$BASELINE" "$WORK/cand2.tsv" 2>&1 || true)"
case "$gout2" in *"NOT MEASURED"*) bad "a measured FAIL is not mislabelled as unmeasured" "no NOT MEASURED" "$gout2";;
                *) ok "a measured FAIL is not mislabelled as unmeasured";; esac
has "a measured FAIL still fails" "$gout2" "FAIL"
# And a measured PASS must pass, or the guard above would be hiding real regressions.
printf 'determinism_thread_invariant\tPASS\n' > "$WORK/cand3.tsv"
gout3="$(bash "$GATE" "$BASELINE" "$WORK/cand3.tsv" 2>&1 || true)"
has "a measured PASS still passes" "$gout3" "PASS"

# Floor on how many assertions must execute. Raise it when adding tests; if it ever reads
# low, an assertion stopped running rather than started failing.
MIN_ASSERTIONS=22
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
