#!/usr/bin/env bash
# archive_run must keep the job's own log, wherever SLURM wrote it.
#
# The log is where every Delta script prints its configuration, so a result archived without
# it cannot be reproduced or compared later. archive_run used to guess the log's name as
# <jobname>_<jobid>.out in the current directory. That silently missed every script whose
# #SBATCH --output names it differently: fit_hg002_defaults (hg002defaults-%j.log),
# fit_overlap_matrix (overlapfit-%j.log) and stage_raw_pairs (rawpairs-%j.log). It now asks
# scontrol for StdOut/StdErr, keeps the guess as a fallback, and warns when it finds nothing.
#
# lib_report.sh is sourced as-is; scontrol and sacct are stubbed as shell functions.
#
# Usage:
#   scripts/delta/tests/test_archive_run.sh            # run the suite
#   scripts/delta/tests/test_archive_run.sh --mutate   # prove the suite is non-vacuous
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB="${LIB:-$HERE/../lib_report.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$LIB" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        if cmp -s "$LIB" "$WORK/mutant.sh"; then
            printf '  ERROR   %-48s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if LIB="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-47s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
scontrol's paths are ignored@    logs="$(_job_log_paths)" || logs=""@    logs=""
the fallback guess is dropped@    [[ -n "$logs" ]] || logs="${base}.out"$'\n'"${base}.err"@    :
a missing log is skipped silently@    if [[ "$copied" -eq 0 && -n "${SLURM_JOB_ID:-}" ]]; then@    if false; then
only StdOut is read@    for key in StdOut StdErr; do@    for key in StdOut; do
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
is()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$2" "$3"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
hasnt() { case "$2" in *"$3"*) bad "$1" "does NOT contain: $3" "$2";; *) ok "$1";; esac; }

export SCRATCH="$WORK/scratch" RESULTS_DIR="$WORK/results"
mkdir -p "$SCRATCH"
# shellcheck disable=SC1090
source "$LIB"
sacct() { return 1; }

# One archive in a fresh directory. $1 job id, $2 job name; scontrol's reply comes from
# $SCONTROL_REPLY, and an empty reply makes it fail, as on a node where it cannot answer.
run_archive() {
    local jid="$1" name="$2"
    ( cd "$WORK/submit" && SLURM_JOB_ID="$jid" SLURM_JOB_NAME="$name" \
        archive_run testkind "$WORK/out" 2>&1 )
}
scontrol() {
    [[ -n "${SCONTROL_REPLY:-}" ]] || return 1
    printf '%s\n' "$SCONTROL_REPLY"
}
mkdir -p "$WORK/submit" "$WORK/out" "$WORK/logs"

echo "=== a custom-named log is found through scontrol ==="
echo "configuration: fragment Normal(400, 90)" > "$WORK/logs/overlapfit-101.log"
SCONTROL_REPLY="JobId=101 JobName=eidolon-overlapfit StdErr=$WORK/logs/overlapfit-101.log StdIn=/dev/null StdOut=$WORK/logs/overlapfit-101.log Power=" \
    out="$(run_archive 101 eidolon-overlapfit)"
is  "the log is archived" "configuration: fragment Normal(400, 90)" \
    "$(cat "$RESULTS_DIR/testkind/job_101/overlapfit-101.log" 2>/dev/null)"
hasnt "no warning when it is found" "$out" "WARNING"

echo "=== separate StdOut and StdErr are both archived ==="
echo "out text" > "$WORK/logs/j-102.out"; echo "err text" > "$WORK/logs/j-102.err"
SCONTROL_REPLY="JobId=102 StdErr=$WORK/logs/j-102.err StdIn=/dev/null StdOut=$WORK/logs/j-102.out" \
    out="$(run_archive 102 j)"
is "stdout archived" "out text" "$(cat "$RESULTS_DIR/testkind/job_102/j-102.out" 2>/dev/null)"
is "stderr archived" "err text" "$(cat "$RESULTS_DIR/testkind/job_102/j-102.err" 2>/dev/null)"

echo "=== scontrol cannot answer: the <jobname>_<jobid> guess still works ==="
echo "fallback text" > "$WORK/submit/eidolon-realism_103.out"
SCONTROL_REPLY="" out="$(run_archive 103 eidolon-realism)"
is "the guessed log is archived" "fallback text" \
   "$(cat "$RESULTS_DIR/testkind/job_103/eidolon-realism_103.out" 2>/dev/null)"
hasnt "no warning when the guess finds it" "$out" "WARNING"

echo "=== nothing found: say so, and still write the manifest ==="
SCONTROL_REPLY="" out="$(run_archive 104 nolog)"
has "the miss is reported"                "$out" "WARNING: job log not found"
has "and names where it looked"           "$out" "nolog_104.out"
is  "the manifest is still written"       "1" \
    "$([[ -s "$RESULTS_DIR/testkind/job_104/run_manifest.tsv" ]] && echo 1 || echo 0)"

echo "=== outside SLURM there is no job log to miss (must-not-fire) ==="
out="$( cd "$WORK/submit" && unset SLURM_JOB_ID && archive_run testkind "$WORK/out" 2>&1 )"
hasnt "no warning for a local run" "$out" "WARNING"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 && "$PASS" -ge 10 ]]
