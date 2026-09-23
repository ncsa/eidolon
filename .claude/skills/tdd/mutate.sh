#!/usr/bin/env bash
# Mutation-test helper: prove a test is non-vacuous by breaking the code it covers.
#
# WHY THIS EXISTS: the mutation step lies silently. A `sed` or `str.replace` whose pattern
# does not match changes nothing and reports nothing, so an UNAPPLIED edit and a SURVIVING
# mutant produce byte-identical output. This script refuses to run a test until it has
# confirmed the file actually changed, and restores the file afterwards so a mutation can
# never reach a commit.
#
#   ./mutate.sh start <file>              snapshot the pristine file
#   ... edit <file> with your normal editing tool ...
#   ./mutate.sh run   <file> [cargo test args...]   verify changed -> test -> restore
#   ./mutate.sh diff  <file>              show the mutation without running anything
#   ./mutate.sh restore <file>            abandon the mutation
#
# Exit status from `run`: 0 = KILLED (the test caught it, the test is real)
#                         1 = SURVIVOR (the test did not catch it — fix the test)
#                         2 = refused: no snapshot, file unchanged, or the filter ran no tests

set -euo pipefail

SNAP_DIR="${TMPDIR:-/tmp}/eidolon-mutate-$(id -u)"

usage() { sed -n '3,20p' "$0" | sed 's/^# \?//'; exit 2; }

snap_path() {
  # Flatten the absolute path into a single snapshot filename.
  local abs
  abs=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
  printf '%s/%s' "$SNAP_DIR" "$(printf '%s' "$abs" | tr '/' '%')"
}

require_file() {
  [ -f "$1" ] || { echo "mutate: no such file: $1" >&2; exit 2; }
}

# Restore content WITHOUT restoring the mtime.
#
# `cp -p` was the original implementation and it is a trap: it sets the file's mtime back to
# when the snapshot was taken, which is EARLIER than the build artifact produced from the
# mutated source. Cargo (and make, and every other mtime-based build system) then decides the
# crate is unchanged and silently reruns the MUTANT. Measured: after a mutation that disabled
# a bounds check, three consecutive `cargo test` runs used the mutated binary and reported a
# genuinely correct test as failing. `touch` forces the rebuild that proves the restore.
restore_file() {
  cp "$1" "$2"
  touch "$2"
}

require_snapshot() {
  local snap=$1
  [ -f "$snap" ] || {
    echo "mutate: no snapshot for this file — run '$0 start <file>' before mutating." >&2
    exit 2
  }
}

cmd=${1:-}; [ -n "$cmd" ] || usage
shift || true

case "$cmd" in
  start)
    file=${1:-}; [ -n "$file" ] || usage
    require_file "$file"
    mkdir -p "$SNAP_DIR"
    snap=$(snap_path "$file")
    cp -p "$file" "$snap"
    echo "mutate: snapshot taken for $file"
    echo "mutate: now apply ONE mutation — flip a comparison, drop a '+ 1', return the"
    echo "        other branch — then: $0 run $file <test_name>"
    ;;

  diff)
    file=${1:-}; [ -n "$file" ] || usage
    require_file "$file"
    snap=$(snap_path "$file"); require_snapshot "$snap"
    diff -u "$snap" "$file" || true
    ;;

  restore)
    file=${1:-}; [ -n "$file" ] || usage
    snap=$(snap_path "$file"); require_snapshot "$snap"
    restore_file "$snap" "$file"
    rm -f "$snap"
    echo "mutate: $file restored"
    ;;

  run)
    file=${1:-}; [ -n "$file" ] || usage
    shift
    require_file "$file"
    snap=$(snap_path "$file"); require_snapshot "$snap"

    # THE GATE. Without this the whole exercise is theater.
    if cmp -s "$snap" "$file"; then
      echo "mutate: FILE IS UNCHANGED — the mutation was never applied." >&2
      echo "mutate: this is the failure mode this script exists to catch. Nothing was run." >&2
      rm -f "$snap"
      exit 2
    fi

    # Restore on ANY exit path, including a Ctrl-C mid-test or a crash in the reporting
    # below. A mutation left in the tree is how a deliberate break becomes an accidental
    # commit, so the trap goes in BEFORE anything that can fail.
    trap 'restore_file "$snap" "$file"; rm -f "$snap"; echo; echo "mutate: $file restored"' EXIT

    # `diff` exits 1 when the files differ, which is the ONLY case that reaches this line.
    # Unguarded under `set -e` + `pipefail` that aborts the script on success — and the
    # resulting exit 1 is indistinguishable from a SURVIVOR verdict. Guard both pipelines.
    echo "mutate: mutation applied:"
    { diff -u "$snap" "$file" || true; } | sed -n '/^@@/,$p' | sed 's/^/    /'
    echo

    out=$(mktemp)
    set +e
    cargo test --workspace --no-fail-fast "$@" 2>&1 | tee "$out"
    rc=${PIPESTATUS[0]}
    set -e

    # A filter that matches NOTHING exits 0, which scores as SURVIVOR — the loudest possible
    # wrong answer from this script. Measured: `mutate.sh run <file> quality_tail_collapse`
    # matched no test NAME (it is a file name; cargo needs `--test` for that), ran zero tests,
    # and reported "nothing caught this" about a mutation three tests would have killed.
    #
    # Sum what actually ran. `test result:` lines report passed/failed per target, so zero of
    # both across every target means the filter selected nothing.
    ran=$(awk '/^test result:/ { p += $4; f += $6 } END { print p + f + 0 }' "$out")
    rm -f "$out"
    if [ "${ran:-0}" -eq 0 ]; then
      echo
      echo "mutate: NO TESTS RAN — the filter matched nothing, so this says nothing." >&2
      echo "mutate: a filter selects test NAMES. For a whole integration-test file use" >&2
      echo "        --test <file_stem>; for a module use a path like mod::tests." >&2
      exit 2
    fi

    echo
    if [ "$rc" -eq 0 ]; then
      echo "mutate: SURVIVOR — the test PASSED against broken code. It is not a test."
      echo "mutate: strengthen the assertion (assert the value, not that something exists),"
      echo "        or you mutated a line the test never reaches."
      exit 1
    fi
    echo "mutate: KILLED (cargo test exit $rc) — the test caught the mutation."
    echo "mutate: record in the PR body which line was mutated and that the test failed."
    exit 0
    ;;

  *) usage ;;
esac
