#!/usr/bin/env bash
# End-to-end tests for realism-panel's read-length band (#672).
#
# WHY THIS EXISTS. `test_realism_panel.sh` asserts that both arms are *passed* the same band.
# That is the drift guard, and it says nothing about whether the band does anything. This
# suite runs the real binary over a real BAM and checks the counts.
#
# WHAT IT IS GUARDING. The panel compared a trimmed real BAM against untrimmed simulated
# reads for every run to date. Measured on job 21830618, over the 61 loci in >= 21 bp
# homopolymers: real reads averaged 89.7 bp against a simulated BAM that was 99.6% exactly
# 151 bp, and 96.8% of the real side's soft clips came from reads of 40-80 bp. Matched on
# read length the same comparison reads 0.61% against 0.26% (P = 0.094); unmatched it reads
# 22.9x. `cand_per_mb` had been reporting the unmatched figure.
#
# KNOWN ANSWER, computed without the code under test. The fixture is built by truncating
# exactly every second record of a committed BAM to 60 bp, so the split is known by
# construction: a band admitting only one length must keep that class and charge the other to
# `len_filtered`, and the two bands must be exact complements of each other.
#
# NOT A CLIP MEASUREMENT. The truncated records carry a plain `60M` CIGAR, so they have no
# soft clips by construction. `clip_pct` moving between bands here shows the metric responds
# to the band; its VALUE is an artifact of the fixture and means nothing about real data.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
SRC="${SRC:-$ROOT/scripts/delta/realism/src/reader.rs}"   # overridden by --mutate
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

FIXTURE="$ROOT/eidolon-core/test_data/HG00096.chr20_1Mb.bam"

# ── mutation mode ────────────────────────────────────────────────────────────
#
# The code under test is Rust, so a mutation has to patch the source and rebuild. The crate
# is small enough that this costs a couple of seconds per mutant. Each mutation is verified
# to have APPLIED before its result is believed: a pattern that does not match produces a
# "survivor" that is really an unmodified build (rule 2).
if [[ "${1:-}" == "--mutate" ]]; then
    if ! command -v cargo >/dev/null 2>&1 || ! command -v samtools >/dev/null 2>&1; then
        echo "SKIP: --mutate needs cargo and samtools"
        exit 0
    fi
    survived=0
    orig="$(cat "$SRC")"
    restore() { printf '%s' "$orig" > "$SRC"; }
    trap 'restore; rm -rf "$WORK"' EXIT
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        restore
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$SRC"
        if [[ "$(cat "$SRC")" == "$orig" ]]; then
            printf '  ERROR    %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if ! (cd "$ROOT" && cargo build --release -p realism-panel) >/dev/null 2>&1; then
            printf '  caught   %-52s (mutant does not compile)\n' "$label"; continue
        fi
        if bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-52s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
the band admits every read@        len >= self.min && len <= self.max@        let _ = len; true
the lower bound is exclusive@        len >= self.min && len <= self.max@        len > self.min && len <= self.max
the upper bound is exclusive@        len >= self.min && len <= self.max@        len >= self.min && len < self.max
excluded reads are not counted@                    filtered[bi] += 1;@                    let _ = bi;
membership is resolved after the band@        let keep = band.contains(aln.query_len());@        let keep = true;
MUTATIONS
    restore
    (cd "$ROOT" && cargo build --release -p realism-panel) >/dev/null 2>&1
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]
    exit $?
fi

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

# ── prerequisites ────────────────────────────────────────────────────────────
#
# Asserted, not assumed. A missing binary would make every comparison below compare empty
# against empty, which passes.
BIN="${BIN:-${CARGO_TARGET_DIR:-$ROOT/target}/release/realism-panel}"
if ! command -v samtools >/dev/null 2>&1; then
    echo "SKIP: samtools not on PATH (conda activate bioinf)"
    exit 0
fi
if [[ ! -x "$BIN" ]]; then
    if command -v cargo >/dev/null 2>&1; then
        (cd "$ROOT" && cargo build --release -p realism-panel) >/dev/null 2>&1
    fi
fi
[[ -x "$BIN" ]] || { echo "SKIP: realism-panel not built at $BIN (cargo build --release)"; exit 0; }
[[ -f "$FIXTURE" ]] || { echo "FATAL: BAM fixture missing at $FIXTURE" >&2; exit 1; }

# ── fixture: exactly half the records truncated to 60 bp ─────────────────────
CTG="$(samtools idxstats "$FIXTURE" | awk '$3>0{print $1; exit}')"
[[ -n "$CTG" ]] || { echo "FATAL: no mapped reads in the fixture" >&2; exit 1; }
read -r LO HI <<<"$(samtools view "$FIXTURE" | awk '{if(NR==1||$4<mn)mn=$4; if($4>mx)mx=$4} END{print mn-1, mx+200}')"
printf '%s\t%s\t%s\n' "$CTG" "$LO" "$HI" > "$WORK/r.bed"

samtools view -h "$FIXTURE" \
  | awk 'BEGIN{OFS="\t"} /^@/{print; next} {n++; if(n%2==0){$10=substr($10,1,60); $11=substr($11,1,60); $6="60M"}; print}' \
  > "$WORK/mixed.sam"
samtools view -b "$WORK/mixed.sam" > "$WORK/mixed.bam" 2>/dev/null
samtools index "$WORK/mixed.bam"

# The known answer, computed from the fixture with samtools rather than from the panel.
N_FULL="$(samtools view -F 0x904 "$WORK/mixed.bam" | awk 'length($10)==100' | wc -l | tr -d ' ')"
N_SHORT="$(samtools view -F 0x904 "$WORK/mixed.bam" | awk 'length($10)==60' | wc -l | tr -d ' ')"
N_TOTAL=$((N_FULL + N_SHORT))
[[ "$N_FULL" -gt 100 && "$N_SHORT" -gt 100 ]] \
  || { echo "FATAL: fixture is not mixed ($N_FULL full, $N_SHORT short)" >&2; exit 1; }
echo "=== fixture: $N_FULL reads of 100 bp, $N_SHORT of 60 bp, on $CTG ==="

# field 5 = reads, 6 = len_filtered, 11 = clip_pct
panel() { "$BIN" --bam "$WORK/mixed.bam" --regions "$WORK/r.bed" --label T "$@" 2>"$WORK/err"; }
col() { panel "${@:2}" | awk -F'\t' -v c="$1" 'NR==2{print $c}'; }

echo "=== an open band excludes nothing ==="
eq "every countable read is measured"        "$(col 5)" "$N_TOTAL"
eq "and none are charged to len_filtered"    "$(col 6)" "0"
eq "the open band prints no stderr notice"   "$(wc -c < "$WORK/err" | tr -d ' ')" "0"

echo "=== a band keeps its own class and charges the rest ==="
eq "band 100-100 keeps the full-length reads" \
   "$(col 5 --min-read-len 100 --max-read-len 100)" "$N_FULL"
eq "band 100-100 charges the short reads to len_filtered" \
   "$(col 6 --min-read-len 100 --max-read-len 100)" "$N_SHORT"
eq "band 60-60 keeps the short reads" \
   "$(col 5 --min-read-len 60 --max-read-len 60)" "$N_SHORT"
eq "band 60-60 charges the full-length reads to len_filtered" \
   "$(col 6 --min-read-len 60 --max-read-len 60)" "$N_FULL"
# The two bands must partition the same read set, or one of them is losing records.
a="$(col 5 --min-read-len 100 --max-read-len 100)"
b="$(col 5 --min-read-len 60 --max-read-len 60)"
eq "the two bands partition the reads exactly" "$((a + b))" "$N_TOTAL"

echo "=== boundaries are inclusive ==="
eq "a band whose lower bound is the read length still admits it" \
   "$(col 5 --min-read-len 100 --max-read-len 200)" "$N_FULL"
eq "a band whose upper bound is the read length still admits it" \
   "$(col 5 --min-read-len 61 --max-read-len 100)" "$N_FULL"
# Must not fire: a band starting one base above 100 admits neither class, so it excludes
# every read in the fixture and the run fails. That failing IS the assertion — if 100 were
# admitted by a band whose lower bound is 101, this would return a count instead.
if panel --min-read-len 101 --max-read-len 200 >/dev/null; then
    bad "a band starting above every read length admits nothing" "failure" "exit 0"
else
    ok "a band starting above every read length admits nothing"
fi
has "and it excluded both classes" "$(cat "$WORK/err")" "excluded all $N_TOTAL"

echo "=== the band's cost is reported ==="
panel --min-read-len 100 --max-read-len 100 >/dev/null
has "stderr names the band"            "$(cat "$WORK/err")" "read-length band 100-100"
has "stderr reports what it kept"      "$(cat "$WORK/err")" "kept $N_FULL of $N_TOTAL"
has "stderr reports the excluded share" "$(cat "$WORK/err")" "% excluded"
has "the TSV header carries len_filtered" \
    "$(panel --min-read-len 100 --max-read-len 100 | head -1)" "len_filtered"

echo "=== a band that excludes everything is a hard failure ==="
# Rule 4: "the band took them all" and "this BAM has no reads here" are different problems,
# and reporting the first as the second sends the operator to the wrong file.
if panel --min-read-len 61 --max-read-len 99 >/dev/null; then
    bad "an all-excluding band exits non-zero" "failure" "exit 0"
else
    ok "an all-excluding band exits non-zero"
fi
has "and says the band was the cause" "$(cat "$WORK/err")" "excluded all"
has "and names the band"              "$(cat "$WORK/err")" "61-99"

echo "=== an inverted band is refused before any measurement ==="
if "$BIN" --bam "$WORK/mixed.bam" --regions "$WORK/r.bed" \
     --min-read-len 100 --max-read-len 60 >/dev/null 2>"$WORK/err2"; then
    bad "an inverted band exits non-zero" "failure" "exit 0"
else
    ok "an inverted band exits non-zero"
fi
has "and says it is inverted" "$(cat "$WORK/err2")" "inverted"

echo "=== the band changes clip-derived metrics — the confound itself ==="
# This is what #672 is about: one BAM, one region, two numbers, differing only by which
# reads were admitted. The VALUES are fixture artifacts (the truncated records carry a bare
# 60M and cannot clip); that they DIFFER is the point.
open_clip="$(col 11)"
band_clip="$(col 11 --min-read-len 100 --max-read-len 100)"
if [[ "$open_clip" != "$band_clip" ]]; then
    ok "clip_pct differs between banded and unbanded ($open_clip vs $band_clip)"
else
    bad "clip_pct differs between banded and unbanded" "two different values" "both $open_clip"
fi

# Floor on how many assertions must execute. If this reads low, an assertion stopped running
# rather than started failing. Raise it when adding tests.
MIN_ASSERTIONS=22
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
