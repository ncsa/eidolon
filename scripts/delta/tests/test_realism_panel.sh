#!/usr/bin/env bash
# Regression tests for the realism panel's reporting in realism_panel.sbatch.
#
# The metrics themselves are unit-tested in scripts/delta/realism (cargo test). What is NOT
# covered there is the summarising awk in this sbatch — the part that turns per-locus rows
# into medians, ranges and ratios. That is exactly the layer where a bug is invisible: a
# wrong median still prints a plausible number, and nobody checks it against the rows.
#
# The awk block is extracted verbatim from the production script, so the two cannot drift.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PIPELINE:-$HERE/../realism_panel.sbatch}"
SUMMARISER="${SUMMARISER:-$HERE/../realism_summarise.awk}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$SUMMARISER" "$WORK/mutant.awk"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.awk"
        if cmp -s "$SUMMARISER" "$WORK/mutant.awk"; then
            printf '  ERROR   %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if SUMMARISER="$WORK/mutant.awk" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-51s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
median is the mean instead@    mid = a[int((n + 1) / 2)]@    mid = 0; for (i = 1; i <= n; i++) mid += a[i]; mid /= n
range reports the wrong ends@    lo = a[1]@    lo = a[n]
gap is inverted@printf " %9.1fx", real_med / sim_med@printf " %9.1fx", sim_med / real_med
a zero denominator prints a number@if (have_real && have_sim && sim_med == 0 && real_med != 0) printf " %10s", "inf"@if (0) printf " %10s", "inf"
values are not sorted before the median@            if (a[j] < a[i]) { t = a[i]; a[i] = a[j]; a[j] = t }@            if (0) { t = a[i]; a[i] = a[j]; a[j] = t }
a sign change still prints a ratio@else if (have_real && have_sim && (real_med < 0 || sim_med < 0))@else if (0)
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]; exit $?
fi

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
hasnt() { case "$2" in *"$3"*) bad "$1" "does NOT contain: $3" "$2";; *) ok "$1";; esac; }

# The same program the job runs — not a copy, not an extraction.
summarise() {  # < panel.tsv
    awk -F'\t' -f "$SUMMARISER"
}

hdr='label\tcontig\tstart\tend\treads\tspan_bp\tcand_bp\tcand_per_mb\timproper_pct\tclip_pct\tmapq0_pct\tdepth_mean\tdepth_vmr\tdepth_excess\tdepth_acf\tins_n\tins_mean\tins_sd\tins_skew\tins_p99'

# Five REAL loci with a deliberately skewed spread — 4 low and 1 high, mirroring the real
# chr22 measurement (VMR 5.51/6.87/7.85/8.88 at four loci, 36.10 at a fifth). The median
# must land in the cluster, not be dragged by the outlier; that difference is the whole
# reason this is a median and not a mean.
{
  printf "$hdr\n"
  # Deliberately NOT in ascending order. A pre-sorted fixture makes the sort a no-op, and
  # mutating it away survives — which it did, first time.
  for v in 36.10 5.51 8.88 6.87 7.85; do
    printf 'REAL\tchr22\t1\t2\t100\t400000\t260\t650.0\t0.0300\t0.0070\t0.0000\t247.00\t%s\t0.1400\t0.800\t50\t550.0\t157.0\t+0.210\t958\n' "$v"
  done
  for i in 1 2 3 4 5; do
    printf 'SIMULATED\tchr22\t1\t2\t100\t400000\t0\t0.0\t0.0000\t0.0000\t0.0400\t247.00\t1.04\t0.0012\t-0.002\t50\t400.0\t89.0\t+0.070\t608\n'
  done
} > "$WORK/panel.tsv"

out="$(awk -F"\t" -f "$SUMMARISER" "$WORK/panel.tsv")"

echo "=== the median resists a single outlier locus ==="
# Median of 5.51 6.87 7.85 8.88 36.10 is 7.85. The MEAN would be 13.04 — which is not a
# value any locus has, and would set a threshold ~66% too high.
has "depth_vmr median is the middle locus, not the mean" "$out" "7.85"
hasnt "the outlier does not become the headline" "$out" "13.04"

echo "=== the range shows how much of a gap is the locus ==="
has "depth_vmr range spans the loci actually measured" "$out" "[5.51-36.1]"

echo "=== a gap against zero is reported as such, not as a number ==="
# SIMULATED cand_per_mb is 0. A ratio would be a division by zero; printing a plausible
# number there is how "no artifacts at all" gets mistaken for "a small gap".
has "candidate breakpoints against zero reads as inf" "$out" "inf"

echo "=== ratios are real over simulated, in that order ==="
# depth_excess: 0.14 real / 0.0012 sim = ~117x. Inverted it would read 0.0x and look fine.
has "depth_excess ratio is ~116.7x" "$out" "116.7x"

echo "=== a ratio across a sign change is reported as a difference ==="
# Autocorrelation goes negative in simulated data (-0.002) and positive in real (+0.8).
# 0.8 / -0.002 is "-400x", which describes nothing. The two differ by 0.802, which does.
has "acf gap is a difference, suffixed d" "$out" "0.802d"
hasnt "and is not a meaningless negative ratio" "$out" "-400.0x"

echo "=== both sides are reported, never just the gap ==="
has "the REAL column is present" "$out" "REAL (median"
has "the SIM column is present" "$out" "SIM (median"

echo "=== metrics with no gap still appear ==="
# depth_mean is matched by construction (the job simulates at the real BAM's depth). It must
# still be printed: a metric that vanishes when it agrees hides the fact that it was checked.
has "depth_mean is reported even though it matches" "$out" "depth_mean"

echo "=== the depth cap says which columns it invalidates ==="
# A capped run is cheap and mostly incomparable. If it reported a gap table with no caveat,
# a smoke number would get quoted as a measurement — which is exactly how the "~8% high"
# figure from a 1.2 kb event on H1N1 ended up in a summary table for two weeks.
cap="$(sed -n '/DEPTH_CAPPED" -eq 1 /,/^fi$/p' "$PIPELINE" | head -40)"
has "the cap names depth_excess as still comparable" "$cap" "depth_excess"
has "the cap names depth_vmr as NOT comparable" "$cap" "NOT COMPARABLE"
has "the cap names the artifact rates it invalidates" "$cap" "cand_per_mb"
has "the cap says what it is for" "$cap" "smoke run"

echo "=== the cap is off by default ==="
# Defaulting to capped would make every run cheap and every number unquotable.
has "MAX_SIM_DEPTH defaults to 0" "$(grep -o 'MAX_SIM_DEPTH:-[0-9]*' "$PIPELINE" | head -1)" "MAX_SIM_DEPTH:-0"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
