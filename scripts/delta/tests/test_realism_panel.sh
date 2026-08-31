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
CONFIGLIB="${CONFIGLIB:-$HERE/../lib_realism_config.sh}"
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
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$CONFIGLIB" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        if cmp -s "$CONFIGLIB" "$WORK/mutant.sh"; then
            printf '  ERROR   %-52s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if CONFIGLIB="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-51s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'CONFIG_MUTATIONS'
a fragment model does not suppress fragment_mean@    if [[ -n "${FRAGMENT_MODEL:-}" ]]; then@    if false; then
model paths are dropped from the config@        [[ -n "$val" ]] && printf '%s: %s\n' "$key" "$val" >> "$out"@        [[ -n "$val" ]] && true
a Normal ceiling is read as unbounded@else ((.Normal.mean + 4 * .Normal.st_dev) | floor) end@else 999999 end
CONFIG_MUTATIONS
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

# ── the simulation config: what the panel actually asks eidolon for ─────────
#
# Running the real function, not grepping the sbatch for a string. The bug being guarded
# against is conditional -- a key emitted when it should not be -- and a grep cannot see a
# condition.
source "$CONFIGLIB"

echo "=== with no models set, the config names eidolon's defaults by omission ==="
(
  unset GC_BIAS_MODEL FRAGMENT_MODEL SEQ_ERROR_MODEL QUALITY_MODEL MUTATION_MODEL GC_NORMALIZE
  write_sim_config "$WORK/c_default.yml" /ref.fa /out seed 8 30 151 400 90
)
cfg="$(cat "$WORK/c_default.yml")"
has   "an untrained run still sets fragment_mean"    "$cfg" "fragment_mean: 400"
has   "an untrained run still sets fragment_st_dev"  "$cfg" "fragment_st_dev: 90"
hasnt "no gc_bias_model key when unset"              "$cfg" "gc_bias_model:"
hasnt "no fragment_model key when unset"             "$cfg" "fragment_model:"
hasnt "no sequence_error_model key when unset"       "$cfg" "sequence_error_model:"

echo "=== a supplied fragment model REPLACES fragment_mean/st_dev, never joins them ==="
# This is the defect: gen_reads/utils/runner.rs prefers explicit mean/st_dev, so emitting
# both silently discards the trained model -- which is what every panel run did, overriding
# eidolon's own shipped empirical distribution with Normal(400, 90).
(
  unset GC_BIAS_MODEL SEQ_ERROR_MODEL QUALITY_MODEL MUTATION_MODEL GC_NORMALIZE
  FRAGMENT_MODEL=/models/frag.json.gz \
    write_sim_config "$WORK/c_frag.yml" /ref.fa /out seed 8 30 151 400 90
)
cfg="$(cat "$WORK/c_frag.yml")"
has   "the trained fragment model is passed through" "$cfg" "fragment_model: /models/frag.json.gz"
hasnt "fragment_mean must NOT also be emitted"       "$cfg" "fragment_mean:"
hasnt "fragment_st_dev must NOT also be emitted"     "$cfg" "fragment_st_dev:"

echo "=== every model knob reaches the config ==="
(
  GC_BIAS_MODEL=/m/gc.json.gz FRAGMENT_MODEL=/m/f.json.gz SEQ_ERROR_MODEL=/m/e.json.gz \
  QUALITY_MODEL=/m/q.json.gz MUTATION_MODEL=/m/mut.json.gz GC_NORMALIZE=true \
    write_sim_config "$WORK/c_all.yml" /ref.fa /out seed 8 30 151 400 90
)
cfg="$(cat "$WORK/c_all.yml")"
has "gc_bias_model"               "$cfg" "gc_bias_model: /m/gc.json.gz"
has "sequence_error_model"        "$cfg" "sequence_error_model: /m/e.json.gz"
has "quality_score_model"         "$cfg" "quality_score_model: /m/q.json.gz"
has "mutation_model"              "$cfg" "mutation_model: /m/mut.json.gz"
has "gc_bias_normalize_coverage"  "$cfg" "gc_bias_normalize_coverage: true"

echo "=== the fragment model's ceiling is read from the model, not assumed ==="
if command -v jq >/dev/null 2>&1; then
  # Known answer: a Discrete model whose largest value is 1094 tops out at 1094.
  printf '{"Discrete":{"distribution":{"values":[300,700,1094],"weights":[0.2,0.6,1.0]}}}' \
    | gzip > "$WORK/disc.json.gz"
  eq_ceiling="$(frag_model_ceiling "$WORK/disc.json.gz")"
  has "a discrete model reports its largest observed length" "$eq_ceiling" "1094"
  # A Normal is unbounded, so the practical ceiling is mean + 4sd = 400 + 360 = 760.
  printf '{"Normal":{"mean":400.0,"st_dev":90.0}}' | gzip > "$WORK/norm.json.gz"
  has "a normal model reports mean + 4sd" "$(frag_model_ceiling "$WORK/norm.json.gz")" "760"
  # Must not fire: an unreadable model must fail rather than invent a ceiling, or the
  # panel would silently skip the MAX_TLEN comparison it exists to make.
  if frag_model_ceiling "$WORK/nope.json.gz" >/dev/null 2>&1; then
    bad "a missing model yields no ceiling" "non-zero exit" "it succeeded"
  else ok "a missing model yields no ceiling"; fi
else
  echo "  SKIP: jq unavailable"
fi

echo "=== the panel refuses a model path that does not exist ==="
# Falling back to a default while reporting a trained model is worse than not running.
guard="$(sed -n '/Refusing rather than silently falling back/,+3p' "$PIPELINE")"
has "the refusal explains why it is fatal" "$guard" "quietly measured defaults"

echo "=== an untrained run says so in the output ==="
prov="$(sed -n '/NO TRAINED MODELS/,+8p' "$PIPELINE")"
has "it names what is being measured"  "$prov" "eidolon AS SHIPPED"
has "it names GC bias as off"          "$prov" "GC bias is off"
has "it points at the fix"             "$prov" "gen-bam-models"

echo "=== the frag ceiling vs MAX_TLEN mismatch is reported ==="
ceil="$(sed -n '/Frag ceiling:/,+14p' "$PIPELINE")"
has "it warns the gap may be the ceiling" "$ceil" "rather than the simulator"
has "it offers the alignment knob"        "$ceil" "ALIGN_TLEN=1"

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
