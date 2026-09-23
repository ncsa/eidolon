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
adapters are emitted as a flat scalar@printf 'adapters:\n  enabled: true\n  preset: %s\n'@printf 'adapters: %s\n' #
adapters are emitted even when unset@    if [[ -n "${ADAPTERS:-}" ]]; then@    if true; then
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

echo "=== adapters: off by default, and a nested map when set ==="
# clip_pct was 19x below real data with this off, and cand_per_mb exactly 0 as a result.
(
  unset GC_BIAS_MODEL FRAGMENT_MODEL SEQ_ERROR_MODEL QUALITY_MODEL MUTATION_MODEL GC_NORMALIZE ADAPTERS
  write_sim_config "$WORK/c_noad.yml" /ref.fa /out seed 8 30 151 400 90
)
hasnt "no adapters key when unset" "$(cat "$WORK/c_noad.yml")" "adapters"
(


  unset GC_BIAS_MODEL FRAGMENT_MODEL SEQ_ERROR_MODEL QUALITY_MODEL MUTATION_MODEL GC_NORMALIZE
  ADAPTERS=truseq write_sim_config "$WORK/c_ad.yml" /ref.fa /out seed 8 30 151 400 90
)
cfg="$(cat "$WORK/c_ad.yml")"
# gen-reads parses `adapters` as a NESTED MAP -- `adapters: truseq` on one line is silently
# ignored (enabled defaults to false), which would leave read-through off while the run
# reported it on. That is the failure this asserts against.
has "adapters is emitted as a mapping"  "$cfg" "adapters:"
has "with enabled: true"                "$cfg" "  enabled: true"
has "and the preset nested under it"    "$cfg" "  preset: truseq"
hasnt "not as a bare scalar"            "$cfg" "adapters: truseq"

echo "=== the panel refuses a preset it cannot pass through ==="
# `custom` needs r1/r2 sequences the job has no way to carry, so accepting it would set
# adapters.enabled with empty sequences -- read-through on, nothing to read through into.
guard="$(sed -n '/ADAPTERS must be truseq or nextera/,+3p' "$PIPELINE")"
has "custom is named as unsupported here" "$guard" "custom"
has "and says why"                        "$guard" "r1/r2 sequences"

echo "=== an adapters-off run says what that costs ==="
prov="$(sed -n '/no adapter read-through/,+1p' "$PIPELINE")"
has "it names the consequence" "$prov" "soft"

echo "=== the job loads its own tools, it does not inherit them ==="
# Job 21674488 died on `bwa-mem2: not found`. The preflight caught it in a second rather
# than at [2/3], but the job should not have needed the submitting shell to have provided
# it: samtools is a module on Delta and bwa-mem2 is in the bioinf conda env.
setup="$(grep -vE '^[[:space:]]*#' "$PIPELINE" | sed -n '1,60p')"
has "it loads the samtools module"  "$setup" "module load samtools"
has "and activates the conda env"   "$setup" "conda_activate bioinf"
# The loading has to happen BEFORE the preflight, or the preflight rejects tools the job
# would have loaded for itself.
load_line="$(grep -n 'conda_activate bioinf' "$PIPELINE" | head -1 | cut -d: -f1)"
pre_line="$(grep -n 'required tool(s) not found' "$PIPELINE" | head -1 | cut -d: -f1)"
[[ -n "$load_line" && -n "$pre_line" && "$load_line" -lt "$pre_line" ]] \
  && ok "tools are loaded before they are checked" \
  || bad "tools are loaded before they are checked" "load < preflight" "$load_line vs $pre_line"

echo "=== regions without reads are dropped BEFORE the expensive work ==="
# realism-panel refuses a read-less region, but in the MEASURE step -- after simulate and
# align. On NA12878 chr22 a peri-centromeric window that the reference calls real sequence
# held exactly zero reads, and the job would have found out an hour in.
pre="$(sed -n '/Does the REAL BAM actually have reads/,+40p' "$PIPELINE")"
has "it counts reads per region up front"   "$pre" 'samtools view -c -F 0x904 "$REAL_BAM"'
has "an empty region is named, not silent"  "$pre" "NO READS in the real BAM"
has "and dropped rather than fatal"         "$pre" "dropped before simulating"
has "falling too low IS fatal"              "$pre" "not a measurement"
has "the floor is a knob"                   "$pre" 'MIN_REGIONS:-'
# Must not fire: with every region covered, nothing is dropped and the BED is untouched.
has "an all-covered run leaves the BED alone" "$pre" 'rm -f "$KEPT"'
# It has to run BEFORE the simulation, or it saves nothing.
pre_line="$(grep -n 'Does the REAL BAM actually have reads' "$PIPELINE" | cut -d: -f1)"
sim_line="$(grep -n '1/3. simulating' "$PIPELINE" | cut -d: -f1)"
[[ -n "$pre_line" && -n "$sim_line" && "$pre_line" -lt "$sim_line" ]] \
  && ok "the precheck runs before the simulation" \
  || bad "the precheck runs before the simulation" "precheck < simulate" "$pre_line vs $sim_line"

echo "=== the alignment reference can differ from the simulation reference ==="
# MAPQ is a statement about COMPETING PLACEMENTS. The HCC1395 normal BAM carries 2779 @SQ
# lines, 2555 of them alt/unplaced/decoy; aligning simulated reads to a 3-contig subset
# gives them nothing to be ambiguous against, so mapq0_pct reads 0 for reasons that have
# nothing to do with the simulator.
sep="$(sed -n '/ALIGN_TO=/,+14p' "$PIPELINE")"
has "the aligner uses ALIGN_TO, not REFERENCE"   "$sep" 'bwa-mem2 mem -t "$THREADS" "$ALIGN_TO"'
has "and the index is built for the same one"    "$sep" 'index_reference_locked "$ALIGN_TO"'
has "a missing ALIGN_REFERENCE is fatal"         "$sep" "ALIGN_REFERENCE not found"
has "defaulting is called out, not silent"       "$sep" "not comparable"
# Must not fire: unset must behave exactly as before, or every existing run changes meaning.
has "unset falls back to the simulation reference" \
    "$(grep -n 'ALIGN_TO=' "$PIPELINE")" 'ALIGN_TO="${ALIGN_REFERENCE:-$REFERENCE}"'
has "the knob defaults to empty" "$(grep -o 'ALIGN_REFERENCE:-[^}]*' "$PIPELINE" | head -1)" "ALIGN_REFERENCE:-"

echo "=== the provenance block names the alignment reference either way ==="
prov="$(sed -n '/align ref/,+3p' "$PIPELINE")"
has "it says when the two differ"    "$prov" "simulated from"
has "and warns when they do not"     "$prov" "NOT"

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

echo "=== INPUT_VCF reaches the config, and is absent when unset ==="
# Drawn germline SVs arrive this way rather than through a model (#685). Supplied variants
# are added to the de novo mutations, so the run keeps its SNPs and small indels.
write_sim_config "$WORK/c_novcf.yml" /ref.fa /out seed 8 30 151 400 90
hasnt "no input_vcf key when unset" "$(cat "$WORK/c_novcf.yml")" "input_vcf"
INPUT_VCF=/v/drawn.vcf write_sim_config "$WORK/c_vcf.yml" /ref.fa /out seed 8 30 151 400 90
has "INPUT_VCF reaches the config" "$(cat "$WORK/c_vcf.yml")" "input_vcf: /v/drawn.vcf"
# Must-not-fire: supplying a VCF must not silently disable SV generation or anything else.
has "setting INPUT_VCF leaves sv_rate_scale alone" \
    "$(cat "$WORK/c_vcf.yml")" "sv_rate_scale: 0.0"

echo "=== SV_RATE_SCALE is overridable, and defaults to off ==="
# cand_per_mb compares a real genome carrying structural variants against a simulation
# that plants none. The knob has to exist before that can be tested (#684); it must also
# stay off by default so every existing panel run is unchanged.
write_sim_config "$WORK/c_sv_default.yml" /ref.fa /out seed 8 30 151 400 90
has "sv_rate_scale defaults to 0.0" \
    "$(cat "$WORK/c_sv_default.yml")" "sv_rate_scale: 0.0"
SV_RATE_SCALE=0.005 write_sim_config "$WORK/c_sv_on.yml" /ref.fa /out seed 8 30 151 400 90
has "SV_RATE_SCALE reaches the config" \
    "$(cat "$WORK/c_sv_on.yml")" "sv_rate_scale: 0.005"
# Exactly one such key. `hasnt "sv_rate_scale: 0.0"` cannot express this -- that string
# is a prefix of "sv_rate_scale: 0.005" and matches its own replacement.
has "the override replaces the default rather than adding a second key" \
    "$(grep -c '^sv_rate_scale:' "$WORK/c_sv_on.yml")" "1"
# Must-not-fire: setting it must not disturb anything else in the config.
if diff -q <(grep -v '^sv_rate_scale:' "$WORK/c_sv_default.yml") \
           <(grep -v '^sv_rate_scale:' "$WORK/c_sv_on.yml") >/dev/null; then
    ok "no other key changes when SV_RATE_SCALE is set"
else
    bad "no other key changes when SV_RATE_SCALE is set" "identical apart from sv_rate_scale" \
        "$(diff <(grep -v '^sv_rate_scale:' "$WORK/c_sv_default.yml") \
                <(grep -v '^sv_rate_scale:' "$WORK/c_sv_on.yml"))"
fi
# The calibration warning has to survive, or someone reads 1.0 as realistic.
if grep -q "UNION across the cohort" "$HERE/../lib_realism_config.sh"; then
    ok "the sites-VCF calibration warning is recorded"
else
    bad "the sites-VCF calibration warning is recorded" "a warning about cohort union" "absent"
fi

echo "=== both arms get the SAME read-length band (#672) ==="
# THE two-component invariant. The panel reports a ratio between two BAMs; if the arms are
# filtered differently the ratio measures the filter. Nothing asserted they had to match
# before, and the arms ran with different read-length distributions for every run to date.
#
# Assert on the shared variable NAME, not on a value: two literal flag lists that happen to
# agree today are exactly the drift this is meant to prevent.
inv_real="$(grep -c -- '--dump-candidates "$OUTDIR/real_candidates.tsv"' "$PIPELINE")"
inv_sim="$(grep -c -- '--dump-candidates "$OUTDIR/sim_candidates.tsv"' "$PIPELINE")"
[[ "$inv_real" == "1" ]] && ok "there is exactly one REAL panel invocation" \
  || bad "there is exactly one REAL panel invocation" "1" "$inv_real"
[[ "$inv_sim" == "1" ]] && ok "there is exactly one SIMULATED panel invocation" \
  || bad "there is exactly one SIMULATED panel invocation" "1" "$inv_sim"

# Each invocation is a backslash-continued command; join continuations before matching so a
# flag on a different physical line still counts as part of its own invocation.
joined="$WORK/pipeline_joined.sh"
sed -e :a -e '/\\$/N; s/\\\n//; ta' "$PIPELINE" > "$joined"
band_uses="$(grep -c 'BAND_ARGS\[@\]' "$joined")"
[[ "$band_uses" == "2" ]] && ok "both panel invocations expand BAND_ARGS" \
  || bad "both panel invocations expand BAND_ARGS" "2" "$band_uses"
# And they must be the two panel invocations, not one of them twice.
real_line="$(grep -n 'dump-candidates "\$OUTDIR/real_candidates.tsv"' "$joined" | cut -d: -f1)"
sim_line="$(grep -n 'dump-candidates "\$OUTDIR/sim_candidates.tsv"' "$joined" | cut -d: -f1)"
if [[ -n "$real_line" ]] && sed -n "${real_line}p" "$joined" | grep -q 'BAND_ARGS\[@\]'; then
    ok "the REAL arm is filtered by the band"
else
    bad "the REAL arm is filtered by the band" "BAND_ARGS on the REAL invocation" "absent"
fi
if [[ -n "$sim_line" ]] && sed -n "${sim_line}p" "$joined" | grep -q 'BAND_ARGS\[@\]'; then
    ok "the SIMULATED arm is filtered by the band"
else
    bad "the SIMULATED arm is filtered by the band" "BAND_ARGS on the SIMULATED invocation" "absent"
fi
# Non-vacuity: the matcher must be able to see an arm that LACKS the band, or the two checks
# above prove nothing about where BAND_ARGS is.
stripped="$WORK/pipeline_nobands.sh"
sed 's/"${BAND_ARGS\[@\]}" //g' "$joined" > "$stripped"
if sed -n "${real_line}p" "$stripped" | grep -q 'BAND_ARGS\[@\]'; then
    bad "the band check can detect a missing band" "no match after stripping" "still matched"
else
    ok "the band check can detect a missing band"
fi

echo "=== the band is built once, and reports what it dropped ==="
# Built from READ_LEN so the two cannot be set independently; a hardcoded 145 would drift the
# moment someone runs at a different read length.
grep -q 'BAND_LO=$(( READ_LEN - READ_LEN_TOL ))' "$PIPELINE" \
  && ok "the band is derived from READ_LEN, not hardcoded" \
  || bad "the band is derived from READ_LEN, not hardcoded" "BAND_LO from READ_LEN" "absent"
# Rule 4: a filter that drops data has to say how much.
grep -q 'len_filtered' "$HERE/../realism/src/main.rs" \
  && ok "the TSV carries a len_filtered column" \
  || bad "the TSV carries a len_filtered column" "len_filtered emitted" "absent"
# MATCH_READ_LEN=0 must be loud, or a historical-mode run reads like a matched one.
# Anchored on the "Read len:" label: a bare 'NOT MATCHED' also matches the depth-cap banner
# further down, and a mutation sweep showed the unanchored version passing with this notice
# deleted.
grep -q 'Read len:.*NOT MATCHED' "$PIPELINE" \
  && ok "an unmatched run says so in its banner" \
  || bad "an unmatched run says so in its banner" "a NOT MATCHED notice" "absent"
# A lower bound below 1 is a configuration error, not a wide-open band.
grep -q 'BAND_LO" -ge 1' "$PIPELINE" \
  && ok "an impossible lower bound is refused" \
  || bad "an impossible lower bound is refused" "a guard on BAND_LO" "absent"

echo "=== the binary is checked for the flags this script passes (#672) ==="
# Job 21840744 generated its reads, aligned BOTH arms, and died at the measurement step on
# `unknown argument --min-read-len`: a current checkout against a $PANEL_BIN built before
# #687. The sbatch checked that the binary EXISTED and never that it understood the flags.
# That is the same two-component drift test_regression_collector.sh pins for the collector.
#
# Derived, not restated. The flag list is extracted from the invocations themselves, so a
# flag added below without being added to PANEL_FLAGS fails here rather than on Delta.
grep -q '^PANEL_FLAGS=' "$PIPELINE" \
  && ok "the script declares the flag set it depends on" \
  || bad "the script declares the flag set it depends on" "a PANEL_FLAGS list" "absent"
declared="$(grep -m1 '^PANEL_FLAGS=' "$PIPELINE" | cut -d'"' -f2)"
# Flags actually handed to the binary: the two panel invocations plus the BAND_ARGS array
# they expand. --bam/--regions/--label are the required arguments and predate any of this.
passed="$(grep -E 'dump-candidates "\$OUTDIR/(real|sim)_candidates.tsv"|^ *BAND_ARGS=\(--' "$joined" \
          | grep -oE '\-\-[a-z][a-z-]+' | sort -u | grep -vE '^--(bam|regions|label)$')"
[[ -n "$passed" ]] && ok "flags were extracted from the invocations" \
  || bad "flags were extracted from the invocations" "a flag list" "nothing matched"
unprobed=""
for f in $passed; do
    case " $declared " in *" $f "*) ;; *) unprobed="$unprobed $f";; esac
done
if [[ -z "$unprobed" ]]; then
    ok "every flag passed to the binary is one the preflight probes"
else
    bad "every flag passed to the binary is one the preflight probes" "no unprobed flags" "unprobed:$unprobed"
fi
# Non-vacuity: the comparison must be able to see an unprobed flag at all.
case " $declared " in
    *" --a-flag-that-is-not-declared "*) bad "the unprobed check can detect one" "no match" "matched";;
    *) ok "the unprobed check can detect one";;
esac

echo "=== and it is checked BEFORE the expensive steps ==="
# Failing after gen-reads and two alignments costs an allocation; failing at job start costs
# a resubmit. The guard's position in the file is the whole value of it.
guard_ln="$(grep -n 'does not accept' "$PIPELINE" | head -1 | cut -d: -f1)"
work_ln="$(grep -n 'gen-reads -c' "$PIPELINE" | head -1 | cut -d: -f1)"
if [[ -n "$guard_ln" && -n "$work_ln" && "$guard_ln" -lt "$work_ln" ]]; then
    ok "the flag preflight runs before gen-reads (line $guard_ln < $work_ln)"
else
    bad "the flag preflight runs before gen-reads" "guard before gen-reads" "guard=$guard_ln work=$work_ln"
fi
has "the failure names the stale-binary cause" "$(cat "$PIPELINE")" "the binary is older than this script"
# Anchored on the linker export, which appears only in this hint. A bare
# "cargo build --release" also occurs in the older "not built" message, and a mutation
# sweep showed the unanchored version passing with the whole hint deleted.
has "and gives the Delta rebuild command"      "$(cat "$PIPELINE")" "CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc cargo build"

# Floor on how many assertions must execute. This file had none, which is how four
# assertions placed inside a `( ... )` subshell -- where PASS/FAIL increments are
# discarded -- ran without changing the count. Raise it when adding tests.
MIN_ASSERTIONS=88
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
