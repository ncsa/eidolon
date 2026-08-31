#!/bin/bash
# SLURM job: measure how far eidolon's reads sit from REAL human sequencing data.
#
# WHY THIS EXISTS. Spot-checked on matched chr22 sequence, real NA12878 against eidolon
# 3.2.0:
#
#     candidate breakpoints / Mb     654.5   vs      0.0
#     improper pairs                 3.15%   vs     0.00%
#     soft clip >= 20 bp             0.67%   vs     0.00%
#     depth excess dispersion       0.1422   vs   0.0012
#
# Zero, not "fewer". Across a full simulated megabase with no planted SVs there is not one
# position where three reads share a clip boundary. A caller tuned on eidolon calibrates
# every clustering threshold and false-positive filter against an empty background, then
# meets ~654 candidates per megabase on real data. That is how a simulator produces a tool
# that works on synthetic data and fails on real data — and it means eidolon-measured
# PRECISION has never predicted real precision. Recall is a separate question: planted
# events are real signal.
#
# SCOPE. Somatic simulation is human — the models are COSMIC/PCAWG-derived and the tissue
# types are human — so realism is measured against real human data and nothing else.
# Germline is species-agnostic, and its correctness is checked a different way: across
# deliberately varied contig shapes, which is where #625 and #607's H1N1 artifacts were
# caught. Those are correctness questions, not realism ones, and this job does not touch
# them.
#
# THIS JOB MEASURES. It does not gate, and it fixes nothing. No metric here has a pass
# threshold, deliberately: locus-to-locus spread in real data is large — VMR across five
# chr22 loci reads 5.51, 6.87, 7.85, 8.88 and 36.10 — so a threshold set from one regime
# would be wrong by 4-6x in either direction. That is the mistake #630 already made once.
# The deliverable is a recorded baseline with its spread; thresholds come after, if at all.
#
# The one hard failure is INABILITY TO MEASURE. A region that names an absent contig, or
# matches no reads, exits non-zero rather than reporting 0.0 artifacts — "no artifacts
# here" and "nothing was measured here" are indistinguishable in a table, and telling them
# apart is the entire point (rule 4).
#
# Prereqs:
#   * REAL_BAM — a real human BAM aligned to $REFERENCE. `bash scripts/delta/stage_hg002.sh`
#     produces one; NA12878 works equally well if you have it staged.
#   * $REFERENCE — the SAME reference the real BAM was aligned to. Comparing across
#     references compares the references.
#   * bwa-mem2 + samtools (bioinf conda env), and a built eidolon + realism-panel.
#
# Usage:
#   # smoke: prove the plumbing cheaply. Capping the depth makes most metrics
#   # incomparable — the job says so, loudly, in its own output.
#   REAL_BAM=... MAX_SIM_DEPTH=30 N_REGIONS=4 sbatch scripts/delta/realism_panel.sbatch
#
#   # real measurement: depths matched, nothing capped.
#   REAL_BAM=$SCRATCH/neat_data/hg002/hg002.bam sbatch scripts/delta/realism_panel.sbatch
#
#   # AS SHIPPED vs TRAINED. With no model set the job measures eidolon out of the box,
#   # which is a real question but not the same one as "what is left after training on
#   # this sample". Build models from the SAME BAM and pass them in:
#   #   eidolon gen-bam-models -c cfg.yml     # frag_length + gc_bias in one BAM walk
#   REAL_BAM=... GC_BIAS_MODEL=$D/gc.json.gz FRAGMENT_MODEL=$D/frag.json.gz \
#     ALIGN_TLEN=1 sbatch scripts/delta/realism_panel.sbatch
#
#   # ADAPTERS=truseq turns on 3' read-through, which is a direct source of soft clips and
#   # also stops short fragments being rejected. clip_pct was 19x below real data without
#   # it, and cand_per_mb was exactly 0 as a consequence.
#   REAL_BAM=... GC_BIAS_MODEL=... FRAGMENT_MODEL=... ADAPTERS=truseq ALIGN_TLEN=1 \
#     sbatch scripts/delta/realism_panel.sbatch
#
#   # ALIGN_TLEN=1 measures both sides over the fragment model's own range. Without it the
#   # panel measures real inserts out to MAX_TLEN while the simulator is capped at the
#   # model's ceiling, and the difference reads as an insert-shape gap that is really the
#   # two filters disagreeing.

#SBATCH --job-name=eidolon-realism
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# Merge stderr into stdout. `exec` with redirections and no command applies them to the
# shell permanently, which is how index_reference_locked silently discarded every
# diagnostic in job 21532276 (#583). Here it is deliberate: one stream, in order, so a
# diagnostic sits next to the step that produced it.
exec 2>&1
echo "[selftest] stderr reaches this log (if this line is missing, diagnostics are being lost)"

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"
source "$REPO_ROOT/scripts/delta/lib_realism_contig.sh"
source "$REPO_ROOT/scripts/delta/lib_realism_config.sh"

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
REAL_BAM="${REAL_BAM:-}"
OUTDIR="${OUTDIR:-$SCRATCH/realism_$SLURM_JOB_ID}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-64}}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-400}"
FRAG_SD="${FRAG_SD:-90}"
N_REGIONS="${N_REGIONS:-10}"
REGION_BP="${REGION_BP:-400000}"
CONTIG="${CONTIG:-}"                 # default: EVERY contig with reads that is also in $REFERENCE
SEED="${SEED:-realism panel}"
# Cap the simulated depth. 0 = match the real BAM exactly, which is what a real measurement
# wants. A smoke run does not: real WGS is often 200-300x, and simulating that across a whole
# chromosome is hours. See the caveat printed below when the cap actually bites.
MAX_SIM_DEPTH="${MAX_SIM_DEPTH:-0}"

# ── Trained models. Unset = eidolon's built-in defaults ─────────────────────
#
# WHY THIS EXISTS. Job 21623108 measured a 16.2x depth_excess gap and a 10.7x insert-skew
# gap and reported them as properties of the simulator. They were properties of a
# CONFIGURATION: gen-reads accepts seven realism inputs and the panel set none of them, so
# the run measured eidolon with GC bias disabled and a hardcoded Gaussian fragment length
# against real human data. That is a legitimate measurement -- it is what a user gets out of
# the box -- but it is not "how close can eidolon get", and the output did not say which of
# the two it was. It does now: every run prints the models it used, defaults included.
#
# Build them from the SAME BAM you are comparing against:
#   eidolon gen-bam-models -c cfg.yml     # frag_length + gc_bias in one BAM walk
#   eidolon gen-seq-error-model -c cfg.yml
GC_BIAS_MODEL="${GC_BIAS_MODEL:-}"
FRAGMENT_MODEL="${FRAGMENT_MODEL:-}"
SEQ_ERROR_MODEL="${SEQ_ERROR_MODEL:-}"
QUALITY_MODEL="${QUALITY_MODEL:-}"
MUTATION_MODEL="${MUTATION_MODEL:-}"
GC_NORMALIZE="${GC_NORMALIZE:-}"     # gc_bias_normalize_coverage; unset = eidolon's default
# Adapter read-through: unset (off) | truseq | nextera. A direct source of soft clips, and
# clip_pct is the metric furthest from real data after training the models -- 19x low on
# job 21636456, which is why cand_per_mb is exactly 0 there (a candidate breakpoint needs
# three reads sharing a clip boundary, and there are not three clipped reads to share one).
ADAPTERS="${ADAPTERS:-}"

# Thresholds. Passed to BOTH sides from one place: a gap measured with different settings
# on each side would be measuring the settings.
MIN_CLIP="${MIN_CLIP:-20}"
MIN_SUPPORT="${MIN_SUPPORT:-3}"
MAX_TLEN="${MAX_TLEN:-2000}"
DEPTH_LAG="${DEPTH_LAG:-500}"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"
PANEL_BIN="$CARGO_TARGET_DIR/release/realism-panel"

[[ -n "$REAL_BAM" ]] || { echo "ERROR: REAL_BAM is required — this job compares against real data, and without it there is nothing to compare to." >&2; exit 1; }
[[ -f "$REAL_BAM" ]]  || { echo "ERROR: REAL_BAM not found: $REAL_BAM" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: REFERENCE not found: $REFERENCE" >&2; exit 1; }
[[ -x "$EIDOLON_BIN" ]] || { echo "ERROR: eidolon not built at $EIDOLON_BIN" >&2; exit 1; }
case "${ADAPTERS:-}" in
    ""|truseq|nextera) ;;
    *) echo "ERROR: ADAPTERS must be truseq or nextera (or unset), got '$ADAPTERS'." >&2
       echo "       gen-reads accepts a 'custom' preset too, but it needs r1/r2 sequences" >&2
       echo "       that this job has no way to pass -- write a config by hand for that." >&2
       exit 1 ;;
esac

for pair in "GC_BIAS_MODEL:$GC_BIAS_MODEL" "FRAGMENT_MODEL:$FRAGMENT_MODEL" \
            "SEQ_ERROR_MODEL:$SEQ_ERROR_MODEL" "QUALITY_MODEL:$QUALITY_MODEL" \
            "MUTATION_MODEL:$MUTATION_MODEL"; do
    name="${pair%%:*}"; path="${pair#*:}"
    [[ -z "$path" || -f "$path" ]] || {
        echo "ERROR: $name is set to '$path', which does not exist." >&2
        echo "       Refusing rather than silently falling back to the default -- a run that" >&2
        echo "       quietly measured defaults while reporting a trained model is worse than" >&2
        echo "       no run at all." >&2
        exit 1; }
done

[[ -x "$PANEL_BIN" ]]   || { echo "ERROR: realism-panel not built at $PANEL_BIN — cargo build --release builds the whole workspace, including it." >&2; exit 1; }

mkdir -p "$OUTDIR"
echo "=== eidolon realism panel: $SLURM_JOB_ID ==="
echo "Real BAM:   $REAL_BAM"
echo "Reference:  $REFERENCE"
echo "Outdir:     $OUTDIR"
echo "Thresholds: min_clip=$MIN_CLIP min_support=$MIN_SUPPORT max_tlen=$MAX_TLEN depth_lag=$DEPTH_LAG"

# ── Provenance: the table below is about a CONFIGURATION, so name it ────────
#
# Rule 4 applied to settings rather than to counts. Every number this job prints is
# conditional on which models were loaded, and a run with everything defaulted looks
# identical in the output to a fully trained one unless it says so.
N_TRAINED=0
echo "Models:"
for pair in "gc_bias:$GC_BIAS_MODEL" "fragment:$FRAGMENT_MODEL" \
            "sequence_error:$SEQ_ERROR_MODEL" "quality_score:$QUALITY_MODEL" \
            "mutation:$MUTATION_MODEL"; do
    name="${pair%%:*}"; path="${pair#*:}"
    if [[ -n "$path" ]]; then
        printf '              %-16s %s\n' "$name" "$path"
        N_TRAINED=$((N_TRAINED + 1))
    else
        printf '              %-16s (eidolon default)\n' "$name"
    fi
done
if [[ -z "$FRAGMENT_MODEL" ]]; then
    printf '              %-16s Normal(%s, %s) -- OVERRIDES the shipped empirical distribution\n' \
           "  fragment via" "$FRAG_MEAN" "$FRAG_SD"
fi
if [[ -n "$ADAPTERS" ]]; then
    printf '              %-16s %s (read-through ON; short fragments are clipped, not rejected)\n' \
           "adapters" "$ADAPTERS"
else
    printf '              %-16s off -- no adapter read-through, so one real source of soft\n' "adapters"
    printf '              %-16s clipping is absent from the simulated side by construction\n' ""
fi
if [[ "$N_TRAINED" -eq 0 ]]; then
    echo ""
    echo "  NO TRAINED MODELS. This measures eidolon AS SHIPPED against real data, which is"
    echo "  what a user gets out of the box -- but it is NOT the gap that remains after"
    echo "  training on this sample. GC bias is off, so depth_vmr/depth_excess/depth_acf"
    echo "  compare a Poisson simulation against real coverage; the fragment length is a"
    echo "  Gaussian, so ins_skew compares a symmetric distribution against a real one."
    echo "  Build models with gen-bam-models and pass GC_BIAS_MODEL/FRAGMENT_MODEL to"
    echo "  measure what is left after that."
fi

# ── The fragment model's ceiling vs what the panel MEASURES ─────────────────
#
# A fragment model has a hard maximum: the builder trims outliers, so on HCC1395 normal it
# came out at 1094 bp while MAX_TLEN defaults to 2000. The simulator then CANNOT emit an
# insert above 1094, the panel happily measures real inserts out to 2000, and the leftover
# difference gets reported as an insert-shape gap. It is not -- it is the two filters
# disagreeing. Measuring your own ceiling and calling it realism is the same error that
# made the first frag-model validation read 99.92% off.
if [[ -n "$FRAGMENT_MODEL" ]]; then
    if command -v jq >/dev/null 2>&1; then
        FRAG_MAX="$(frag_model_ceiling "$FRAGMENT_MODEL" || true)"
        if [[ -n "${FRAG_MAX:-}" && "$FRAG_MAX" != "null" ]]; then
            echo ""
            echo "Frag ceiling: the model tops out at ${FRAG_MAX} bp; panel MAX_TLEN is ${MAX_TLEN} bp"
            if [[ "$FRAG_MAX" -lt "$MAX_TLEN" ]]; then
                if [[ "${ALIGN_TLEN:-0}" != "0" ]]; then
                    MAX_TLEN="$FRAG_MAX"
                    echo "  ALIGN_TLEN set: lowering MAX_TLEN to ${MAX_TLEN} so both sides are measured"
                    echo "  over the same range. The real side now excludes inserts above ${MAX_TLEN},"
                    echo "  which the simulator could not have produced anyway."
                else
                    echo "  WARNING: the simulator cannot produce an insert above ${FRAG_MAX}, but the"
                    echo "  panel measures real inserts up to ${MAX_TLEN}. Part of any ins_sd / ins_p99 /"
                    echo "  ins_skew gap below is that ceiling rather than the simulator's behaviour."
                    echo "  Set ALIGN_TLEN=1 to measure both sides over the model's own range instead."
                fi
            fi
        fi
    else
        echo "  NOTE: jq unavailable, so the fragment model's ceiling was not checked against"
        echo "        MAX_TLEN=${MAX_TLEN}. If the model tops out lower, the insert columns will"
        echo "        overstate the gap."
    fi
fi

# ── Depth of the real BAM, so the simulation can match it ────────────────────
#
# Matching depth is not cosmetic. VMR scales with the mean for multiplicative coverage
# noise, so comparing a 30x simulation against a 250x real BAM inflates the apparent gap by
# ~8x. `depth_excess` in the panel output is depth-independent and is the number to quote
# across datasets, but matching anyway keeps every OTHER metric comparable too — clip rates
# and candidate-breakpoint counts both rise with depth.
# The contig must carry reads in the real BAM AND exist in the reference at the SAME
# length. select_contig lives in lib_realism_contig.sh so its tests exercise this exact
# code; the reasoning for each half of the check is documented there.
[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

samtools idxstats "$REAL_BAM" > "$OUTDIR/real.idxstats" 2>/dev/null \
    || { echo "ERROR: samtools idxstats failed on $REAL_BAM — is it indexed?" >&2; exit 1; }

# Capture to a file, not a pipe: `while read` on the right of a pipeline runs in a
# subshell and loses everything it accumulated, and $(...) under `set -e` swallows the
# non-zero status this function uses to refuse an unmeasurable pairing.
select_contigs "$REAL_BAM" "${REFERENCE}.fai" "$OUTDIR/real.idxstats" "${CONTIG:-}" \
    > "$OUTDIR/contigs.tsv" || exit 1
N_CONTIGS="$(wc -l < "$OUTDIR/contigs.tsv")"
echo "Contigs:    $N_CONTIGS usable (name, length, agreeing between BAM and reference)"
while IFS=$'\t' read -r c l; do printf '              %-10s %s bp\n' "$c" "$l"; done < "$OUTDIR/contigs.tsv"
# ── Regions: evenly spaced, across every usable contig ──────────────────────
#
# Many loci, never one. Real dispersion varies 4-6x between loci on the same chromosome, so
# a single region would set any future envelope wrong in whichever direction it landed.
#
# Spread across ALL usable contigs rather than piling onto the first. The simulation covers
# the whole reference regardless of what gets measured, so on a three-chromosome reference
# measuring one contig pays for three and reads one — and the extra loci come from
# different GC and repeat contexts, which is the spread the baseline is supposed to record.
REGIONS="$OUTDIR/regions.bed"
: > "$REGIONS"
MARGIN=2000000

# place_regions lives in lib_realism_contig.sh so the distribution arithmetic is tested
# rather than trusted. It reports contigs too short to host a region instead of dropping
# them in silence, and refuses to place zero.
# Where the reference has no sequence. Without this every contig's first window lands at
# exactly $MARGIN, which on any acrocentric chromosome (13, 14, 15, 21, 22) is inside the
# p-arm N block — job 21622644 died on chr21:2000000-2400000 for precisely this reason.
NMASK="$OUTDIR/nmask.bed"
n_mask "$REFERENCE" "$OUTDIR/contigs.tsv" > "$NMASK" || exit 1
N_MASKED="$(wc -l < "$NMASK")"
if [[ "$N_MASKED" -gt 0 ]]; then
    echo "N blocks:   $N_MASKED region(s) of the reference are unsequenced and will be avoided"
    awk -F'\t' '{ printf "              %-10s %s-%s (%.1f Mb)\n", $1, $2, $3, ($3-$2)/1e6 }' "$NMASK"
fi

PLACE="$(place_regions "$OUTDIR/contigs.tsv" "$N_REGIONS" "$REGION_BP" "$MARGIN" "$REGIONS" "$NMASK")" || exit 1
N_PLACED="${PLACE%%$'\t'*}"
N_ELIGIBLE="${PLACE##*$'\t'}"
echo "Regions:    $N_PLACED x $REGION_BP bp across $N_ELIGIBLE contig(s)"
# Rule 4: report the denominator. Asking for 10 and placing 6 is a different measurement,
# and the median of 6 loci prints just as confidently as the median of 10. place_regions
# has already named each contig that fell short on stderr; this is the total.
[[ "$N_PLACED" -eq "$N_REGIONS" ]] \
    || echo "  NOTE: asked for $N_REGIONS regions, placed $N_PLACED — see the notes above"

# ── Depth of the real BAM, measured where we know there is sequence ──────────
#
# This probe used to be a fixed $FIRST_CONTIG:1000000-1200000, which is the same bug as the
# one above: on an acrocentric contig that window is entirely N, the probe reads 0x, and the
# job aborts with "nothing to match" on a BAM that is perfectly fine. Probing inside the
# first PLACED region reuses the N check instead of repeating the assumption.
PROBE="$(head -1 "$REGIONS" | awk -F'\t' '{ printf "%s:%d-%d", $1, $2+1, ($2+200000 < $3 ? $2+200000 : $3) }')"
REAL_DEPTH="$(samtools depth -a -r "$PROBE" "$REAL_BAM" \
              | awk '{s+=$3; n++} END{ if(n) printf "%.0f", s/n; else print 0 }')"
[[ "${REAL_DEPTH:-0}" -gt 0 ]] || {
    echo "ERROR: measured 0x depth in the real BAM over $PROBE — nothing to match." >&2
    echo "       That window is real sequence and inside a region carrying reads, so this" >&2
    echo "       is about the BAM rather than about where it was probed." >&2
    exit 1; }
SIM_DEPTH="$REAL_DEPTH"
DEPTH_CAPPED=0
if [[ "$MAX_SIM_DEPTH" -gt 0 && "$REAL_DEPTH" -gt "$MAX_SIM_DEPTH" ]]; then
    SIM_DEPTH="$MAX_SIM_DEPTH"
    DEPTH_CAPPED=1
fi
echo "Probe:      $PROBE"
echo "Real depth: ${REAL_DEPTH}x"
if [[ "$DEPTH_CAPPED" -eq 1 ]]; then
    echo "Sim depth:  ${SIM_DEPTH}x  ** CAPPED by MAX_SIM_DEPTH **"
    echo ""
    echo "  DEPTHS DO NOT MATCH, so most metrics below are NOT comparable between the two"
    echo "  sides. Coverage noise is largely multiplicative, so VMR scales with the mean;"
    echo "  clip rates and candidate-breakpoint counts rise with depth as well. At ${REAL_DEPTH}x"
    echo "  against ${SIM_DEPTH}x those columns compare depth, not realism."
    echo ""
    echo "  STILL COMPARABLE:  depth_excess (Poisson floor removed, depth-normalised),"
    echo "                     ins_skew, and the shape of ins_sd/ins_p99 relative to ins_mean."
    echo "  NOT COMPARABLE:    depth_vmr, cand_per_mb, clip_pct, improper_pct, mapq0_pct."
    echo ""
    echo "  Use a cap for a smoke run that proves the plumbing. Unset it for a number anyone"
    echo "  is going to quote."
    echo ""
else
    echo "Sim depth:  ${SIM_DEPTH}x (matched)"
fi


# ── Simulate at matched depth, no variants ───────────────────────────────────
#
# sv_rate_scale 0 and no input VCF: every artifact the panel finds on this side is one the
# simulator produced on its own, which is exactly the background a caller's false-positive
# filters get calibrated against.
SIMDIR="$OUTDIR/sim"
mkdir -p "$SIMDIR"
write_sim_config "$SIMDIR/config.yml" "$REFERENCE" "$SIMDIR" "$SEED" "$THREADS" \
                 "$SIM_DEPTH" "$READ_LEN" "$FRAG_MEAN" "$FRAG_SD"

echo "[1/3] simulating ${SIM_DEPTH}x with no variants..."
"$EIDOLON_BIN" gen-reads -c "$SIMDIR/config.yml"

echo "[2/3] aligning the simulated reads..."
index_reference_locked "$REFERENCE"
bwa-mem2 mem -t "$THREADS" "$REFERENCE" "$SIMDIR/sim_r1.fastq.gz" "$SIMDIR/sim_r2.fastq.gz" \
  2> "$SIMDIR/bwa.log" | samtools sort -@ 8 -o "$SIMDIR/sim.bam" -
samtools index "$SIMDIR/sim.bam"

# ── Measure both sides, identical settings ───────────────────────────────────
echo "[3/3] measuring..."
# --dump-candidates on the REAL side only. cand_per_mb is the headline gap and a count
# cannot say what produced it: clip boundaries agreed on by many reads at a fixed locus come
# from a shared cause -- a structural variant, a repeat boundary, an assembly error -- and
# scattered read-level artifacts never cluster at all. Adapter read-through brought clip_pct
# from 19x low to exactly 1.0x on job 21638729 and cand_per_mb stayed at zero, which is that
# distinction measured. Classify the real sites before trying to reproduce them.
"$PANEL_BIN" --bam "$REAL_BAM" --regions "$REGIONS" --label REAL \
    --min-clip "$MIN_CLIP" --min-support "$MIN_SUPPORT" \
    --max-tlen "$MAX_TLEN" --depth-lag "$DEPTH_LAG" \
    --dump-candidates "$OUTDIR/real_candidates.tsv" > "$OUTDIR/real.tsv"
"$PANEL_BIN" --bam "$SIMDIR/sim.bam" --regions "$REGIONS" --label SIMULATED \
    --min-clip "$MIN_CLIP" --min-support "$MIN_SUPPORT" \
    --max-tlen "$MAX_TLEN" --depth-lag "$DEPTH_LAG" \
    --dump-candidates "$OUTDIR/sim_candidates.tsv" > "$OUTDIR/sim.tsv"

cat "$OUTDIR/real.tsv" > "$OUTDIR/panel.tsv"
tail -n +2 "$OUTDIR/sim.tsv" >> "$OUTDIR/panel.tsv"

# ── Report ───────────────────────────────────────────────────────────────────
#
# Per-metric: median across loci, and the min-max range. The RANGE is not decoration — it
# is what says whether a gap is a property of the data or of the locus someone happened to
# pick, and it is what any future threshold has to be derived from.
echo ""
echo "════════════════════════════════════════════════════════════════"
if [[ "$DEPTH_CAPPED" -eq 1 ]]; then
    echo "REALISM PANEL — $N_PLACED loci across $N_ELIGIBLE contig(s) — REAL ${REAL_DEPTH}x vs SIM ${SIM_DEPTH}x (CAPPED)"
else
    echo "REALISM PANEL — $N_PLACED loci across $N_ELIGIBLE contig(s), ${REAL_DEPTH}x both sides"
fi
echo ""
awk -F'\t' -f "$REPO_ROOT/scripts/delta/realism_summarise.awk" "$OUTDIR/panel.tsv"
echo ""
# What the candidate sites look like, not just how many. Bands are descriptive, not a
# verdict: they say which explanation to go looking for.
if [[ -s "$OUTDIR/real_candidates.tsv" ]]; then
    echo ""
    echo "Candidate sites, REAL side — where the cand_per_mb gap comes from:"
    awk -F'\t' 'NR>1 {
            n++
            if ($8 >= 0.5) rep++            # half the overlapping reads unmappable
            else if ($7 >= 0.30) sv++       # a third of the depth agrees on one boundary
            else scat++
            sup += $4; if ($4 > mx) mx = $4
        }
        END {
            if (!n) { print "              none"; exit }
            printf "              %d sites, median support %.0f reads (max %d)\n", n, sup/n, mx
            printf "              %5.1f%%  in repeat-like context (>=50%% of local depth is MAPQ 0)\n", rep*100/n
            printf "              %5.1f%%  strongly supported in mappable sequence (>=30%% of depth agrees)\n", sv*100/n
            printf "              %5.1f%%  weakly supported, mappable — read-level noise\n", scat*100/n
        }' "$OUTDIR/real_candidates.tsv"
    echo "              full list: $OUTDIR/real_candidates.tsv"
fi

echo "Per-locus rows: $OUTDIR/panel.tsv"
echo ""
if [[ "$DEPTH_CAPPED" -eq 1 ]]; then
    echo "DEPTHS WERE NOT MATCHED (${REAL_DEPTH}x vs ${SIM_DEPTH}x). Only depth_excess and the"
    echo "insert-shape columns are comparable; the rest compare depth. See the note above."
    echo ""
fi
echo "This job MEASURES. Nothing above is a pass or a fail. A ratio is a gap to decide"
echo "about, and the [min-max] beside each median is how much of it is the locus rather"
echo "than the data — real dispersion varies several-fold between loci on one chromosome."
echo "════════════════════════════════════════════════════════════════"

archive_run "realism" "$OUTDIR" "panel.tsv" "real_candidates.tsv" || true
