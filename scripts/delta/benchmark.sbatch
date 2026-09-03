#!/bin/bash
# SLURM job: NEAT 4 (Python) vs eidolon (Rust) throughput benchmark on Delta.
#
# Two phases (see config below):
#   A. Head-to-head eidolon vs NEAT 4 at 1 thread across genome sizes -> the
#      "eidolon is faster + uses less memory" claim (wall-clock + peak RSS).
#   B. eidolon-only thread-scaling sweep on a MULTI-CONTIG genome -> eidolon's
#      parallel speedup curve (NEAT's 1-thread time from A is the baseline).
#
# This replaces the old full (genomes x 8 threads x reps x 2 tools) cross
# product, which redundantly swept threads on single-contig genomes (where
# eidolon can't parallelize) and ran NEAT at every thread count — that overran
# the walltime. Results land in $OUTDIR/median.tsv (one row per genome/threads/
# tool) and per_rep.tsv (raw).
#
# Usage:
#   sbatch scripts/delta/benchmark.sbatch
# Overrides (no edits needed):
#   SCALING_GENOME=chrset:/path/multi.fa SCALING_THREAD_MODES="1 8 64" \
#   REPS=2 sbatch scripts/delta/benchmark.sbatch

#SBATCH --job-name=eidolon-benchmark
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
# Lean A/B design: ~3 genomes x 2 tools @ 1 thread + ~4 eidolon scaling points.
# Far smaller than the old cross-product; 6h is ample (NEAT at 1 thread is the
# slowest single run, ~13 min/rep on chr22 at 10x).
#SBATCH --time=06:00:00
# Exclusive node: this is the headline eidolon-vs-NEAT4 performance comparison, so
# we want the whole node to ourselves — co-located jobs add noise to wall-clock
# and memory-bandwidth measurements. (Justified only for the benchmark; the
# other jobs are I/O- or single-contig-bound and don't need it.)
#SBATCH --exclusive
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Paths (edit these) ──────────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with EIDOLON_REPO.
REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution
DATA_DIR="${DATA_DIR:-$SCRATCH/neat_data}"        # reference FASTAs live here
OUTDIR="${OUTDIR:-$SCRATCH/benchmark_$SLURM_JOB_ID}"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"
CONDA_ENV_NAME="neat4"

# PHASE A — head-to-head eidolon vs NEAT 4 at a single thread count (HEADTOHEAD_THREADS)
# across genome sizes. This is the core "eidolon is faster + leaner" claim. We do
# NOT sweep threads with both tools: most of these genomes are few-contig, and
# eidolon parallelizes per-contig (sub-contig chunking off by default), so a thread
# sweep on a single-contig genome (e.g. chr22) is pure redundancy — every mode
# returns the same number — and NEAT-at-every-thread-count is what blew the
# 8h walltime. Missing files are skipped.
# Overridable as a space-separated "label:path" list, because a caller that needs only
# some of these has no other way to say so. The regression gate reads FOUR numbers out of
# this job — ecoli and chr22, eidolon, 1 thread — and discards the rest, so it was paying
# ~3 h for the GRCh38 arm and ~45 min for the NEAT arm to obtain ~3 min of gated compute,
# and timing out (job 21766280) doing it.
GENOMES=(${GENOMES:-
    "ecoli:$DATA_DIR/ecoli.fa"
    "yeast:$DATA_DIR/yeast.fa"
    "chr22:$DATA_DIR/chr22.fa"
    "human:$DATA_DIR/GRCh38.fa"
})
HEADTOHEAD_THREADS="${HEADTOHEAD_THREADS:-1}"

# PHASE B — eidolon thread-scaling curve. Only meaningful on a MULTI-CONTIG genome
# (eidolon's unit of parallelism is the contig), so chr22 alone can't show scaling.
# eidolon only: NEAT 4's single-thread time from Phase A is the reference baseline;
# sweeping NEAT across threads is the expensive, low-value part we're cutting.
SCALING_GENOME="${SCALING_GENOME:-yeast:$DATA_DIR/yeast.fa}"
SCALING_THREAD_MODES=(${SCALING_THREAD_MODES:-1 4 16 64 128})

# Simulation parameters — must be identical between NEAT 4 and eidolon for
# a fair comparison.
COVERAGE=10
READ_LEN=151
FRAG_MEAN=350
FRAG_SD=50
REPS="${REPS:-3}"   # repetitions per config; median is reported.
# Was a plain `REPS=3`, which silently ignored the environment while line 20 of
# this header documented `REPS=2` as a supported override. A documented knob that
# does nothing is worse than no knob: it sends a caller looking for the real cause
# somewhere else.

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
setup_conda
conda_activate "$CONDA_ENV_NAME"
NEAT_BIN="$(conda run -n "$CONDA_ENV_NAME" which neat 2>/dev/null || echo '')"

NCORES=$SLURM_CPUS_PER_TASK
mkdir -p "$OUTDIR"

echo "=== eidolon benchmark: $SLURM_JOB_ID ==="
echo "Cores available: $NCORES"
echo "eidolon:  $EIDOLON_BIN"
echo "NEAT 4: ${NEAT_BIN:-NOT FOUND}"
echo "Output: $OUTDIR"
echo "Genomes: ${GENOMES[*]}"

PERREP="$OUTDIR/per_rep.tsv"
MEDIAN="$OUTDIR/median.tsv"
printf "genome\tsize_mb\tthreads\ttool\trep\twall_s\tpeak_rss_mb\texit\n" > "$PERREP"
printf "genome\tsize_mb\tthreads\ttool\twall_s_med\tpeak_rss_mb_med\tn_ok\n" > "$MEDIAN"

# ── Helpers ─────────────────────────────────────────────────────────────
parse_time() {
    awk '
        /Elapsed \(wall clock\)/ {
            n=split($NF,a,":");
            if(n==3) wall=a[1]*3600+a[2]*60+a[3];
            else if(n==2) wall=a[1]*60+a[2];
            else wall=a[1];
        }
        /Maximum resident set size/ { rss=$NF/1024 }
        END { printf "%.2f %.1f", wall, rss }
    ' "$1"
}

median_of() {
    printf "%s\n" "$@" | sort -n | awk '{v[NR]=$1}
        END{ if(NR==0){print "NA"; exit}
             if(NR%2){print v[(NR+1)/2]}
             else{printf "%.2f",(v[NR/2]+v[NR/2+1])/2} }'
}

# Size in whole MB, or empty if it cannot be determined.
#
# The `-L` is load-bearing, and job 21766280 is why. `$SCRATCH/neat_data/GRCh38.fa` is a
# relative SYMLINK to `GRCh38.chr.fa`, and `stat -c%s` without `-L` reports the size of the
# LINK -- 13 bytes, the character length of that target path -- not the 3,144,230,986-byte
# genome. 13/1048576 rounds to 0.0 MB, so the NEAT size cap compared `(( 0 > 200 ))`, found
# it false, and ran the ~12 h/rep human arm the cap exists to prevent, burning the whole
# 6 h wall clock for zero perf rows.
#
# The asymmetry is what hid it: `[[ -f ]]` DOES follow symlinks, so the file passed every
# existence check and eidolon read it correctly for 77 minutes. Only the size was wrong.
#
# An unknown size is still returned as EMPTY rather than 0, so it can never reach a numeric
# comparison and fail open a second way. Callers must handle the empty case explicitly.
# Second arg is the printf format, defaulting to the whole-MB integer the size cap needs
# for its `(( ))` comparison. run_one passes %.1f: the per_rep.tsv `size_mb` column has
# always carried one decimal (4.5, 11.7, 49.4), and collapsing it to `4`/`49` would quietly
# change a reported artifact to suit an internal caller.
genome_size_mb() {
    local bytes fmt="${2:-%.0f}"
    bytes="$(stat -Lc%s "$1" 2>/dev/null)" || return 0
    [[ "$bytes" =~ ^[0-9]+$ ]] || return 0
    awk -v b="$bytes" -v f="$fmt" 'BEGIN{printf f, b/1048576}'
}

run_one() {
    local label="$1" fasta="$2" nthreads="$3" tool="$4" cfg="$5" run_dir="$6"
    local size_mb walls=() rsss=() ok=0

    size_mb="$(genome_size_mb "$fasta" '%.1f')"
    if [[ -z "$size_mb" ]]; then
        # Reported, not swallowed: a 0.0 in the TSV is indistinguishable from a tiny
        # genome, and that ambiguity is what hid the cap failure for a whole 6 h job.
        echo "  WARNING: cannot determine size of $fasta — recording size_mb=NA" >&2
        size_mb="NA"
    fi

    for rep in $(seq 1 "$REPS"); do
        local tf="$OUTDIR/${tool}_${label}_t${nthreads}_r${rep}.time"
        local log="$OUTDIR/${tool}_${label}_t${nthreads}_r${rep}.log"
        rm -rf "$run_dir"; mkdir -p "$run_dir"
        echo "  [$tool] $label threads=$nthreads rep=$rep/$REPS"

        # Capture the tool's real exit code without tripping `set -e`: the
        # `|| ec=$?` keeps a failed run from aborting the whole sweep, so one
        # bad config is recorded (exit!=0) rather than killing the job.
        local ec=0
        if [[ "$tool" == "neat" ]]; then
            # Match eidolon's minimal logging (WARNING) for an equal-to-equal comparison —
            # NEAT defaults to INFO. NEAT has no per-read logging (its hot-loop debug
            # lines are commented out), so this barely changes its time; it just makes
            # both tools rigorously symmetric on log verbosity.
            /usr/bin/time -v -o "$tf" "$NEAT_BIN" --log-level WARNING read-simulator -c "$cfg" \
                -o "$run_dir" -p bench > "$log" 2>&1 || ec=$?
        else
            /usr/bin/time -v -o "$tf" "$EIDOLON_BIN" --log-level warn gen-reads -c "$cfg" > "$log" 2>&1 || ec=$?
        fi

        local parsed wall rss
        parsed=$(parse_time "$tf"); wall=${parsed% *}; rss=${parsed#* }
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$label" "$size_mb" "$nthreads" "$tool" "$rep" "$wall" "$rss" "$ec" >> "$PERREP"
        echo "    wall=${wall}s rss=${rss}MB exit=$ec"
        [[ "$ec" -eq 0 ]] && { walls+=("$wall"); rsss+=("$rss"); ok=$((ok+1)); }
    done

    local wmed rmed
    wmed=$(median_of "${walls[@]:-}")
    rmed=$(median_of "${rsss[@]:-}")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$label" "$size_mb" "$nthreads" "$tool" "$wmed" "$rmed" "$ok" >> "$MEDIAN"
}

# run_eidolon / run_neat <label> <fasta> <nthreads>: write the tool's config and
# time it. Factored out so Phase A and Phase B share one definition.
run_eidolon() {
    local label="$1" fasta="$2" nthreads="$3"
    local out="$OUTDIR/eidolon_${label}_t${nthreads}" cfg="$OUTDIR/eidolon_${label}_t${nthreads}.yml"
    cat > "$cfg" <<YAML
reference: $fasta
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: false
produce_bam: false
overwrite_output: true
rng_seed: delta benchmark seed
num_threads: $nthreads
output_dir: $out
output_filename: bench
YAML
    run_one "$label" "$fasta" "$nthreads" "eidolon" "$cfg" "$out"
}
run_neat() {
    local label="$1" fasta="$2" nthreads="$3"
    [[ -n "${NEAT_BIN:-}" ]] || { echo "  NEAT binary not found — skipping NEAT for $label"; return; }
    # NEAT 4 (Python) wall scales with genome size: ~12 min/rep on chr22 (50 Mb), so a
    # full human genome (~3.1 Gb) is ~12 h/rep — a 3-rep sweep can never finish and
    # times out the whole job (job 19891600 died here). Cap the NEAT arm by genome
    # size; eidolon still runs on the large genome for its own numbers, and that eidolon-
    # only result (eidolon does human in ~1 h; NEAT can't) is itself the comparison.
    local neat_cap_mb="${NEAT_MAX_GENOME_MB:-200}"
    local sz_mb; sz_mb="$(genome_size_mb "$fasta")"
    # FAIL CLOSED. An unknown size skips the NEAT arm rather than risking it. Getting this
    # backwards cost job 21766280 its entire 6 h wall clock and produced no perf rows.
    if [[ -z "$sz_mb" ]]; then
        echo "  [neat] SKIP $label: cannot determine genome size, so the ${neat_cap_mb} MB cap"
        echo "         cannot be enforced. Skipping rather than risking an unbounded NEAT arm."
        return
    fi
    # A cap of 0 disables the NEAT arm outright — how the regression suite asks for
    # eidolon-only perf rows without editing this file.
    if (( neat_cap_mb <= 0 )); then
        echo "  [neat] SKIP $label (NEAT_MAX_GENOME_MB=${neat_cap_mb}): NEAT arm disabled"
        return
    fi
    if (( sz_mb > neat_cap_mb )); then
        echo "  [neat] SKIP $label (${sz_mb} MB > NEAT_MAX_GENOME_MB=${neat_cap_mb} MB): NEAT 4 is intractable at this size"
        return
    fi
    local cfg="$OUTDIR/neat_${label}_t${nthreads}.yml"
    cat > "$cfg" <<YAML
reference: $fasta
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_bam: false
produce_vcf: false
rng_seed: 42
threads: $nthreads
overwrite_output: true
YAML
    run_one "$label" "$fasta" "$nthreads" "neat" "$cfg" "$OUTDIR/neat_${label}_t${nthreads}"
}

# ── Phase A: head-to-head (eidolon vs NEAT) at HEADTOHEAD_THREADS ───────────
echo
echo "=== Phase A: head-to-head @ ${HEADTOHEAD_THREADS} thread(s) ==="
for entry in "${GENOMES[@]}"; do
    label="${entry%%:*}"; fasta="${entry#*:}"
    [[ -f "$fasta" ]] || { echo "MISSING genome: $fasta — skipping"; continue; }
    run_eidolon "$label" "$fasta" "$HEADTOHEAD_THREADS"
    run_neat  "$label" "$fasta" "$HEADTOHEAD_THREADS"
done

# ── Phase B: eidolon thread-scaling on a multi-contig genome ───────────────
sc_label="${SCALING_GENOME%%:*}"; sc_fasta="${SCALING_GENOME#*:}"
if [[ -f "$sc_fasta" ]]; then
    echo
    echo "=== Phase B: eidolon thread-scaling on $sc_label (multi-contig) ==="
    ncontigs=$(grep -c '^>' "$sc_fasta" 2>/dev/null || echo '?')
    echo "    ($sc_label has $ncontigs contigs; eidolon parallelizes per contig)"
    for nthreads in "${SCALING_THREAD_MODES[@]}"; do
        # Skip the threads value already measured in Phase A (no need to repeat).
        [[ "$nthreads" == "$HEADTOHEAD_THREADS" ]] && continue
        run_eidolon "$sc_label" "$sc_fasta" "$nthreads"
    done
else
    echo "Phase B skipped: scaling genome not staged ($sc_fasta)."
fi

echo
echo "=== MEDIAN RESULTS (REPS=$REPS) ==="
column -t -s$'\t' "$MEDIAN"
echo
echo "Full per-rep data: $PERREP"

# Persist report artifacts off scratch + record provenance/resources.
archive_run benchmark "$OUTDIR" median.tsv per_rep.tsv