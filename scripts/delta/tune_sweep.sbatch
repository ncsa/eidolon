#!/bin/bash
# SLURM job: eidolon HPC tuning-isolation sweep on Delta.
#
# eidolon's per-contig parallelism REGRESSES with threads on Delta's dual-socket
# EPYC (yeast 30x: 1thr 77s -> 4thr 110s) while it SCALES on a single-socket box
# (1thr 77s -> 2thr 57s). The regression is platform-specific, so we isolate the
# cause here — WITHOUT rebuilding eidolon — by crossing thread count with two
# levers:
#   * NUMA placement:   plain  vs  numactl --interleave=all
#   * allocator:        glibc  vs  mimalloc  vs  jemalloc  (via LD_PRELOAD)
# If an allocator swap removes the regression -> malloc lock contention (the hot
# loop's per-fragment get_subseq/reverse_complement allocations). If numactl
# removes it -> NUMA memory placement. This tells us whether the fix is code
# (reduce allocations) or launch-config (numactl/allocator) before we invest in
# either. Results -> $OUTDIR/tune.tsv (variant, threads, wall_s, rss_mb).
#
# Prereqs: setup.sh; a multi-contig genome staged (yeast). mimalloc/jemalloc
# from conda-forge (auto-detected, or set MIMALLOC_SO / JEMALLOC_SO).
#
# Usage:
#   sbatch scripts/delta/tune_sweep.sbatch
#   GENOME=yeast:$SCRATCH/neat_data/yeast.fa THREADS="1 2 4 8 16 64" REPS=2 \
#     sbatch scripts/delta/tune_sweep.sbatch

#SBATCH --job-name=eidolon-tune
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=06:00:00
#SBATCH --exclusive
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

GENOME="${GENOME:-yeast:$SCRATCH/neat_data/yeast.fa}"
GLABEL="${GENOME%%:*}"; GFASTA="${GENOME#*:}"
OUTDIR="${OUTDIR:-$SCRATCH/tune_${SLURM_JOB_ID}}"
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
THREADS="${THREADS:-1 2 4 8 16 32 64}"
REPS="${REPS:-2}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

source "$HOME/.cargo/env"
setup_conda
conda_activate bioinf   # convenient env that has the allocator libs (or install there)

mkdir -p "$OUTDIR"
[[ -f "$GFASTA" ]] || { echo "Genome not staged: $GFASTA" >&2; exit 1; }

# ── Locate optional allocators (no-op variant if absent) ─────────────────
# Always return 0: a non-zero return here would propagate through the bare
# `MIMALLOC_SO=$(find_so ...)` assignment and trip `set -e` (silently exiting
# the whole job before any output) whenever the allocator isn't installed.
find_so() { local p; for p in "$@"; do [[ -e "$p" ]] && { printf '%s\n' "$p"; return 0; }; done; return 0; }
# Look in the ACTIVE bioinf env's lib ($CONDA_PREFIX, set by conda_activate
# above). Do NOT use `conda info --base`/envs/bioinf — on Delta user envs live
# under ~/.conda/envs, not under the module base prefix, so that path is wrong
# and silently found no allocators (only baseline/numa variants ran). Override
# with CONDA_LIB / MIMALLOC_SO / JEMALLOC_SO if needed.
CONDA_LIB="${CONDA_LIB:-${CONDA_PREFIX:-}/lib}"
MIMALLOC_SO="${MIMALLOC_SO:-$(find_so "$CONDA_LIB"/libmimalloc.so* 2>/dev/null)}"
JEMALLOC_SO="${JEMALLOC_SO:-$(find_so "$CONDA_LIB"/libjemalloc.so* 2>/dev/null)}"
HAVE_NUMACTL=0; command -v numactl >/dev/null 2>&1 && HAVE_NUMACTL=1

echo "=== eidolon tuning sweep: $SLURM_JOB_ID ==="
echo "Genome:    $GLABEL ($GFASTA, $(grep -c '^>' "$GFASTA") contigs)  ${COVERAGE}x"
echo "Threads:   $THREADS   Reps: $REPS"
echo "mimalloc:  ${MIMALLOC_SO:-<not found>}"
echo "jemalloc:  ${JEMALLOC_SO:-<not found>}"
echo "numactl:   $([[ $HAVE_NUMACTL == 1 ]] && echo yes || echo no)"

TSV="$OUTDIR/tune.tsv"
printf "variant\tthreads\trep\twall_s\trss_mb\texit\n" > "$TSV"

parse_time() { awk '
    /Elapsed \(wall clock\)/ { n=split($NF,a,":"); if(n==3)w=a[1]*3600+a[2]*60+a[3]; else if(n==2)w=a[1]*60+a[2]; else w=a[1] }
    /Maximum resident set size/ { r=$NF/1024 }
    END { printf "%.2f %.1f", w, r }' "$1"; }

# run_variant <label> <numactl_flag:0|1> <preload_so|"">  <threads>
run_variant() {
    local label="$1" use_numa="$2" preload="$3" nt="$4"
    local run_dir="$OUTDIR/${label}_t${nt}" cfg="$OUTDIR/${label}_t${nt}.yml"
    cat > "$cfg" <<YAML
reference: $GFASTA
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: 2
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
produce_fastq: true
produce_vcf: false
produce_bam: false
overwrite_output: true
rng_seed: tune sweep
num_threads: $nt
output_dir: $run_dir
output_basename: b
YAML
    local pre=(); [[ "$use_numa" == "1" && "$HAVE_NUMACTL" == "1" ]] && pre=(numactl --interleave=all)
    local envp=();   [[ -n "$preload" ]] && envp=(env "LD_PRELOAD=$preload")
    local rep tf ec parsed
    for rep in $(seq 1 "$REPS"); do
        rm -rf "$run_dir"; mkdir -p "$run_dir"
        tf="$OUTDIR/${label}_t${nt}_r${rep}.time"
        ec=0
        /usr/bin/time -v -o "$tf" "${pre[@]}" "${envp[@]}" "$EIDOLON_BIN" --log-level warn gen-reads -c "$cfg" \
            >"$run_dir/log" 2>&1 || ec=$?
        parsed=$(parse_time "$tf")
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$nt" "$rep" "${parsed% *}" "${parsed#* }" "$ec" >> "$TSV"
        echo "  [$label] t=$nt rep=$rep -> wall=${parsed% *}s rss=${parsed#* }MB exit=$ec"
    done
}

# Variants: only those whose tool is present are exercised.
for nt in $THREADS; do
    run_variant "baseline"            0 ""             "$nt"
    [[ $HAVE_NUMACTL == 1 ]]    && run_variant "numa_interleave" 1 ""             "$nt"
    [[ -n "$MIMALLOC_SO" ]]     && run_variant "mimalloc"        0 "$MIMALLOC_SO" "$nt"
    [[ -n "$JEMALLOC_SO" ]]     && run_variant "jemalloc"        0 "$JEMALLOC_SO" "$nt"
    [[ -n "$MIMALLOC_SO" && $HAVE_NUMACTL == 1 ]] && run_variant "mimalloc_numa" 1 "$MIMALLOC_SO" "$nt"
done

echo
echo "=== TUNING SWEEP RESULTS ($GLABEL ${COVERAGE}x) ==="
column -t -s "$(printf '\t')" "$TSV"
echo
echo "Mean wall per (variant,threads) [ok runs only]:"
awk -F'\t' 'NR>1 && $6==0 {k=sprintf("%-16s t%03d",$1,$2); s[k]+=$4; n[k]++}
            END{for(k in s) printf "  %s  %.1fs\n", k, s[k]/n[k]}' "$TSV" | sort

archive_run tune "$OUTDIR" tune.tsv
