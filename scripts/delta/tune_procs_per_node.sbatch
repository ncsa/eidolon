#!/bin/bash
# SLURM job: how many single-thread eidolon processes saturate ONE Delta node?
#
# The whole-genome strategy is region-sharding (genome_array.sbatch): many
# single-thread eidolon processes, because sharding parallelizes sub-contig and packs
# far denser than threading one process (which caps at the genome's contig count).
# The knob is how many shards to pack per node. On the fixed build (v1.19.1, #340)
# packing is near-linear — a real 25 Mb window holds ~99% per-process efficiency to
# 64/node, dropping to ~75% at 128 (hyperthread oversubscription); the earlier
# "efficient zone <=8" was an artifact of per-process trace-log I/O, since removed.
#
# This launches K identical single-thread eidolon processes CONCURRENTLY on one
# exclusive node for K in PROC_COUNTS, and measures aggregate throughput
# (K / makespan) and per-process efficiency. The knee (where efficiency finally
# falls off) is the optimal shards/node — currently ~64 (one per physical core).
#
# Per-process work is one full GENOME at COVERAGE (a stand-in for a shard);
# outputs go to node-local /tmp (740 GB SSD) so the measurement reflects
# CPU/memory bandwidth, not Lustre I/O.
#
# Prereqs: setup.sh (eidolon built); GENOME staged.
# Usage:
#   sbatch scripts/delta/tune_procs_per_node.sbatch
#   GENOME=$SCRATCH/neat_data/ecoli.fa COVERAGE=30 PROC_COUNTS="1 2 4 8 16 32 64 128" \
#     sbatch scripts/delta/tune_procs_per_node.sbatch

#SBATCH --job-name=eidolon-ppn
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --exclusive
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

GENOME="${GENOME:-$SCRATCH/neat_data/ecoli.fa}"     # per-process workload (a shard stand-in)
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
PROC_COUNTS="${PROC_COUNTS:-1 2 4 8 16 32 64 128}"
OUTDIR="${OUTDIR:-$SCRATCH/ppn_$SLURM_JOB_ID}"
WORK="${WORK_LOCAL:-/tmp/ppn_$SLURM_JOB_ID}"        # node-local SSD for per-proc output
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

source "$HOME/.cargo/env" 2>/dev/null || true
mkdir -p "$OUTDIR" "$WORK"
[[ -f "$GENOME" ]] || { echo "genome not staged: $GENOME" >&2; exit 1; }

echo "=== procs-per-node sweep: $SLURM_JOB_ID ==="
echo "Genome: $GENOME  ${COVERAGE}x   K values: $PROC_COUNTS"
echo "Per-proc output on node-local: $WORK"

TSV="$OUTDIR/ppn.tsv"
printf "procs\tmakespan_s\tproc_median_s\tthroughput_per_s\tefficiency_pct\n" > "$TSV"

# median of integer/float args via sort (no gawk asort dependency)
median() { printf '%s\n' "$@" | sort -n | awk '{v[NR]=$1} END{ if(NR==0){print 0;exit}
    print (NR%2)?v[(NR+1)/2]:(v[NR/2]+v[NR/2+1])/2 }'; }

t1=""   # single-process makespan baseline (for efficiency)
for K in $PROC_COUNTS; do
    rm -rf "$WORK"/*
    echo "── K=$K concurrent single-thread processes ──"
    start=$SECONDS
    for i in $(seq 1 "$K"); do
        d="$WORK/p$i"; mkdir -p "$d"
        cat > "$d/sim.yml" <<YAML
reference: $GENOME
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
output_dir: $d
output_basename: s
rng_seed: ppn $K $i
num_threads: 1
YAML
        # /usr/bin/time per process → per-proc elapsed; run in background.
        ( /usr/bin/time -f '%e' -o "$d/t" "$EIDOLON_BIN" --log-level warn gen-reads -c "$d/sim.yml" >/dev/null 2>"$d/log" ) &
    done
    wait
    makespan=$(( SECONDS - start ))
    [[ -z "$t1" ]] && t1="$makespan"

    walls=(); for i in $(seq 1 "$K"); do
        w=$(cat "$WORK/p$i/t" 2>/dev/null || echo 0); walls+=("$w")
    done
    pmed=$(median "${walls[@]}")
    # throughput = jobs completed per second; efficiency = K=1 makespan / this makespan.
    thru=$(awk -v k="$K" -v m="$makespan" 'BEGIN{printf "%.4f", (m>0)? k/m : 0}')
    eff=$(awk -v t="$t1" -v m="$makespan" 'BEGIN{printf "%.1f", (m>0)? 100*t/m : 0}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$K" "$makespan" "$pmed" "$thru" "$eff" >> "$TSV"
    echo "   makespan=${makespan}s  proc_median=${pmed}s  throughput=${thru}/s  efficiency=${eff}%"
done

echo
echo "=== PROCS-PER-NODE RESULTS ($(basename "$GENOME") ${COVERAGE}x) ==="
column -t -s "$(printf '\t')" "$TSV"
echo
# Optimal K = max throughput (the bandwidth-saturation knee).
best=$(awk -F'\t' 'NR>1 && $4>mx{mx=$4; k=$1} END{print k}' "$TSV")
echo "Peak aggregate throughput at K=$best processes/node."
echo "→ genome_array.sbatch: set --cpus-per-task=$(( 128 / (best>0?best:1) )) (≈$best shards/node), or throttle the array with %$best per node-equivalent."

rm -rf "$WORK"
archive_run ppn "$OUTDIR" ppn.tsv
