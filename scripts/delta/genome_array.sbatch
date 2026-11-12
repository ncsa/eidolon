#!/bin/bash
# SLURM array job: sharded whole-genome simulation on Delta — EXCLUSIVE-node packing.
#
# Each array task owns ONE whole node (--exclusive) and runs a BATCH of
# SHARDS_PER_NODE windows on it concurrently, one single-threaded eidolon process
# per window. This is the key to predictable performance on Delta:
#
#   Use the SHARED `cpu` partition and a per-shard core slice lands our process
#   alongside other tenants' jobs; they steal memory bandwidth and make per-window
#   time (and thus the whole-genome wall-clock) unpredictable. So we take whole
#   nodes with --exclusive: the only contention is our own SHARDS_PER_NODE
#   processes, which is controllable and reproducible.
#
#   Packing density: the procs-per-node sweep on the FIXED build (v1.19.1, #340 —
#   the earlier "efficient zone <=8" number was an artifact of per-process trace-log
#   I/O, since corrected) shows K independent single-thread processes pack
#   NEAR-LINEARLY — ~99% per-process efficiency to 64/node on a real 25 Mb window,
#   dropping to ~75% at 128 (hyperthread oversubscription). So the default is 64
#   (one per physical core): whole human genome = 306 windows -> ~5 node-tasks.
#   Use 128 for the fewest nodes (~75% eff), or lower K if a genome's largest contig
#   is big and RAM is tight (peak RSS/node ~= K * largest-contig-loaded).
#
# Reads/variants come out in absolute genome coordinates (target_bed is a
# generation-time filter), so merge_shards.sh concatenates the per-shard outputs
# into a whole-genome dataset.
#
# Workflow:
#   samtools faidx REF.fa
#   N=$(bash scripts/delta/make_shard_beds.sh REF.fa.fai 25000000 $SCRATCH/shards)
#   K=64; T=$(( (N + K - 1) / K ))         # node-tasks = ceil(N / SHARDS_PER_NODE)
#   SHARD_DIR=$SCRATCH/shards OUTROOT=$SCRATCH/wg_$(date +%s) REFERENCE=REF.fa \
#     SHARDS_PER_NODE=$K \
#     sbatch --array=0-$((T-1))%20 scripts/delta/genome_array.sbatch
#   # then, after the array finishes:
#   SHARD_OUTROOT=<OUTROOT> bash scripts/delta/merge_shards.sh
#
# For replicated / swept runs, drive it with run_genome_reps.sh instead.
#
# The `%20` throttle caps concurrently-running NODES (each task = 1 exclusive
# node), tune to your node availability. Straggler-proofing: --requeue is set, so a
# task stuck on a sick node is killed and rescheduled elsewhere — finished shards in
# the batch are skipped on the retry (idempotent), only the unfinished ones rerun.
# On the fixed build a batch finishes in minutes (windows run in ~1-2 min even
# packed 64-deep), so the 1 h --time is generous headroom, not a tight ceiling.

#SBATCH --job-name=eidolon-wgshard
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --requeue
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# NOTE: not `-e` — we manage per-shard failures explicitly so one bad window
# doesn't abort the whole node-batch before the others finish.
set -uo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
SHARD_DIR="${SHARD_DIR:-$SCRATCH/shards}"
OUTROOT="${OUTROOT:-$SCRATCH/wg_$SLURM_ARRAY_JOB_ID}"
SHARDS_PER_NODE="${SHARDS_PER_NODE:-64}"   # windows packed per exclusive node (~99% eff on the fixed build; see header)
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
PLOIDY="${PLOIDY:-2}"
SV_RATE_SCALE="${SV_RATE_SCALE:-0}"        # germline-style by default; >0 for SVs
SEED_ROOT="${SEED_ROOT:-eidolon wg shard}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

source "$HOME/.cargo/env" 2>/dev/null || true

TASK="${SLURM_ARRAY_TASK_ID:?this script must be submitted with --array}"
[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }

# This node-task owns shard ids [FIRST, LAST].
FIRST=$(( TASK * SHARDS_PER_NODE ))
LAST=$(( FIRST + SHARDS_PER_NODE - 1 ))
echo "=== wg node-task $TASK : shards $FIRST..$LAST on $(hostname) (exclusive, ${SHARDS_PER_NODE}/node) ==="
echo "ref=$REFERENCE  cov=${COVERAGE}x  sv_rate_scale=$SV_RATE_SCALE  outroot=$OUTROOT"

# Simulate one window. Single-threaded (fastest per-process on Delta); per-shard
# distinct seed so the shards are independent draws; target_bed restricts
# generation to this window. Idempotent: skips a window already complete (so a
# requeued task only redoes the unfinished windows in its batch).
run_one_shard() {
    local id="$1"
    local bed; bed="$(printf '%s/shard_%04d.bed' "$SHARD_DIR" "$id")"
    [[ -f "$bed" ]] || return 0                 # past the end of the shard list
    local outdir="$OUTROOT/shard_$(printf '%04d' "$id")"
    mkdir -p "$outdir"
    if [[ -f "$outdir/s_r1.fastq.gz" && -f "$outdir/s.vcf.gz" ]]; then
        echo "  shard $id already done — skipping."; return 0
    fi
    cat > "$outdir/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: $PLOIDY
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $outdir
output_basename: s
rng_seed: $SEED_ROOT $id
sv_rate_scale: $SV_RATE_SCALE
target_bed: $bed
num_threads: 1
YAML
    echo "  shard $id start: $(cat "$bed")"
    /usr/bin/time -v -o "$outdir/time.txt" "$EIDOLON_BIN" --log-level warn gen-reads -c "$outdir/sim.yml" \
        > "$outdir/run.log" 2>&1
}

# Launch the whole batch concurrently on the exclusive node, then reap.
pids=(); ids=()
for (( id=FIRST; id<=LAST; id++ )); do
    run_one_shard "$id" &
    pids+=("$!"); ids+=("$id")
done

rc=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        e="$(grep -m1 'Elapsed' "$OUTROOT/shard_$(printf '%04d' "${ids[$i]}")/time.txt" 2>/dev/null || true)"
        echo "  shard ${ids[$i]} complete: $e"
    else
        echo "  shard ${ids[$i]} FAILED (see shard_$(printf '%04d' "${ids[$i]}")/run.log)" >&2
        rc=1
    fi
done

# Non-zero exit lets SLURM --requeue retry this node-batch on another node;
# the idempotent skip above means only the unfinished windows rerun.
exit $rc
