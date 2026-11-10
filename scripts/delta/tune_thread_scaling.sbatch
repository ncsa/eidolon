#!/bin/bash
# SLURM job: does eidolon's in-process (per-contig) threading SCALE on the fixed
# (v1.19.1) build? — the counterpart to tune_procs_per_node.sbatch.
#
# WHY THIS EXISTS
#   §3.6 concluded in-process threading *regressed* (yeast 1thr 218s -> 16thr 337s)
#   and pivoted the whole-genome strategy to multi-process sharding. But every §3.6
#   scaling number was measured under the old default log-level `trace`, whose
#   per-base debug I/O — emitted PER THREAD — turned added threads into log-write
#   contention. The fixed build's benchmark (§3.2 Phase B) showed yeast threading
#   now scales ~3.9x (1thr 20.7s -> 16thr 5.3s). Open question: does that recovery
#   hold on a genome BIGGER than cache-friendly yeast? eidolon parallelizes per
#   contig, so this needs a MULTI-CONTIG genome (soy: 20 chr + hundreds of
#   scaffolds, ~1 Gb) — chr22 (single contig) can't show it.
#
#   Measures ONE eidolon process at each thread count. Output goes to node-local /tmp
#   (740 GB SSD) so the number reflects CPU/memory bandwidth, not Lustre I/O, and
#   so a big genome's FASTQ never touches the (project-shared) scratch quota.
#
# Prereqs: setup.sh (FIXED eidolon built); a MULTI-CONTIG GENOME staged.
# Usage:
#   GENOME=$SCRATCH/neat_data/soy/soy.fa COVERAGE=10 THREADS="1 2 4 8 16 32 64" \
#     sbatch scripts/delta/tune_thread_scaling.sbatch

#SBATCH --job-name=eidolon-threadscale
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=06:00:00
#SBATCH --exclusive
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

GENOME="${GENOME:-$SCRATCH/neat_data/yeast.fa}"     # MUST be multi-contig (per-contig parallelism)
COVERAGE="${COVERAGE:-10}"
READ_LEN="${READ_LEN:-151}"
THREADS="${THREADS:-1 2 4 8 16 32 64}"
REPS="${REPS:-2}"
OUTDIR="${OUTDIR:-$SCRATCH/threadscale_$SLURM_JOB_ID}"   # small artifacts (TSV) only
WORK="${WORK_LOCAL:-/tmp/threadscale_$SLURM_JOB_ID}"     # node-local SSD for FASTQ
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

source "$HOME/.cargo/env" 2>/dev/null || true
mkdir -p "$OUTDIR" "$WORK"
[[ -f "$GENOME" ]] || { echo "genome not staged: $GENOME" >&2; exit 1; }
ncontigs=$(grep -c '^>' "$GENOME" 2>/dev/null || echo '?')
echo "=== thread-scaling sweep: $SLURM_JOB_ID ==="
echo "Genome: $GENOME ($ncontigs contigs)  ${COVERAGE}x   threads: $THREADS   reps: $REPS"
[[ "$ncontigs" == "1" ]] && echo "WARNING: single-contig genome — eidolon parallelizes per contig, so threads won't help." >&2
echo "Per-thread output on node-local: $WORK"

TSV="$OUTDIR/threadscale.tsv"
printf "threads\twall_med_s\trss_med_mb\tspeedup\tefficiency_pct\n" > "$TSV"

median() { printf '%s\n' "$@" | sort -n | awk '{v[NR]=$1} END{ if(NR==0){print 0;exit}
    if(NR%2)print v[(NR+1)/2]; else printf "%.2f",(v[NR/2]+v[NR/2+1])/2 }'; }
parse_time() { awk '
    /Elapsed \(wall clock\)/ { n=split($NF,a,":"); if(n==3)w=a[1]*3600+a[2]*60+a[3]; else if(n==2)w=a[1]*60+a[2]; else w=a[1] }
    /Maximum resident set size/ { r=$NF/1024 }
    END { printf "%.2f %.1f", w, r }' "$1"; }

t1med=""   # single-thread median wall, baseline for speedup/efficiency
for T in $THREADS; do
    echo "── threads=$T ──"
    walls=(); rsss=()
    for rep in $(seq 1 "$REPS"); do
        d="$WORK/t${T}_r${rep}"; rm -rf "$d"; mkdir -p "$d"
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
rng_seed: threadscale $T $rep
num_threads: $T
YAML
        tf="$d/time"
        /usr/bin/time -v -o "$tf" "$EIDOLON_BIN" --log-level warn gen-reads -c "$d/sim.yml" >"$d/log" 2>&1
        parsed=$(parse_time "$tf"); walls+=("${parsed% *}"); rsss+=("${parsed#* }")
        echo "   rep=$rep wall=${parsed% *}s rss=${parsed#* }MB"
        rm -rf "$d"/*.fastq.gz    # reclaim node-local space between reps
    done
    wmed=$(median "${walls[@]}"); rmed=$(median "${rsss[@]}")
    [[ -z "$t1med" ]] && t1med="$wmed"
    spd=$(awk -v a="$t1med" -v b="$wmed" 'BEGIN{printf "%.2f", (b>0)? a/b : 0}')
    eff=$(awk -v s="$spd" -v t="$T" 'BEGIN{printf "%.1f", (t>0)? 100*s/t : 0}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$T" "$wmed" "$rmed" "$spd" "$eff" >> "$TSV"
    echo "   -> median wall=${wmed}s  speedup=${spd}x  efficiency=${eff}%"
done

echo
echo "=== THREAD-SCALING RESULTS ($(basename "$GENOME") $ncontigs contigs, ${COVERAGE}x) ==="
column -t -s "$(printf '\t')" "$TSV"
echo
best=$(awk -F'\t' 'NR>1 && $4+0>mx{mx=$4; t=$1} END{printf "%s (%.2fx)", t, mx}' "$TSV")
echo "Peak speedup at threads=$best."
echo "Compare to §3.6 (pre-fix, trace logging): threading REGRESSED. If this scales,"
echo "the regression was per-thread log I/O — reweigh threads-vs-shards for whole-genome."

rm -rf "$WORK"
archive_run threadscale "$OUTDIR" threadscale.tsv
