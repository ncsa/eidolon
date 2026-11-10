#!/usr/bin/env bash
# Replicated timing of gen-reads' final-FASTQ stage at genome scale: shuffle off, the shuffle as
# merged in #875, and the #876 rewrite (#876).
#
# WHY: one run put the shuffle at 2 h 12 m against 6 m 46 s for plain concatenation, but that
# comparison mixed three things. The shuffled run had TMPDIR on Lustre scratch and the baseline
# did not, and the run fell on a night of filesystem work. This holds TMPDIR constant (scratch,
# for every task) and runs the three conditions side by side, one of each per round, so each
# round sees the same filesystem.
#
# History: run 23620433 (2026-11-07) measured an earlier #876 draft, a single scatter log per
# mate: shuffle-and-write took 66-68 min against v1's 43 min for the whole stage, from ~650k
# small random reads on Lustre. The branch now uses sequential per-bucket files instead.
#
# DESIGN: a 9-task array, %2. Task i runs condition i % 3:
#   0  off  shuffle_fastq: false, $BIN_DEV   (concatenation)
#   1  v1   shuffle on,           $BIN_DEV   (#875: a file per bucket per mate, reopened per append)
#   2  new  shuffle on,           $BIN_876   (the #876 branch under test)
# %2, not %3: three at once ran out of scratch quota (tasks 4 and 5 of 23620433). v1 writes its
# bucket files uncompressed (~164 GB of temp at this scale) on top of the output.
# Each task regenerates the same reads from the same config (CONFIG_DIR/sim.yml, the 21x CHM13
# simulation), so v1 and v2 must produce byte-identical FASTQs; the md5 column checks that at
# scale. Outputs are deleted after each task, so disk stays bounded (~300 GB at peak per round).
#
# USAGE (two builds: develop into the usual target dir, the #876 branch into another):
#   BIN_876=$SCRATCH/cargo-target/eidolon-876/release/eidolon \
#     sbatch scripts/delta/shuffle_bench.sbatch
# then, when all nine are done (each task writes its own row, so concurrent tasks never race):
#   (printf 'task\tcond\tnode\tstage_s\twall_s\tmaxrss_gb\tr1_bytes\tr2_bytes\tr1_md5\n'; \
#     cat $SCRATCH/shuffle_bench/result_*.tsv) | column -t
# and check the exit state of every task, not just the table:
#   sacct -j <id> --format=JobID,State,ExitCode,Elapsed,MaxRSS

#SBATCH --job-name=eidolon-shufbench
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --array=0-8%2
#SBATCH --output=shufbench-%A_%a.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"

CONFIG_DIR="${CONFIG_DIR:-$SCRATCH/raw_pairs/chm13_shuffled_sim}"
BIN_DEV="${BIN_DEV:-$SCRATCH/cargo-target/eidolon/release/eidolon}"
BIN_876="${BIN_876:?set BIN_876 to the #876 branch build}"
OUT_ROOT="${OUT_ROOT:-$SCRATCH/shuffle_bench}"
TASK="${SLURM_ARRAY_TASK_ID:?run as an array job}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

case $((TASK % 3)) in
    0) cond=off; bin="$BIN_DEV"; shuffle=false ;;
    1) cond=v1;  bin="$BIN_DEV"; shuffle=true ;;
    2) cond=new; bin="$BIN_876"; shuffle=true ;;
    *) echo "FATAL: unreachable condition for task $TASK" >&2; exit 1 ;;
esac
for f in "$CONFIG_DIR/sim.yml" "$bin"; do
    [[ -s "$f" ]] || { echo "FATAL: missing $f" >&2; exit 1; }
done
work="$OUT_ROOT/task_${TASK}_$cond"
tmp="$OUT_ROOT/tmp_${TASK}"
rm -rf "$work" "$tmp"
mkdir -p "$work" "$tmp"

echo "=== shuffle_bench task $TASK: $cond on $(hostname) ==="
echo "binary: $("$bin" --version)"
# The same config, its outputs redirected, shuffle set explicitly either way.
grep -v -E '^(output_dir|output_filename|output_prefix|num_threads|shuffle_fastq):' "$CONFIG_DIR/sim.yml" > "$work/c.yml"
printf 'output_dir: %s\noutput_prefix: o\nnum_threads: %s\nshuffle_fastq: %s\n' \
    "$work" "$THREADS" "$shuffle" >> "$work/c.yml"

rc=0
TMPDIR="$tmp" /usr/bin/time -f '%e %M' -o "$work/time.txt" \
    "$bin" gen-reads -c "$work/c.yml" > "$work/run.log" 2>&1 || rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FATAL: gen-reads exited $rc. Last lines:" >&2
    tail -20 "$work/run.log" >&2
    exit 1
fi

# Seconds from "Read generation complete" to the first "Successfully wrote": the final-FASTQ stage.
clean=$(sed 's/\x1b\[[0-9;]*m//g' "$work/run.log")
stage=$(awk '
    function secs(t,  a) { split(t, a, ":"); return a[1] * 3600 + a[2] * 60 + a[3] }
    /Read generation complete/ && !s { s = secs($1) }
    /Successfully wrote fastq/ && s && !e { e = secs($1) }
    END { if (!s || !e) exit 1; d = e - s; if (d < 0) d += 86400; print d }' <<< "$clean") \
    || { echo "FATAL: could not find the stage boundaries in run.log" >&2; exit 1; }
grep -E 'shuffle (scatter|and write)' <<< "$clean" || true
read -r wall maxrss_kb < "$work/time.txt"
size1=$(stat -c %s "$work/o_r1.fastq.gz")
size2=$(stat -c %s "$work/o_r2.fastq.gz")
md5=$(md5sum "$work/o_r1.fastq.gz" | cut -c1-12)

printf '%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\n' "$TASK" "$cond" "$(hostname -s)" "$stage" "$wall" \
    "$(awk -v k="$maxrss_kb" 'BEGIN { print k / 1048576 }')" "$size1" "$size2" "$md5" > "$OUT_ROOT/result_${TASK}.tsv"
echo "stage ${stage}s, wall ${wall}s, peak $((maxrss_kb / 1048576)) GB, r1 $size1 bytes, md5 $md5"

archive_run "shuffle_bench/task_${TASK}_$cond" "$work" c.yml run.log time.txt
rm -rf "$work" "$tmp"
echo "DONE."
