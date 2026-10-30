#!/usr/bin/env bash
# Simulate CHM13 reads matched to a real CHM13 library, for the telomere comparison (#861).
#
# The real library is staged and aligned by stage_raw_pairs.sbatch (TAG=chm13_novaseq). This
# fits the fragment model from that BAM, because the read taper at a chromosome end is set by
# fragment length, then simulates at the real library's depth and read length. It writes the
# reads where stage_raw_pairs.sbatch looks for local FASTQs, so the simulated reads are aligned
# by exactly the same bwa-mem2 | fixmate | sort | markdup path as the real ones:
#
#   id=$(sbatch --parsable scripts/delta/telomere_sim.sbatch)
#   id2=$(TAG=chm13_sim REFERENCE=$SCRATCH/neat_data/chm13/hs1.fa \
#     sbatch --parsable --dependency=afterok:$id \
#     scripts/delta/stage_raw_pairs.sbatch)
#   sbatch --dependency=afterok:$id2 scripts/delta/telomere_ends.sbatch
#
# NO chrY. CHM13 is a 46,XX cell line, and hs1 carries HG002's chrY. Simulating from it would
# put Y reads, including the Y copies of the pseudoautosomal telomeres, in the simulated library
# and nowhere in the real one. So the reads are drawn from hs1 minus chrY, and both libraries
# are aligned to the full hs1, as the real one already was.
#
# Quality and error models are the shipped defaults (HG002). Only the fragment model is fitted
# from this library. Depth comparisons are normalized by interior control windows downstream,
# so the depth match here need not be exact.

#SBATCH --job-name=eidolon-telosim
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=telosim-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chm13/hs1.fa}"
REAL_BAM="${REAL_BAM:-$SCRATCH/raw_pairs/chm13_novaseq/chm13_novaseq.bam}"
WORK="${WORK:-$SCRATCH/raw_pairs/chm13_sim}"
EXCLUDE="${EXCLUDE:-chrY}"
READ_LEN="${READ_LEN:-151}"
EIDOLON="${EIDOLON:-$SCRATCH/cargo-target/eidolon/release/eidolon}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

module load samtools/1.22-cce19.0.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH" >&2; exit 1; }
for f in "$REFERENCE" "$REFERENCE.fai" "$REAL_BAM" "$EIDOLON"; do
    [[ -s "$f" ]] || { echo "FATAL: missing input: $f" >&2; exit 1; }
done
flagstat="$(dirname "$REAL_BAM")/flagstat.txt"
[[ -s "$flagstat" ]] || { echo "FATAL: no $flagstat; stage_raw_pairs.sbatch writes it" >&2; exit 1; }
mkdir -p "$WORK"
if [[ -s "$WORK/r1.fastq.gz" || -s "$WORK/r2.fastq.gz" ]]; then
    echo "FATAL: $WORK already holds reads. Delete r1/r2.fastq.gz (and any BAM) to resimulate." >&2
    exit 1
fi

echo "=== telomere_sim (#861) ==="
echo "repo:    $REPO ($(git -C "$REPO" rev-parse --short HEAD))"
echo "binary:  $("$EIDOLON" --version)"

# ── reference without the contigs the sample does not have ──────────────────────
sim_ref="$WORK/sim_reference.fa"
if [[ ! -s "$sim_ref.fai" ]]; then
    keep=$(awk -v ex=",$EXCLUDE," 'index(ex, "," $1 ",") == 0 { print $1 }' "$REFERENCE.fai")
    # shellcheck disable=SC2086
    samtools faidx "$REFERENCE" $keep > "$sim_ref.part"
    mv "$sim_ref.part" "$sim_ref"
    samtools faidx "$sim_ref"
fi
dropped=$(( $(wc -l < "$REFERENCE.fai") - $(wc -l < "$sim_ref.fai") ))
genome=$(awk '{ s += $2 } END { print s }' "$sim_ref.fai")
echo "simulating from $(wc -l < "$sim_ref.fai") contigs, $genome bp ($dropped excluded: $EXCLUDE)"
[[ "$dropped" -eq $(awk -F, '{ print NF }' <<< "$EXCLUDE") ]] \
    || { echo "FATAL: expected to drop $EXCLUDE, dropped $dropped contig(s)" >&2; exit 1; }

# ── depth matched to the real library ───────────────────────────────────────────
primary=$(awk '$NF == "primary" && $4 == "primary" { print $1 }' "$flagstat")
[[ "$primary" =~ ^[0-9]+$ && "$primary" -gt 0 ]] || { echo "FATAL: no primary count in $flagstat" >&2; exit 1; }
# gen-reads takes an integer coverage (per-base read depth, both mates counted), so this rounds;
# telomere_ends.sbatch normalizes each library by its own interior depth.
real_cov=$(awk -v n="$primary" -v l="$READ_LEN" -v g="$genome" 'BEGIN { printf "%.2f", n * l / g }')
COVERAGE="${COVERAGE:-$(awk -v c="$real_cov" 'BEGIN { printf "%.0f", c }')}"
echo "real library: $primary primary reads = ${real_cov}x at ${READ_LEN} bp; simulating ${COVERAGE}x"

# ── fragment model from the real BAM, same settings as the shipped one ──────────
frag="$WORK/frag_chm13.json.gz"
if [[ ! -s "$frag" ]]; then
    printf 'input_file: %s\noutput_file: %s\noverwrite_output: true\ndistribution: discrete\nmin_reads: 100\n' \
        "$REAL_BAM" "$frag" > "$WORK/frag.yml"
    "$EIDOLON" gen-frag-length-model -c "$WORK/frag.yml" > "$WORK/frag.log" 2>&1 \
        || { echo "FATAL: gen-frag-length-model failed:" >&2; tail -20 "$WORK/frag.log" >&2; exit 1; }
fi
bash "$REPO/scripts/delta/validate_frag_model.sh" "$REAL_BAM" "$frag" 2>&1 | tee "$WORK/frag_validate.log"

# ── simulate ────────────────────────────────────────────────────────────────────
cat > "$WORK/sim.yml" <<EOF
reference: $sim_ref
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_model: $frag
mutation_rate: 0
produce_fastq: true
produce_bam: false
produce_vcf: false
overwrite_output: true
output_dir: $WORK
output_filename: sim
rng_seed: chm13 telomere 861
num_threads: $THREADS
EOF
echo "[sim] gen-reads at ${COVERAGE}x ..."
"$EIDOLON" gen-reads -c "$WORK/sim.yml" > "$WORK/sim.log" 2>&1 \
    || { echo "FATAL: gen-reads failed:" >&2; tail -20 "$WORK/sim.log" >&2; exit 1; }
mv "$WORK/sim_r1.fastq.gz" "$WORK/r1.fastq.gz"
mv "$WORK/sim_r2.fastq.gz" "$WORK/r2.fastq.gz"
n1=$(zcat "$WORK/r1.fastq.gz" | wc -l)
n2=$(zcat "$WORK/r2.fastq.gz" | wc -l)
[[ "$n1" -eq "$n2" && "$n1" -gt 0 ]] || { echo "FATAL: r1 $n1 lines, r2 $n2 lines" >&2; exit 1; }
echo "simulated pairs: $((n1 / 4)) (real: $((primary / 2)))"

archive_run "telomere_861/sim" "$WORK" sim.yml sim.log frag.yml frag.log frag_validate.log
echo "DONE. Next: stage_raw_pairs.sbatch with TAG=chm13_sim (see header)."
