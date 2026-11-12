#!/bin/bash
# SLURM job: capture a trusted SIMULATOR baseline for eidolon on Delta.
#
# Unlike germline_e2e.sbatch (which measures *caller* recall/precision), this
# job measures the *simulator itself* and records a machine-readable
# `baseline.json` to diff future versions against. It is the safety net for
# feature work that touches the read/variant path — most immediately the
# polyploidy/allele-dosage work (issue #266), whose acceptance criterion is
# "diploid output stays byte-identical".
#
# What it captures (diploid, fixed seed):
#   1. DETERMINISM   — content-identical FASTQ + truth VCF across two thread
#                      counts. Catches the class of bug where output is correct
#                      at 1 thread but races/duplicates above 1 (e.g. the
#                      sub-contig-chunking R2 temp-file collision).
#   2. VARIANT CENSUS — SNP / indel / SV counts + Ts/Tv from the truth VCF.
#   3. COVERAGE      — mean/sd depth and breadth from the golden BAM.
#   4. HET ALLELE FRACTION — observed alt fraction at heterozygous SNP sites
#                      from the golden BAM. Diploid should center ~0.5; this is
#                      the EXACT quantity the dosage work changes, so this run
#                      is the diploid reference #266/#270 validate against.
#
# Uses eidolon's own golden BAM (produce_bam: true) so the measurement reflects
# the simulator's intent, not aligner artifacts. No bwa-mem2/GATK needed.
#
# Prerequisites: run scripts/delta/setup.sh once first.
#
# Usage:
#   sbatch scripts/delta/baseline_capture.sbatch
#   REFERENCE=$SCRATCH/neat_data/chr22.fa COVERAGE=30 sbatch scripts/delta/baseline_capture.sbatch

#SBATCH --job-name=eidolon-baseline
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Parameters ───────────────────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with EIDOLON_REPO.
REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

# chr22 is the right size for a baseline: big enough to be meaningful, small
# enough to run cheaply and repeat per release.
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/baseline_$SLURM_JOB_ID}"

COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
PLOIDY="${PLOIDY:-2}"
# Seed string — determinism + the baseline depend on this being constant across
# runs/versions. When run as a job array (`sbatch --array=1-3`) for triplicate
# replication, the array task id is appended so each replicate uses a distinct
# seed (independent draw) while staying constant within a run; SLURM schedules
# the tasks across (likely different) nodes. Plain submit → no suffix.
SEED="${SEED:-eidolon delta baseline v1}${SLURM_ARRAY_TASK_ID:+ rep$SLURM_ARRAY_TASK_ID}"
# Two thread counts for the determinism check. Output must be content-identical.
THREADS_A=1
THREADS_B="${THREADS_B:-8}"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

# ── Environment ──────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf   # bcftools (mpileup/stats/query)

mkdir -p "$OUTDIR"
DIR_A="$OUTDIR/run_t${THREADS_A}"
DIR_B="$OUTDIR/run_t${THREADS_B}"
mkdir -p "$DIR_A" "$DIR_B"
PREFIX="sample"
BASELINE_JSON="$OUTDIR/baseline.json"

# Prefer the binary's own --version; fall back to the workspace Cargo.toml +
# git describe (always available in the repo) for older binaries without the flag.
EIDOLON_VERSION="$("$EIDOLON_BIN" --version 2>/dev/null || awk -F"'" '/^version/{print $2; exit}' "$REPO_ROOT/Cargo.toml")"
EIDOLON_GIT="$(git -C "$REPO_ROOT" describe --tags --always --dirty 2>/dev/null || echo unknown)"

echo "=== eidolon baseline capture: $SLURM_JOB_ID ==="
echo "eidolon:     $EIDOLON_VERSION"
echo "Reference: $REFERENCE"
echo "Ploidy:    $PLOIDY   Coverage: ${COVERAGE}x   Seed: '$SEED'"
echo "Threads:   A=$THREADS_A  B=$THREADS_B"
echo "Output:    $OUTDIR"

[[ -f "$REFERENCE" ]] || { echo "Reference not found: $REFERENCE" >&2; exit 1; }
[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

# ── Helper: write a gen-reads config ──────────────────────────────────────
# $1 = output dir, $2 = num_threads, $3 = produce_bam (true/false)
write_cfg() {
    local dir="$1" threads="$2" bam="$3"
    cat > "$dir/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: $PLOIDY
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: $bam
overwrite_output: true
output_dir: $dir
output_basename: $PREFIX
rng_seed: $SEED
num_threads: $threads
YAML
}

# ── Stage 1: two simulation runs (same seed, different thread counts) ──────
# Run A also emits the golden BAM (used for coverage + allele-fraction).
if [[ -f "$DIR_A/${PREFIX}_r1.fastq.gz" && -f "$DIR_A/${PREFIX}.bam" ]]; then
    echo "[1/5] Run A outputs present — skipping."
else
    echo "[1/5] Simulation run A (threads=$THREADS_A, +BAM)..."
    write_cfg "$DIR_A" "$THREADS_A" true
    "$EIDOLON_BIN" --log-level warn gen-reads -c "$DIR_A/sim.yml"
fi
if [[ -f "$DIR_B/${PREFIX}_r1.fastq.gz" ]]; then
    echo "[1/5] Run B outputs present — skipping."
else
    echo "[1/5] Simulation run B (threads=$THREADS_B)..."
    write_cfg "$DIR_B" "$THREADS_B" false
    "$EIDOLON_BIN" --log-level warn gen-reads -c "$DIR_B/sim.yml"
fi

TRUTH_VCF="$DIR_A/${PREFIX}.vcf.gz"
GOLDEN_BAM="$DIR_A/${PREFIX}.bam"

# ── Stage 2: determinism (content-identical across thread counts) ─────────
# Compare DECOMPRESSED content so a differing gzip mtime/header can't cause a
# false mismatch. For the VCF, ignore '##' header lines (they may carry a
# timestamp); the #CHROM line + records are what must match.
echo "[2/5] Determinism check (A vs B)..."
md5_fastq_a=$(zcat "$DIR_A/${PREFIX}_r1.fastq.gz" "$DIR_A/${PREFIX}_r2.fastq.gz" | md5sum | awk '{print $1}')
md5_fastq_b=$(zcat "$DIR_B/${PREFIX}_r1.fastq.gz" "$DIR_B/${PREFIX}_r2.fastq.gz" | md5sum | awk '{print $1}')
md5_vcf_a=$(zcat "$TRUTH_VCF"            | grep -v '^##' | md5sum | awk '{print $1}')
md5_vcf_b=$(zcat "$DIR_B/${PREFIX}.vcf.gz" | grep -v '^##' | md5sum | awk '{print $1}')

det_fastq=FAIL; [[ "$md5_fastq_a" == "$md5_fastq_b" ]] && det_fastq=PASS
det_vcf=FAIL;   [[ "$md5_vcf_a"   == "$md5_vcf_b"   ]] && det_vcf=PASS
echo "      FASTQ across threads: $det_fastq ($md5_fastq_a vs $md5_fastq_b)"
echo "      VCF   across threads: $det_vcf ($md5_vcf_a vs $md5_vcf_b)"
n_reads=$(( $(zcat "$DIR_A/${PREFIX}_r1.fastq.gz" | wc -l) / 4 ))

# ── Stage 3: variant census ────────────────────────────────────────────────
echo "[3/5] Variant census..."
[[ -f "${TRUTH_VCF}.tbi" ]] || bcftools index -f -t "$TRUTH_VCF"
stats="$OUTDIR/bcftools_stats.txt"
bcftools stats "$TRUTH_VCF" > "$stats"
n_snp=$(awk -F'\t' '/^SN/ && /number of SNPs:/        {print $NF}' "$stats")
n_indel=$(awk -F'\t' '/^SN/ && /number of indels:/    {print $NF}' "$stats")
ts_tv=$(awk -F'\t' '/^TSTV/ {print $5; exit}' "$stats")
# Symbolic SVs aren't counted by bcftools stats as SNP/indel; count from ALT.
n_sv=$(zcat "$TRUTH_VCF" | grep -v '^#' | grep -cE 'SVTYPE=|<DEL>|<DUP>|<INV>|<INS>|<BND>|<CNV>' || true)
echo "      SNPs=$n_snp  indels=$n_indel  SVs=$n_sv  Ts/Tv=$ts_tv"

# ── Stage 4: coverage from the golden BAM ──────────────────────────────────
echo "[4/5] Coverage..."
if [[ ! -f "${GOLDEN_BAM}.bai" ]]; then
    # eidolon emits a coordinate-sorted golden BAM, but sort defensively.
    if ! samtools index "$GOLDEN_BAM" 2>/dev/null; then
        samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$DIR_A/${PREFIX}.sorted.bam" "$GOLDEN_BAM"
        GOLDEN_BAM="$DIR_A/${PREFIX}.sorted.bam"
        samtools index "$GOLDEN_BAM"
    fi
fi
# Mean/sd over COVERED positions (depth>0) + breadth, so reference N-gaps
# (which eidolon leaves uncovered) don't drag the mean down.
read -r cov_mean cov_sd cov_breadth < <(
    samtools depth -a "$GOLDEN_BAM" | awk '
        { tot++; if ($3>0){ c++; s+=$3; ss+=$3*$3 } }
        END{ if(c==0){print "0 0 0"; exit}
             m=s/c; v=ss/c-m*m; if(v<0)v=0;
             printf "%.2f %.2f %.4f\n", m, sqrt(v), c/tot }'
)
echo "      mean(covered)=${cov_mean}x  sd=${cov_sd}  breadth=${cov_breadth}"

# ── Stage 5: het-SNP allele fraction ───────────────────────────────────────
# The diploid reference for the dosage work: alt fraction at het SNP sites.
echo "[5/5] Het-SNP allele fraction..."
het_sites="$OUTDIR/het_snps.vcf.gz"
bcftools view -v snps -g het -Oz -o "$het_sites" "$TRUTH_VCF"
bcftools index -f -t "$het_sites"
n_het=$(bcftools view -H "$het_sites" | wc -l)

af_hist="$OUTDIR/het_af_histogram.tsv"
af_mean="NA"
if [[ "$n_het" -gt 0 ]]; then
    # Pile up the golden BAM at het SNP sites ONCE, pull per-site AD (ref,alt),
    # compute alt/(ref+alt), bin into deciles (-> histogram file via stderr) and
    # report the mean (-> stdout, captured into af_mean).
    af_mean=$(
        bcftools mpileup -f "$REFERENCE" -R "$het_sites" -a FORMAT/AD -Ou "$GOLDEN_BAM" 2>/dev/null \
            | bcftools query -f '[%AD]\n' \
            | awk -F',' '
                NF>=2 && ($1+$2)>0 { af=$2/($1+$2); sum+=af; n++; b=int(af*10); if(b>9)b=9; h[b]++ }
                END{
                    for(i=0;i<10;i++) printf "%.1f-%.1f\t%d\n", i/10, (i+1)/10, h[i]+0 > "/dev/stderr"
                    if(n>0) printf "%.4f", sum/n; else printf "NA"
                }' 2> "$af_hist"
    )
fi
echo "      het SNPs=$n_het  mean alt-AF=$af_mean  (expect ~0.5 for diploid)"
[[ -f "$af_hist" ]] && { echo "      AF histogram:"; sed 's/^/        /' "$af_hist"; }

# ── Emit baseline.json ─────────────────────────────────────────────────────
cat > "$BASELINE_JSON" <<JSON
{
  "eidolon_version": "$EIDOLON_VERSION",
  "job_id": "$SLURM_JOB_ID",
  "reference": "$(basename "$REFERENCE")",
  "params": { "ploidy": $PLOIDY, "coverage": $COVERAGE, "read_len": $READ_LEN,
              "frag_mean": $FRAG_MEAN, "frag_sd": $FRAG_SD, "seed": "$SEED" },
  "determinism": {
    "threads_compared": [$THREADS_A, $THREADS_B],
    "fastq_identical": "$det_fastq",
    "vcf_identical": "$det_vcf",
    "fastq_md5": "$md5_fastq_a",
    "vcf_records_md5": "$md5_vcf_a"
  },
  "reads": { "pairs": $n_reads },
  "variant_census": { "snps": ${n_snp:-0}, "indels": ${n_indel:-0}, "svs": ${n_sv:-0}, "ts_tv": "${ts_tv:-NA}" },
  "coverage": { "mean_covered": $cov_mean, "sd": $cov_sd, "breadth": $cov_breadth },
  "het_snp_allele_fraction": { "n_sites": $n_het, "mean_alt_af": "$af_mean", "histogram_file": "$(basename "$af_hist")" }
}
JSON

echo
echo "════════════════════════════════════════════════════════════════"
echo "Baseline captured — Job $SLURM_JOB_ID"
echo "  $BASELINE_JSON"
echo "════════════════════════════════════════════════════════════════"
cat "$BASELINE_JSON"
echo
if [[ "$det_fastq" != "PASS" || "$det_vcf" != "PASS" ]]; then
    echo "WARNING: determinism check did NOT pass — investigate before relying on this baseline." >&2
fi
echo "Commit baseline.json alongside the release tag; diff future runs against it."

# Persist report artifacts off scratch + record provenance/resources.
archive_run baseline "$OUTDIR" baseline.json het_af_histogram.tsv bcftools_stats.txt
