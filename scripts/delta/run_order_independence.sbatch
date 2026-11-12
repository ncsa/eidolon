#!/bin/bash
# SLURM job: reproducibility + determinism validation (report §3.10).
#
# Confirms eidolon's output is reproducible and invariant to things that MUST NOT
# change a correct simulation, while being explicit about what legitimately does
# change. All runs use ONE genome (default yeast: multi-contig + fast) at a FIXED
# seed; only the variable under test changes between runs. Each run emits a truth
# VCF (reads off — variants are the signal); we compare header-less, fully-sorted
# variant bodies by md5 (set-equality), and for the name check compare the
# (POS,REF,ALT) multiset (ignoring the legitimately-relabeled CHROM).
#
# Checks (each PASS/FAIL printed in a parseable RESULTS block at the end):
#   A determinism_1thread      same config+seed, run twice, 1 thread   -> identical
#   B determinism_Nthread      same config+seed, run twice, N threads  -> identical
#   C thread_invariant         1-thread run vs N-thread run            -> identical set + count
#   D contig_name_invariant    rename every contig (order+seq fixed)   -> identical (POS,REF,ALT) + count
#   E shard_order_independent  same per-shard seeds, two merge orders  -> identical set
#   F shard_disjoint           no variant emitted by >1 shard window   -> 0 duplicates
# Informational (NOT pass/fail — differs by design):
#   - sharded union vs monolithic (shards use distinct per-(contig,chunk) seeds)
#
# NOTE ON CONTIG ORDER: reproducibility means "same inputs -> same output" (ACM /
# The Turing Way). A reordered FASTA is a *different input file* (different checksum,
# a different artifact), so producing a different realization is correct behavior, not
# a defect — see report §3.10. We therefore do NOT test contig-order invariance; doing
# so would assert a non-invariance and read as a spurious failure. Check D pins the real
# boundary instead: eidolon's realization depends on contig ORDER + SEQUENCE, not NAME
# (the RNG derives one stream per contig from its index, not its name).
#
# Why this matters: validates that HPC region-sharding (genome_array.sbatch +
# merge_shards.sh) is correct and reproducible, and that threading never perturbs
# results. Fast (~few min on yeast).
#
# Prereqs: setup.sh (eidolon build + bioinf env w/ bcftools, samtools). Submit from repo root.
# Usage:
#   sbatch scripts/delta/run_order_independence.sbatch
#   REFERENCE=$SCRATCH/neat_data/yeast.fa COVERAGE=30 MT_THREADS=8 NSHARDS=4 \
#     sbatch scripts/delta/run_order_independence.sbatch

#SBATCH --job-name=eidolon-orderindep
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -uo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # setup_conda/conda_activate + $SCRATCH

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/yeast.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/orderindep_$SLURM_JOB_ID}"
SEED="${SEED:-eidolon order independence}"   # one fixed seed for every monolithic run
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
MT_THREADS="${MT_THREADS:-8}"              # threads for the multi-thread / thread-invariance checks
NSHARDS="${NSHARDS:-4}"                     # target shard count for the sharding checks
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0 2>/dev/null || true
setup_conda
conda_activate bioinf   # bcftools, samtools, seqkit

[[ -x "$EIDOLON_BIN" ]] || { echo "eidolon binary not found: $EIDOLON_BIN (run setup.sh)" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "reference not found: $REFERENCE" >&2; exit 1; }
mkdir -p "$OUTDIR"
[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

echo "=== reproducibility / determinism check: $SLURM_JOB_ID ==="
echo "ref=$REFERENCE  seed='$SEED'  cov=${COVERAGE}x  mt_threads=$MT_THREADS  nshards=$NSHARDS"
echo "out=$OUTDIR"

# ── one eidolon run -> $1/sample.vcf.gz ────────────────────────────────────
# run_eidolon <outdir> <reference> <seed> <threads> [target_bed]
run_eidolon() {
    local wd="$1" ref="$2" seed="$3" threads="$4" bed="${5:-}"
    mkdir -p "$wd"
    {
        echo "reference: $ref"
        echo "read_len: $READ_LEN"
        echo "coverage: $COVERAGE"
        echo "ploidy: 2"
        echo "paired_ended: true"
        echo "fragment_mean: 350"
        echo "fragment_st_dev: 50"
        echo "produce_fastq: false"   # variants only — the determinism signal
        echo "produce_vcf: true"
        echo "produce_bam: false"
        echo "overwrite_output: true"
        echo "output_dir: $wd"
        echo "output_basename: sample"
        echo "rng_seed: $seed"
        echo "sv_rate_scale: 0"
        [[ -n "$bed" ]] && echo "target_bed: $bed"
        echo "num_threads: $threads"
    } > "$wd/sim.yml"
    "$EIDOLON_BIN" --log-level warn gen-reads -c "$wd/sim.yml" > "$wd/run.log" 2>&1
    [[ -f "$wd/sample.vcf.gz" ]] || { echo "  no VCF from $wd (see run.log)" >&2; return 1; }
}

# header-less, fully-sorted body -> stable md5 (set-equality across runs)
md5of() { bcftools view -H "$1" 2>/dev/null | LC_ALL=C sort | md5sum | cut -d' ' -f1; }
nvar()  { bcftools view -H "$1" 2>/dev/null | wc -l; }
# (POS,REF,ALT) multiset md5 — ignores CHROM, for the contig-name check (rename
# relabels CHROM but must not move any variant).
posrefalt() { bcftools view -H "$1" 2>/dev/null | awk '{print $2"\t"$4"\t"$5}' | LC_ALL=C sort | md5sum | cut -d' ' -f1; }

# ── A/B: determinism, 1 thread and N threads ─────────────────────────────
echo; echo "[A] determinism @ 1 thread (run twice)..."
run_eidolon "$OUTDIR/a1" "$REFERENCE" "$SEED" 1
run_eidolon "$OUTDIR/a2" "$REFERENCE" "$SEED" 1
echo "[B] determinism @ $MT_THREADS threads (run twice)..."
run_eidolon "$OUTDIR/b1" "$REFERENCE" "$SEED" "$MT_THREADS"
run_eidolon "$OUTDIR/b2" "$REFERENCE" "$SEED" "$MT_THREADS"
# [C] thread_invariant is A1 vs B1 (single-thread vs multi-thread) — no extra run needed.

# ── D: contig-name invariance ────────────────────────────────────────────
# Rename every contig (prepend "renamed_" to each header name; order + sequence
# untouched). Output relabels CHROM, but the (POS,REF,ALT) realization must be
# identical, since the RNG keys on contig index, not name.
echo "[D] contig-name invariance (rename every contig, order preserved)..."
NAME_REF="$OUTDIR/ref_renamed.fa"
sed 's/^>/>renamed_/' "$REFERENCE" > "$NAME_REF"
samtools faidx "$NAME_REF"
run_eidolon "$OUTDIR/d1" "$NAME_REF" "$SEED" 1

# ── E/F: shard-order independence + disjointness ─────────────────────────
echo "[E/F] sharding (~$NSHARDS windows, SAME per-shard seeds, two merge orders)..."
GENOME_BP=$(awk '{s+=$2} END{print s}' "${REFERENCE}.fai")
SHARD_BP=$(( (GENOME_BP + NSHARDS - 1) / NSHARDS ))
SHARD_DIR="$OUTDIR/shards"
N=$(bash "$REPO_ROOT/scripts/delta/make_shard_beds.sh" "${REFERENCE}.fai" "$SHARD_BP" "$SHARD_DIR")
echo "  made $N shard bed(s)"
SHVCF=()
for (( id=0; id<N; id++ )); do
    bed="$(printf '%s/shard_%04d.bed' "$SHARD_DIR" "$id")"
    sd="$OUTDIR/shard_$(printf '%04d' "$id")"
    run_eidolon "$sd" "$REFERENCE" "$SEED shard $id" 1 "$bed"   # fixed per-shard seed
    v="$sd/sample.vcf.gz"; bcftools index -f -t "$v" 2>/dev/null || true
    SHVCF+=("$v")
done
# merge in natural order and in reversed order; both must give the same set
REVVCF=(); for (( i=${#SHVCF[@]}-1; i>=0; i-- )); do REVVCF+=("${SHVCF[$i]}"); done
bcftools concat -a -Ou "${SHVCF[@]}"  2>/dev/null | bcftools sort -Oz -o "$OUTDIR/merge_fwd.vcf.gz" 2>/dev/null
bcftools concat -a -Ou "${REVVCF[@]}" 2>/dev/null | bcftools sort -Oz -o "$OUTDIR/merge_rev.vcf.gz" 2>/dev/null
# disjointness: a CHROM:POS:REF:ALT key emitted by >1 shard window = double-count
DUP=$(for v in "${SHVCF[@]}"; do bcftools view -H "$v" | awk '{print $1":"$2":"$4":"$5}'; done \
        | LC_ALL=C sort | uniq -d | wc -l)

# ── verdicts ─────────────────────────────────────────────────────────────
A1=$(md5of "$OUTDIR/a1/sample.vcf.gz"); A2=$(md5of "$OUTDIR/a2/sample.vcf.gz")
B1=$(md5of "$OUTDIR/b1/sample.vcf.gz"); B2=$(md5of "$OUTDIR/b2/sample.vcf.gz")
MF=$(md5of "$OUTDIR/merge_fwd.vcf.gz"); MR=$(md5of "$OUTDIR/merge_rev.vcf.gz")
# contig-name: compare (POS,REF,ALT) multiset + count (CHROM legitimately relabels)
BASE_PRA=$(posrefalt "$OUTDIR/a1/sample.vcf.gz"); NAME_PRA=$(posrefalt "$OUTDIR/d1/sample.vcf.gz")
BASE_N=$(nvar "$OUTDIR/a1/sample.vcf.gz");        NAME_N=$(nvar "$OUTDIR/d1/sample.vcf.gz")
verdict() { [[ "$2" == "$3" ]] && echo "PASS" || echo "FAIL"; }
name_verdict() { { [[ "$BASE_PRA" == "$NAME_PRA" && "$BASE_N" == "$NAME_N" ]]; } && echo "PASS" || echo "FAIL"; }

echo
echo "════════════════════════════════════════════════════════════════"
echo "REPRODUCIBILITY RESULTS — $SLURM_JOB_ID ($(basename "$REFERENCE"))"
echo "════════════════════════════════════════════════════════════════"
printf '  %-28s %s\n' "determinism_1thread:"       "$(verdict x "$A1" "$A2")  ($A1 vs $A2)"
printf '  %-28s %s\n' "determinism_${MT_THREADS}thread:" "$(verdict x "$B1" "$B2")  ($B1 vs $B2)"
printf '  %-28s %s\n' "thread_invariant:"          "$(verdict x "$A1" "$B1")  (1thr $A1 vs ${MT_THREADS}thr $B1)"
printf '  %-28s %s\n' "contig_name_invariant:"     "$(name_verdict)  (POS,REF,ALT $BASE_PRA vs $NAME_PRA; n $BASE_N vs $NAME_N)"
printf '  %-28s %s\n' "shard_order_independent:"   "$(verdict x "$MF" "$MR")  (fwd $MF vs rev $MR)"
printf '  %-28s %s\n' "shard_disjoint:"            "$([[ "$DUP" -eq 0 ]] && echo PASS || echo FAIL)  ($DUP duplicate keys)"
echo "  ── informational (differences here are expected, not failures) ──"
printf '  %-28s %s\n' "variant_counts:"            "base=$(nvar "$OUTDIR/a1/sample.vcf.gz") renamed=$(nvar "$OUTDIR/d1/sample.vcf.gz") shard_union=$(nvar "$OUTDIR/merge_fwd.vcf.gz")"
echo "════════════════════════════════════════════════════════════════"
