#!/bin/bash
# SLURM job: SCN Phase 2 real-data validation for per-variant allele_fraction (#398).
#
# End-to-end test that gen-reads REPRODUCES the real pool's allele-frequency spectrum:
#   real pool VCF (per-site AF from AD) --input_vcf--> gen-reads (honors allele_fraction)
#   --> golden VCF measured AF  ≈  real pool AF   (Joao's question, realism epic #311)
#
# Steps:
#   1. Derive per-site AF for the pool: add INFO/AF from FORMAT/AD on the staged pool VCF,
#      keep SNVs only (AF from AD is cleanest for SNVs; see docs/scn_phase2_af_design.md §7).
#   2. gen-reads with that VCF as input_vcf, mutation_rate=0 and sv_rate_scale=0 so ONLY the
#      supplied sites are placed (nothing de novo confounds the per-site AF comparison), using
#      the Phase-1 models so reads/coverage are realistic.
#   3. Compare golden-VCF AF vs pool AF at shared sites with `eidolon compare-af`.
#
# The feature is automatic: gen-reads honors a per-variant AF whenever the input VCF carries
# INFO/AF (or FORMAT/AD). There is no CLI flag to enable it.
#
# Prereqs: stage_scn.sh (REF + BAM + VCF) and model_builders already run for this strain;
# eidolon built from the feat/variant-allele-fraction branch (setup.sh); the staged pool VCF must
# carry FORMAT/AD (bcftools call -a AD — stage_scn.sh's default). Pass MODELS=<build>/models.
#
# Usage:
#   MODELS=$SCRATCH/modelbuild_<JOB>/models sbatch scripts/delta/run_scn_af_validation.sh
#   MODELS=… ACC=GCA_040805705.1 SRR=SRR27329600 COV=60 sbatch scripts/delta/run_scn_af_validation.sh

#SBATCH --job-name=eidolon-scnaf
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

D="${DATA_DIR:-$SCRATCH/neat_data/scn}"
ACC="${ACC:-GCA_040805935.1}"
SRR="${SRR:-SRR27329602}"
REF="${REFERENCE:-$D/${ACC}.fa}"
POOL_VCF="${POOL_VCF:-$D/${SRR}.vcf.gz}"      # staged short-variant VCF
BAM="${INPUT_BAM:-$D/${SRR}.bam}"             # staged pool BAM (fallback AD source)
MODELS="${MODELS:?set MODELS=<model_builders output>/models}"
COV="${COV:-50}"                               # match the real pool depth so AF noise matches
MIN_DEPTH="${MIN_DEPTH:-20}"                    # gate low-depth sites in the comparison
OUTDIR="${OUTDIR:-$SCRATCH/scn_af_${SLURM_JOB_ID:-manual}}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="${EIDOLON_BIN:-$CARGO_TARGET_DIR/release/eidolon}"

source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
module load bcftools/1.22 2>/dev/null || module load bcftools 2>/dev/null || true

MUT="$MODELS/mut_model.json.gz"
[[ -s "$REF" ]]      || { echo "reference not staged: $REF (run stage_scn.sh)" >&2; exit 1; }
[[ -s "$POOL_VCF" ]] || { echo "pool VCF not found: $POOL_VCF (run stage_scn.sh)" >&2; exit 1; }
[[ -s "$MUT" ]]      || { echo "mut_model not found: $MUT (run model_builders)" >&2; exit 1; }
[[ -x "$EIDOLON_BIN" ]] || { echo "eidolon binary not found: $EIDOLON_BIN (run setup.sh on the feat branch)" >&2; exit 1; }

echo "=== banner: ACC=$ACC SRR=$SRR ref=$REF pool_vcf=$POOL_VCF cov=$COV ==="

mkdir -p "$OUTDIR"

# ── Step 1: per-site AF from AD, SNVs only ───────────────────────────────────
# The pooled allele fraction is the READ fraction FORMAT/AD (alt/(ref+alt)), NOT the
# GT-derived AF. Do NOT run `bcftools +fill-tags -t AF` here — that tag is computed from
# GT, so on a diploid-called pool it collapses every het to 0.5 and every hom to 1.0
# (the exact bug that produced a piled-at-0.5 "truth" and r=0.61 in job 20270457). Instead
# pass the AD-carrying VCF straight through with NO INFO/AF: both gen-reads (from_file) and
# `eidolon compare-af` then derive the fraction from FORMAT/AD, which is the real pool spectrum.
# stage_scn.sh emits AD (mpileup -a FORMAT/AD); VCFs staged before that fix lack it, so fall
# back to re-deriving AD from the BAM at just the SNV sites (fast, -R restricted). POOL_AF is
# BOTH the gen-reads input and the truth, so the comparison stays self-consistent either way.
# mpileup max per-site depth; its default of 250 would silently downsample a
# deep pool. Explicit so any cap is deliberate.
MPILEUP_MAX_DEPTH="${MPILEUP_MAX_DEPTH:-2000}"
POOL_AF="$OUTDIR/pool.af.vcf.gz"
if bcftools view -h "$POOL_VCF" | grep -q '##FORMAT=<ID=AD,'; then
    echo "pool VCF carries FORMAT/AD — using AD-based fraction directly"
    bcftools view -v snps "$POOL_VCF" -Oz -o "$POOL_AF"
elif [[ -s "$BAM" ]]; then
    echo "pool VCF lacks FORMAT/AD — re-deriving AD from BAM at the SNV sites"
    SITES="$OUTDIR/snv_sites.vcf.gz"
    bcftools view -v snps "$POOL_VCF" -Oz -o "$SITES"
    bcftools index -t "$SITES"
    # No `bcftools call` here (#450). `call -m` makes a diploid ML genotype call and a
    # hom-ref call discards the uncalled allele TOGETHER WITH ITS READ COUNT, so the
    # low-AF end of a pool spectrum — the part this harness exists to validate — was
    # being destroyed before the comparison saw it. Worse than the subclonal case,
    # which at least passed -C alleles. Measuring a fraction needs no genotype call;
    # `eidolon compare-af` reads the per-allele AD directly.
    bcftools mpileup -a FORMAT/AD -d "$MPILEUP_MAX_DEPTH" -f "$REF" -R "$SITES" \
        "$BAM" -Oz -o "$POOL_AF" 2>/dev/null
else
    echo "ABORT: $POOL_VCF has no FORMAT/AD and no BAM at $BAM to re-derive it." >&2
    echo "       Re-stage with the current stage_scn.sh (emits FORMAT/AD), or set INPUT_BAM." >&2
    exit 1
fi
bcftools index -t "$POOL_AF"
naf=$(bcftools view -H "$POOL_AF" | wc -l)
echo "pool AD-SNV sites: $naf"
[[ "$naf" -gt 0 ]] || { echo "ABORT: 0 SNV sites — check the pool VCF / BAM." >&2; exit 1; }
# Sanity: the AD-derived truth must NOT be piled at 0.5/1.0 (that = GT-derived, the bug above).
# Report the spread so a green run can't hide a collapsed spectrum.
echo "truth AF spread (from FORMAT/AD, should span deciles, not pile at 0.5):"
bcftools query -f '[%AD]\n' "$POOL_AF" \
  | awk -F',' '{r=$1;a=$2;t=r+a; if(t>0){printf "%.3f\n", a/t}}' \
  | awk '{b=int($1*10); if(b>9)b=9; c[b]++} END{for(i=0;i<10;i++)printf "  [%.1f,%.1f) %d\n", i/10,(i+1)/10,c[i]+0}'

# ── Step 2: gen-reads honoring the input AF (nothing de novo) ─────────────────
cat > "$OUTDIR/af.yml" <<YML
reference: $REF
read_len: 151
coverage: $COV
ploidy: 2
paired_ended: true
mutation_rate: 0.0
sv_rate_scale: 0.0
input_vcf: $POOL_AF
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $OUTDIR
output_filename: scn_af
rng_seed: scn af replay
num_threads: $THREADS
mutation_model: $MUT
YML
[[ -s "$MODELS/seq_error.json.gz" ]]   && echo "sequence_error_model: $MODELS/seq_error.json.gz" >> "$OUTDIR/af.yml"
[[ -s "$MODELS/frag_length.json.gz" ]] && echo "fragment_model: $MODELS/frag_length.json.gz"     >> "$OUTDIR/af.yml"
[[ -s "$MODELS/gc_bias.json.gz" ]]     && echo "gc_bias_model: $MODELS/gc_bias.json.gz"           >> "$OUTDIR/af.yml"

echo "=== gen-reads (honor input AF, mut_rate=0, cov=$COV) ==="
"$EIDOLON_BIN" --log-level info gen-reads -c "$OUTDIR/af.yml"

GOLDEN="$OUTDIR/scn_af.vcf.gz"
[[ -s "$GOLDEN" ]] || { echo "no golden VCF at $GOLDEN — did produce_vcf run?" >&2; exit 1; }
ngold=$(bcftools view -H "$GOLDEN" | wc -l)
echo "golden VCF records: $ngold"

# ── Step 3: compare per-site AF (real vs simulated) ──────────────────────────
echo
echo "════════════════════════════════════════════════════════════════"
echo "SCN Phase 2 AF reproduction — real pool vs simulated golden VCF"
echo "════════════════════════════════════════════════════════════════"
"$EIDOLON_BIN" compare-af \
    --truth "$POOL_AF" --sim "$GOLDEN" --min-depth "$MIN_DEPTH"
echo "════════════════════════════════════════════════════════════════"
echo "Inputs: pool AF=$POOL_AF  golden=$GOLDEN  (cov=$COV, min-depth=$MIN_DEPTH)"
echo "Interpretation: r>=0.95 + low per-decile MAE = the AF spectrum was reproduced."
