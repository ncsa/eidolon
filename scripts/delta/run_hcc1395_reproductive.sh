#!/bin/bash
# SLURM job: REPRODUCTIVE somatic-VAF validation on REAL cancer data (SEQC2 HCC1395, #405).
#
# The soy validation proved eidolon reproduces a SYNTHETIC subclonal VAF spectrum. This
# proves it on a REAL tumor's EMPIRICAL spectrum: take HCC1395's high-confidence somatic
# SNVs with their observed VAF, replay them through gen-cancer-reads, align the merged
# reads, and confirm the observed VAF tracks NEAT_VAF (= the input VAF).
#
#   HCC1395 somatic SNVs + observed VAF  --somatic_vcf-->  gen-cancer-reads (reproductive)
#     --bwa-mem2--> merged BAM --> observed VAF  ≈  NEAT_VAF (= the real tumor VAF)
#
# Unlike the synthetic run, the input is a real tumor's full, continuous VAF distribution
# (real subclonal structure), so Pearson r is meaningful here too — this is the reproductive
# analogue of the #398 pooled-AF validation, on cancer data.
#
# WHAT THIS DOES (one allocation)
#   1. Carve a single-contig GRCh38 reference (fast index/align).
#   2. Build the somatic-VAF input: HCC1395 SNVs on that contig, each carrying an OBSERVED
#      per-variant VAF — from the truth's own INFO/AF if present, else derived from the REAL
#      tumor BAM by forced-allele genotyping (no low-VAF dropout).
#   3. Hand off to run_subclonal_vaf_validation.sh (SOMATIC_VCF mode) via `bash` (its #SBATCH
#      lines are inert when sourced this way), so the one validated engine does the sim +
#      align + compare + fidelity gate.
#
# PREREQS
#   * stage_hcc1395.sh already run: $HCC1395_DIR/somatic.vcf.gz + $HCC1395_DIR/tumor.bam
#     (REGION must include the contig used here; default chr1).
#   * GRCh38.fa (chr-prefixed) staged; eidolon built (setup.sh); bwa-mem2 (conda/module).
#   * ADJUST module/conda lines to Delta's current names (as in the engine script).
#
# USAGE
#   sbatch scripts/delta/run_hcc1395_reproductive.sh
#   REGION=chr1 PURITY=0.9 COV=200 sbatch scripts/delta/run_hcc1395_reproductive.sh

#SBATCH --job-name=eidolon-hcc1395repro
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

D="${HCC1395_DIR:-$SCRATCH/neat_data/hcc1395}"
REF_FULL="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"   # full GRCh38 the tumor BAM was aligned to
CTG="${REGION:-chr1}"                                   # single contig; must match stage_hcc1395 + chr-prefix
PURITY="${PURITY:-0.9}"                                 # HCC1395 ~pure; high purity minimises VAF>purity clamping
COV="${COV:-200}"
OUTDIR="${OUTDIR:-$SCRATCH/hcc1395_repro_${SLURM_JOB_ID:-manual}}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"

source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
module load bcftools/1.22 2>/dev/null || module load bcftools 2>/dev/null || true

SOMATIC="$D/somatic.vcf.gz"
TUMOR_BAM="$D/tumor.bam"
[[ -s "$REF_FULL" ]]  || { echo "GRCh38 not staged: $REF_FULL" >&2; exit 1; }
[[ -s "$SOMATIC" ]]   || { echo "HCC1395 somatic truth not found: $SOMATIC (run stage_hcc1395.sh)" >&2; exit 1; }
[[ -s "$TUMOR_BAM" ]] || { echo "HCC1395 tumor BAM not found: $TUMOR_BAM (run stage_hcc1395.sh)" >&2; exit 1; }
mkdir -p "$OUTDIR"
echo "=== banner: HCC1395 reproductive  contig=$CTG  purity=$PURITY  cov=$COV ==="

# ── 1. single-contig reference ───────────────────────────────────────────────
CREF="$OUTDIR/${CTG}.fa"
if [[ ! -s "$CREF" ]]; then
  echo "=== carving $CTG from $REF_FULL ==="
  samtools faidx "$REF_FULL" "$CTG" > "$CREF"
fi
samtools faidx "$CREF"
[[ -s "$CREF.fai" ]] || { echo "failed to index $CREF (bad contig name '$CTG'? check chr-prefix)" >&2; exit 1; }

# ── 2. somatic SNVs on this contig, each with an OBSERVED per-variant VAF ─────
# Prefer the truth's own AF/AD; otherwise derive the observed VAF from the REAL tumor
# BAM (forced-allele genotyping so low-VAF somatic sites are not dropped). Either way
# eidolon's from_file reads the fraction (INFO/AF or FORMAT/AD) into allele_fraction.
SITES="$OUTDIR/hcc1395_${CTG}_snv.vcf.gz"
bcftools view -r "$CTG" -v snps "$SOMATIC" -Oz -o "$SITES"
bcftools index -t "$SITES"
nsnv=$(bcftools view -H "$SITES" | wc -l)
echo "HCC1395 somatic SNVs on $CTG: $nsnv"
[[ "$nsnv" -gt 0 ]] || { echo "ABORT: 0 somatic SNVs on $CTG — was stage_hcc1395 run for this contig?" >&2; exit 1; }

SOM_AF="$OUTDIR/hcc1395_${CTG}_somatic_af.vcf.gz"
if bcftools view -h "$SITES" | grep -q '##INFO=<ID=AF,' || bcftools view -h "$SITES" | grep -q '##FORMAT=<ID=AD,'; then
  echo "somatic truth carries AF/AD — using its observed VAF directly"
  cp "$SITES" "$SOM_AF"
else
  echo "somatic truth lacks AF/AD — deriving observed VAF from the real tumor BAM ($CTG)"
  ALLELES="$OUTDIR/alleles.tsv.gz"
  bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$SITES" | bgzip > "$ALLELES"
  tabix -s1 -b2 -e2 "$ALLELES"
  # -C alleles -T forces the known alt allele so AD is measured at every site (incl. low
  # VAF). mpileup uses the FULL ref the BAM was aligned to. pipefail-safe: no head.
  bcftools mpileup -a FORMAT/AD -f "$REF_FULL" -R "$SITES" "$TUMOR_BAM" -Ou 2>/dev/null \
    | bcftools call -m -C alleles -T "$ALLELES" -Oz -o "$SOM_AF" 2>/dev/null
fi
bcftools index -t "$SOM_AF"
naf=$(bcftools view -H "$SOM_AF" | wc -l)
echo "somatic SNVs with an observed VAF: $naf"
[[ "$naf" -gt 0 ]] || { echo "ABORT: 0 SNVs with a usable VAF — check the tumor BAM coverage on $CTG." >&2; exit 1; }
# Read-the-artifact: show the real tumor's VAF spectrum we're about to reproduce.
echo "HCC1395 observed VAF spread (real tumor spectrum, from AF/AD):"
bcftools query -f '[%AD]\n' "$SOM_AF" 2>/dev/null \
  | awk -F',' '{r=$1;a=$2;t=r+a; if(t>0){v=a/t; b=int(v*10); if(b>9)b=9; c[b]++}}
      END{for(i=0;i<10;i++)printf "  [%.1f,%.1f) %d\n",i/10,(i+1)/10,c[i]+0}' \
  || bcftools query -f '%INFO/AF\n' "$SOM_AF" \
     | awk '{b=int($1*10); if(b>9)b=9; c[b]++} END{for(i=0;i<10;i++)printf "  [%.1f,%.1f) %d\n",i/10,(i+1)/10,c[i]+0}'

# ── 3. hand off to the validated engine (reproductive mode) ──────────────────
# `bash` (not sbatch) so its #SBATCH lines are inert and it runs in THIS allocation.
echo
echo "=== handing off to run_subclonal_vaf_validation.sh (SOMATIC_VCF mode) ==="
REFERENCE="$CREF" \
SOMATIC_VCF="$SOM_AF" \
PURITY="$PURITY" \
COV="$COV" \
OUTDIR="$OUTDIR/validate" \
CARGO_TARGET_DIR="$CARGO_TARGET_DIR" \
EIDOLON_REPO="$REPO_ROOT" \
  bash "$REPO_ROOT/scripts/delta/run_subclonal_vaf_validation.sh"
