#!/bin/bash
# SLURM job: real-data validation of subclonal / somatic VAF reproduction (#405).
#
# The unit + H1N1 e2e tests prove gen-cancer-reads STAMPS the composed fraction
# (purity x dosage x CCF) and records it as INFO/NEAT_VAF. What they cannot prove —
# because it needs a real aligner at real coverage — is that the emitted READS
# actually pile up to that fraction once aligned. This job closes that loop:
#
#   architecture (inline subclones OR a real somatic VCF)
#     --gen-cancer-reads--> merged FASTQ + merged_truth.vcf.gz (carries NEAT_VAF)
#     --bwa-mem2--> merged BAM
#     --mpileup at the somatic sites--> OBSERVED VAF
#     compare OBSERVED VAF  vs  NEAT_VAF (the intended observed VAF)
#
# A high correlation + low MAE means the mixed-sample reads reproduce the planted
# VAF spectrum, i.e. NEAT_VAF is the ground truth a caller will measure. This is
# the generative composition (#405) and, with SOMATIC_VCF, the reproductive replay.
#
# WHY NEAT_VAF, NOT FORMAT/AF: FORMAT/AF in the truth is measured per-pass
# (tumor-only) = dosage x CCF; the observed sample VAF after tumor/normal mixing is
# purity x that = NEAT_VAF. We map NEAT_VAF -> INFO/AF on the truth side so the
# existing scn_af_compare.py (truth INFO/AF vs sim FORMAT/AD) does the comparison.
#
# PREREQS
#   * eidolon built from develop (>= #405 merged) via setup.sh -> $EIDOLON_BIN
#   * a REFERENCE FASTA (bwa-mem2-indexable) — the ONLY dataset you must stage
#   * a tumor model — defaults to the BUNDLED tools/cosmic_v104_pancancer_model.json.gz
#     (committed in the repo; no model_builders run needed). Override with TUMOR_MODEL=
#     or MODELS=<model_builders output>/models.
#   * conda `bioinf` env providing bwa-mem2 (as in cancer_pipeline.sbatch)
#   * ADJUST the `module load` / conda lines below to Delta's current names.
#
# USAGE
#   # generative (self-contained — just a reference + the bundled COSMIC model):
#   REFERENCE=$SCRATCH/neat_data/soy/soy.fa \
#     sbatch scripts/delta/run_subclonal_vaf_validation.sh
#   # reproductive (replay a real somatic VCF; INFO/AF or FORMAT/AD = observed VAF):
#   REFERENCE=... SOMATIC_VCF=$SCRATCH/data/tumor_somatic.vcf.gz PURITY=0.6 \
#     sbatch scripts/delta/run_subclonal_vaf_validation.sh

#SBATCH --job-name=eidolon-subvaf
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

# ── Inputs ───────────────────────────────────────────────────────────────────
REF="${REFERENCE:?set REFERENCE=<bwa-mem2-indexable FASTA>}"
# Tumor model shapes WHICH somatic variants are placed (not their VAF, which is what
# we validate), so the BUNDLED COSMIC model works out of the box — no model_builders
# staging needed. Override with a freshly-built one via TUMOR_MODEL=<mut_model.json.gz>
# or MODELS=<model_builders output>/models.
if [[ -n "${TUMOR_MODEL:-}" ]]; then
    :
elif [[ -n "${MODELS:-}" ]]; then
    TUMOR_MODEL="$MODELS/mut_model.json.gz"
else
    TUMOR_MODEL="$REPO_ROOT/tools/cosmic_v104_pancancer_model.json.gz"
fi
SOMATIC_VCF="${SOMATIC_VCF:-}"                 # set → reproductive replay; unset → generative
PURITY="${PURITY:-0.7}"
COV="${COV:-80}"                               # total (merged) depth; keep high for tight VAF
MIN_DEPTH="${MIN_DEPTH:-25}"                   # gate low-depth sites in the comparison
R_MIN="${R_MIN:-0.90}"                         # PASS gate: Pearson r floor
MAE_MAX="${MAE_MAX:-0.05}"                     # PASS gate: mean abs error ceiling
# Inline subclonal architecture used when SOMATIC_VCF is unset (clonal + 2 subclones):
SUBCLONES="${SUBCLONES:-1.0:0.5,0.5:0.3,0.2:0.2}"   # ccf:weight,ccf:weight,...
OUTDIR="${OUTDIR:-$SCRATCH/subvaf_${SLURM_JOB_ID:-manual}}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="${EIDOLON_BIN:-$CARGO_TARGET_DIR/release/eidolon}"

source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
module load bcftools/1.22 2>/dev/null || module load bcftools 2>/dev/null || true
conda_activate bioinf   # bwa-mem2 (not a module; see cancer_pipeline.sbatch)

[[ -s "$REF" ]]          || { echo "reference not found: $REF" >&2; exit 1; }
[[ -s "$TUMOR_MODEL" ]]  || { echo "tumor model not found: $TUMOR_MODEL (bundled default is tools/cosmic_v104_pancancer_model.json.gz)" >&2; exit 1; }
[[ -x "$EIDOLON_BIN" ]]  || { echo "eidolon not built: $EIDOLON_BIN (setup.sh on develop)" >&2; exit 1; }
command -v bwa-mem2 >/dev/null || { echo "bwa-mem2 not on PATH (conda bioinf?)" >&2; exit 1; }
[[ -z "$SOMATIC_VCF" || -s "$SOMATIC_VCF" ]] || { echo "SOMATIC_VCF set but missing: $SOMATIC_VCF" >&2; exit 1; }

mkdir -p "$OUTDIR"
if [[ -n "$SOMATIC_VCF" ]]; then MODE="reproductive (somatic_vcf)"; else MODE="generative (subclones=$SUBCLONES)"; fi
echo "=== banner: mode=$MODE  ref=$REF  tumor_model=$TUMOR_MODEL  purity=$PURITY  cov=$COV ==="

# ── Step 1: cancer YAML (subclonal architecture on the tumor pass) ───────────
YML="$OUTDIR/cancer.yml"
{
  echo "reference: $REF"
  echo "output_dir: $OUTDIR"
  echo "output_prefix: subvaf"
  echo "total_coverage: $COV"
  echo "purity: $PURITY"
  echo "read_len: 151"
  echo "paired_ended: true"
  echo "fragment_mean: 350"
  echo "fragment_st_dev: 50"
  echo "tumor_model: $TUMOR_MODEL"
  echo "overwrite_output: true"
  echo "rng_seed: subvaf-${SLURM_JOB_ID:-manual}"
  if [[ -n "$SOMATIC_VCF" ]]; then
    # Pure replay: no de-novo somatic, so every somatic site traces to the input VCF.
    echo "tumor_mutation_rate: 0.0"
    echo "somatic_vcf: $SOMATIC_VCF"
  else
    echo "tumor_mutation_rate: 0.00001"
    echo "subclones:"
    IFS=',' read -ra CLONES <<< "$SUBCLONES"
    for c in "${CLONES[@]}"; do
      echo "  - {ccf: ${c%%:*}, weight: ${c##*:}}"
    done
  fi
} > "$YML"
echo "--- cancer.yml ---"; cat "$YML"

echo "=== gen-cancer-reads ==="
"$EIDOLON_BIN" --log-level info gen-cancer-reads -c "$YML"

TRUTH="$OUTDIR/subvaf_merged_truth.vcf.gz"
R1="$OUTDIR/subvaf_merged_r1.fastq.gz"
R2="$OUTDIR/subvaf_merged_r2.fastq.gz"
for f in "$TRUTH" "$R1" "$R2"; do [[ -s "$f" ]] || { echo "missing output: $f" >&2; exit 1; }; done

# ── Step 2: somatic sites, with NEAT_VAF surfaced as INFO/AF (the truth) ──────
# scn_af_compare.py reads INFO/AF; NEAT_VAF is the intended OBSERVED VAF.
SITES="$OUTDIR/somatic_sites.vcf.gz"
{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=AF,Number=A,Type=Float,Description="intended observed VAF (from NEAT_VAF)">'
  echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
  bcftools view -H -i 'INFO/NEAT_ORIGIN="somatic"' "$TRUTH" \
    | awk -F'\t' 'match($8,/NEAT_VAF=[0-9.]+/){v=substr($8,RSTART+9,RLENGTH-9);
        print $1"\t"$2"\t.\t"$4"\t"$5"\t.\tPASS\tAF="v}'
} | bgzip > "$SITES"
bcftools index -t "$SITES"
nsom=$(bcftools view -H "$SITES" | wc -l)
echo "somatic sites carrying NEAT_VAF: $nsom"
[[ "$nsom" -gt 0 ]] || { echo "ABORT: 0 somatic NEAT_VAF sites — is this build >= #405?" >&2; exit 1; }
# Read-the-artifact guard: NEAT_VAF must span a spectrum, not pile at one value.
echo "NEAT_VAF spread (should span deciles for a subclonal architecture):"
bcftools query -f '%INFO/AF\n' "$SITES" \
  | awk '{b=int($1*10); if(b>9)b=9; c[b]++} END{for(i=0;i<10;i++)printf "  [%.1f,%.1f) %d\n",i/10,(i+1)/10,c[i]+0}'

# Tab-delimited allele list for forced genotyping (-C alleles): CHROM POS REF,ALT.
ALLELES="$OUTDIR/alleles.tsv.gz"
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$SITES" | bgzip > "$ALLELES"
tabix -s1 -b2 -e2 "$ALLELES"

# ── Step 3: align the MERGED (mixed) reads — the sample a caller sees ─────────
if [[ ! -s "${REF}.bwt.2bit.64" ]]; then
  echo "=== bwa-mem2 index (one-time) ==="; bwa-mem2 index "$REF" 2>&1 | tail -3
fi
MERGED_BAM="$OUTDIR/merged.bam"
echo "=== bwa-mem2 mem + sort (merged tumor+normal reads) ==="
bwa-mem2 mem -t "$THREADS" -R '@RG\tID:subvaf\tSM:merged\tPL:ILLUMINA' "$REF" "$R1" "$R2" \
    2>"$OUTDIR/bwa.log" \
  | samtools sort -@ "$THREADS" -o "$MERGED_BAM"
samtools index "$MERGED_BAM"

# ── Step 4: OBSERVED VAF at the somatic sites (forced genotyping) ─────────────
# -C alleles -T $ALLELES forces the KNOWN alt allele so AD is measured even for
# low-VAF subclonal sites that a de-novo caller's LoD would drop — exactly the
# sites #405 is about. FORMAT/AD then carries the observed alt fraction.
SIM="$OUTDIR/observed.vcf.gz"
echo "=== mpileup + forced-allele genotyping at the somatic sites ==="
bcftools mpileup -a FORMAT/AD -f "$REF" -R "$SITES" "$MERGED_BAM" -Ou 2>/dev/null \
  | bcftools call -m -C alleles -T "$ALLELES" -Oz -o "$SIM" 2>/dev/null
bcftools index -t "$SIM"
nsim=$(bcftools view -H "$SIM" | wc -l)
echo "sites genotyped in merged BAM: $nsim"
[[ "$nsim" -gt 0 ]] || { echo "ABORT: 0 sites genotyped — check alignment / -C alleles." >&2; exit 1; }

# ── Step 5: compare intended NEAT_VAF (truth) vs observed merged-BAM VAF ──────
echo
echo "════════════════════════════════════════════════════════════════"
echo "#405 subclonal VAF reproduction — intended NEAT_VAF vs observed merged VAF"
echo "════════════════════════════════════════════════════════════════"
CMP="$OUTDIR/compare.txt"
python3 "$REPO_ROOT/scripts/delta/scn_af_compare.py" \
    --truth "$SITES" --sim "$SIM" --min-depth "$MIN_DEPTH" | tee "$CMP"
echo "════════════════════════════════════════════════════════════════"

# ── PASS/FAIL gate (don't trust a green run — parse the actual numbers) ───────
# Parse the summary stat lines only: "n=.. Pearson r=.. Spearman rho=.." and
# "MAE=.. RMSE=..". Anchor with ^ so the per-decile "  [..) .. MAE=.." lines and the
# "  target: r>=.." hint line can't leak into the captured values.
r=$(awk -F'Pearson r=' '/^n=.*Pearson r=/{split($2,a," "); print a[1]}' "$CMP")
n=$(awk -F'n=' '/^n=.*Pearson r=/{split($2,a," "); print a[1]}' "$CMP")
mae=$(awk -F'MAE=' '/^MAE=/{split($2,a," "); print a[1]}' "$CMP")
echo "gate: n=$n r=$r (>=${R_MIN})  MAE=$mae (<=${MAE_MAX})"
verdict=FAIL
if [[ -n "$r" && -n "$mae" ]] \
   && awk -v r="$r" -v rm="$R_MIN" 'BEGIN{exit !(r>=rm)}' \
   && awk -v m="$mae" -v mm="$MAE_MAX" 'BEGIN{exit !(m<=mm)}'; then
  verdict=PASS
fi
echo "VERDICT: $verdict  (mode: $MODE)"
echo "  observed VAF should track NEAT_VAF = purity x dosage x CCF (reproductive: the input VAF)."

archive_run subvaf "$OUTDIR" "$CMP" "$YML" "$OUTDIR/bwa.log" || true
[[ "$verdict" == PASS ]]
