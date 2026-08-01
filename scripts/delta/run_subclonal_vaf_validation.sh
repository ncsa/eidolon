#!/bin/bash
# SLURM job: real-data validation of subclonal / somatic VAF reproduction (#405).
#
# The unit + H1N1 e2e tests prove gen-cancer-reads STAMPS the composed fraction
# (purity x dosage x CCF) and records it as INFO/EIDOLON_VAF. What they cannot prove —
# because it needs a real aligner at real coverage — is that the emitted READS
# actually pile up to that fraction once aligned. This job closes that loop:
#
#   architecture (inline subclones OR a real somatic VCF)
#     --gen-cancer-reads--> merged FASTQ + merged_truth.vcf.gz (carries EIDOLON_VAF)
#     --bwa-mem2--> merged BAM
#     --mpileup at the somatic sites--> OBSERVED VAF
#     compare OBSERVED VAF  vs  EIDOLON_VAF (the intended observed VAF)
#
# A high correlation + low MAE means the mixed-sample reads reproduce the planted
# VAF spectrum, i.e. EIDOLON_VAF is the ground truth a caller will measure. This is
# the generative composition (#405) and, with SOMATIC_VCF, the reproductive replay.
#
# WHY EIDOLON_VAF, NOT FORMAT/AF: FORMAT/AF in the truth is measured per-pass
# (tumor-only) = dosage x CCF; the observed sample VAF after tumor/normal mixing is
# purity x that = EIDOLON_VAF. We map EIDOLON_VAF -> INFO/AF on the truth side so the
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
#SBATCH --mem=64G
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
COV="${COV:-150}"                              # total (merged) depth; high so low-CCF sites
                                               # don't drop out and per-site VAF noise is tight
MIN_DEPTH="${MIN_DEPTH:-25}"                   # gate low-depth sites in the comparison
# mpileup's max per-site depth. Its default of 250 would silently downsample a
# 200x+ run; set explicitly so any cap is deliberate and visible.
MPILEUP_MAX_DEPTH="${MPILEUP_MAX_DEPTH:-2000}"
# PASS is gated on FIDELITY metrics that hold for a discrete subclonal architecture:
#   |mean(observed-intended)| = bias (systematic composition error) and MAE (per-site
#   accuracy vs the binomial noise floor). Pearson r is ADVISORY only — for tight,
#   discrete EIDOLON_VAF clusters at moderate depth, per-site noise depresses r even when
#   the reproduction is unbiased and accurate (it fits a CONTINUOUS spectrum like #398,
#   not clusters). Raise COV and/or widen the CCF architecture to make r meaningful.
BIAS_MAX="${BIAS_MAX:-0.02}"                   # PASS gate: |mean(observed-intended)|
MAE_MAX="${MAE_MAX:-0.05}"                     # PASS gate: mean abs error ceiling
R_MIN="${R_MIN:-0.90}"                         # advisory only (warned, not gated)
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
# bwa-mem2 comes from the `bioinf` conda env OR a module. A batch shell needs the
# conda profile sourced before `conda activate`, else you get the (non-fatal here)
# "Run 'conda init'" error. Tolerate failure — we verify bwa-mem2 on PATH below.
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
fi
conda_activate bioinf 2>/dev/null || true
module load bwa-mem2 2>/dev/null || true   # in case it's a module, not conda

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
  # Only the MERGED reads are aligned; drop the per-pass FASTQ copies to ~third the
  # FASTQ footprint (matters at genome scale, where disk is the limiter).
  echo "keep_per_pass: false"
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

# ── Step 2: somatic sites, with EIDOLON_VAF surfaced as INFO/AF (the truth) ──────
# scn_af_compare.py reads INFO/AF; EIDOLON_VAF is the intended OBSERVED VAF.
SITES="$OUTDIR/somatic_sites.vcf.gz"
{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=AF,Number=A,Type=Float,Description="intended observed VAF (from EIDOLON_VAF)">'
  echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
  bcftools view -H -i 'INFO/EIDOLON_ORIGIN="somatic"' "$TRUTH" \
    | awk -F'\t' 'match($8,/EIDOLON_VAF=[0-9.]+/){v=substr($8,RSTART,RLENGTH); sub(/^EIDOLON_VAF=/,"",v);
        print $1"\t"$2"\t.\t"$4"\t"$5"\t.\tPASS\tAF="v}'
} | bgzip > "$SITES"
bcftools index -t "$SITES"
nsom=$(bcftools view -H "$SITES" | wc -l)
echo "somatic sites carrying EIDOLON_VAF: $nsom"
[[ "$nsom" -gt 0 ]] || { echo "ABORT: 0 somatic EIDOLON_VAF sites — is this build >= #405?" >&2; exit 1; }
# The record count above is necessary but not sufficient, and this is the exact site
# where that bit: it passed while every value was the malformed string `AF=AF=0.3000`,
# because counting records says nothing about their content. bcftools would not have
# complained either — it silently converts a type-mismatched value to `.`.
validate_artifacts "$EIDOLON_BIN" "somatic sites" "$SITES" || exit 1
# Read-the-artifact guard: EIDOLON_VAF must span a spectrum, not pile at one value.
echo "EIDOLON_VAF spread (should span deciles for a subclonal architecture):"
bcftools query -f '%INFO/AF\n' "$SITES" \
  | awk '{b=int($1*10); if(b>9)b=9; c[b]++} END{for(i=0;i<10;i++)printf "  [%.1f,%.1f) %d\n",i/10,(i+1)/10,c[i]+0}'


# ── Step 3: align the MERGED (mixed) reads — the sample a caller sees ─────────
if [[ ! -s "${REF}.bwt.2bit.64" ]]; then
  echo "=== bwa-mem2 index (one-time) ==="; bwa-mem2 index "$REF" 2>&1 | tail -3
fi
MERGED_BAM="$OUTDIR/merged.bam"
SORT_TMP="$OUTDIR/sort_tmp"; mkdir -p "$SORT_TMP"
echo "=== bwa-mem2 mem + sort (merged tumor+normal reads) ==="
echo "  scratch free before align: $(df -h "$OUTDIR" | awk 'NR==2{print $4" free ("$5" used)"}')"
# Explicit -T keeps sort spill files on this scratch dir (a node-local $TMPDIR can be
# too small or vanish). -m 1G caps per-thread sort RAM: with -@ N that is N GB, which
# must sit alongside bwa-mem2's index+buffers under --mem (a chr-scale reference's
# index is several GB — 2G/thread OOM-killed a chr1 job at --mem=48G). pipefail turns a
# bwa-mem2 crash into a hard failure here rather than a downstream "0 sites". A large
# reference at high COV can also exhaust scratch — that surfaces as a samtools write
# error, so report space + the bwa log on failure.
if ! bwa-mem2 mem -t "$THREADS" -R '@RG\tID:subvaf\tSM:merged\tPL:ILLUMINA' "$REF" "$R1" "$R2" \
      2>"$OUTDIR/bwa.log" \
    | samtools sort -@ "$THREADS" -m 1G -T "$SORT_TMP/st" -o "$MERGED_BAM"; then
  echo "ALIGN/SORT FAILED. scratch: $(df -h "$OUTDIR" | awk 'NR==2{print $4" free, "$5" used"}')" >&2
  echo "--- last 20 lines of bwa.log ---" >&2; tail -20 "$OUTDIR/bwa.log" >&2
  echo "  If this is a disk/quota issue, re-run with a single-scaffold REFERENCE and/or lower COV." >&2
  exit 1
fi
samtools index "$MERGED_BAM"

# ── Step 4: OBSERVED VAF at the somatic sites (mpileup AD, no genotype call) ──
# Read FORMAT/AD straight from mpileup. NO `bcftools call` (#450).
#
# The previous version piped mpileup into `bcftools call -m -C alleles`, on the
# reasoning that -C alleles "forces the KNOWN alt allele so AD is measured even for
# low-VAF subclonal sites that a de-novo caller's LoD would drop". -C alleles does
# constrain which alleles are considered, but `call -m` still makes a diploid ML
# GENOTYPE call — and a hom-ref call discards the uncalled allele together with its
# read count. Measured in isolation at 7% VAF:
#
#   mpileup alone          REF=G ALT=T,<*>   AD=93,7,0     <- 7/100 = 0.07, correct
#   mpileup | call -m      REF=G ALT=.       AD=93         <- allele and count gone
#
# Coverage cannot rescue it: identical at 100x, 300x and 600x, because the decision
# is driven by the allele FRACTION against a diploid model — no diploid genotype
# predicts 7%. So the step meant to preserve low-VAF sites was destroying exactly
# them, and 160 of 567 planted sites (the entire lowest CCF cluster) never reached
# the comparison. Measuring an allele fraction needs no genotype call at all.
#
# scn_af_compare.py selects this allele's own AD element from the ALT list, so
# mpileup's `<*>` non-ref placeholder beside the real base is handled.
SIM="$OUTDIR/observed.vcf.gz"
echo "=== mpileup at the somatic sites (AD only, no genotype calling) ==="
# -d: mpileup's default max depth is 250, which would silently downsample a 200x+
# run. Set it explicitly so the cap is visible and generous.
bcftools mpileup -a FORMAT/AD -d "$MPILEUP_MAX_DEPTH" -f "$REF" -R "$SITES" \
    "$MERGED_BAM" -Oz -o "$SIM" 2>/dev/null
bcftools index -t "$SIM"
nsim=$(bcftools view -H "$SIM" | wc -l)
echo "sites piled up in merged BAM: $nsim"
[[ "$nsim" -gt 0 ]] || { echo "ABORT: 0 sites piled up — check alignment / -R sites." >&2; exit 1; }
# Read-the-artifact guard: a site with coverage must carry a per-allele AD. If AD is
# absent the comparison would silently find nothing to join.
nad=$(bcftools query -f '%CHROM\t%POS\t[%AD]\n' "$SIM" 2>/dev/null | awk -F'\t' '$3!="" && $3!="."' | wc -l)
echo "sites carrying FORMAT/AD: $nad"
[[ "$nad" -gt 0 ]] || { echo "ABORT: no site carries FORMAT/AD — mpileup -a FORMAT/AD failed." >&2; exit 1; }

# ── Step 5: compare intended EIDOLON_VAF (truth) vs observed merged-BAM VAF ──────
echo
echo "════════════════════════════════════════════════════════════════"
echo "#405 subclonal VAF reproduction — intended EIDOLON_VAF vs observed merged VAF"
echo "════════════════════════════════════════════════════════════════"
CMP="$OUTDIR/compare.txt"
# scn_af_compare.py exits non-zero when too much of the planted set went unscored.
# Capture that instead of letting `set -e` + pipefail abort here, which would skip both
# the verdict and archive_run — the artifacts are most worth keeping when it fails.
# `$(cmd; echo $?)` cannot be used for this: the inherited `set -e` kills the subshell
# before echo runs. `|| rc=$?` is the form that works.
CMP_RC=0
python3 "$REPO_ROOT/scripts/delta/scn_af_compare.py" \
    --truth "$SITES" --sim "$SIM" --min-depth "$MIN_DEPTH" | tee "$CMP" || CMP_RC=$?
echo "════════════════════════════════════════════════════════════════"

# ── PASS/FAIL gate (don't trust a green run — parse the actual numbers) ───────
# Parse only the summary stat lines: "n=.. Pearson r=.. Spearman rho=.." and
# "MAE=.. RMSE=.. mean(sim-truth)=..". Anchor with ^ so the per-decile "  [..) MAE=.."
# and "  target: r>=.." lines can't leak into the captured values.
r=$(awk -F'Pearson r=' '/^n=.*Pearson r=/{split($2,a," "); print a[1]}' "$CMP")
n=$(awk -F'n=' '/^n=.*Pearson r=/{split($2,a," "); print a[1]}' "$CMP")
mae=$(awk -F'MAE=' '/^MAE=/{split($2,a," "); print a[1]}' "$CMP")
bias=$(awk -F'mean\\(sim-truth\\)=' '/^MAE=.*mean\(sim-truth\)=/{split($2,a," "); print a[1]}' "$CMP")
absbias=$(awk -v b="$bias" 'BEGIN{b=b+0; print (b<0)?-b:b}')
echo "gate: n=$n  |bias|=$absbias (<=${BIAS_MAX})  MAE=$mae (<=${MAE_MAX})  [advisory r=$r vs ${R_MIN}]"
# Advisory: a low r on a discrete architecture is expected, not a failure — flag it.
if [[ -n "$r" ]] && awk -v r="$r" -v rm="$R_MIN" 'BEGIN{exit !(r<rm)}'; then
  echo "  NOTE: Pearson r=$r < ${R_MIN} — expected for tight discrete EIDOLON_VAF clusters at" >&2
  echo "        this depth; raise COV and/or widen SUBCLONES for a Pearson-meaningful run." >&2
fi
verdict=FAIL
# Coverage first: bias and MAE computed over a subset that dropped a whole VAF stratum
# look BETTER than an honest full-coverage run, so they can never be the only gate.
# That is #450 exactly — it reported clean numbers and PASS while excluding the sites
# it existed to test.
if [[ "$CMP_RC" -ne 0 ]]; then
  echo "  coverage gate FAILED (scn_af_compare.py rc=$CMP_RC) — see the per-decile" >&2
  echo "  planted/scored table above. bias and MAE are not usable." >&2
elif [[ -n "$bias" && -n "$mae" ]] \
   && awk -v b="$absbias" -v bm="$BIAS_MAX" 'BEGIN{exit !(b<=bm)}' \
   && awk -v m="$mae" -v mm="$MAE_MAX" 'BEGIN{exit !(m<=mm)}'; then
  verdict=PASS
fi
echo "VERDICT: $verdict  (mode: $MODE)"
echo "  Fidelity: unbiased (|bias|<=$BIAS_MAX) + accurate (MAE<=$MAE_MAX vs the binomial"
echo "  noise floor) means the merged reads reproduce EIDOLON_VAF = purity x dosage x CCF"
echo "  (reproductive: the input VAF). r is a spread-dependent summary, not the gate."

archive_run subvaf "$OUTDIR" "$CMP" "$YML" "$OUTDIR/bwa.log" || true
[[ "$verdict" == PASS ]]
