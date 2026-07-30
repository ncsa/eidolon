#!/bin/bash
# SLURM job: structural-variant (SV) validation pipeline on Delta.
#
# WHY THIS EXISTS: gen-cancer-reads emits structural variants (DEL/DUP/INV/INS/BND
# and CNVs) into the tumor truth VCF, but every validation so far ran through
# Mutect2 + som.py, which score SNVs/indels ONLY. The SV machinery has therefore
# never been checked against a real SV caller. This pipeline closes that gap: it
# simulates an SV-rich tumor/normal pair, calls somatic SVs with Manta, and scores
# recall/precision against eidolon's own SV truth with truvari — the same
# "align + call + compare" loop that surfaced the strand bugs (#285/#287) for SNVs.
#
# Stages:
#   1. Simulate tumor/normal FASTQ + SV-tagged truth VCF (cancer_simulate.sh,
#      SV_RATE_SCALE > 0).
#   2-3. Align normal + merged-tumor to reference with BWA-MEM2; sort + index.
#   3b. PRE-FLIGHT: inspect eidolon's SV truth records (SVTYPE/SVLEN/END/symbolic
#      ALT) and extract a somatic-SV-only truth set; abort early if it's empty or
#      malformed so we don't waste a Manta run. (First-run sanity for VCF
#      compatibility — adjust the SV filter below to match the schema it prints.)
#   4. Call somatic SVs with Manta (tumor/normal mode).
#   5. Score Manta calls vs the SV truth with truvari (overall + per SV type).
#
# Each stage is idempotent — existing outputs are skipped, so a re-run resumes.
#
# PREREQUISITES (one-time; Manta needs python2, truvari python3 — separate envs):
#   bash scripts/delta/setup.sh              # reference build + bioinf env
#   conda create -y -n manta   -c bioconda -c conda-forge manta
#   conda create -y -n truvari -c bioconda -c conda-forge truvari
#   (override env names with MANTA_ENV / TRUVARI_ENV if you named them differently)
#
# Usage:
#   sbatch scripts/delta/sv_pipeline.sbatch
#   REFERENCE=/path/ref.fa SV_RATE_SCALE=2.0 TUMOR_MODEL=/path/model.json.gz \
#     sbatch scripts/delta/sv_pipeline.sbatch
#   # per-tissue SV spectrum check (BRCA DUP-dominant / skin BND / lung DEL):
#   TUMOR_MODEL=$REPO_ROOT/tools/cosmic_pancancer_sv_BRCA.json.gz sbatch ...

#SBATCH --job-name=eidolon-sv
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Parameters ───────────────────────────────────────────────────────────
REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/sv_$SLURM_JOB_ID}"

TOTAL_COVERAGE="${TOTAL_COVERAGE:-30}"
PURITY="${PURITY:-0.6}"
READ_LEN="${READ_LEN:-151}"
PAIRED="${PAIRED:-true}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
# SVs are the whole point — must be > 0. Scale up to enrich the truth set.
SV_RATE_SCALE="${SV_RATE_SCALE:-1.0}"

# Tumor model: pan-cancer SV model by default; swap a per-tissue SV model
# (cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz) to check tissue-specific spectra.
TUMOR_MODEL="${TUMOR_MODEL:-$REPO_ROOT/tools/cosmic_v104_pancancer_model.json.gz}"

MANTA_ENV="${MANTA_ENV:-manta}"
DELLY_ENV="${DELLY_ENV:-delly}"
TRUVARI_ENV="${TRUVARI_ENV:-truvari}"
# Run Delly as a 2nd somatic SV caller alongside Manta (cross-validation; the
# BND-recovery question). Set RUN_DELLY=0 to skip.
RUN_DELLY="${RUN_DELLY:-1}"
# Run GATK somatic CNV to validate the CNV machinery through a real copy-number
# caller (eidolon's CNVs were previously only depth-validated). No new env (gatk4 in
# bioinf). Set RUN_CNV=0 to skip.
RUN_CNV="${RUN_CNV:-1}"
# eidolon SVs reach ~1 Mb; truvari's default --sizemax is 50 kb, which would silently
# EXCLUDE most large eidolon SVs from scoring. Raise it to cover the full range.
TRUVARI_SIZEMAX="${TRUVARI_SIZEMAX:-5000000}"

THREADS=$SLURM_CPUS_PER_TASK
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"
PREFIX="$OUTDIR/sample"

if awk "BEGIN{exit !($SV_RATE_SCALE <= 0)}"; then
    echo "ERROR: SV_RATE_SCALE must be > 0 for the SV pipeline (got $SV_RATE_SCALE)." >&2
    exit 1
fi

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda

mkdir -p "$OUTDIR"
REF_NAME="$(basename "$REFERENCE")"

echo "=== eidolon SV validation pipeline: $SLURM_JOB_ID ==="
echo "Reference:  $REFERENCE"
echo "Output:     $OUTDIR"
echo "Coverage:   ${TOTAL_COVERAGE}x  Purity: $PURITY  SV scale: $SV_RATE_SCALE"
echo "Tumor model:$TUMOR_MODEL"
echo "Threads:    $THREADS"

# ── Stage 1: Simulate (SV-rich tumor/normal) ─────────────────────────────
conda_activate bioinf   # bwa-mem2, bcftools (samtools via module)
MERGED_R1="${PREFIX}_merged_r1.fastq.gz"
NORMAL_R1="${PREFIX}_normal_r1.fastq.gz"
TRUTH_VCF="${PREFIX}_merged_truth.vcf.gz"

if [[ -f "$MERGED_R1" && -f "$NORMAL_R1" && -f "$TRUTH_VCF" ]]; then
    echo "[1/5] Simulation outputs already present — skipping."
else
    echo "[1/5] Simulating SV-rich tumor/normal reads..."
    PAIRED_FLAG=""
    [[ "$PAIRED" == "true" ]] && PAIRED_FLAG="--paired-ended --fragment-mean $FRAG_MEAN --fragment-st-dev $FRAG_SD"
    RNG_SEED_ROOT="${RNG_SEED_ROOT:-sv-validate}${SLURM_ARRAY_TASK_ID:+-rep$SLURM_ARRAY_TASK_ID}"
    bash "$REPO_ROOT/tools/cancer_simulate.sh" \
        --reference "$REFERENCE" \
        --output-prefix "$PREFIX" \
        --total-coverage "$TOTAL_COVERAGE" \
        --purity "$PURITY" \
        --read-len "$READ_LEN" \
        --tumor-model "$TUMOR_MODEL" \
        --sv-rate-scale "$SV_RATE_SCALE" \
        --rng-seed "$RNG_SEED_ROOT" \
        --eidolon-bin "$EIDOLON_BIN" \
        $PAIRED_FLAG
    echo "[1/5] Simulation complete."
fi

# ── Stage 2 & 3: Index reference, align ──────────────────────────────────
NORMAL_BAM="$OUTDIR/normal.bam"
TUMOR_BAM="$OUTDIR/tumor.bam"

index_and_align() {
    local in_fq="$1" out_bam="$2" sample="$3"
    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then
        echo "  $sample BAM already present — skipping."; return 0
    fi
    # Guard on -s, not -f: a build killed mid-write leaves a 0-byte .bwt.2bit.64
    # that -f treats as present, poisoning every later run (bwa-mem2: "Unexpected
    # end of file"). -s rebuilds on empty/missing.
    if [[ ! -s "${REFERENCE}.bwt.2bit.64" ]] && [[ ! -s "${REFERENCE}.bwt" ]]; then
        echo "  Indexing reference for BWA-MEM2 (one-time)..."
        bwa-mem2 index "$REFERENCE" 2>&1 | tail -3
    fi
    [[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
    echo "  Aligning $sample..."
    local rg="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    if [[ "$PAIRED" == "true" ]]; then
        local in_r2="${in_fq/_r1./_r2.}"
        bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$in_fq" "$in_r2" 2>"$OUTDIR/${sample}_bwa.log" \
            | samtools sort -@ "$THREADS" -o "$out_bam"
    else
        bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$in_fq" 2>"$OUTDIR/${sample}_bwa.log" \
            | samtools sort -@ "$THREADS" -o "$out_bam"
    fi
    samtools index "$out_bam"
    echo "  $sample BAM ready: $out_bam"
}

echo "[2/5] Aligning normal reads..."
index_and_align "$NORMAL_R1" "$NORMAL_BAM" "NORMAL"
echo "[3/5] Aligning tumor (merged) reads..."
index_and_align "${PREFIX}_merged_r1.fastq.gz" "$TUMOR_BAM" "TUMOR"

# ── Stage 3b: PRE-FLIGHT — inspect + extract the somatic SV truth ────────
# eidolon's SVs are symbolic records (ALT like <DEL>/<DUP>/<INV>, INFO SVTYPE/
# SVLEN/END). Print a sample so we can confirm the schema, then extract the
# somatic SV subset truvari will score. EIDOLON_ORIGIN=somatic drops germline
# carry-through; the symbolic-ALT test keeps only SVs (not SNVs/indels).
# Normalize pre-v2.1.0 NEAT_* tokens if this truth VCF predates the rename (e.g.
# re-scoring an archived artifact). Idempotent and a no-op on a current file, so
# the unconditional call is safe. `|| rc=$?` (not `; echo $?`) so the non-zero
# --check status is *tested* — under set -e an untested non-zero would abort here.
migrate_rc=0
bash "$REPO_ROOT/tools/migrate_legacy_tokens.sh" --check "$TRUTH_VCF" >/dev/null 2>&1 \
    || migrate_rc=$?
if [[ "$migrate_rc" -eq 10 ]]; then
    MIGRATED_TRUTH="$OUTDIR/truth_migrated.vcf.gz"
    bash "$REPO_ROOT/tools/migrate_legacy_tokens.sh" "$TRUTH_VCF" "$MIGRATED_TRUTH"
    TRUTH_VCF="$MIGRATED_TRUTH"
fi

TRUTH_SV="$OUTDIR/truth_sv.vcf.gz"
echo "[3b/5] Inspecting eidolon SV truth..."
echo "  --- SV INFO header lines ---"
bcftools view -h "$TRUTH_VCF" | grep -iE 'INFO=<ID=(SVTYPE|SVLEN|END|EIDOLON_ORIGIN)' || \
    echo "  (no SVTYPE/SVLEN/END/EIDOLON_ORIGIN INFO headers found — schema may differ!)"
echo "  --- first symbolic-ALT (SV) records ---"
bcftools view -H "$TRUTH_VCF" | awk '$5 ~ /^</' | head -5 || true
echo "  --- SVTYPE counts in truth (all origins) ---"
bcftools view -H "$TRUTH_VCF" | grep -oE 'SVTYPE=[A-Z]+' | sort | uniq -c || true

if [[ -f "$TRUTH_SV" ]]; then
    echo "  somatic SV truth already extracted — skipping."
else
    # Permissive SV selector: symbolic ALT OR an SVTYPE INFO tag, AND somatic.
    # If the pre-flight above shows a different schema, adjust this -i expression.
    bcftools view -i '(ALT[0]~"<" || INFO/SVTYPE!=".") && INFO/EIDOLON_ORIGIN="somatic"' \
        -O z -o "$TRUTH_SV" "$TRUTH_VCF"
    bcftools sort -O z -o "$TRUTH_SV" "$TRUTH_SV" 2>/dev/null || true
    bcftools index -f -t "$TRUTH_SV"
fi
SV_TRUTH_N=$(bcftools view -H "$TRUTH_SV" 2>/dev/null | wc -l)
echo "  somatic SV truth records: $SV_TRUTH_N"
if [[ "$SV_TRUTH_N" -eq 0 ]]; then
    echo "ERROR: no somatic SV records extracted. Check the pre-flight schema dump" >&2
    echo "       above and adjust the SV selector. Raising SV_RATE_SCALE also helps." >&2
    exit 1
fi

# Pre-generate per-type truth subsets NOW, while bcftools (bioinf) is active, so
# the truvari scoring stage needs only truvari on its PATH (no bcftools in that
# env). Empty types are removed so the scoring loop skips them.
for svt in DEL DUP INV INS BND CNV; do
    typed="$OUTDIR/truth_sv_${svt}.vcf.gz"
    bcftools view -i "INFO/SVTYPE=\"$svt\"" -O z -o "$typed" "$TRUTH_SV" 2>/dev/null || continue
    bcftools index -f -t "$typed" 2>/dev/null || true
    n=$(bcftools view -H "$typed" 2>/dev/null | wc -l)
    echo "  truth $svt: $n"
    [[ "$n" -eq 0 ]] && rm -f "$typed" "${typed}.tbi"
done

# ── Stage 4: Manta somatic SV calling ────────────────────────────────────
MANTA_DIR="$OUTDIR/manta"
MANTA_SV="$MANTA_DIR/results/variants/somaticSV.vcf.gz"
if [[ -f "$MANTA_SV" ]]; then
    echo "[4/5] Manta somaticSV already present — skipping."
else
    echo "[4/5] Running Manta (tumor/normal somatic SV)..."
    conda_activate "$MANTA_ENV"
    rm -rf "$MANTA_DIR"
    # Manta wants a bgzipped+faidx'd reference and indexed BAMs (we have both).
    configManta.py \
        --normalBam "$NORMAL_BAM" \
        --tumorBam  "$TUMOR_BAM" \
        --referenceFasta "$REFERENCE" \
        --runDir "$MANTA_DIR" 2>&1 | tail -5
    "$MANTA_DIR/runWorkflow.py" -j "$THREADS" 2>&1 | tail -10
    [[ -f "$MANTA_SV" ]] || { echo "ERROR: Manta produced no somaticSV.vcf.gz" >&2; exit 1; }
    # Count without bcftools — the manta env doesn't have it; zcat/grep are env-independent.
    echo "[4/5] Manta complete: $(zcat "$MANTA_SV" | grep -vc '^#' || true) somatic SV calls."
fi

# Diagnostic: compare PRE-filter candidates vs POST-filter somatic calls. This
# separates "reads carry no SV signal" (simulator/coverage issue) from "signal
# present but somatic-classified away" (support/purity issue) when recall is low.
CAND_SV="$MANTA_DIR/results/variants/candidateSV.vcf.gz"
SOM_N=$(zcat "$MANTA_SV" | grep -vc '^#' || true)
CAND_N=0; [[ -f "$CAND_SV" ]] && CAND_N=$(zcat "$CAND_SV" | grep -vc '^#' || true)
echo "  Manta SV signal: $CAND_N candidate(s) (pre-filter) -> $SOM_N somatic (post-filter)"
if [[ "$CAND_N" -eq 0 ]]; then
    echo "  NOTE: 0 candidates -> reads present NO detectable SV signal to the aligner"
    echo "        (coverage too low, too few truth SVs, or an eidolon SV-read-emission gap)."
elif [[ "$SOM_N" -eq 0 ]]; then
    echo "  NOTE: candidates exist but 0 somatic -> signal IS there; Manta's somatic"
    echo "        classifier dropped them (low support at this coverage/purity)."
fi

# ── Stage 4b: Delly somatic SV calling (2nd caller vs Manta) ─────────────
# Cross-validation: an independent caller confirms eidolon's SVs aren't recovered
# only by Manta — and is the direct test of the BND-recovery question (does Delly
# recover the translocations Manta re-represents/misses?). RUN_DELLY=0 to skip.
DELLY_SV=""
if [[ "$RUN_DELLY" == "1" ]]; then
    DELLY_SV="$OUTDIR/delly.somatic.vcf.gz"
    if [[ -f "$DELLY_SV" ]]; then
        echo "[4b] Delly somatic SV already present — skipping."
    else
        echo "[4b] Running Delly (tumor/normal somatic SV)..."
        conda_activate "$DELLY_ENV"
        # delly call takes both BAMs (indexed) + a faidx'd reference (we have both).
        delly call -g "$REFERENCE" -o "$OUTDIR/delly.bcf" "$TUMOR_BAM" "$NORMAL_BAM" 2>&1 | tail -3
        # Somatic filter needs a sample->{tumor,control} TSV; our @RG SM tags are
        # TUMOR / NORMAL (set in index_and_align).
        printf 'TUMOR\ttumor\nNORMAL\tcontrol\n' > "$OUTDIR/delly_samples.tsv"
        delly filter -f somatic -s "$OUTDIR/delly_samples.tsv" \
            -o "$OUTDIR/delly.somatic.bcf" "$OUTDIR/delly.bcf" 2>&1 | tail -3
        # Convert BCF->VCF.gz for truvari (bcftools lives in bioinf, not the delly env).
        conda_activate bioinf
        bcftools view "$OUTDIR/delly.somatic.bcf" -O z -o "$DELLY_SV"
        bcftools index -f -t "$DELLY_SV"
        echo "[4b] Delly complete: $(zcat "$DELLY_SV" | grep -vc '^#' || true) somatic SV calls."
    fi
fi

# ── Stage 4c: GATK somatic CNV (validate the CNV machinery via a caller) ─
# eidolon's CNVs were only depth-validated; here a real copy-number caller must
# recover them. GATK4 (bioinf) — no new env. RUN_CNV=0 to skip.
CNV_TRUTH="$OUTDIR/truth_sv_CNV.vcf.gz"
if [[ "$RUN_CNV" == "1" && -f "$CNV_TRUTH" ]]; then
    conda_activate bioinf
    # GATK CNV is a SECONDARY validation — never let it abort the run (truvari
    # scoring of Manta/Delly must still happen). Guarded; full GATK output -> cnv_gatk.log.
    run_gatk_cnv() {
        local log="$OUTDIR/cnv_gatk.log"; : > "$log"
        [[ -f "${REFERENCE%.fa}.dict" ]] || gatk CreateSequenceDictionary -R "$REFERENCE" >>"$log" 2>&1 || return 1
        gatk PreprocessIntervals -R "$REFERENCE" --bin-length 1000 --padding 0 \
            --interval-merging-rule OVERLAPPING_ONLY -O "$OUTDIR/cnv.interval_list" >>"$log" 2>&1 || return 1
        gatk CollectReadCounts -I "$OUTDIR/tumor.bam" -L "$OUTDIR/cnv.interval_list" \
            --interval-merging-rule OVERLAPPING_ONLY -O "$OUTDIR/tumor.counts.hdf5" >>"$log" 2>&1 || return 1
        # No panel-of-normals: CreateReadCountPanelOfNormals is a Spark tool and
        # Spark 3.5 fails on the env's Java 25 (Subject.getSubject removed). With a
        # single matched normal a PoN is overkill anyway — DenoiseReadCounts without
        # one standardizes the tumor against its own median, which still surfaces the
        # copy-ratio deviations at eidolon's CNV loci.
        gatk DenoiseReadCounts -I "$OUTDIR/tumor.counts.hdf5" \
            --standardized-copy-ratios "$OUTDIR/tumor.stdCR.tsv" \
            --denoised-copy-ratios "$OUTDIR/tumor.denoisedCR.tsv" >>"$log" 2>&1 || return 1
        gatk ModelSegments --denoised-copy-ratios "$OUTDIR/tumor.denoisedCR.tsv" \
            -O "$OUTDIR" --output-prefix cnv >>"$log" 2>&1 || return 1
        gatk CallCopyRatioSegments -I "$OUTDIR/cnv.cr.seg" -O "$OUTDIR/cnv_called.seg" >>"$log" 2>&1 || return 1
    }
    if [[ -f "$OUTDIR/cnv_called.seg" ]]; then
        echo "[4c] GATK CNV segments already present — skipping call."
    elif run_gatk_cnv; then
        echo "[4c] GATK CNV calling complete."
    else
        echo "[4c] WARNING: GATK CNV failed (see $OUTDIR/cnv_gatk.log) — continuing without CNV validation." >&2
    fi
    # Score: does a called gain/loss segment overlap each truth CNV with the right
    # direction? (CN>2 -> gain '+', CN<2 -> loss '-'.) Manta/truvari can't score CNV,
    # so this is the CNV feature's caller-level validation.
    if [[ -f "$OUTDIR/cnv_called.seg" ]]; then
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/CN\n' "$CNV_TRUTH" \
            | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,($4>2?"+":($4<2?"-":"0"))}' > "$OUTDIR/cnv_truth.bed"
        awk '!/^@/ && $1!="CONTIG" && ($6=="+"||$6=="-"){print $1,$2,$3,$6}' \
            "$OUTDIR/cnv_called.seg" > "$OUTDIR/cnv_called.bed"
        awk 'NR==FNR{n++;c[n]=$1;st[n]=$2;e[n]=$3;d[n]=$4;next}
             {for(i=1;i<=n;i++) if($1==c[i] && $2<=e[i] && $3>=st[i] && $4==d[i]) hit[i]=1}
             END{tp=0; for(i=1;i<=n;i++) tp+=hit[i];
                 printf "[4c] CNV recovery: %d/%d truth CNVs called with matching direction\n", tp+0, n+0}' \
             "$OUTDIR/cnv_truth.bed" "$OUTDIR/cnv_called.bed"
    fi
elif [[ "$RUN_CNV" == "1" ]]; then
    echo "[4c] no CNV truth (truth_sv_CNV.vcf.gz) — skipping CNV validation (raise SV_RATE_SCALE / coverage)."
fi

# ── Stage 5: Score each caller with truvari (overall + per SV type) ──────
echo "[5/5] Scoring SVs with truvari..."
conda_activate "$TRUVARI_ENV"

run_truvari() {  # <truth.vcf.gz> <out-label> <comp.vcf.gz>
    # Separate `local`s on purpose: a combined `local a=$1 b=$2 c=...$b` expands
    # $b before it's assigned, which trips `set -u` (unbound variable).
    local truth="$1"
    local label="$2"
    local comp="$3"
    local out="$OUTDIR/truvari_$label"
    rm -rf "$out"   # truvari requires the output dir not pre-exist
    truvari bench -b "$truth" -c "$comp" -f "$REFERENCE" -o "$out" \
        --passonly --sizemax "$TRUVARI_SIZEMAX" 2>&1 | tail -3 || true
    if [[ -f "$out/summary.json" ]]; then
        python3 - "$out/summary.json" "$label" <<'PY'
import json, sys
s = json.load(open(sys.argv[1])); name = sys.argv[2]
# truvari writes null (None) for metrics when a subset has no comparisons;
# `or 0` coerces None -> 0 (s.get(k, 0) only defaults a MISSING key, not a null one).
print(f"  {name:<18} recall={(s.get('recall') or 0):.3f}  precision={(s.get('precision') or 0):.3f}  "
      f"f1={(s.get('f1') or 0):.3f}  TP={s.get('TP-base') or 0}  FN={s.get('FN') or 0}  FP={s.get('FP') or 0}")
PY
    else
        echo "  $label: truvari produced no summary.json (see $out)"
    fi
}

# Score one caller's VCF against the truth — overall + each per-type subset.
# Outputs land under truvari_<caller>_<overall|TYPE>/, so Manta vs Delly is a
# direct side-by-side (esp. for BND).
score_caller() {  # <caller> <comp.vcf.gz>
    local caller="$1" comp="$2"
    [[ -f "$comp" ]] || { echo "  ($caller: no calls to score)"; return 0; }
    echo "  === $caller ==="
    run_truvari "$TRUTH_SV" "${caller}_overall" "$comp"
    for typed in "$OUTDIR"/truth_sv_*.vcf.gz; do
        [[ -f "$typed" ]] || continue
        local svt; svt=$(basename "$typed" .vcf.gz); svt=${svt#truth_sv_}
        run_truvari "$typed" "${caller}_${svt}" "$comp"
    done
}

score_caller manta "$MANTA_SV"
[[ -n "$DELLY_SV" ]] && score_caller delly "$DELLY_SV"

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

════════════════════════════════════════════════════════════════
SV validation pipeline complete — Job $SLURM_JOB_ID

Inputs:
  Reference:    $REFERENCE
  Coverage:     ${TOTAL_COVERAGE}x  Purity: $PURITY  SV scale: $SV_RATE_SCALE
  Tumor model:  $TUMOR_MODEL

Key outputs in $OUTDIR:
  SV truth:     truth_sv.vcf.gz            ($SV_TRUTH_N somatic SVs)
  Manta calls:  manta/results/variants/somaticSV.vcf.gz
  Delly calls:  delly.somatic.vcf.gz       (RUN_DELLY=1)
  Scores:       truvari_<caller>_<overall|TYPE>/summary.json   (manta + delly)

Interpretation: compare the two callers side-by-side. A type eidolon emitted heavily
that BOTH callers miss = candidate SV-generation issue; one caller recovering what
the other misses (especially BND) = a caller limitation, not a simulator defect.
Low overall with high per-type recall = truvari SVTYPE/SVLEN representation
mismatch (inspect the pre-flight schema dump).
EOF

# Reclaim scratch once scoring is done: FASTQ is the bulk (~GBs/run at 60x) and
# isn't needed afterward — BAMs, Manta calls, the truth VCFs, and the truvari
# summaries are kept, so re-scoring and re-calling still work. This keeps a
# parameter sweep (per-tissue / purity) from exhausting the shared project quota.
# Set PRUNE_FASTQ=0 to retain the reads.
if [[ "${PRUNE_FASTQ:-1}" == "1" ]]; then
    echo "[prune] removing FASTQ intermediates (set PRUNE_FASTQ=0 to keep)..."
    rm -f "$OUTDIR"/sample_*_r1.fastq.gz "$OUTDIR"/sample_*_r2.fastq.gz
fi

# Persist report artifacts off scratch + record provenance/resources.
# (manta overall always exists; delly summaries remain under $OUTDIR/truvari_delly_*)
archive_run sv "$OUTDIR" truvari_manta_overall/summary.json
