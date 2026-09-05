#!/bin/bash
# SLURM job: end-to-end cancer pipeline on Delta.
#
# Stages:
#   1. Simulate tumor/normal FASTQ + origin-tagged truth VCF via
#      tools/cancer_simulate.sh (wraps eidolon gen-reads twice).
#   2. Align normal + merged-tumor FASTQs to reference with BWA-MEM2.
#   3. Sort, index BAMs with SAMtools.
#   4. Call somatic variants with GATK Mutect2 (tumor/normal mode).
#   5. Filter calls with GATK FilterMutectCalls.
#   6. Score filtered calls against the eidolon truth VCF using hap.py
#      som.py inside an Apptainer container (Docker not available on HPC).
#
# Each stage is idempotent — outputs that already exist are skipped,
# so re-running after a failure resumes from the last successful stage.
#
# Prerequisites:
#   - Run scripts/delta/setup.sh once before submitting this job.
#   - Stage your reference FASTA into $SCRATCH (see REFERENCE below).
#   - Modules: bwa-mem2, samtools, gatk, bcftools must be available.
#     Adjust `module load` lines to match Delta's current module names.
#
# Usage:
#   sbatch scripts/delta/cancer_pipeline.sbatch
#
# To override key parameters without editing:
#   REFERENCE=/path/to/ref.fa PURITY=0.4 sbatch scripts/delta/cancer_pipeline.sbatch

#SBATCH --job-name=eidolon-cancer
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

# ── Parameters (edit these) ─────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with EIDOLON_REPO.
REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

# Reference genome. For a smoke test use chr22; for a full-scale run swap
# in GRCh38. Must be uncompressed (.fa, not .fa.gz) for BWA-MEM2 indexing.
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"

# Output root — all intermediate files and final scores go here.
OUTDIR="${OUTDIR:-$SCRATCH/cancer_$SLURM_JOB_ID}"

# Simulation parameters.
TOTAL_COVERAGE="${TOTAL_COVERAGE:-30}"
PURITY="${PURITY:-0.6}"          # tumor cell fraction
# Tumor somatic per-base mutation rate. Empty → cancer_simulate.sh's default
# (1e-5, a typical solid tumor). Set to sweep tumor burden, or "model" to defer
# to the tumor model's fitted (corpus-aggregated) rate.
TUMOR_MUTATION_RATE="${TUMOR_MUTATION_RATE:-}"
READ_LEN="${READ_LEN:-151}"
PAIRED="${PAIRED:-true}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
SV_RATE_SCALE="${SV_RATE_SCALE:-1.0}"   # >0 enables SV simulation in tumor pass
# PRUNE=1 deletes the bulky FASTQ + BAM intermediates after a SUCCESSFUL score,
# keeping only the small result artifacts (scored*.stats.csv, VCFs, logs). Use for
# large sweeps so per-point output doesn't exhaust the shared scratch project quota.
PRUNE="${PRUNE:-0}"

# Tumor mutation model (COSMIC pan-cancer). Bundled in the repo.
TUMOR_MODEL="$REPO_ROOT/tools/cosmic_v104_pancancer_model.json.gz"
# Optional. Set to compare error models against an otherwise identical run -- the
# somatic_snv_recall question in #680 needs exactly that A/B.
SEQ_ERROR_MODEL="${SEQ_ERROR_MODEL:-}"
SEQ_ERR_FLAG=""
# `if`, not `[[ ]] && x=y`. The && form survives at top level but returns 1 when the test
# is false, so it aborts under `set -e` the moment anyone moves it inside a function --
# the footgun cancer_simulate.sh already documents at its own config writer.
if [[ -n "$SEQ_ERROR_MODEL" ]]; then
    SEQ_ERR_FLAG="--seq-error-model $SEQ_ERROR_MODEL"
fi

# Path to the hap.py Apptainer SIF (built by setup.sh).
HAPPY_SIF="${HAPPY_SIF:-$SCRATCH/sif/happy_v0.3.12.sif}"

# Secondary somatic SNV/indel callers for cross-validation vs Mutect2. VarScan2 is
# the working cross-check; Strelka2 2.9.10 SIGSEGVs on Delta's stack so it's OFF by
# default (set RUN_STRELKA=1 to attempt it). Both are non-fatal — a failure never
# aborts the primary Mutect2 result.
RUN_VARSCAN="${RUN_VARSCAN:-1}"
VARSCAN_ENV="${VARSCAN_ENV:-varscan}"
RUN_STRELKA="${RUN_STRELKA:-0}"
STRELKA_ENV="${STRELKA_ENV:-strelka}"

# Default the SLURM vars so a direct (non-sbatch) run doesn't abort under `set -u`.
: "${SLURM_JOB_ID:=manual}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"
PREFIX="$OUTDIR/sample"

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf   # provides bwa-mem2, gatk4, bcftools (not available as modules)

mkdir -p "$OUTDIR"
# Refuse to run concurrently in an OUTDIR another job already owns —
# interleaved writes yield results computed from corrupted intermediates.
lock_outdir "$OUTDIR" || exit 1
REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"

echo "=== eidolon cancer pipeline: $SLURM_JOB_ID ==="
echo "Reference:  $REFERENCE"
echo "Output:     $OUTDIR"
echo "Coverage:   ${TOTAL_COVERAGE}x  Purity: $PURITY"
echo "Threads:    $THREADS"
echo "SV scale:   $SV_RATE_SCALE"

# ── Stage 1: Simulate ────────────────────────────────────────────────────
MERGED_R1="${PREFIX}_merged_r1.fastq.gz"
NORMAL_R1="${PREFIX}_normal_r1.fastq.gz"
TRUTH_VCF="${PREFIX}_merged_truth.vcf.gz"

# Refuse to score fresh caller output against a truth VCF written by a different
# build of eidolon (see lib_report.sh::check_outdir_version).
check_outdir_version "$EIDOLON_BIN" "$OUTDIR" || exit 1

SIM_REGENERATED=0
if [[ "$FORCE_RESIM" == "0" && -f "$MERGED_R1" && -f "$NORMAL_R1" && -f "$TRUTH_VCF" ]]; then
    echo "[1/6] Simulation outputs already present — skipping."
else
    echo "[1/6] Simulating tumor/normal reads..."
    PAIRED_FLAG=""
    [[ "$PAIRED" == "true" ]] && PAIRED_FLAG="--paired-ended --fragment-mean $FRAG_MEAN --fragment-st-dev $FRAG_SD"
    TMR_FLAG=""
    [[ -n "$TUMOR_MUTATION_RATE" ]] && TMR_FLAG="--tumor-mutation-rate $TUMOR_MUTATION_RATE"
    # Under `sbatch --array=1-3` (triplicate), each replicate uses a distinct
    # seed root (independent draw); SLURM spreads tasks across nodes. Plain
    # submit uses the default root.
    RNG_SEED_ROOT="${RNG_SEED_ROOT:-cancer-simulate}${SLURM_ARRAY_TASK_ID:+-rep$SLURM_ARRAY_TASK_ID}"
    bash "$REPO_ROOT/tools/cancer_simulate.sh" \
        --reference "$REFERENCE" \
        --output-prefix "$PREFIX" \
        --total-coverage "$TOTAL_COVERAGE" \
        --purity "$PURITY" \
        --read-len "$READ_LEN" \
        --tumor-model "$TUMOR_MODEL" \
        $SEQ_ERR_FLAG \
        --sv-rate-scale "$SV_RATE_SCALE" \
        --rng-seed "$RNG_SEED_ROOT" \
        --eidolon-bin "$EIDOLON_BIN" \
        $TMR_FLAG \
        $PAIRED_FLAG
    echo "[1/6] Simulation complete."
    SIM_REGENERATED=1
    stamp_outdir_version "$OUTDIR"
fi

# ── Stage 2 & 3: Index reference, align FASTQs ──────────────────────────
NORMAL_BAM="$OUTDIR/normal.bam"
TUMOR_BAM="$OUTDIR/tumor.bam"

# If the simulation was (re)generated this run, any pre-existing downstream artifacts
# (BAMs, caller VCFs, scoring truth + scores) were built from a DIFFERENT truth VCF and
# are now stale. The idempotent stages below would otherwise skip on those and score the
# fresh truth against stale calls (a re-run with a fixed OUTDIR + PRUNE=1 hit exactly this).
# Invalidate them so they rebuild consistently against the new simulation.
if [[ "$SIM_REGENERATED" == "1" ]]; then
    echo "[1/6] Simulation regenerated — invalidating stale downstream outputs."
    rm -f "$NORMAL_BAM" "${NORMAL_BAM}.bai" "$TUMOR_BAM" "${TUMOR_BAM}.bai" \
          "$OUTDIR"/mutect2.*.vcf.gz "$OUTDIR"/mutect2.*.vcf.gz.tbi "$OUTDIR"/mutect2.*.vcf.gz.stats \
          "$OUTDIR"/varscan.somatic.vcf.gz* "$OUTDIR"/strelka.somatic.vcf.gz* \
          "$OUTDIR/scoring_truth.vcf.gz" "$OUTDIR/scoring_truth.vcf.gz.tbi" \
          "$OUTDIR"/scored*.csv "$OUTDIR"/scored*.json "$OUTDIR"/scored*.vcf.gz* 2>/dev/null || true
fi

index_and_align() {
    local in_fq="$1" out_bam="$2" sample="$3"

    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then
        echo "  $sample BAM already present — skipping."
        return 0
    fi

    # Index reference once, serialized so concurrent jobs can't race the shared
    # index. The -s (not -f) guard against a 0-byte .bwt.2bit.64 from a build
    # killed mid-write lives in index_reference_locked (lib_report.sh).
    index_reference_locked "$REFERENCE"

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

echo "[2/6] Aligning normal reads..."
index_and_align "$NORMAL_R1" "$NORMAL_BAM" "NORMAL"

echo "[3/6] Aligning tumor (merged) reads..."
# Merged FASTQ is the purity-mixed tumor biopsy.
TUMOR_MERGED_R1="${PREFIX}_merged_r1.fastq.gz"
index_and_align "$TUMOR_MERGED_R1" "$TUMOR_BAM" "TUMOR"

# ── Stage 4: Mutect2 ─────────────────────────────────────────────────────
RAW_VCF="$OUTDIR/mutect2.unfiltered.vcf.gz"
if [[ -f "$RAW_VCF" ]]; then
    echo "[4/6] Mutect2 raw VCF already present — skipping."
else
    echo "[4/6] Running Mutect2 (tumor/normal mode)..."
    # Sequence dictionary required by GATK.
    REF_DICT="${REFERENCE%.fa}.dict"
    [[ -f "$REF_DICT" ]] || gatk CreateSequenceDictionary -R "$REFERENCE"

    gatk Mutect2 \
        -R "$REFERENCE" \
        -I "$TUMOR_BAM"  --tumor-sample  TUMOR \
        -I "$NORMAL_BAM" --normal-sample NORMAL \
        --native-pair-hmm-threads "$THREADS" \
        -O "$RAW_VCF" 2>&1 | tail -20
    echo "[4/6] Mutect2 complete."
fi

# ── Stage 5: FilterMutectCalls ───────────────────────────────────────────
FILTERED_VCF="$OUTDIR/mutect2.filtered.vcf.gz"
if [[ -f "$FILTERED_VCF" ]]; then
    echo "[5/6] FilterMutectCalls output already present — skipping."
else
    echo "[5/6] Filtering Mutect2 calls..."
    gatk FilterMutectCalls \
        -R "$REFERENCE" \
        -V "$RAW_VCF" \
        -O "$FILTERED_VCF"
    echo "[5/6] Filtering complete."
fi

# ── Stage 5b: Strelka2 somatic SNV/indel (2nd caller vs Mutect2) ─────────
# Cross-validation: a second, independently-modeled somatic caller confirms
# eidolon's somatic SNV/indel data isn't recovered only by Mutect2. Reuses the same
# tumor/normal BAMs. RUN_STRELKA=0 to skip.
STRELKA_VCF=""
if [[ "$RUN_STRELKA" == "1" ]]; then
    STRELKA_DIR="$OUTDIR/strelka"
    STRELKA_VCF="$OUTDIR/strelka.somatic.vcf.gz"
    if [[ -f "$STRELKA_VCF" ]]; then
        echo "[5b] Strelka2 somatic VCF already present — skipping."
    else
        echo "[5b] Running Strelka2 (tumor/normal somatic SNV/indel)..."
        # Strelka2 is a SECONDARY cross-check and must NEVER abort the run: the
        # bioconda strelka 2.9.10 binary SIGSEGVs on some HPC glibc/CPU combos.
        # Run it guarded — on failure, warn and continue so the primary Mutect2
        # result is still scored (Stage 6).
        run_strelka() {
            conda_activate "$STRELKA_ENV"
            rm -rf "$STRELKA_DIR"
            configureStrelkaSomaticWorkflow.py \
                --normalBam "$NORMAL_BAM" --tumorBam "$TUMOR_BAM" \
                --referenceFasta "$REFERENCE" --runDir "$STRELKA_DIR" \
                > "$OUTDIR/strelka_config.log" 2>&1 || return 1
            "$STRELKA_DIR/runWorkflow.py" -m local -j "$THREADS" \
                > "$OUTDIR/strelka_run.log" 2>&1 || return 1
        }
        if run_strelka; then
            conda_activate bioinf   # Strelka emits SNVs + indels separately; merge.
            bcftools concat -a \
                "$STRELKA_DIR/results/variants/somatic.snvs.vcf.gz" \
                "$STRELKA_DIR/results/variants/somatic.indels.vcf.gz" \
                | bcftools sort -Oz -o "$STRELKA_VCF"
            bcftools index -f -t "$STRELKA_VCF"
            echo "[5b] Strelka2 complete: $(zcat "$STRELKA_VCF" | grep -vc '^#' || true) somatic calls."
        else
            echo "[5b] WARNING: Strelka2 failed (see $OUTDIR/strelka_run.log) — continuing without it; Mutect2 unaffected." >&2
            STRELKA_VCF=""
        fi
    fi
fi

# ── Stage 5c: VarScan2 somatic SNV/indel (robust 2nd caller vs Mutect2) ──
# The working SNV/indel cross-check (Strelka2 is broken on Delta). Combined
# normal+tumor mpileup -> varscan somatic -> processSomatic (high-confidence
# somatic) -> merge SNP+indel -> scored vs the same somatic truth. Guarded:
# never aborts the run. RUN_VARSCAN=0 to skip.
VARSCAN_VCF=""
if [[ "$RUN_VARSCAN" == "1" ]]; then
    VARSCAN_VCF="$OUTDIR/varscan.somatic.vcf.gz"
    if [[ -f "$VARSCAN_VCF" ]]; then
        echo "[5c] VarScan2 somatic VCF already present — skipping."
    else
        echo "[5c] Running VarScan2 (tumor/normal somatic SNV/indel)..."
        # Everything (calling AND bcftools post-processing) is inside the guard so
        # a VarScan-VCF quirk can't abort the run. VarScan VCFs lack ##contig lines,
        # which makes bcftools sort fail on the undefined contig — add them from the
        # reference .fai via `reheader --fai` before sorting.
        run_varscan() {
            local vlog="$OUTDIR/varscan_run.log"
            conda_activate "$VARSCAN_ENV"   # varscan (samtools comes from the module)
            # Combined mpileup, NORMAL first then TUMOR (--mpileup 1 expects this order).
            samtools mpileup -B -q 1 -f "$REFERENCE" "$NORMAL_BAM" "$TUMOR_BAM" \
                > "$OUTDIR/varscan.mpileup" 2>"$OUTDIR/varscan_mpileup.log" || return 1
            varscan somatic "$OUTDIR/varscan.mpileup" "$OUTDIR/varscan" \
                --mpileup 1 --output-vcf 1 > "$vlog" 2>&1 || return 1
            varscan processSomatic "$OUTDIR/varscan.snp.vcf"   >> "$vlog" 2>&1 || return 1
            varscan processSomatic "$OUTDIR/varscan.indel.vcf" >> "$vlog" 2>&1 || return 1
            conda_activate bioinf   # bcftools
            local v
            for v in snp indel; do
                bcftools reheader --fai "${REFERENCE}.fai" "$OUTDIR/varscan.$v.Somatic.hc.vcf" 2>>"$vlog" \
                    | bcftools sort -Oz -o "$OUTDIR/varscan.$v.somatic.vcf.gz" 2>>"$vlog" || return 1
                bcftools index -f -t "$OUTDIR/varscan.$v.somatic.vcf.gz" 2>>"$vlog" || return 1
            done
            bcftools concat -a "$OUTDIR/varscan.snp.somatic.vcf.gz" "$OUTDIR/varscan.indel.somatic.vcf.gz" 2>>"$vlog" \
                | bcftools sort -Oz -o "$VARSCAN_VCF" 2>>"$vlog" || return 1
            bcftools index -f -t "$VARSCAN_VCF" 2>>"$vlog" || return 1
            rm -f "$OUTDIR/varscan.mpileup"   # large intermediate — reclaim space
        }
        if run_varscan; then
            echo "[5c] VarScan2 complete: $(zcat "$VARSCAN_VCF" | grep -vc '^#' || true) somatic calls."
        else
            echo "[5c] WARNING: VarScan2 failed (see $OUTDIR/varscan_run.log) — continuing without it; Mutect2 unaffected." >&2
            VARSCAN_VCF=""
        fi
    fi
fi

# ── Stage 6: Score each caller with hap.py som.py ───────────────────────
# Somatic-only truth (EIDOLON_ORIGIN=somatic); germline carry-through (shared) is
# correctly ignored by T/N callers, so scoring against it would be meaningless.
# Strelka (5b) may leave us in its conda env (especially on its failure path);
# scoring needs bioinf (bcftools) + apptainer, so re-activate it here.
conda_activate bioinf
# Normalize pre-v2.1.0 NEAT_* tokens if this truth VCF predates the rename (e.g.
# re-scoring an archived artifact). The converter is idempotent and a no-op on a
# current file, so calling it unconditionally is safe. It matters because bcftools
# rejects an undefined INFO tag in the -i expression below at parse time.
# `|| rc=$?` (not `; echo $?`) so the non-zero --check status is *tested* — under
# set -e an untested non-zero would abort before we could read it. Only 0 and 10 are
# meaningful: any other status means the pre-flight itself failed (missing input, no
# bcftools, bad converter path) and must NOT be mistaken for "already current", or the
# auto-convert silently no-ops and the run dies later on a raw undefined-tag error.
# stderr is deliberately not suppressed so the reason is visible.
migrate_rc=0
bash "$REPO_ROOT/tools/migrate_legacy_tokens.sh" --check "$TRUTH_VCF" >/dev/null \
    || migrate_rc=$?
case "$migrate_rc" in
    0) ;;   # already current
    10) MIGRATED_TRUTH="$OUTDIR/truth_migrated.vcf.gz"
        bash "$REPO_ROOT/tools/migrate_legacy_tokens.sh" "$TRUTH_VCF" "$MIGRATED_TRUTH"
        TRUTH_VCF="$MIGRATED_TRUTH" ;;
    *)  echo "ERROR: legacy-token pre-flight failed (rc=$migrate_rc) on $TRUTH_VCF" >&2
        echo "  Not proceeding: an unconverted truth VCF would fail the filter below." >&2
        exit 1 ;;
esac

SCORE_TRUTH="$OUTDIR/scoring_truth.vcf.gz"
if [[ ! -f "$SCORE_TRUTH" ]]; then
    echo "[6/6] Filtering truth to somatic records..."
    bcftools view -i 'INFO/EIDOLON_ORIGIN="somatic"' -O z -o "$SCORE_TRUTH" "$TRUTH_VCF"
    bcftools index -f -t "$SCORE_TRUTH"
fi
somatic_count=$(bcftools view -H "$SCORE_TRUTH" | wc -l)
echo "      Somatic truth records: $somatic_count"
if [[ "$somatic_count" -eq 0 ]]; then
    echo "ERROR: no somatic truth records after filtering with:" >&2
    echo "         INFO/EIDOLON_ORIGIN=\"somatic\"" >&2
    echo "  Scoring against an empty truth would report every caller as 0 recall with" >&2
    echo "  no indication why, so this is fatal. Check INFO/EIDOLON_ORIGIN is present" >&2
    echo "  in $TRUTH_VCF (pre-v3.0.0 truth VCFs use NEAT_ORIGIN — convert with" >&2
    echo "  tools/migrate_legacy_tokens.sh), or that the tumor pass planted somatic" >&2
    echo "  variants at all (raise TUMOR_MUTATION_RATE / coverage)." >&2
    exit 1
fi
[[ -f "$HAPPY_SIF" ]] || { echo "hap.py SIF not found: $HAPPY_SIF — run setup.sh first" >&2; exit 1; }

# hap.py/som.py hardcodes UCSC chr-prefixing of input VCFs (perl s/^([0-9XYM])/chr/),
# then bcftools-norm's against the reference. With an Ensembl-named reference (1/2/3)
# the renamed VCFs (chr1) no longer match the reference (1) and som.py dies
# ("sequence chr1 not found"). If the reference isn't chr-prefixed, build a
# chr-prefixed copy for scoring so the renamed VCFs match. (chr22.fa etc. pass
# through unchanged.) Alignment/Mutect2 already used $REFERENCE as-is — unaffected.
SCORE_REF_DIR="$REF_DIR"; SCORE_REF_NAME="$REF_NAME"
# Only build the chr-prefixed copy when som.py will ACTUALLY rename the VCF contigs.
# Its perl is `s/^([0-9XYM])/chr/`, so it renames 1/2/X/Y/MT but leaves anything else
# alone — including GenBank accession names (CM000851.4, soybean). Keying off "not
# chr-prefixed" lumped those in with Ensembl-numeric and prefixed the REFERENCE while
# som.py left the VCF unprefixed, so scoring died with a contig mismatch
# ("Error running BCFTOOLS ... vcfcheck") on a reference that needed no help at all.
if head -1 "${REFERENCE}.fai" 2>/dev/null | grep -qE '^[0-9XYM]'; then
    SCORE_REF="${REFERENCE%.fa}.chr.fa"
    if [[ ! -f "$SCORE_REF" ]]; then
        echo "[6/6] Reference is Ensembl-named; building chr-prefixed copy for som.py..."
        sed '/^>/ s/^>/>chr/' "$REFERENCE" > "$SCORE_REF"
        samtools faidx "$SCORE_REF"
    fi
    SCORE_REF_DIR="$(dirname "$SCORE_REF")"; SCORE_REF_NAME="$(basename "$SCORE_REF")"
fi

# som_score <query.vcf.gz> <out-label>. -N (--normalize-all) runs bcftools norm on
# BOTH truth and query — essential for indels (representation/left-align mismatch
# otherwise double-counts FN+FP). HGREF resolves som.py's reference lookup.
som_score() {
    local query="$1" label="$2"
    [[ -f "$query" ]] || { echo "  ($label: no query VCF to score)"; return 0; }
    if [[ -f "$OUTDIR/${label}.stats.csv" ]]; then
        echo "[6/6] $label score already present — skipping."; return 0
    fi
    echo "[6/6] Scoring $label with som.py (Apptainer)..."
    apptainer exec \
        --env "HGREF=/refs/$SCORE_REF_NAME" \
        --bind "$SCORE_REF_DIR:/refs:ro" \
        --bind "$OUTDIR:/work" \
        "$HAPPY_SIF" \
        /opt/hap.py/bin/som.py \
            "/work/$(basename "$SCORE_TRUTH")" \
            "/work/$(basename "$query")" \
            -r "/refs/$SCORE_REF_NAME" \
            -o "/work/$label" \
            -N \
            --feature-table generic
}

som_score "$FILTERED_VCF" scored                                    # Mutect2
[[ -n "$VARSCAN_VCF" ]] && som_score "$VARSCAN_VCF" scored_varscan   # VarScan2 (2nd caller)
[[ -n "$STRELKA_VCF" ]] && som_score "$STRELKA_VCF" scored_strelka   # Strelka2 (opt-in)
echo "[6/6] Scoring complete."

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

════════════════════════════════════════════════════════════════
Cancer pipeline complete — Job $SLURM_JOB_ID

Inputs:
  Reference:      $REFERENCE
  Coverage:       ${TOTAL_COVERAGE}x  Purity: $PURITY
  SV rate scale:  $SV_RATE_SCALE

Key outputs in $OUTDIR:
  Simulation:     sample_merged_{r1,r2}.fastq.gz
  Truth VCF:      sample_merged_truth.vcf.gz  (EIDOLON_ORIGIN-tagged)
  Somatic truth:  scoring_truth.vcf.gz         (somatic only)
  Mutect2 calls:  mutect2.filtered.vcf.gz   -> scored.stats.csv
  VarScan2 calls: varscan.somatic.vcf.gz    -> scored_varscan.stats.csv  (RUN_VARSCAN=1)
  Strelka2 calls: strelka.somatic.vcf.gz    -> scored_strelka.stats.csv  (RUN_STRELKA=1, opt-in)

Top-line results (Mutect2):
EOF
[[ -f "$OUTDIR/scored.stats.csv" ]] && column -t -s',' "$OUTDIR/scored.stats.csv" || echo "  (Mutect2 score not found)"
if [[ -n "$VARSCAN_VCF" && -f "$OUTDIR/scored_varscan.stats.csv" ]]; then
    echo; echo "Top-line results (VarScan2):"
    column -t -s',' "$OUTDIR/scored_varscan.stats.csv"
fi
if [[ -n "$STRELKA_VCF" && -f "$OUTDIR/scored_strelka.stats.csv" ]]; then
    echo; echo "Top-line results (Strelka2):"
    column -t -s',' "$OUTDIR/scored_strelka.stats.csv"
fi

# Persist report artifacts off scratch + record provenance/resources.
archive_run cancer "$OUTDIR" scored.stats.csv

# Optional cleanup of bulky intermediates once scored + archived, so large sweeps
# don't exhaust the shared scratch project quota. Only prunes after a SUCCESSFUL
# score (scored.stats.csv present) so a failed run keeps its inputs for resume.
if [[ "$PRUNE" == "1" && -f "$OUTDIR/scored.stats.csv" ]]; then
    echo "[prune] PRUNE=1 — removing FASTQ + BAM intermediates (keeping scores, VCFs, logs)..."
    rm -f "$OUTDIR"/*.fastq.gz \
          "$OUTDIR"/normal.bam "$OUTDIR"/normal.bam.bai \
          "$OUTDIR"/tumor.bam  "$OUTDIR"/tumor.bam.bai
fi
