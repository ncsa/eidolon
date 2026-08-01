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
#   truvari now installs into `bioinf` alongside bcftools (solve verified clean), so
#   scoring no longer switches envs — that split existed only because Manta needs
#   python2, and it forced per-type subsets to be pre-generated in an earlier stage
#   since bcftools was absent during scoring.
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
TRUVARI_ENV="${TRUVARI_ENV:-bioinf}"   # truvari lives in bioinf; see prereqs
# Run Delly as a 2nd somatic SV caller alongside Manta (cross-validation; the
# BND-recovery question). Set RUN_DELLY=0 to skip.
RUN_DELLY="${RUN_DELLY:-1}"
# Run GATK somatic CNV to validate the CNV machinery through a real copy-number
# caller (eidolon's CNVs were previously only depth-validated). No new env (gatk4 in
# bioinf). Set RUN_CNV=0 to skip.
RUN_CNV="${RUN_CNV:-1}"
# eidolon SVs reach ~1 Mb; truvari's default --sizemax is 50 kb, which would silently
# EXCLUDE most large eidolon SVs from scoring. Raise it to cover the full range.
# truvari DISCARDS any record larger than this before comparison. The previous 5Mb
# default silently removed every SV above it — and sv_max_length_fraction (0.25)
# permits 14Mb events on a 56Mb scaffold, while BND mate positions are drawn across
# the whole contig, so junction spans reach tens of Mb. Default is now effectively
# unlimited; what it would exclude is logged before scoring.
TRUVARI_SIZEMAX="${TRUVARI_SIZEMAX:-100000000}"
# Passed to truvari as --bnddist: max per-breakpoint offset when matching BND-to-BND.
# Generous against short-read breakpoint precision, tight against inter-junction
# spacing at realistic SV rates.
BND_DIST="${BND_DIST:-500}"

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
# Refuse to run concurrently in an OUTDIR another job already owns —
# interleaved writes yield results computed from corrupted intermediates.
lock_outdir "$OUTDIR" || exit 1
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

# Refuse to score fresh caller output against a truth VCF written by a different
# build of eidolon (see lib_report.sh::check_outdir_version).
check_outdir_version "$EIDOLON_BIN" "$OUTDIR" || exit 1

if [[ "$FORCE_RESIM" == "0" && -f "$MERGED_R1" && -f "$NORMAL_R1" && -f "$TRUTH_VCF" ]]; then
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
    stamp_outdir_version "$OUTDIR"
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
    # Serialized: concurrent jobs must not race the shared index (lib_report.sh).
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
# Only 0 and 10 are meaningful: any other status means the pre-flight itself failed
# (missing input, no bcftools, bad converter path) and must NOT be mistaken for
# "already current", or the auto-convert silently no-ops and the run dies later on a
# raw undefined-tag error. stderr is deliberately not suppressed.
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

# The aggregate must cover only what truvari can match. A BND has no SVLEN and
# END==POS, and CNV is a copy-ratio call — truvari matches neither, so including
# them made every aggregate a guaranteed partial miss (#457). They are scored
# separately: BND by breakpoint proximity, CNV by the stage-4c direction check.
# What would --sizemax silently discard? Truvari drops oversized records BEFORE
# comparison, so an unreported ceiling turns real events into guaranteed misses.
n_over=$(bcftools query -f '%INFO/SVLEN\n' "$TRUTH_SV" 2>/dev/null \
         | awk -v m="$TRUVARI_SIZEMAX" '{v=($1<0?-$1:$1)+0; if (v>m) n++} END{print n+0}')
echo "  truth records exceeding --sizemax=$TRUVARI_SIZEMAX: $n_over"
[[ "$n_over" -gt 0 ]] && echo "    NOTE: those would be discarded before comparison — raise TRUVARI_SIZEMAX." >&2

TRUTH_SCOREABLE="$OUTDIR/truth_sv_scoreable.vcf.gz"
bcftools view -i 'INFO/SVTYPE!="BND" && INFO/SVTYPE!="CNV"' \
    -O z -o "$TRUTH_SCOREABLE" "$TRUTH_SV" 2>/dev/null || true
if [[ -f "$TRUTH_SCOREABLE" ]]; then
    bcftools index -f -t "$TRUTH_SCOREABLE" 2>/dev/null || true
    n_score=$(bcftools view -H "$TRUTH_SCOREABLE" | wc -l)
    n_bnd=$(bcftools view -H "$TRUTH_SV" -i 'INFO/SVTYPE="BND"' 2>/dev/null | wc -l)
    n_cnv=$(bcftools view -H "$TRUTH_SV" -i 'INFO/SVTYPE="CNV"' 2>/dev/null | wc -l)
    echo "  truvari-scoreable truth: $n_score record(s)"
    echo "    excluded from the aggregate: $n_bnd BND record(s) (= $((n_bnd / 2)) junction(s)), $n_cnv CNV"
    echo "    -> BND scored by proximity below; CNV by the stage-4c direction check."
    [[ "$n_score" -eq 0 ]] && rm -f "$TRUTH_SCOREABLE" "${TRUTH_SCOREABLE}.tbi"
fi

# ── Stage 4: Manta somatic SV calling ────────────────────────────────────
# Calibration/coverage flags. Declared here because stage 4a can set CAL_FAIL, which
# is well before the scoring stage that consumes it (set -u).
CAL_FAIL=0
SCORE_COVERAGE_FAIL=0

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

# ── Stage 4a: normalise Manta's inversion representation ─────────────────
# Manta stopped emitting <INV> around v1.5 and reports inversions as BND PAIRS
# instead, shipping libexec/convertInversion.py to convert them back. Delta runs
# 1.6.0, so `manta_INV recall=0.000` is DOCUMENTED TOOL BEHAVIOUR, not a miss.
# Measured in job 20684366: Manta placed 18 breakends within 500bp of 5 of the 6
# truth inversions (4 per inversion = 2 junctions x 2 records, 2 for a 1588bp one)
# while emitting zero INV records. Those 18 were then counted as manta_BND false
# positives, so BOTH numbers were representation artifacts:
#     manta_INV recall 0.000  -> 5/6 detected
#     manta_BND precision 0.438 -> every breakend accounted for
#
# The caller side is normalised with the CALLER'S OWN converter rather than by
# deriving a breakend representation of our truth. build_bnd_spans was correct code
# built on a representation we chose, and its negative control killed it; Manta
# ships the canonical mapping, so use it.
#
# Falsifiable expectation, checked below: eidolon's 7 junctions are LONE
# head-to-head adjacencies with no tail-to-tail partner, so they must NOT convert
# (14 BND records survive); the 18 properly-bracketed inversion breakends must.
convert_manta_inversions() {  # <in.vcf.gz> <out.vcf.gz>
    # rc 0 = converted, 2 = converter unavailable (older Manta), 1 = error.
    local src="$1" dst="$2"
    local samtools_bin; samtools_bin="$(command -v samtools || true)"
    [[ -n "$samtools_bin" ]] || { echo "  ERROR: samtools not on PATH for convertInversion" >&2; return 1; }

    conda_activate "$MANTA_ENV" || return 1
    local conv; conv="$(dirname "$(command -v configManta.py 2>/dev/null || echo /nonexistent)")/../libexec/convertInversion.py"
    if [[ ! -f "$conv" ]]; then
        echo "  Manta inversion conversion: convertInversion.py not found — this Manta" >&2
        echo "    predates the BND-pair inversion representation, so no conversion needed." >&2
        conda_activate bioinf || true
        return 2
    fi

    local plain="$OUTDIR/.manta_conv_in.vcf" out="$OUTDIR/.manta_conv_out.vcf" err="$OUTDIR/.manta_conv_err"
    zcat "$src" > "$plain"
    local rc=0
    # usage: convertInversion.py <samtools path> <reference fasta> <vcf file>
    python "$conv" "$samtools_bin" "$REFERENCE" "$plain" > "$out" 2> "$err" || rc=$?
    conda_activate bioinf || true
    if [[ "$rc" -ne 0 ]]; then
        echo "  ERROR: convertInversion.py failed (rc=$rc):" >&2
        sed 's/^/      /' "$err" >&2
        rm -f "$plain" "$out" "$err"; return 1
    fi

    # Verify the conversion did what it claims, rather than trusting a zero exit.
    local bnd_in bnd_out inv_out
    # grep -c prints 0 and exits 1 on no match; `|| true` keeps that from aborting
    # under `set -e` while still yielding the count.
    bnd_in=$(grep -c 'SVTYPE=BND' "$plain" || true)
    bnd_out=$(grep -c 'SVTYPE=BND' "$out" || true)
    inv_out=$(grep -c 'SVTYPE=INV' "$out" || true)
    if [[ "$inv_out" -eq 0 ]]; then
        echo "  ERROR: convertInversion.py produced NO <INV> records from $bnd_in breakend(s)." >&2
        echo "    Either no inversion-shaped BND pairs were present, or the conversion" >&2
        echo "    silently no-opped. Not scoring against an unconverted file." >&2
        rm -f "$plain" "$out" "$err"; return 1
    fi
    if [[ "$bnd_out" -ge "$bnd_in" ]]; then
        echo "  ERROR: BND count did not drop ($bnd_in -> $bnd_out) yet $inv_out INV" >&2
        echo "    record(s) appeared — the output is internally inconsistent." >&2
        rm -f "$plain" "$out" "$err"; return 1
    fi

    bcftools sort -O z -o "$dst" "$out" 2>/dev/null || { rm -f "$plain" "$out" "$err"; return 1; }
    bcftools index -f -t "$dst" 2>/dev/null || true
    echo "  Manta inversion conversion: $bnd_in BND -> $bnd_out BND + $inv_out INV"
    echo "    (Manta reports inversions as breakend pairs since ~v1.5; converted with"
    echo "     its own libexec/convertInversion.py so INV truth is scored like-for-like)"
    rm -f "$plain" "$out" "$err"
    return 0
}

MANTA_SV_RAW="$MANTA_SV"
if [[ -f "$MANTA_SV" ]]; then
    MANTA_SV_CONV="$OUTDIR/manta_somaticSV_invconv.vcf.gz"
    conv_rc=0
    convert_manta_inversions "$MANTA_SV" "$MANTA_SV_CONV" || conv_rc=$?
    case "$conv_rc" in
        0) MANTA_SV="$MANTA_SV_CONV" ;;   # score against the normalised calls
        2) ;;                              # no converter: nothing to normalise
        *) echo "ERROR: Manta inversion conversion failed — scoring against the raw" >&2
           echo "  calls would report manta_INV=0.000 for inversions Manta DID find," >&2
           echo "  and count those breakends as manta_BND false positives." >&2
           CAL_FAIL=1 ;;
    esac
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

# Pick the tumor sample from a caller VCF. truvari otherwise takes the FIRST sample,
# which for Manta's somaticSV is NORMAL — so every somatic comparison was reading the
# normal sample's genotypes.
comp_sample() {  # <comp.vcf.gz>
    local s
    s=$(bcftools query -l "$1" 2>/dev/null | grep -ixm1 -E 'tumor|tumour' || true)
    [[ -n "$s" ]] || s=$(bcftools query -l "$1" 2>/dev/null | tail -1)
    printf '%s' "$s"
}

run_truvari() {  # <truth.vcf.gz> <out-label> <comp.vcf.gz> [extra truvari args...]
    # Separate `local`s on purpose: a combined `local a=$1 b=$2 c=...$b` expands
    # $b before it's assigned, which trips `set -u` (unbound variable).
    local truth="$1"
    local label="$2"
    local comp="$3"
    shift 3
    local out="$OUTDIR/truvari_$label"
    local csample; csample=$(comp_sample "$comp")
    rm -rf "$out"   # truvari requires the output dir not pre-exist
    truvari bench -b "$truth" -c "$comp" -f "$REFERENCE" -o "$out" \
        --passonly --sizemax "$TRUVARI_SIZEMAX" \
        ${csample:+--cSample "$csample"} "$@" 2>&1 | tail -3 || true
    if [[ -f "$out/summary.json" ]]; then
        python3 - "$out/summary.json" "$label" <<'PY'
import json, sys
s = json.load(open(sys.argv[1])); name = sys.argv[2]
# truvari writes null (None) for metrics when a subset has no comparisons;
# `or 0` coerces None -> 0 (s.get(k, 0) only defaults a MISSING key, not a null one).
print(f"  {name:<18} recall={(s.get('recall') or 0):.3f}  precision={(s.get('precision') or 0):.3f}  "
      f"f1={(s.get('f1') or 0):.3f}  TP={s.get('TP-base') or 0}  FN={s.get('FN') or 0}  FP={s.get('FP') or 0}")
PY
        check_denominator "$truth" "$out" "$label" || true
    else
        echo "  $label: truvari produced no summary.json (see $out)"
        SCORE_COVERAGE_FAIL=$((SCORE_COVERAGE_FAIL + 1))
    fi
}

# Coverage check: truvari's TP-base+FN is its POST-filter base count. If that is less than
# the number of records in the truth we handed it, records were discarded before comparison
# and the reported recall was computed over a SUBSET — a wrong denominator, silently.
#
# This checks the OUTCOME rather than enumerating causes, so it generalises past any one
# filter. An audit of the archived runs found 4 such cases in 152 per-type checks (a DEL in
# job_19758751 reported 0.800 where the full denominator gives 0.727; an INV in
# job_20636298 reported 0.800 vs 0.667) and nothing in the output said so.
#
# It also ATTRIBUTES the loss: every base record truvari considered lands in tp-base.vcf.gz
# or fn.vcf.gz, so truth-minus-those is exactly the discarded set, printed with SVTYPE/SVLEN
# so the responsible filter is obvious. eidolon's own floor is MIN_SV_LENGTH_BP=50, equal to
# truvari's --sizemin default, so a size-driven drop here means the LARGE end (--sizemax).
check_denominator() {  # <truth.vcf.gz> <truvari_out_dir> <label>
    local truth="$1" out="$2" label="$3"
    local n_truth base
    n_truth=$(bcftools view -H "$truth" 2>/dev/null | wc -l)
    base=$(python3 -c "import json;d=json.load(open('$out/summary.json'));print((d.get('TP-base') or 0)+(d.get('FN') or 0))" 2>/dev/null || echo 0)
    [[ "$n_truth" -eq "$base" ]] && return 0

    echo "  ERROR: $label scored only $base of $n_truth truth record(s); $((n_truth - base)) discarded" >&2
    echo "    before comparison, so the recall above used the WRONG DENOMINATOR." >&2
    if [[ "$base" -gt 0 ]]; then
        python3 -c "
import json
d = json.load(open('$out/summary.json')); tp = d.get('TP-base') or 0
print(f'    reported {tp}/$base={tp/$base:.3f}  ->  over the full truth {tp}/$n_truth={tp/$n_truth:.3f}')" >&2
    fi

    # Attribute it: every base record truvari considered lands in tp-base or fn, so
    # truth-minus-those is exactly what was discarded. awk-join (not join(1)) to avoid
    # depending on key collation matching between two independently sorted streams.
    local considered="$OUTDIR/.considered_$label" f found=0
    : > "$considered"
    for f in "$out/tp-base.vcf.gz" "$out/fn.vcf.gz"; do
        if [[ -f "$f" ]]; then
            bcftools query -f '%CHROM\t%POS\n' "$f" >> "$considered" 2>/dev/null && found=1
        fi
    done
    if [[ "$found" -eq 1 ]]; then
        echo "    discarded records (CHROM POS SVTYPE SVLEN); --sizemax=$TRUVARI_SIZEMAX:" >&2
        bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$truth" 2>/dev/null \
          | awk -F'\t' -v seen="$considered" '
              BEGIN { while ((getline line < seen) > 0) { split(line, k, "\t"); c[k[1] "~" k[2]] = 1 } }
              !((($1 "~" $2) in c)) { n++; if (n <= 10) printf "      %s %s %s %s\n", $1, $2, $3, $4 }
              END { if (n > 10) printf "      ... and %d more\n", n - 10 }' >&2
    else
        echo "    (no tp-base/fn VCF in $out — cannot attribute; inspect it directly)" >&2
    fi
    rm -f "$considered"

    SCORE_COVERAGE_FAIL=$((SCORE_COVERAGE_FAIL + 1))
    return 1
}

# ── BND junction spans: a comparable representation ──────────────────────
# A BND pair at (A,B) and a caller's <DUP:TANDEM>/<DEL>/<INV> spanning A..B describe the
# SAME event, but truvari matches variants as INTERVALS and these are not
# interval-equivalent: the BND record is a point (END==POS) while the caller's is a span
# of tens of Mb. No similarity threshold or --typeignore can bridge that — verified on
# job 20636298, where Manta called a truth junction as <DUP:TANDEM> at POS 1bp away with
# END exactly on the mate, and truvari still reported TP=0 with --typeignore --pctsize 0
# --pctseq 0 --sizemax 100000000.
#
# So convert the truth: emit ONE record per junction, POS=min, END=max. That is a derived
# COMPARISON artifact, not an output format — END on a real BND must equal POS. Only
# intra-contig junctions can become spans; inter-contig ones are counted and skipped.
build_bnd_spans() {  # <bnd_truth.vcf.gz> <out.vcf.gz>
    # Exit status is three-way on purpose, because "no spans" has two very different
    # meanings and the caller must not treat them alike:
    #   0 = spans written
    #   2 = nothing to write because every junction is inter-contig (legitimate)
    #   1 = anything else (a defect; the previous version reported this as case 2)
    local src="$1" dst="$2"
    [[ -f "$src" ]] || { echo "  ERROR: BND spans: no BND truth at $src" >&2; return 1; }
    local body="$OUTDIR/.bnd_span_body" stats="$OUTDIR/.bnd_span_stats"
    local hdr="$OUTDIR/.bnd_span_hdr" err="$OUTDIR/.bnd_span_err"
    rm -f "$body" "$stats" "$hdr" "$err"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$src" \
      | awk -F'\t' -v stats="$stats" 'BEGIN{OFS="\t"}
          { alt=$4; gsub(/[][]/, "|", alt); n=split(alt, p, "|"); loc=""
            for (i=1;i<=n;i++) if (p[i] ~ /:/) { loc=p[i]; break }
            if (loc == "") { bad++; next }
            split(loc, m, ":"); mc=m[1]; mp=m[2]+0
            if (mc != $1) { inter++; next }
            lo = ($2+0 < mp) ? $2+0 : mp; hi = ($2+0 < mp) ? mp : $2+0
            key = $1 ":" lo "-" hi
            if (key in seen) { dup++; next }
            seen[key]=1
            print $1, lo, ".", $3, "<BND_SPAN>", 37, "PASS", \
                  "SVTYPE=BND;END=" hi ";SVLEN=" (hi-lo), "GT", "0/1" }
          END{ print length(seen), dup+0, inter+0, bad+0 > stats
               printf("  BND spans: %d junction(s) from %d record(s) "\
                      "(inter-contig skipped=%d, unparsable ALT=%d)\n",
                      length(seen), length(seen)+dup+inter+bad, inter+0, bad+0) > "/dev/stderr" }' \
      > "$body"

    local n_junc dup inter bad
    if ! read -r n_junc dup inter bad < "$stats" 2>/dev/null; then
        echo "  ERROR: BND spans: derivation produced no statistics — awk did not run." >&2
        rm -f "$body" "$stats"; return 1
    fi
    rm -f "$stats"
    if [[ "$n_junc" -eq 0 ]]; then
        rm -f "$body"
        if [[ "$bad" -eq 0 && "$inter" -gt 0 ]]; then
            echo "  BND spans: none derived — all $inter junction record(s) are inter-contig," >&2
            echo "    which cannot be expressed as a span. Expected on a multi-contig reference." >&2
            return 2
        fi
        echo "  ERROR: BND spans: 0 junctions derived from a non-empty BND truth" >&2
        echo "    (unparsable ALT=$bad, inter-contig=$inter). The ALT parser or the truth is wrong." >&2
        return 1
    fi

    # Inherit the source header, then backfill ONLY what is missing. The body writes
    # SVTYPE/END/SVLEN and a GT column; an undeclared END is parsed as Type=String and
    # kills bcftools outright, which the old `2>/dev/null` hid while still returning 0.
    bcftools view -h "$src" | grep -vE '^#CHROM' > "$hdr"
    local decl id kind
    for decl in '##ALT=<ID=BND_SPAN,Description="Derived: intra-contig BND junction as a span (comparison only)">' \
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' \
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">' \
                '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' \
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'; do
        id="${decl#*<ID=}"; id="${id%%,*}"
        kind="${decl%%=*}"
        grep -q "^${kind}=<ID=${id}," "$hdr" || echo "$decl" >> "$hdr"
    done

    local rc=0
    { cat "$hdr"
      bcftools view -h "$src" | grep -E '^#CHROM'
      cat "$body"; } | bcftools sort -O z -o "$dst" 2> "$err" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        echo "  ERROR: BND spans: bcftools sort failed (rc=$rc):" >&2
        sed 's/^/      /' "$err" >&2
        rm -f "$body" "$hdr" "$err" "$dst"; return 1
    fi
    bcftools index -f -t "$dst" 2>/dev/null || true

    # Coverage of its own input: the file must hold every junction awk counted.
    local n_written; n_written=$(bcftools view -H "$dst" 2>/dev/null | wc -l)
    if [[ "$n_written" -ne "$n_junc" ]]; then
        echo "  ERROR: BND spans: derived $n_junc junction(s) but wrote $n_written record(s)." >&2
        sed 's/^/      /' "$err" >&2
        rm -f "$body" "$hdr" "$err"; return 1
    fi
    echo "  BND spans: $n_written junction span(s) -> $(basename "$dst")"
    rm -f "$body" "$hdr" "$err"
    return 0
}

# ── Rule 0: calibrate the scorer before believing it ────────────────────
# A scorer must reproduce a PERFECT score against itself in the configuration the real
# run uses. If truth-vs-truth is not 1.000, the configuration is broken and nothing it
# reports about a caller means anything. This check, run once, would have immediately
# localised the BND problem that instead took a week and five wrong explanations:
# `truvari bench -b truth -c truth --sizemax 5000000` cannot match a 29Mb junction to
# itself, because sizemax discards it before comparison.
selftest_truvari() {  # <truth.vcf.gz> <label> [extra args...]
    local truth="$1" label="$2"; shift 2
    local out="$OUTDIR/truvari_selftest_$label"
    rm -rf "$out"
    truvari bench -b "$truth" -c "$truth" -f "$REFERENCE" -o "$out" \
        --passonly --sizemax "$TRUVARI_SIZEMAX" "$@" >/dev/null 2>&1 || true
    local r=0
    [[ -f "$out/summary.json" ]] && r=$(python3 -c "import json,sys;print(json.load(open('$out/summary.json')).get('recall') or 0)")
    if awk -v r="$r" 'BEGIN{exit !(r > 0.999)}'; then
        # recall=1.000 is necessary but NOT sufficient: a filter that discards a record
        # removes it from base AND comp alike, so the subset still matches itself
        # perfectly. Only the denominator check sees that.
        check_denominator "$truth" "$out" "selftest_$label" || return 1
        echo "  selftest $label: PASS (truth vs truth recall=$r, all $(bcftools view -H "$truth" | wc -l) records scored)"
        return 0
    fi
    echo "ERROR: scorer selftest FAILED for $label — truth vs truth recall=$r, expected 1.000." >&2
    echo "  The scoring configuration cannot match the truth to ITSELF, so any recall it" >&2
    echo "  reports against a caller is meaningless. Usual causes: --sizemax discarding" >&2
    echo "  records before comparison, or a size/sequence similarity threshold that a" >&2
    echo "  point-like record cannot satisfy. See $out." >&2
    return 1
}

# Negative control: the same truth shifted far enough that nothing should match. Guards
# against a configuration so permissive it matches anything.
decoy_truvari() {  # <truth.vcf.gz> <label> [extra args...]
    local truth="$1" label="$2"; shift 2
    local decoy="$OUTDIR/.decoy_${label}.vcf.gz" out="$OUTDIR/truvari_decoy_$label"
    { bcftools view -h "$truth"
      bcftools view -H "$truth" | awk -F'\t' 'BEGIN{OFS="\t"}{$2=$2+2000000; print}'; } \
      | bcftools sort -O z -o "$decoy" 2>/dev/null || { echo "  decoy $label: could not build (skipped)"; return 0; }
    bcftools index -f -t "$decoy" 2>/dev/null || true
    rm -rf "$out"
    truvari bench -b "$truth" -c "$decoy" -f "$REFERENCE" -o "$out" \
        --passonly --sizemax "$TRUVARI_SIZEMAX" "$@" >/dev/null 2>&1 || true
    local r=0
    [[ -f "$out/summary.json" ]] && r=$(python3 -c "import json,sys;print(json.load(open('$out/summary.json')).get('recall') or 0)")
    if awk -v r="$r" 'BEGIN{exit !(r < 0.001)}'; then
        echo "  decoy    $label: PASS (shifted truth recall=$r, expected 0)"
    else
        echo "  WARNING: decoy control for $label matched (recall=$r) — the configuration" >&2
        echo "  is matching on position far more loosely than intended. See $out." >&2
    fi
    return 0
}

# Score one caller's VCF against the truth — aggregate + per-type + BND (#457).
#
# What this fixes, all of which produced numbers that were not measurements:
#
#  1. Per-type runs now filter the QUERY as well as the truth. Passing the full caller
#     VCF against a type-filtered truth counted every record of another type as a false
#     positive, giving arithmetically impossible output — manta_BND FP=39 against
#     manta_overall FP=20. A subset cannot have more false positives than the whole.
#  2. BND is scored TWO ways, because callers represent junctions two ways: BND-to-BND
#     via --bnddist for callers that emit breakends, and span-to-span against the derived
#     junction spans for callers that re-represent a junction as <DUP:TANDEM>/<DEL>/<INV>.
#     Reported separately and labelled; the span comparison is what recovers Manta.
#  3. CNV stays out of truvari entirely — it is a copy-ratio call, scored by the stage-4c
#     direction check instead. The aggregate covers types truvari can match, and the log
#     states what was excluded and where it went.
#  4. Every scorer is calibrated first (selftest + decoy). See selftest_truvari.
score_caller() {  # <caller> <comp.vcf.gz>
    local caller="$1" comp="$2"
    [[ -f "$comp" ]] || { echo "  ($caller: no calls to score)"; return 0; }
    echo "  === $caller (cSample=$(comp_sample "$comp")) ==="

    if [[ -f "$TRUTH_SCOREABLE" ]]; then
        # The aggregate truth excludes BND and CNV, so the query must exclude them too.
        # Otherwise every BND/CNV call the caller made is a false positive against a
        # truth that never contained one — job 20675480 reported manta_overall FP=24
        # from a 24-record truth, with Manta's BND calls making up the bulk. `-e`
        # (exclude) rather than `-i`, so a record with no SVTYPE at all is kept instead
        # of silently dropped.
        local qs="$OUTDIR/${caller}_calls_scoreable.vcf.gz"
        if ! bcftools view -e 'INFO/SVTYPE="BND" || INFO/SVTYPE="CNV"' \
                -O z -o "$qs" "$comp" 2>/dev/null; then
            echo "  ERROR: $caller: could not build the aggregate query subset —" >&2
            echo "    skipping the aggregate rather than scoring an unfiltered query." >&2
            SCORE_COVERAGE_FAIL=$((SCORE_COVERAGE_FAIL + 1))
        else
            bcftools index -f -t "$qs" 2>/dev/null || true
            local n_all n_sub
            n_all=$(bcftools view -H "$comp" 2>/dev/null | wc -l)
            n_sub=$(bcftools view -H "$qs" 2>/dev/null | wc -l)
            echo "  (aggregate query: $n_sub of $n_all $caller call(s); $((n_all - n_sub)) BND/CNV call(s) excluded, scored separately)"
            run_truvari "$TRUTH_SCOREABLE" "${caller}_overall" "$qs"
        fi
    else
        echo "  ($caller: no truvari-scoreable truth types — aggregate skipped)"
    fi

    for typed in "$OUTDIR"/truth_sv_*.vcf.gz; do
        [[ -f "$typed" ]] || continue
        local svt; svt=$(basename "$typed" .vcf.gz); svt=${svt#truth_sv_}
        case "$svt" in
            BND|BNDspan|CNV|scoreable) continue ;;   # handled separately below / 4c
        esac
        local q="$OUTDIR/${caller}_calls_${svt}.vcf.gz"
        bcftools view -i "INFO/SVTYPE=\"$svt\"" -O z -o "$q" "$comp" 2>/dev/null || true
        if [[ -f "$q" ]] && [[ "$(bcftools view -H "$q" 2>/dev/null | wc -l)" -gt 0 ]]; then
            bcftools index -f -t "$q" 2>/dev/null || true
            run_truvari "$typed" "${caller}_${svt}" "$q"
        else
            rm -f "$q" "${q}.tbi"
            # Recall is 0 by definition; say so rather than emit a precision figure
            # computed from an empty query.
            echo "  ${caller}_${svt} recall=0.000  (caller emitted no $svt records)"
        fi
    done

    # BND, representation 1: breakend-to-breakend. The query is the caller's BND
    # records ONLY — comparing breakends to breakends. Passing the whole call set here
    # made every non-BND call a false positive: job 20675480 reported manta_BND FP=43
    # against manta_overall FP=24, and a subset cannot have more false positives than
    # the whole. That was #457's defect surviving in the one run it never covered.
    if [[ -f "$OUTDIR/truth_sv_BND.vcf.gz" ]]; then
        local qb="$OUTDIR/${caller}_calls_BND.vcf.gz"
        bcftools view -i 'INFO/SVTYPE="BND"' -O z -o "$qb" "$comp" 2>/dev/null || true
        if [[ -f "$qb" ]] && [[ "$(bcftools view -H "$qb" 2>/dev/null | wc -l)" -gt 0 ]]; then
            bcftools index -f -t "$qb" 2>/dev/null || true
            run_truvari "$OUTDIR/truth_sv_BND.vcf.gz" "${caller}_BND" "$qb" \
                --bnddist "$BND_DIST"
        else
            rm -f "$qb" "${qb}.tbi"
            echo "  ${caller}_BND recall=0.000  (caller emitted no BND records)"
        fi
    fi
    # BND, representation 2: junction spans vs whatever type the caller used. -t
    # (--typeignore) is required: the caller legitimately calls a junction a DUP.
    #
    # The query is deliberately UNFILTERED here — the whole point is that the caller may
    # have used any type for the junction, so it cannot be narrowed by type. That makes
    # recall the measurement and FP/precision uninterpretable (every unrelated call is
    # an unmatched comp record), so say so rather than let the number be quoted.
    if [[ -f "$OUTDIR/truth_sv_BNDspan.vcf.gz" && "$BNDSPAN_VALID" -eq 1 ]]; then
        run_truvari "$OUTDIR/truth_sv_BNDspan.vcf.gz" "${caller}_BNDspan" "$comp" -t
        echo "    ^ BNDspan: recall is the measurement; precision/FP/f1 are NOT — the query"
        echo "      is intentionally unfiltered because the caller may use any SVTYPE here."
    elif [[ -f "$OUTDIR/truth_sv_BNDspan.vcf.gz" ]]; then
        echo "  ${caller}_BNDspan: not scored (representation bridge failed calibration)"
    fi
}

# ── Rule 0: calibrate every scorer before believing any of it ────────────
# Derive comparable junction spans from the BND pairs (see build_bnd_spans). This runs
# HERE, not in stage 3b, because bash resolves function names at call time: the call used
# to sit ~250 lines above the definition, so it was `command not found` (127) on every
# run and the `||` branch reported it as "no intra-contig junctions". The BNDspan scoring
# path — the one that recovers a caller which re-represents a junction — therefore never
# executed. Handle every status explicitly; an unexpected rc must not read as success.
# Initialised near the top (before stage 4a, which can set it).
BND_SPAN_RC=0
# Set to 0 by probe_bndspan.sh when truvari cannot distinguish a correct junction call
# from a displaced one at this scale — see the probe for why that is not tunable.
BNDSPAN_VALID=1
if [[ -f "$OUTDIR/truth_sv_BND.vcf.gz" ]]; then
    build_bnd_spans "$OUTDIR/truth_sv_BND.vcf.gz" "$OUTDIR/truth_sv_BNDspan.vcf.gz" \
        || BND_SPAN_RC=$?
    case "$BND_SPAN_RC" in
        0) ;;                                   # spans written
        2) echo "  (BNDspan scoring skipped: no intra-contig junctions to span)" ;;
        *) echo "ERROR: BND span derivation failed (rc=$BND_SPAN_RC) — BNDspan scoring" >&2
           echo "  cannot run, so a caller that re-represents junctions will score 0" >&2
           echo "  recall for reasons that are ours, not the caller's." >&2
           CAL_FAIL=1 ;;
    esac
fi

echo "[5/5] Calibrating scorers (truth-vs-truth must be perfect; decoy must be zero)..."
if [[ -f "$TRUTH_SCOREABLE" ]]; then
    selftest_truvari "$TRUTH_SCOREABLE" scoreable || CAL_FAIL=1
    decoy_truvari    "$TRUTH_SCOREABLE" scoreable
fi
if [[ -f "$OUTDIR/truth_sv_BND.vcf.gz" ]]; then
    selftest_truvari "$OUTDIR/truth_sv_BND.vcf.gz" BND --bnddist "$BND_DIST" || CAL_FAIL=1
fi
if [[ -f "$OUTDIR/truth_sv_BNDspan.vcf.gz" ]]; then
    selftest_truvari "$OUTDIR/truth_sv_BNDspan.vcf.gz" BNDspan -t || CAL_FAIL=1
    # The selftest above compares spans to THEMSELVES, which matches trivially. It says
    # nothing about the premise BNDspan actually rests on: that truvari -t matches a
    # derived <BND_SPAN> to the <DUP:TANDEM> a caller emits for the same junction. If
    # that does not hold, BNDspan recall measures our comparison, not the caller — the
    # exact error this line of work began with. Known-answer, both directions.
    "$REPO_ROOT/scripts/delta/probe_bndspan.sh" \
        --spans "$OUTDIR/truth_sv_BNDspan.vcf.gz" \
        --reference "$REFERENCE" --work "$OUTDIR/.bndspan_probe" || BNDSPAN_VALID=0
    # A failed bridge disables BNDspan ONLY. It says nothing about DEL/DUP/INV scoring,
    # so failing the whole calibration here would throw away sound measurements.
    if [[ "$BNDSPAN_VALID" -eq 0 ]]; then
        echo "  BNDspan scoring DISABLED — the probe above shows this configuration cannot" >&2
        echo "  tell a correct junction call from a displaced one, so any recall it produced" >&2
        echo "  would not be a measurement. DEL/DUP/INV and BND-to-BND are unaffected." >&2
    fi
fi
# Skip scoring rather than exit: the calibration artifacts and the truth VCFs are worth
# archiving precisely BECAUSE the scorer misbehaved, and the exit below runs after
# archive_run. Scoring itself is pointless with a scorer that cannot score a known answer.
if [[ "$CAL_FAIL" -ne 0 ]]; then
    echo "ERROR: at least one scorer cannot match the truth to itself — refusing to" >&2
    echo "  report caller results computed with a configuration that demonstrably" >&2
    echo "  cannot score a perfect answer. Skipping caller scoring; see the selftest" >&2
    echo "  output above and the archived truvari_selftest_* directories." >&2
else
    score_caller manta "$MANTA_SV"
    [[ -n "$DELLY_SV" ]] && score_caller delly "$DELLY_SV"
fi

# ── Summary ──────────────────────────────────────────────────────────────
# Never print an unqualified "complete" over a run whose numbers are not usable — a
# summary banner that reads as success regardless of outcome is how a bad run gets
# believed. The detail follows in the deferred failure gate below.
if [[ "$CAL_FAIL" -ne 0 ]]; then
    RUN_STATUS="FAILED (scorer calibration) — no caller results"
elif [[ "$SCORE_COVERAGE_FAIL" -gt 0 ]]; then
    RUN_STATUS="FAILED ($SCORE_COVERAGE_FAIL comparison(s) scored a SUBSET of the truth)"
else
    RUN_STATUS="complete"
fi
cat <<EOF

════════════════════════════════════════════════════════════════
SV validation pipeline $RUN_STATUS — Job $SLURM_JOB_ID

Inputs:
  Reference:    $REFERENCE
  Coverage:     ${TOTAL_COVERAGE}x  Purity: $PURITY  SV scale: $SV_RATE_SCALE
  Tumor model:  $TUMOR_MODEL

Key outputs in $OUTDIR:
  SV truth:     truth_sv.vcf.gz            ($SV_TRUTH_N somatic SVs)
  Manta calls:  manta/results/variants/somaticSV.vcf.gz   (raw: $MANTA_SV_RAW)
  Manta scored: $MANTA_SV
                (differs from raw when inversions were converted from BND pairs)
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

# ── Deferred failure gate ────────────────────────────────────────────────
# Deliberately last. Every discrepancy is reported inline next to the comparison it
# affects (so the scale is visible rather than just the first case), the artifacts are
# archived (so a failure is diagnosable after the fact), and only then does the job
# exit non-zero — because a run that published recalls over a subset of the truth, or
# scored with an uncalibrated scorer, must not look successful to whoever reads the
# SLURM exit status.
if [[ "$CAL_FAIL" -ne 0 || "$SCORE_COVERAGE_FAIL" -gt 0 ]]; then
    echo "" >&2
    [[ "$CAL_FAIL" -ne 0 ]] && \
        echo "FAILED: scorer calibration (Rule 0) — no caller results were produced." >&2
    if [[ "$SCORE_COVERAGE_FAIL" -gt 0 ]]; then
        echo "FAILED: $SCORE_COVERAGE_FAIL comparison(s) did not score the whole truth set." >&2
        echo "  Those recalls are over a subset and must not be reported as-is. Raise" >&2
        echo "  TRUVARI_SIZEMAX (currently $TRUVARI_SIZEMAX) or relax the size filters," >&2
        echo "  then re-score. Artifacts are archived." >&2
    fi
    exit 1
fi
