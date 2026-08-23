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
# 128, not 64: --mem=240G on a 256 GB node already reserves the WHOLE node, so asking
# for half its cores left 64 idle at zero scheduling benefit. Measured on job 20719077
# (soy whole assembly, full pipeline): 2h54m wall, 14.13% CPU efficiency at 64 cores —
# most time is in serial phases, so this speeds up bwa-mem2 and the callers rather than
# the whole run. Do NOT lower --mem to schedule more easily: that same job peaked at
# 91.78 GB on a 1 Gb genome, and GRCh38 is ~3x, so 240G is closer to necessary than
# padded. (Job 20745149's 39.86 GB is NOT a counter-example — that run re-scored existing
# artifacts and never simulated or aligned.)
#SBATCH --cpus-per-task=128
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
# Say WHICH build produced these numbers, in the log, every time. This pipeline never
# rebuilds — it runs whatever is at $EIDOLON_BIN — and on a fresh OUTDIR there is no
# stamp to compare against, so a stale binary runs silently and its results are
# indistinguishable from a current one. Job 20682989 reproduced job 20675480's caller
# numbers exactly and the log could not answer "was the fix even in this build?".
echo "eidolon binary: $EIDOLON_BIN"
echo "eidolon version: ${EIDOLON_VERSION:-unknown}"
echo "repo HEAD:      $(git -C "$REPO_ROOT" describe --tags --always --dirty 2>/dev/null || echo unknown)"
# A binary older than the checkout is almost always a forgotten rebuild.
#
# The SEMVER half of this check cannot see a stale binary within a version, which is every
# build in this repo (#513). `develop` carries 116 commits since v3.0.0 under one Cargo
# version, and that version reads 3.1.0 against main's 3.0.0 with no release behind it — it
# was bumped partly to make this very check fire. So the commit is compared too, from the
# stamp `eidolon/build.rs` puts in `--version` ("3.1.0+abc1234", "+abc1234-dirty", or
# "+unknown"). Job 20682989 reproduced job 20675480's numbers exactly and the log could not
# answer "was the fix even in this build?"; the semver alone still could not.
check_binary_provenance() {  # <eidolon --version output> <repo root>
    local ver="$1" root="$2" rc=0
    # dirty MUST be initialised, not merely declared: `local dirty` leaves it unset, and
    # "${repo_sha}${dirty}" then aborts the whole job under `set -u` on a CLEAN checkout —
    # i.e. the normal case. Caught by the suite before it reached Delta.
    local repo_ver repo_sha want
    local dirty=""

    repo_ver="$(awk -F\' '/^version = /{print $2; exit}' "$root/Cargo.toml" 2>/dev/null)"
    if [[ -n "$repo_ver" && "$ver" != *"$repo_ver"* ]]; then
        echo "ERROR: binary reports '${ver:-unknown}' but the checkout is $repo_ver." >&2
        echo "  The binary is stale — rebuild before running, or the results describe code" >&2
        echo "  that is not the code you think you are testing." >&2
        rc=1
    fi

    repo_sha="$(git -C "$root" rev-parse --short=7 HEAD 2>/dev/null)"
    # -uno: the pipeline's own SLURM .out/.err files land in the submit dir, which IS the
    # repo, so counting untracked files would report every run as dirty.
    [[ -n "$(git -C "$root" status --porcelain --untracked-files=no 2>/dev/null)" ]] && dirty="-dirty"
    want="${repo_sha}${dirty}"

    if [[ -z "$repo_sha" ]]; then
        echo "WARNING: cannot read the checkout's git commit; build provenance unverified." >&2
    elif [[ "$ver" == *"+unknown"* ]]; then
        # "Cannot verify" must not read as "verified" — the #509 lesson, one guard over.
        echo "ERROR: the binary carries no git stamp (+unknown), so it cannot be matched to" >&2
        echo "  this checkout. Rebuild where git is available." >&2
        rc=1
    elif [[ "$ver" != *"+$want"* ]]; then
        echo "ERROR: the binary was built from a DIFFERENT COMMIT than this checkout — the" >&2
        echo "  forgotten rebuild the version check above cannot see." >&2
        echo "    binary:   $ver" >&2
        echo "    checkout: $want" >&2
        echo "  Run 'cargo build --release' before scoring." >&2
        rc=1
    fi

    [[ "$rc" -eq 0 ]] && return 0
    echo "  Set ALLOW_STALE_BIN=1 to override (the results then describe an unknown build)." >&2
    [[ "${ALLOW_STALE_BIN:-0}" == "1" ]]
}
check_binary_provenance "${EIDOLON_VERSION:-}" "$REPO_ROOT" || exit 1

# Disk state, in every log, before any work. Job 20904141 rep 5 died ~3h in when the
# project quota filled at the FASTQ merge, having already spent simulation and alignment —
# and the log gave no way to see how close the preceding replicates had come. Printed RAW
# and never parsed: `df` misreports inside a project-quota'd directory (it shows the quota,
# not the mount), the quota table and df have disagreed for /scratch, and a guard built on
# fragile parsing is worse than no guard.
#
# Peak footprint of ONE replicate, ANCHORED on the only full measurement available: ~214 GB for
# GRCh38 (3.1 Gbp) at 30x — FASTQ 116 GB (merged 29+29, tumor 17+17, normal 12+12, itemised from
# job 20884022) plus BAMs 98 GB (the figure prune_bams itself reported across campaign 20925151),
# coexisting because the prune only runs at the end. An earlier "113 GB" figure counted the FASTQ
# alone and was half of reality.
#
# SCALED, not constant. Both FASTQ and BAM grow with genome size x coverage, so a fixed 214 GB
# demanded 61x more than a chr22 replicate actually needs (~3.5 GB) and would have refused every
# smoke run once /scratch passed 336 GB used — breaking the smoke-first workflow exactly when
# disk is tight and a cheap run matters most. Job 21035422 passed only because 343 GB happened to
# be free.
scaled_replicate_peak_gb() {  # <reference_bp> <coverage>; echoes whole GB
    awk -v bp="$1" -v cov="$2" 'BEGIN {
        # Cannot scale -> assume the worst rather than wave a run through. "Unknown" must not
        # read as "small", which is the same rule as the calibration controls.
        if (bp <= 0 || cov <= 0) { print 214; exit }
        gb = 214.0 * (bp / 3100000000.0) * (cov / 30.0)
        # Floor: below this the reference index, Manta workspace and temp files dominate and the
        # read-proportional model understates the real need.
        if (gb < 5) gb = 5
        printf "%d", (gb == int(gb) ? gb : int(gb) + 1)
    }'
}

# Reference size in bases. Prefers the .fai (exact); falls back to the FASTA byte count, which is
# a fine proxy for a capacity gate (~60 bases per 61 bytes with newlines) and — unlike the .fai —
# always exists, since the index is not built until stage 2.
reference_bp() {  # <reference.fa>
    local ref="$1" bytes
    if [[ -s "$ref.fai" ]]; then
        awk '{s += $2} END {print s + 0}' "$ref.fai"
        return
    fi
    bytes=$(stat -c %s "$ref" 2>/dev/null || echo 0)
    awk -v b="$bytes" 'BEGIN {printf "%d", b * 0.98}'
}

REPLICATE_PEAK_GB="${REPLICATE_PEAK_GB:-$(scaled_replicate_peak_gb \
    "$(reference_bp "$REFERENCE")" "$TOTAL_COVERAGE")}"

# Extract used/limit in KB from `lfs quota -p` output. Separated from the check so the
# decision logic is testable without Lustre: the parsing is the fragile half, and the real
# output is not guessable — an over-quota row marks values with a trailing '*' and adds a
# grace column, which naive field-splitting silently mangles.
#   Filesystem   used  bquota  blimit  bgrace       files iquota ilimit igrace
#     /scratch 549.9G*   500G    550G  6d6h28m43s  21516* 850000 935000      -
# Called WITHOUT -h so the values are plain KB and no unit parsing is needed.
parse_lfs_quota_kb() {  # <lfs-quota-output> <filesystem-path>; echoes "<used_kb> <limit_kb>"
    awk -v fs="$2" '
        $1 == fs {
            used = $2; hard = $4
            gsub(/\*/, "", used); gsub(/\*/, "", hard)
            # blimit (hard) is what actually stops a write; 0 means "no hard limit set",
            # in which case fall back to bquota (soft).
            if (hard + 0 == 0) { soft = $3; gsub(/\*/, "", soft); hard = soft }
            if (used ~ /^[0-9]+$/ && hard ~ /^[0-9]+$/) { print used, hard; found = 1 }
        }
        END { exit !found }
    ' <<< "$1"
}

# Refuse to START a replicate that cannot fit. Task 4 of array 20904141 burned 8 h 45 m
# against a full filesystem to produce an 8 KB output directory, and tasks 5-8 then failed
# in turn; the log's only clue was `cat: write error: Disk quota exceeded` on stderr. The
# footprint is known and the quota is one command, so this is answerable in a second, up
# front. Advisory when the quota cannot be read — a check that cannot run must say so, but
# must not invent a reason to stop a run that would have been fine.
echo "── disk pre-flight (one replicate needs ~${REPLICATE_PEAK_GB} GB peak) ──"
_scratch_root="$(dirname "$OUTDIR")"
_projid="$(lfs project -d "$_scratch_root" 2>/dev/null | awk '{print $1}')" || true
if [[ -z "${_projid:-}" || "$_projid" == "0" ]]; then
    echo "  WARNING: no Lustre project id on $_scratch_root — free space UNVERIFIED." >&2
else
    _fsmount="$(df -P "$_scratch_root" 2>/dev/null | awk 'NR==2{print $6}')"
    lfs quota -h -p "$_projid" "$_fsmount" 2>&1 || true
    _raw="$(lfs quota -p "$_projid" "$_fsmount" 2>/dev/null)" || true
    if _kb="$(parse_lfs_quota_kb "${_raw:-}" "$_fsmount")"; then
        _free_gb=$(( ( ${_kb% *} > ${_kb#* } ? 0 : ${_kb#* } - ${_kb% *} ) / 1048576 ))
        echo "  free before start: ${_free_gb} GB (need ~${REPLICATE_PEAK_GB} GB)"
        if [[ "$_free_gb" -lt "$REPLICATE_PEAK_GB" ]]; then
            echo "ERROR: only ${_free_gb} GB free; one replicate needs ~${REPLICATE_PEAK_GB} GB." >&2
            echo "  Refusing to start rather than fail hours in with a partial output dir" >&2
            echo "  that itself holds ~207 GB and starves the next replicate (job 20884022)." >&2
            echo "  Largest consumers:" >&2
            du -sh "$_scratch_root"/* 2>/dev/null | sort -rh | head -5 >&2
            echo "  A FAILED replicate keeps its FASTQ, so clear those directories first." >&2
            echo "  Set SKIP_DISK_PREFLIGHT=1 to override." >&2
            [[ "${SKIP_DISK_PREFLIGHT:-0}" == "1" ]] || exit 1
        fi
    else
        echo "  WARNING: could not parse the quota output — free space UNVERIFIED." >&2
    fi
fi

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

# Gate the truth BEFORE stage 3b inspects it or any caller is scored against it. The
# truth VCF is the denominator of every number this pipeline reports, so a malformed one
# does not produce an error — it produces a plausible wrong answer.
validate_artifacts "$EIDOLON_BIN" "simulated truth" \
    "$TRUTH_VCF" "${PREFIX}_normal.vcf.gz" "${PREFIX}_tumor.vcf.gz" || exit 1

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

# The aligned BAMs are what Manta, Delly and GATK read. `samtools quickcheck` would
# catch a truncated one; nothing in this pipeline was running it.
validate_artifacts "$EIDOLON_BIN" "aligned BAMs" "$NORMAL_BAM" "$TUMOR_BAM" || exit 1

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
# env). Empty types are removed so the scoring loop skips them — and RECORDED, because a
# skipped type is an unmeasured stratum that has to be reported rather than dropped.
# The canonical SVTYPE list, declared ONCE. The scoring loop iterates this same array
# rather than globbing truth_sv_*.vcf.gz, so a derived artifact can never be mistaken
# for an SVTYPE — see SV_TYPES_SCORED_IN_LOOP below.
SV_TYPES=(DEL DUP INV INS BND CNV)
# Types the generic scoring loop handles. BND is scored by breakpoint proximity and CNV
# by the stage-4c direction check, so both are excluded here and handled explicitly.
# Derived from SV_TYPES so a newly added type is scored rather than silently dropped.
SV_TYPES_SCORED_IN_LOOP=()
for svt in "${SV_TYPES[@]}"; do
    case "$svt" in BND|CNV) continue ;; esac
    SV_TYPES_SCORED_IN_LOOP+=("$svt")
done
# A type the model ASKS FOR but that no record realises is an UNMEASURED STRATUM, not a
# result. The zero-count file is deleted here and the scoring loop then skips it on
# `[[ -f "$typed" ]] || continue`, so before this change the run printed a bare
# `truth INS: 0` and never mentioned INS again — and the aggregate recall silently meant
# "of the types we happened to draw". Job 20853511 reported manta_overall=0.861 over
# DEL/DUP/INV only, with INS absent from 123 records.
#
# NOT a hard failure, deliberately: the pancancer model draws Ins at p=0.0279, so ~6% of
# runs legitimately plant none and failing those would be a false alarm. What is not
# acceptable is the silence, so the stratum is named here and again in the summary. Chasing
# a real generation defect needs replicates, not one run — see docs/sv_support_matrix.md,
# which records symbolic <INS> as a silent no-op from input_vcf (#500) while asserting the
# de novo path synthesises sequence. That claim is untested and this is the run that would
# have surfaced it.
# Sets TYPES_UNMEASURED (global) to every SV_TYPES entry with zero truth records. A
# function rather than an inline loop so it is reachable from the test suite: the
# behaviour that matters is "an empty stratum is named", and asserting that by grepping
# the source would only check that the code looks right.
split_truth_by_type() {  # <truth_sv.vcf.gz> <outdir>
    local truth="$1" dir="$2" svt typed n
    TYPES_UNMEASURED=()
    for svt in "${SV_TYPES[@]}"; do
        typed="$dir/truth_sv_${svt}.vcf.gz"
        bcftools view -i "INFO/SVTYPE=\"$svt\"" -O z -o "$typed" "$truth" 2>/dev/null || continue
        bcftools index -f -t "$typed" 2>/dev/null || true
        n=$(bcftools view -H "$typed" 2>/dev/null | wc -l)
        echo "  truth $svt: $n"
        if [[ "$n" -eq 0 ]]; then
            rm -f "$typed" "${typed}.tbi"
            TYPES_UNMEASURED+=("$svt")
            echo "    NOT MEASURED: 0 $svt planted, so no $svt recall exists for any caller." >&2
            echo "      This is a hole in the run's coverage, not a passing stratum." >&2
        fi
    done
    if [[ "${#TYPES_UNMEASURED[@]}" -gt 0 ]]; then
        echo "  SV types with an EMPTY denominator: ${TYPES_UNMEASURED[*]} (${#TYPES_UNMEASURED[@]} of ${#SV_TYPES[@]})"
        echo "    Every recall below is over the remaining types only."
    fi
}
TYPES_UNMEASURED=()
split_truth_by_type "$TRUTH_SV" "$OUTDIR"

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
# Strata whose ENTIRE query was filtered out before comparison (#511). Not fatal — a caller
# emitting only low-confidence calls in one stratum is normal — but it must reach the banner,
# because TP=0/FP=0 is indistinguishable from "the caller emitted nothing".
QUERY_FILTERED=()

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
        0)  MANTA_SV="$MANTA_SV_CONV"     # score against the normalised calls
            # A residual the raw/converted split CANNOT remove, stated here so the INV
            # precision figure is not misread. An INVERSION-ORIENTED junction (t]p] or
            # [p[t) has the same shape as an inversion, so convertInversion.py turns the
            # caller's breakend pair for OUR junction into an <INV> record, which then
            # scores as an INV false positive against a truth that never contained it.
            #
            # This used to assert that ALL our junctions were head-to-head and print
            # n_records/2 as the count. That was hardcoded from when geometry was
            # hardcoded, and v3.1.0 falsified it: job 20719077 measured 226 direct /
            # 128 head-to-head / 112 tail-to-tail across 466 junctions. The number is
            # now MEASURED, because a printed claim nobody checks is how the harness
            # ends up lying in exactly the way it exists to prevent.
            n_inv_or=$(bcftools query -f '%ALT\n' "$OUTDIR/truth_sv_BND.vcf.gz" 2>/dev/null \
                | grep -cE '^([ACGTNacgtn]+\]|\[)' || true)
            if [[ "${n_inv_or:-0}" -gt 0 ]]; then
                echo "    NOTE: $((n_inv_or / 2)) of the truth's junction(s) are inversion-oriented"
                echo "      (t]p] or [p[t) — a caller's breakends for those convert to <INV> and"
                echo "      count as INV false positives. Expect INV precision to be depressed by"
                echo "      up to that many; it is our representation, not a caller result."
                echo "      Direct junctions (t[p[ / ]p]t) do NOT convert; they are called"
                echo "      <DEL>/<DUP:TANDEM>, which is the BND-recall ceiling noted in stage 5."
            fi
            ;;
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
printf 'label\tfilter\tTP\tFN\tFP\trecall\tprecision\tf1\n' > "$OUTDIR/scorer_recall.tsv"

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
    local all_calls="${TRUVARI_ALL_CALLS:-0}"
    local filter_name="PASS-only"
    local suffix=""
    local filter_args=(--passonly)
    if [[ "$all_calls" == "1" ]]; then
        filter_name="all-calls"
        suffix="_allcalls"
        filter_args=()
    fi
    local out="$OUTDIR/truvari_${label}${suffix}"
    local csample; csample=$(comp_sample "$comp")
    rm -rf "$out"   # truvari requires the output dir not pre-exist
    truvari bench -b "$truth" -c "$comp" -f "$REFERENCE" -o "$out" \
        "${filter_args[@]}" --sizemax "$TRUVARI_SIZEMAX" \
        ${csample:+--cSample "$csample"} "$@" 2>&1 | tail -3 || true
    if [[ -f "$out/summary.json" ]]; then
        python3 - "$out/summary.json" "${label}${suffix}" <<'PY'
import json, sys
s = json.load(open(sys.argv[1])); name = sys.argv[2]
# truvari writes null (None) for metrics when a subset has no comparisons;
# `or 0` coerces None -> 0 (s.get(k, 0) only defaults a MISSING key, not a null one).
print(f"  {name:<18} recall={(s.get('recall') or 0):.3f}  precision={(s.get('precision') or 0):.3f}  "
      f"f1={(s.get('f1') or 0):.3f}  TP={s.get('TP-base') or 0}  FN={s.get('FN') or 0}  FP={s.get('FP') or 0}")
PY
        python3 - "$out/summary.json" "$label" "$filter_name" >> "$OUTDIR/scorer_recall.tsv" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))
def n(key):
    return s.get(key) or 0
def f(key):
    value = s.get(key)
    return "" if value is None else f"{value:.6f}"
print("\t".join([sys.argv[2], sys.argv[3], str(n("TP-base")), str(n("FN")),
                   str(n("FP")), f("recall"), f("precision"), f("f1")]))
PY
        check_denominator "$truth" "$out" "${label}${suffix}" || true
        check_query_denominator "$comp" "$out" "${label}${suffix}" || QUERY_FILTERED+=("${label}${suffix}")
    else
        echo "  $label: truvari produced no summary.json (see $out)"
        SCORE_COVERAGE_FAIL=$((SCORE_COVERAGE_FAIL + 1))
    fi
}

# Caller results are reported twice: the existing PASS-only score remains the headline, while
# the all-calls score attributes matches hidden solely by the pipeline's filter. The paired
# summaries and scorer_recall.tsv make the distinction durable instead of leaving it implicit
# in the query-coverage warning (#541).
run_truvari_both() {  # <truth.vcf.gz> <out-label> <comp.vcf.gz> [extra truvari args...]
    run_truvari "$@"
    TRUVARI_ALL_CALLS=1 run_truvari "$@"
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
# Coverage of the QUERY side, which nothing checked (#511). `check_denominator` below
# reconciles the TRUTH (n_truth vs TP-base+FN) and hard-fails when records were discarded
# before comparison. Nothing did the equivalent for the comp VCF, so every `precision` this
# pipeline reports was computed over a denominator no one verified — the same defect on the
# other axis.
#
# The motivating case, campaign 20925151: `manta_INS` reported TP=0 FP=0, which reads as "the
# caller emitted nothing". It emitted three INS calls, one of them a TRUE detection
# (chr7:3198935 SVLEN 61 against truth chr7:3198936 SVLEN 61 — one base off, exact size). All
# three were non-PASS, so `--passonly` removed them BEFORE comparison and the stratum scored
# 0/0. The reported 0.000 was a statement about the FILTER, not the caller, and nothing said so.
#
# Split into a numbers-only verdict plus thin plumbing, so the decision logic is testable
# without truvari or a BAM — the same seam as parse_lfs_quota_kb.
#
# LOUD, NOT FATAL, deliberately. A caller emitting only low-confidence calls in one stratum is
# normal behaviour, not a harness bug; failing the run would let it abort an 8-hour campaign.
# It surfaces in the closing banner instead, so it cannot be missed while reading the numbers.
report_query_coverage() {  # <label> <n_query> <n_scored> <n_nonpass>; rc 1 = wholly filtered
    local label="$1" n_query="$2" scored="$3" nonpass="$4"

    [[ "$n_query" -gt 0 ]] || return 0            # nothing was handed over; not our business
    if [[ "$scored" -gt "$n_query" ]]; then
        # More classified than handed over should be impossible. Report rather than compute a
        # negative "dropped" count and present it as a fact.
        echo "  WARNING: $label scored $scored of $n_query query record(s) — MORE than were" >&2
        echo "    supplied, so this arithmetic is not understood. Do not quote its precision." >&2
        return 0
    fi
    [[ "$scored" -eq "$n_query" ]] && return 0    # fully accounted for

    local dropped=$((n_query - scored))
    local why="size/other filters"
    [[ "$nonpass" -gt 0 ]] && why="--passonly ($nonpass non-PASS)"

    if [[ "$scored" -eq 0 ]]; then
        # The dangerous case: TP=0 FP=0 is indistinguishable from "the caller emitted nothing".
        echo "ERROR: $label scored NONE of its $n_query query record(s) — all removed before" >&2
        echo "  comparison by $why. recall/precision for this stratum therefore describe the" >&2
        echo "  FILTER, not the caller, and must not be reported as a caller result (#511)." >&2
        [[ "$nonpass" -gt 0 ]] && \
            echo "  --passonly is THIS PIPELINE'S choice, not a property of the caller." >&2
        return 1
    fi
    echo "    ^ query coverage: $scored of $n_query record(s) scored, $dropped dropped by $why"
    return 0
}

check_query_denominator() {  # <comp.vcf.gz> <truvari_out_dir> <label>
    local comp="$1" out="$2" label="$3"
    local n_query nonpass kv tp fp
    n_query=$(bcftools view -H "$comp" 2>/dev/null | wc -l)
    # Read TP-comp and FP together, and FAIL if either key is absent rather than defaulting to
    # zero — an absent key would make every record look dropped and manufacture the exact alarm
    # this function exists to raise. #509's lesson, one guard over.
    kv=$(python3 - "$out/summary.json" <<'PYQ' 2>/dev/null
import json, sys
d = json.load(open(sys.argv[1]))
for k in ("TP-comp", "FP"):
    if k not in d:
        sys.exit(3)
print(int(d["TP-comp"] or 0), int(d["FP"] or 0))
PYQ
    ) || {
        echo "  ($label: query coverage UNVERIFIED — summary.json lacks TP-comp/FP)" >&2
        return 0
    }
    read -r tp fp <<< "$kv"
    nonpass=$(bcftools view -H -i 'FILTER!="PASS" && FILTER!="."' "$comp" 2>/dev/null | wc -l)
    report_query_coverage "$label" "$n_query" "$((tp + fp))" "$nonpass"
}

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

# ── BND geometry: does the truth describe all four adjacencies, reciprocally? ──
# Before v3.1.0 every emitted junction was head-to-head (t]p]). A caller therefore only
# ever saw inversion-like adjacencies, and "BND recall=1.000" was a measurement of one
# quarter of the problem that read like the whole of it. Geometry is now sampled across
# all four VCF 4.2 §5.4 forms, so this asserts two things the old harness never did:
#
#   1. COVERAGE OF ITS OWN INPUT — every ALT parses into one of the four forms, and the
#      per-form counts are reported. A geometry distribution over an unknown denominator
#      is not a result; if a form is missing, that must be visible rather than inferred
#      from a recall number that moved.
#   2. RECIPROCITY — each junction's two records form a spec-valid pair. This is the
#      #451 defect stated as an assertion: the truth described a rearrangement the reads
#      did not carry, and nothing checked the two sides agreed, so it looked like a
#      caller failure for two months (v1.13.1 2026-06-02 -> v3.1.0 2026-08-01).
#
# Pairings are the spec's own worked examples: t[p[ <-> ]p]t (direct), t]p] <-> t]p]
# (head-to-head), [p[t <-> [p[t (tail-to-tail).
check_bnd_geometry() {  # <bnd_truth.vcf.gz>
    local src="$1"
    [[ -f "$src" ]] || { echo "  ERROR: BND geometry: no BND truth at $src" >&2; return 1; }
    bcftools query -f '%CHROM\t%POS\t%ALT\n' "$src" 2>/dev/null \
      | awk -F'\t' '
          # Which of the four forms is this ALT? "" means unparsable, which is a defect
          # and is counted separately rather than being silently dropped.
          function form(a) {
              if (a ~ /^[ACGTNacgtn]+\[/) return "t[p["   # case 1  DeletionLike
              if (a ~ /^[ACGTNacgtn]+\]/) return "t]p]"   # case 2  HeadToHead
              if (a ~ /^\[/)              return "[p[t"   # case 3  TailToTail
              if (a ~ /^\]/)              return "]p]t"   # case 4  DuplicationLike
              return ""
          }
          function mate(a,   b,n,p,i) {
              b=a; gsub(/[][]/,"|",b); n=split(b,p,"|")
              for (i=1;i<=n;i++) if (p[i] ~ /:/) return p[i]
              return ""
          }
          # The spec pairing. Note it is NOT "same form": a direct junction is described
          # by two DIFFERENT forms, which is exactly the case a naive check would miss.
          function mates_ok(a,b) {
              if (a=="t[p[") return b=="]p]t"
              if (a=="]p]t") return b=="t[p["
              if (a=="t]p]") return b=="t]p]"
              if (a=="[p[t") return b=="[p[t"
              return 0
          }
          { n++
            f=form($3); m=mate($3)
            if (f=="" || m=="") { bad++; badex=(badex==""?$3:badex); next }
            cnt[f]++; parsed++
            self=$1":"$2; g[self]=f; tgt[self]=m; ord[++k]=self }
          END{
            printf("  BND geometry: %d record(s), %d parsed into a form, %d unparsable\n",
                   n+0, parsed+0, bad+0)
            split("t[p[ t]p] [p[t ]p]t", forms, " ")
            split("direct/deletion-like head-to-head tail-to-tail direct/duplication-like",
                  names, " ")
            for (i=1;i<=4;i++)
                printf("    %-5s %-22s %d\n", forms[i], names[i], cnt[forms[i]]+0)
            # A truth that is entirely one form is internally CONSISTENT, so every check
            # below passes on it — that is exactly what pre-v3.1.0 output looked like.
            # Flagged, not failed: fitted weights may legitimately be skewed, so a hard
            # threshold here would be a flaky assertion about a random draw. The point is
            # that a silent regression to one geometry must not be invisible.
            nd=0; for (i=1;i<=4;i++) if (cnt[forms[i]]>0) nd++
            if (nd==1 && parsed>=24)
                printf("    NOTE: all %d junction record(s) share ONE geometry. That is the\n"\
                       "      pre-v3.1.0 hardcoded behaviour; under uniform weights it is\n"\
                       "      vanishingly unlikely at this count. Check bnd_geometry_weights\n"\
                       "      in the model before reading any BND number below.\n",
                       parsed+0)
            # Reciprocity, checked per record against the record its ALT points at.
            for (i=1;i<=k;i++) {
                s=ord[i]; t=tgt[s]
                if (!(t in g)) { unpaired++
                    if (uex=="") uex=s" -> "t" (no record there)"
                    continue }
                if (!mates_ok(g[s], g[t])) { mismatch++
                    if (mex=="") mex=s" "g[s]" <-> "t" "g[t] }
            }
            printf("    reciprocity: %d unpaired, %d mispaired\n", unpaired+0, mismatch+0)
            if (bad)      printf("    first unparsable ALT: %s\n", badex)
            if (unpaired) printf("    first unpaired:       %s\n", uex)
            if (mismatch) printf("    first mispaired:      %s\n", mex)
            # A record whose mate is absent or contradicts it describes a rearrangement
            # the reads cannot carry. That is not a caller result to be interpreted, it
            # is our artifact being wrong, so it is a hard failure.
            exit (bad+unpaired+mismatch > 0) ? 1 : 0
          }'
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

# Read one numeric metric out of a truvari summary.json. Prints the value and returns 0;
# returns 1 if the summary is MISSING or unparseable.
#
# WHY THE DISTINCTION IS THE WHOLE POINT: both calibration controls used to open with
# `local r=0` and then overwrite it only `[[ -f "$out/summary.json" ]]`. A truvari that
# crashed therefore left r=0, which the selftest reads as "recall 0 → FAIL" (safe) but the
# decoy reads as "recall 0 → PASS" — a negative control that passes *because it never ran*.
# `truvari … >/dev/null 2>&1 || true` then discarded the only evidence of why.
# Job 20853511 printed no decoy verdict at all with an empty truvari_decoy_scoreable/,
# which is what sent us looking. An absent measurement must never be spellable as a
# passing one; that is the same "assert content, not existence" rule as everywhere else.
truvari_metric() {  # <summary.json> <key>
    [[ -f "$1" ]] || return 1
    python3 - "$1" "$2" <<'PY' || return 1
import json, sys
v = json.load(open(sys.argv[1])).get(sys.argv[2])
print(float(v) if v is not None else 0.0)
PY
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
    local err="$OUTDIR/.selftest_${label}.stderr"
    rm -rf "$out"
    truvari bench -b "$truth" -c "$truth" -f "$REFERENCE" -o "$out" \
        --passonly --sizemax "$TRUVARI_SIZEMAX" "$@" >"$err" 2>&1 || true
    # "Did not run" is reported separately from "ran and scored 0". Both fail, but only
    # one of them is a scoring-threshold problem, and conflating them cost a day.
    local r
    if ! r=$(truvari_metric "$out/summary.json" recall); then
        echo "ERROR: scorer selftest for $label DID NOT RUN — no readable" >&2
        echo "  $out/summary.json. truvari's own output follows:" >&2
        sed -n '1,20p' "$err" >&2 2>/dev/null || true
        return 1
    fi
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
    # Displace each record by MORE THAN ITS OWN LENGTH. A flat shift does not decoy a
    # variant larger than the shift: truvari matches sized variants on reciprocal overlap,
    # so a 19.6 Mb deletion moved 2 Mb still overlaps 90% and matching it is CORRECT.
    # Measured on job 20824321, where the old flat 2 Mb shift gave decoy recall=0.202 with
    # TP=17 — and all 17 were larger than the shift (min 2,048,336 bp vs a 2,000,000 bp
    # shift, median 4.3 Mb, max 19.6 Mb). Not one match was smaller than the displacement,
    # so the control was reporting "the shift is too small", not "the scorer is too loose",
    # while its warning text claimed the latter.
    # A breakend has no SVLEN and END==POS, so it falls back to the 2 Mb floor, which
    # displaces a point event completely.
    { bcftools view -h "$truth"
      bcftools view -H "$truth" | awk -F'\t' 'BEGIN{OFS="\t"}
        { len=0
          if (match($8, /SVLEN=-?[0-9]+/)) { len=substr($8, RSTART+6, RLENGTH-6); if (len<0) len=-len }
          shift_by = 2000000; if (2*len > shift_by) shift_by = 2*len
          $2 = $2 + shift_by
          # END must move with POS or the record describes a span it no longer occupies.
          if (match($8, /END=[0-9]+/)) {
              e = substr($8, RSTART+4, RLENGTH-4) + shift_by
              $8 = substr($8, 1, RSTART-1) "END=" e substr($8, RSTART+RLENGTH)
          }
          print }'; } \
      | bcftools sort -O z -o "$decoy" 2>/dev/null || { echo "  decoy $label: could not build (skipped)"; return 0; }
    bcftools index -f -t "$decoy" 2>/dev/null || true
    rm -rf "$out"
    local err="$OUTDIR/.decoy_${label}.stderr"
    truvari bench -b "$truth" -c "$decoy" -f "$REFERENCE" -o "$out" \
        --passonly --sizemax "$TRUVARI_SIZEMAX" "$@" >"$err" 2>&1 || true
    local r
    if ! r=$(truvari_metric "$out/summary.json" recall); then
        echo "ERROR: decoy control for $label DID NOT RUN — no readable" >&2
        echo "  $out/summary.json. This is NOT a pass: the decoy is the only evidence that" >&2
        echo "  a recall figure means detection rather than a matching window loose enough" >&2
        echo "  to accept anything, so every number from this run would be uncalibrated." >&2
        echo "  truvari's own output follows:" >&2
        sed -n '1,20p' "$err" >&2 2>/dev/null || true
        return 1
    fi
    if awk -v r="$r" 'BEGIN{exit !(r < 0.001)}'; then
        echo "  decoy    $label: PASS (shifted truth recall=$r, expected 0)"
        return 0
    fi
    echo "ERROR: decoy control for $label matched (recall=$r). Displacement scales with" >&2
    echo "  SVLEN, so this is not an undersized shift — the scorer is matching on position" >&2
    echo "  more loosely than intended, or records ran off the end of their contig." >&2
    echo "  Inspect $out/tp-base.vcf.gz." >&2
    return 1
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
# <caller> <comp.vcf.gz> [raw_comp.vcf.gz]
#
# `comp` is what most comparisons score against — for Manta that is the
# inversion-normalised file, so <INV> truth is scored like-for-like. `raw_comp` is the
# PRE-conversion call set and exists for one reason: convertInversion.py consumes the
# breakend records it converts, so scoring BND against the converted file compares the
# BND truth against a file with no breakends in it.
#
# Job 20699004 did exactly that. The conversion reported "32 BND -> 0 BND + 16 INV" and
# manta_BND then read recall=0.000 "(caller emitted no BND records)" — down from 1.000
# the run before. That was a harness artifact, not a caller result.
#
# So: BND scores against RAW, everything else against CONVERTED. The two are answers to
# different questions and are NOT additive — the same underlying call can appear as a
# breakend pair in raw and as one <INV> in converted.
score_caller() {  # <caller> <comp.vcf.gz> [raw_comp.vcf.gz]
    local caller="$1" comp="$2" raw_comp="${3:-$2}"
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
            run_truvari_both "$TRUTH_SCOREABLE" "${caller}_overall" "$qs"
        fi
    else
        echo "  ($caller: no truvari-scoreable truth types — aggregate skipped)"
    fi

    # Iterate the CANONICAL type list, never a glob of whatever files happen to exist.
    #
    # This used to be `for typed in "$OUTDIR"/truth_sv_*.vcf.gz` with a hand-maintained
    # exclusion list, which meant any newly derived truth artifact was silently promoted
    # to an SVTYPE. Adding truth_sv_BNDinv.vcf.gz did exactly that in job 20745149:
    #     manta_BNDinv recall=0.000  (caller emitted no BNDinv records)
    #     manta_BNDinv       recall=0.738 ... TP=354  FN=126
    # Two contradictory lines for one label, and for Delly only the false one appeared.
    # The bug was not the missing exclusion, it was that the loop's input was the
    # directory contents; adding to an exclusion list would leave the next artifact to
    # find the same trap. SV_TYPES_SCORED_IN_LOOP is derived from SV_TYPES so a new type
    # is picked up automatically while a new derived FILE cannot inject a fake one.
    local svt
    for svt in "${SV_TYPES_SCORED_IN_LOOP[@]}"; do
        local typed="$OUTDIR/truth_sv_${svt}.vcf.gz"
        [[ -f "$typed" ]] || continue
        local q="$OUTDIR/${caller}_calls_${svt}.vcf.gz"
        bcftools view -i "INFO/SVTYPE=\"$svt\"" -O z -o "$q" "$comp" 2>/dev/null || true
        if [[ -f "$q" ]] && [[ "$(bcftools view -H "$q" 2>/dev/null | wc -l)" -gt 0 ]]; then
            bcftools index -f -t "$q" 2>/dev/null || true
            run_truvari_both "$typed" "${caller}_${svt}" "$q"
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
        # From the RAW calls: the converted file has had its breakends consumed.
        bcftools view -i 'INFO/SVTYPE="BND"' -O z -o "$qb" "$raw_comp" 2>/dev/null || true
        if [[ -f "$qb" ]] && [[ "$(bcftools view -H "$qb" 2>/dev/null | wc -l)" -gt 0 ]]; then
            bcftools index -f -t "$qb" 2>/dev/null || true
            run_truvari_both "$OUTDIR/truth_sv_BND.vcf.gz" "${caller}_BND" "$qb" \
                --bnddist "$BND_DIST"
            # RETIRED: the BND/BNDinv split (#507). It existed because the de novo
            # sampler hardcoded the mate to the anchor's own contig, and within one
            # chromosome a direct adjacency genuinely IS a deletion or tandem
            # duplication — so a caller emitted <DEL>/<DUP:TANDEM> and all 452
            # direct-form records of job 20719077 were unmatchable as BND by
            # construction, flooring the aggregate at 0.380.
            #
            # BND is now inter-chromosomal only. Across two chromosomes there is no
            # deletion or duplication interpretation, so a caller emits BND for all four
            # geometries and the direct-form records match fine. Measured on job
            # 20824321 the two collapsed (BND 0.935 vs BNDinv 0.929), and pooled over
            # campaign 20925151 they agree to within noise (manta 0.868 vs 0.876, delly
            # 0.417 vs 0.412). `${caller}_BND` is simply the BND number now.
            #
            # It was also actively misleading: BNDinv restricted the TRUTH to the
            # inversion-oriented records while leaving the QUERY as every BND call, so
            # each unmatched call became a false positive and precision read 0.245
            # against BND's 0.547 on an identical query set (#512). The surrounding text
            # told the reader to quote that row.
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
        run_truvari_both "$OUTDIR/truth_sv_BNDspan.vcf.gz" "${caller}_BNDspan" "$raw_comp" -t
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
    # Runs BEFORE the spans are derived: if the junctions are not reciprocal, every
    # number downstream describes our own broken artifact rather than any caller.
    # `|| rc=$?` and not `$(cmd; echo $?)` — the latter loses the status under set -e.
    BND_GEOM_RC=0
    check_bnd_geometry "$OUTDIR/truth_sv_BND.vcf.gz" || BND_GEOM_RC=$?
    if [[ "$BND_GEOM_RC" -ne 0 ]]; then
        echo "ERROR: the BND truth is not internally consistent (see above). Scoring it" >&2
        echo "  would measure our emission defect and report it as caller performance." >&2
        CAL_FAIL=1
    fi
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
    decoy_truvari    "$TRUTH_SCOREABLE" scoreable || CAL_FAIL=1
fi
if [[ -f "$OUTDIR/truth_sv_BND.vcf.gz" ]]; then
    selftest_truvari "$OUTDIR/truth_sv_BND.vcf.gz" BND --bnddist "$BND_DIST" || CAL_FAIL=1
    # BND had a positive control and NO negative one, which is the half that bounds
    # whether a recall figure means detection or coincidence. `--bnddist` sets a
    # proximity window; without a decoy nothing measures how often a *displaced*
    # breakend still lands inside it. Cheap, and this is exactly where the decoy works
    # best: a breakend is a point (END==POS), so any shift displaces it fully.
    decoy_truvari    "$OUTDIR/truth_sv_BND.vcf.gz" BND --bnddist "$BND_DIST" || CAL_FAIL=1
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
    echo "ERROR: at least one scorer failed calibration — either it cannot match the" >&2
    echo "  truth to itself (positive control) or it matched a deliberately displaced" >&2
    echo "  copy of it (negative control), and a control that DID NOT RUN counts as a" >&2
    echo "  failure, not a pass. Refusing to report caller results from a configuration" >&2
    echo "  not demonstrated to discriminate. Skipping caller scoring; see the output" >&2
    echo "  above and the archived truvari_selftest_*/truvari_decoy_* directories." >&2
else
    score_caller manta "$MANTA_SV" "$MANTA_SV_RAW"
    # Delly emits no BND at all and is not converted, so raw and scored are the same file.
    [[ -n "$DELLY_SV" ]] && score_caller delly "$DELLY_SV"
fi

# Count, per probe, how many input SEQUENCES contain it or its reverse complement. Probes come
# from a TSV (chrom, pos, inserted-length, probe); sequences arrive on stdin, one per line.
# Split out from verify_planted_ins so the matching is testable without samtools or a BAM.
count_probe_hits() {  # <probes.tsv>  < sequences
    awk -v pf="$1" '
        function rc(s,   i, o, c) {
            o = ""
            for (i = length(s); i >= 1; i--) {
                c = substr(s, i, 1)
                o = o (c == "A" ? "T" : c == "T" ? "A" : c == "C" ? "G" : c == "G" ? "C" : c)
            }
            return o
        }
        BEGIN {
            FS = OFS = "\t"
            while ((getline line < pf) > 0) {
                split(line, f, "\t")
                if (f[4] == "") continue
                n++; chrom[n] = f[1]; pos[n] = f[2]; len[n] = f[3]
                fwd[n] = f[4]; rev[n] = rc(f[4]); hits[n] = 0
            }
        }
        { for (i = 1; i <= n; i++) if (index($0, fwd[i]) || index($0, rev[i])) hits[i]++ }
        END { for (i = 1; i <= n; i++) print chrom[i], pos[i], len[i], hits[i] }
    '
}

# READ-LEVEL verification of what we PLANTED, before anything is pruned.
#
# WHY THIS EXISTS: every campaign to date produced BAMs, scored VCF-against-VCF, and deleted
# the BAMs without ever asking whether the reads contained what the truth VCF declared. #516
# survived all of them that way — insertions above ~a read length are spliced in but only their
# first ~100bp is sequenced, so campaign 20925151 reported Manta INS recall 1/22 and the
# missing evidence was ours. A caller recall of 0 is uninterpretable until you know the
# evidence was there to find.
#
# INS is what gets checked because the truth VCF carries the inserted sequence LITERALLY, so a
# probe is derivable with no extra input and no assumptions. Reads carrying novel sequence are
# either soft-clipped near the locus or have no reference home at all, so both the locus
# neighbourhood and the unmapped pile are searched. Set VERIFY_INS_READS=0 to skip (the
# unmapped scan streams the whole BAM — minutes, against a multi-hour run).
verify_planted_ins() {  # <truth_sv.vcf.gz> <tumor.bam> <outdir>
    local truth="$1" bam="$2" dir="$3"
    [[ "${VERIFY_INS_READS:-1}" == "1" ]] || { echo "  (INS read verification skipped)"; return 0; }
    if [[ ! -f "$bam" ]]; then
        echo "ERROR: $bam is missing — INS read verification was not evaluated" >&2
        return 1
    fi

    local probes="$dir/.ins_probes.tsv"
    # ALT is anchor+novel for a literal insertion, so the inserted bases are ALT past REF's
    # length. Each probe is a 30-mer of that novel sequence: it cannot come from the
    # reference (~4^-30 by chance). Probe the INTERIOR, never the head — a partially
    # realized insertion has a perfectly good head, which is how #516 hid.
    #
    # Five probes, not one. A single midpoint 30-mer is a fragile test of a real property:
    # measured on a 5 Mb de novo run, interior coverage of a heterozygous insertion is
    # ~half the flank's, so ordinary dropout leaves ~3% of insertions with no read over any
    # one chosen base while 68-96% of their length is covered. That is sparse coverage, not
    # a missing insertion, and failing the run on it would cry wolf every time. Requiring
    # only that SOME interior probe is seen keeps the #516 signature (no interior at all)
    # a hard failure while tolerating a thin subclonal insertion.
    if ! bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'INFO/SVTYPE="INS"' "$truth" 2>/dev/null \
      | awk -F'\t' 'BEGIN{OFS="\t"} {
            ins = substr($4, length($3) + 1); L = length(ins)
            if (L < 30) next
            delete seen
            for (k = 1; k <= 5; k++) {
                s = int(L * k / 6) - 14
                if (s < 1) s = 1
                if (s > L - 29) s = L - 29
                if (s in seen) continue
                seen[s] = 1
                print $1, $2, L, substr(ins, s, 30)
            }
        }' > "$probes"; then
        echo "ERROR: could not extract literal INS probes from $truth — verification was not evaluated" >&2
        return 1
    fi

    local n_probes; n_probes=$(wc -l < "$probes" 2>/dev/null || echo 0)
    local n_ins; n_ins=$(cut -f1,2 "$probes" 2>/dev/null | sort -u | wc -l)
    if [[ "$n_probes" -gt 0 && "${n_ins:-0}" -le 0 ]]; then
        echo "ERROR: $n_probes probe(s) collapsed to no insertion — verification was not evaluated" >&2
        return 1
    fi
    if [[ "$n_probes" -eq 0 ]]; then
        echo "  INS read check: no INS record carries >=30bp of novel sequence — nothing to verify"
        return 0
    fi

    local seqs="$dir/.ins_seqs.txt"
    : > "$seqs"
    # One fetch per insertion, not per probe — the five probes of an insertion share a locus.
    # Materialized first, not piped in via <(...): a process substitution that dies is
    # invisible to set -e and would silently run zero iterations.
    local loci
    if ! loci="$(cut -f1,2 "$probes" | sort -u)"; then
        echo "ERROR: could not enumerate INS loci — verification was not evaluated" >&2
        rm -f "$seqs"
        return 1
    fi
    while IFS=$'\t' read -r c p; do
        [[ -n "$c" ]] || continue
        local lo=$(( p > 1000 ? p - 1000 : 1 ))
        if ! samtools view "$bam" "$c:$lo-$((p + 1000))" 2>/dev/null | cut -f10 >> "$seqs"; then
            echo "ERROR: samtools could not read the INS locus $c:$p — verification was not evaluated" >&2
            rm -f "$seqs"
            return 1
        fi
    done <<< "$loci"
    # A read lying wholly inside a large insertion has no reference home, so it lands here.
    if ! samtools view -@ 4 -f 4 "$bam" 2>/dev/null | cut -f10 >> "$seqs"; then
        echo "ERROR: samtools could not read unmapped tumor reads — verification was not evaluated" >&2
        rm -f "$seqs"
        return 1
    fi

    echo "  INS read-level support (5 interior 30-mers per inserted sequence):"
    local supported=0 unsupported=0
    local hit_rows
    if ! hit_rows="$(count_probe_hits "$probes" < "$seqs")"; then
        echo "ERROR: probe matching failed — INS read verification was not evaluated" >&2
        rm -f "$seqs"
        return 1
    fi

    local result_rows=0
    if [[ -n "$hit_rows" ]]; then
        result_rows="$(printf '%s\n' "$hit_rows" | awk -F'\t' '$4 ~ /^[0-9]+$/ { n++ } END { print n + 0 }')"
    fi
    if [[ "$result_rows" -ne "$n_probes" ]]; then
        echo "ERROR: probe matcher evaluated $result_rows of $n_probes probes — INS read verification was not evaluated" >&2
        rm -f "$seqs"
        return 1
    fi

    # Collapse the per-probe rows to one verdict per insertion, keyed on chrom+pos.
    local per_ins
    per_ins="$(printf '%s\n' "$hit_rows" | awk -F'\t' 'BEGIN{OFS="\t"}
        {
            k = $1 SUBSEP $2
            if (!(k in seen)) { seen[k] = 1; order[++m] = k; chrom[k] = $1; pos[k] = $2; len[k] = $3 }
            np[k]++
            if ($4 > 0) { nh[k]++; if ($4 > best[k]) best[k] = $4 }
        }
        END { for (i = 1; i <= m; i++) { k = order[i]
            print chrom[k], pos[k], len[k], nh[k] + 0, np[k], best[k] + 0 } }')"

    while IFS=$'\t' read -r c p l nh np best; do
        [[ -n "$c" ]] || continue
        if [[ "${nh:-0}" -gt 0 ]]; then
            supported=$((supported + 1))
            printf '    %s:%s  %s bp  %s/%s probes  %s read(s)\n' "$c" "$p" "$l" "$nh" "$np" "$best"
        else
            unsupported=$((unsupported + 1))
            printf '    %s:%s  %s bp  NO READ SUPPORT (0/%s interior probes)\n' "$c" "$p" "$l" "$np" >&2
        fi
    done <<< "$per_ins"

    if [[ $((supported + unsupported)) -ne "$n_ins" ]]; then
        echo "ERROR: probe result accounting covered $((supported + unsupported)) of $n_ins insertion(s)" >&2
        rm -f "$seqs"
        return 1
    fi

    echo "  INS read support: $supported of $n_ins planted insertion(s) present in the reads"
    INS_UNSUPPORTED="$unsupported"

    # Rule 4: a metric over an unknown — or tiny — denominator is not a result. The COSMIC
    # somatic model plants very few insertions on a single chromosome (job 21378484: FOUR, on
    # chr22 at SV_RATE_SCALE=30), and "4 of 4, 100%" reads exactly like "validated" while
    # resting on four draws. It is not a hard failure — the run's DEL/DUP/INV results are
    # unaffected and the check still catches a real 0/5 — but it must never be quoted as
    # evidence that insertions work at scale.
    local min_denom="${INS_MIN_DENOM:-20}"
    if [[ "$n_ins" -lt "$min_denom" ]]; then
        INS_DENOM_THIN="$n_ins"
        echo "  DENOMINATOR TOO SMALL: $n_ins insertion(s) is below $min_denom, so the percentage"
        echo "    above cannot support a conclusion about insertion realization at scale. A pass"
        echo "    here means these $n_ins worked, nothing more. Raise SV_RATE_SCALE, use a larger"
        echo "    reference, or include germline insertions to get a real denominator."
    fi
    if [[ "$unsupported" -gt 0 ]]; then
        # Sparse interior coverage costs a few insertions on any real run; a #516-class
        # regression costs nearly all of them. Only the latter invalidates INS recall.
        local min_pct="${INS_SUPPORT_MIN_PCT:-95}"
        local pct=$(( supported * 100 / n_ins ))
        echo "  INS read support: ${pct}% (floor ${min_pct}%, set INS_SUPPORT_MIN_PCT to change)"
        if [[ "$pct" -lt "$min_pct" ]]; then
            echo "ERROR: only ${pct}% of planted insertion(s) have interior read evidence for" >&2
            echo "  their declared length, below the ${min_pct}% floor. Any INS recall reported" >&2
            echo "  above is a statement about our own read generation, not about the caller —" >&2
            echo "  see #516. DEL/DUP/INV/BND/CNV are unaffected." >&2
        else
            echo "  ($unsupported thinly-covered insertion(s) above the floor — reported, not fatal)"
            INS_UNSUPPORTED=0
        fi
    fi
    rm -f "$seqs"
}
INS_UNSUPPORTED=0
INS_VERIFY_FAILED=0
INS_DENOM_THIN=0
if [[ -f "$OUTDIR/truth_sv_INS.vcf.gz" ]]; then
    if ! verify_planted_ins "$OUTDIR/truth_sv_INS.vcf.gz" "$TUMOR_BAM" "$OUTDIR"; then
        INS_VERIFY_FAILED=1
    fi
fi

# ── Summary ──────────────────────────────────────────────────────────────
# Never print an unqualified "complete" over a run whose numbers are not usable — a
# summary banner that reads as success regardless of outcome is how a bad run gets
# believed. The detail follows in the deferred failure gate below.
if [[ "$CAL_FAIL" -ne 0 ]]; then
    RUN_STATUS="FAILED (scorer calibration) — no caller results"
elif [[ "$SCORE_COVERAGE_FAIL" -gt 0 ]]; then
    RUN_STATUS="FAILED ($SCORE_COVERAGE_FAIL comparison(s) scored a SUBSET of the truth)"
elif [[ "${INS_VERIFY_FAILED:-0}" -ne 0 ]]; then
    RUN_STATUS="FAILED (INS read verification was not evaluated)"
elif [[ "${INS_UNSUPPORTED:-0}" -gt 0 ]]; then
    # Deliberately not a hard failure while #516 is open: every campaign planting insertions
    # above ~a read length would fail, which would block DEL/DUP/INV work that is unaffected.
    # It is a distinct, visible status instead — the numbers are usable EXCEPT for INS.
    RUN_STATUS="complete, but INS INVALID ($INS_UNSUPPORTED insertion(s) absent from the reads — #516)"
elif [[ "${INS_DENOM_THIN:-0}" -gt 0 ]]; then
    # Passing on a handful of insertions is not evidence that insertions work. Say so in the
    # banner, where it cannot be skipped, rather than only in the INS block above.
    RUN_STATUS="complete; INS read support passed on only $INS_DENOM_THIN insertion(s) — too few to conclude"
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
                Aggregate/DEL/DUP/INV are scored against this (inversion-normalised)
                file; BND and BNDspan are scored against the RAW calls above, because
                convertInversion.py consumes the breakend records it converts. The two
                are separate questions and their TP counts are not additive.
  Delly calls:  delly.somatic.vcf.gz       (RUN_DELLY=1)
  Scores:       truvari_<caller>_<overall|TYPE>/summary.json   (manta + delly)

Interpretation: compare the two callers side-by-side. A type eidolon emitted heavily
that BOTH callers miss = candidate SV-generation issue; one caller recovering what
the other misses (especially BND) = a caller limitation, not a simulator defect.
Low overall with high per-type recall = truvari SVTYPE/SVLEN representation
mismatch (inspect the pre-flight schema dump).
EOF

# Repeated here on purpose. Stage 3b prints this ~900 lines earlier in the log, and the
# banner above is what actually gets read and pasted — an empty denominator that only
# appears upstream is an empty denominator nobody sees.
if [[ "${#QUERY_FILTERED[@]}" -gt 0 ]]; then
    cat <<EOF
QUERY WHOLLY FILTERED: ${#QUERY_FILTERED[@]} stratum/strata scored NONE of the caller's records
because every one was removed before comparison: ${QUERY_FILTERED[*]}
  Their recall/precision describe the FILTER, not the caller, and must not be quoted as a
  caller result. --passonly is this pipeline's choice; a caller emitting only low-confidence
  calls in a stratum is normal behaviour. To see what was actually emitted:
    bcftools query -f '%CHROM\t%POS\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\n' <caller vcf>
  Campaign 20925151 is the case history: manta_INS read TP=0 FP=0 while Manta had emitted
  three INS calls, one of them a true detection, all non-PASS. See #511.
EOF
fi
if [[ "${#TYPES_UNMEASURED[@]}" -gt 0 ]]; then
    cat <<EOF
COVERAGE HOLE: ${#TYPES_UNMEASURED[@]} of ${#SV_TYPES[@]} SV type(s) had NO truth records and were
not measured at all: ${TYPES_UNMEASURED[*]}
  Every recall above is over the remaining types. There is no caller result for these,
  neither pass nor fail. An empty stratum has TWO causes that look identical here, and
  they are diagnosed in opposite directions:
    (a) nothing was PLANTED  -> the draw is stochastic; replicate before concluding.
    (b) nothing was SELECTABLE -> records exist but this harness cannot see them,
        because the SV selector needs a symbolic ALT or an INFO/SVTYPE tag.
  Check (b) FIRST — it is cheap and it is what actually happened with INS, where the
  truth VCF held insertions as untagged literal records for every run before v3.1.x
  and 'truth INS: 0' was a statement about the selector, not the simulator:
    bcftools view -H -i 'INFO/EIDOLON_ORIGIN="somatic"' <merged truth> \\
      | awk 'length(\$5)-length(\$4) >= 50 && \$5 !~ /</ {n++} END {print n+0" untagged"}'
  A non-zero count there means the stratum is a selector gap, and no number of
  replicates will fix it.
EOF
fi

# Reclaim scratch once scoring is done: FASTQ is the bulk (~GBs/run at 60x) and
# isn't needed afterward — BAMs, Manta calls, the truth VCFs, and the truvari
# summaries are kept, so re-scoring and re-calling still work. This keeps a
# parameter sweep (per-tissue / purity) from exhausting the shared project quota.
# Set PRUNE_FASTQ=0 to retain the reads.
if [[ "${PRUNE_FASTQ:-1}" == "1" ]]; then
    echo "[prune] removing FASTQ intermediates (set PRUNE_FASTQ=0 to keep)..."
    rm -f "$OUTDIR"/sample_*_r1.fastq.gz "$OUTDIR"/sample_*_r2.fastq.gz
fi

# The BAMs are the rest of the bulk — ~101 GB of a ~113 GB replicate (30 GB normal +
# 71 GB tumor, measured on job 20892075). Pruning them here, INSIDE the pipeline, rather
# than from a dependent cleanup job.
#
# WHY INLINE: the external strip job failed four different ways in one week — a script
# path under $HOME that did not exist, a glob pointing at /work/nvme while runs landed on
# /scratch, a missing --account, and a script whose creation was never confirmed. Each
# time it exited quietly and nothing downstream noticed, because a cleanup job that does
# nothing looks exactly like one that had nothing to do. An `--array=4-8%1` campaign then
# accumulated ~100 GB per finished replicate until rep 5 hit the project quota and the
# array died. A separate job is a separate thing to forget; this cannot be forgotten.
#
# GUARD: only when scoring actually produced a summary. A run that failed calibration, or
# died before stage 5, keeps its BAMs — that is exactly when someone needs to re-call.
# Note this makes the replicate re-SCOREABLE but not re-CALLABLE; the FASTQ prune above
# already removed re-alignment, and CLAUDE.md's guidance for a campaign is to clear
# finished replicates down to truvari_*/summary.json, which is all aggregate_sv_reps.sh
# reads. Set PRUNE_BAM=0 to keep them (do that for a single diagnostic run, not a sweep).
prune_bams() {  # <outdir>
    local dir="$1" bam_kb
    [[ "${PRUNE_BAM:-1}" == "1" ]] || { echo "[prune] keeping BAMs (PRUNE_BAM=0)"; return 0; }
    if [[ "${INS_VERIFY_FAILED:-0}" -ne 0 ]]; then
        echo "[prune] KEEPING BAMs: INS read verification failed before evaluating every probe" >&2
        return 0
    fi
    # MUST NOT FIRE without a completed score. Deleting the BAMs of a run that died before
    # stage 5 destroys the only artifact that makes it re-callable, and the FASTQ is already
    # gone by this point, so there is no way back.
    if [[ ! -f "$dir/truvari_manta_overall/summary.json" ]]; then
        echo "[prune] KEEPING BAMs: no truvari_manta_overall/summary.json, so scoring did" >&2
        echo "  not complete and this is a run someone will need to re-call." >&2
        return 0
    fi
    bam_kb=$(du -sk "$dir"/normal.bam "$dir"/tumor.bam 2>/dev/null | awk '{s+=$1} END{print s+0}')
    rm -f "$dir"/normal.bam "$dir"/normal.bam.bai "$dir"/tumor.bam "$dir"/tumor.bam.bai
    echo "[prune] removed aligned BAMs, reclaimed $((bam_kb / 1048576)) GB" \
         "(PRUNE_BAM=0 keeps them; replicate stays re-scoreable, not re-callable)"
}
prune_bams "$OUTDIR"

# Persist report artifacts off scratch + record provenance/resources.
# (manta overall always exists; delly summaries remain under $OUTDIR/truvari_delly_*)
archive_run sv "$OUTDIR" truvari_manta_overall/summary.json scorer_recall.tsv

# ── Deferred failure gate ────────────────────────────────────────────────
# Deliberately last. Every discrepancy is reported inline next to the comparison it
# affects (so the scale is visible rather than just the first case), the artifacts are
# archived (so a failure is diagnosable after the fact), and only then does the job
# exit non-zero — because a run that published recalls over a subset of the truth, or
# scored with an uncalibrated scorer, must not look successful to whoever reads the
# SLURM exit status.
if [[ "$CAL_FAIL" -ne 0 || "$SCORE_COVERAGE_FAIL" -gt 0 || "${INS_VERIFY_FAILED:-0}" -ne 0 ]]; then
    echo "" >&2
    [[ "$CAL_FAIL" -ne 0 ]] && \
        echo "FAILED: scorer calibration (Rule 0) — no caller results were produced." >&2
    if [[ "$SCORE_COVERAGE_FAIL" -gt 0 ]]; then
        echo "FAILED: $SCORE_COVERAGE_FAIL comparison(s) did not score the whole truth set." >&2
        echo "  Those recalls are over a subset and must not be reported as-is. Raise" >&2
        echo "  TRUVARI_SIZEMAX (currently $TRUVARI_SIZEMAX) or relax the size filters," >&2
        echo "  then re-score. Artifacts are archived." >&2
    fi
    [[ "${INS_VERIFY_FAILED:-0}" -eq 0 ]] || \
        echo "FAILED: INS read verification did not evaluate every probe; BAMs were kept for diagnosis." >&2
    exit 1
fi
