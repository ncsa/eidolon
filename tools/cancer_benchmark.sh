#!/usr/bin/env bash
# Score a somatic-variant caller against an eidolon cancer-MVP truth VCF.
#
# Pipeline:
#   1. Index reference (bwa index + samtools faidx) if not already done.
#   2. Align the normal FASTQ → normal.bam (sorted, dedup'd, indexed).
#   3. Align the tumor (merged) FASTQ → tumor.bam (sorted, dedup'd, indexed).
#   4. Run GATK Mutect2 in tumor/normal mode against the indexed BAMs.
#   5. Filter the raw calls with GATK FilterMutectCalls.
#   6. Score the filtered VCF against the eidolon-emitted truth VCF using
#      Illumina's som.py inside the GA4GH-recommended hap.py image.
#
# Every external tool runs in a Docker (or Podman-as-Docker) container with
# the image tags pinned, so the benchmark is reproducible across machines.
# Steps are idempotent: outputs that already exist on disk are skipped, so
# re-running the script after a change to one step doesn't redo the whole
# pipeline.
#
# Expected inputs are the outputs of `tools/cancer_simulate.sh`:
#   - normal-fastq → <prefix>_normal_r1.fastq.gz  (pure-normal sample)
#   - tumor-fastq  → <prefix>_merged_r1.fastq.gz  (purity-mixed biopsy)
#   - truth-vcf    → <prefix>_merged_truth.vcf.gz (origin-tagged truth set)
#
# Usage:
#   tools/cancer_benchmark.sh \
#       --reference ~/code/data/chr22.fa \
#       --normal-fastq smoketest_normal_r1.fastq.gz \
#       --tumor-fastq  smoketest_merged_r1.fastq.gz \
#       --truth-vcf    smoketest_merged_truth.vcf.gz \
#       --output-dir   benchmark_out
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Pinned image versions ───────────────────────────────────────────────
# Pinned so a fresh checkout months from now reproduces the same scoring.
BWA_IMG="quay.io/biocontainers/bwa:0.7.18--he4a0461_0"
SAMTOOLS_IMG="quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
BCFTOOLS_IMG="quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
GATK_IMG="docker.io/broadinstitute/gatk:4.5.0.0"
HAPPY_IMG="docker.io/jmcdani20/hap.py:v0.3.12"

# ── Defaults ────────────────────────────────────────────────────────────
REFERENCE=""
NORMAL_FASTQ=""
TUMOR_FASTQ=""
TRUTH_VCF=""
OUTPUT_DIR="benchmark_out"
THREADS="4"
DOCKER="docker"
# Default truth filter: keep only records whose origin is "somatic" — those
# are the only variants Mutect2 in tumor/normal mode is expected to call.
# Without this filter, the truth VCF's ~99% germline-carry-through records
# (EIDOLON_ORIGIN=shared) inflate the FN count and tank apparent recall even
# when the somatic caller is doing fine. Set to empty (--truth-filter '')
# to skip filtering and score against the full truth as-is.
TRUTH_FILTER='INFO/EIDOLON_ORIGIN="somatic"'

usage() {
    cat <<'EOF'
Usage: cancer_benchmark.sh [options]

Required:
  --reference      Decompressed reference FASTA (.fa, not .fa.gz). Must be
                   on a path the Docker daemon can mount. If a .fai is
                   missing it will be built; same for the BWA index.
  --normal-fastq   Single-end FASTQ from the pure-normal pass of
                   cancer_simulate.sh (typically <prefix>_normal_r1.fastq.gz)
  --tumor-fastq    Single-end FASTQ representing the purity-mixed tumor
                   biopsy (typically <prefix>_merged_r1.fastq.gz, the
                   normal + tumor concatenation)
  --truth-vcf      The origin-tagged truth VCF emitted by the merge step
                   of cancer_simulate.sh (typically <prefix>_merged_truth.vcf.gz)

Output:
  --output-dir     Directory for intermediate BAMs and final scoring
                   outputs (default: benchmark_out, created if absent)

Resources:
  --threads        Threads passed to bwa / samtools / Mutect2 (default: 4)

Truth filtering:
  --truth-filter   bcftools view -i expression applied to the truth VCF
                   before scoring. Default: 'INFO/EIDOLON_ORIGIN="somatic"'
                   — keeps only the records Mutect2 in tumor/normal mode
                   is expected to call (germline-carried-through records
                   would otherwise be counted as FN). Pass an empty string
                   to skip filtering and score against the full truth.

Reproducibility:
  --docker         Container runtime to use (default: docker; podman in
                   compatibility mode also works)

  -h, --help       Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)     REFERENCE="$2"; shift 2 ;;
        --normal-fastq)  NORMAL_FASTQ="$2"; shift 2 ;;
        --tumor-fastq)   TUMOR_FASTQ="$2"; shift 2 ;;
        --truth-vcf)     TRUTH_VCF="$2"; shift 2 ;;
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --truth-filter)  TRUTH_FILTER="$2"; shift 2 ;;
        --docker)        DOCKER="$2"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# ── Validation ──────────────────────────────────────────────────────────
for var in REFERENCE NORMAL_FASTQ TUMOR_FASTQ TRUTH_VCF; do
    if [[ -z "${!var}" ]]; then
        echo "Missing required --$(echo "$var" | tr '[:upper:]_' '[:lower:]-')" >&2
        usage >&2
        exit 2
    fi
done

for f in "$REFERENCE" "$NORMAL_FASTQ" "$TUMOR_FASTQ" "$TRUTH_VCF"; do
    [[ -f "$f" ]] || { echo "File not found: $f" >&2; exit 2; }
done

# Reference must be uncompressed for BWA to index. If the caller passed a
# .fa.gz, point them at the decompression command rather than silently
# decompressing — it's a sizable artifact and they should choose where it
# lives.
if [[ "$REFERENCE" == *.gz ]]; then
    echo "Reference must be uncompressed for BWA indexing." >&2
    echo "  zcat \"$REFERENCE\" > \"\${REFERENCE%.gz}\"" >&2
    echo "  ...then re-run with --reference \"\${REFERENCE%.gz}\"" >&2
    exit 2
fi

if ! command -v "$DOCKER" >/dev/null 2>&1; then
    echo "$DOCKER not found on PATH; install it or pass --docker <runtime>" >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"

# Resolve to absolute paths so Docker volume mounts work regardless of
# the user's cwd. Reference and inputs may live in different directories;
# we mount each parent and refer to files via the in-container path.
REFERENCE="$(realpath "$REFERENCE")"
NORMAL_FASTQ="$(realpath "$NORMAL_FASTQ")"
TUMOR_FASTQ="$(realpath "$TUMOR_FASTQ")"
TRUTH_VCF="$(realpath "$TRUTH_VCF")"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"

REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"
NORMAL_DIR="$(dirname "$NORMAL_FASTQ")"
NORMAL_NAME="$(basename "$NORMAL_FASTQ")"
TUMOR_DIR="$(dirname "$TUMOR_FASTQ")"
TUMOR_NAME="$(basename "$TUMOR_FASTQ")"
TRUTH_DIR="$(dirname "$TRUTH_VCF")"
TRUTH_NAME="$(basename "$TRUTH_VCF")"

NORMAL_BAM="$OUTPUT_DIR/normal.bam"
TUMOR_BAM="$OUTPUT_DIR/tumor.bam"
RAW_VCF="$OUTPUT_DIR/mutect2.unfiltered.vcf.gz"
FILTERED_VCF="$OUTPUT_DIR/mutect2.filtered.vcf.gz"
SCORE_PREFIX="$OUTPUT_DIR/scored"

# ── Helpers ─────────────────────────────────────────────────────────────
# Run a single command inside a container with the reference dir mounted
# read-only at /refs and the output dir mounted read-write at /work.
run_in() {
    local image="$1"; shift
    "$DOCKER" run --rm \
        -v "$REF_DIR:/refs" \
        -v "$OUTPUT_DIR:/work" \
        -v "$NORMAL_DIR:/normal_in:ro" \
        -v "$TUMOR_DIR:/tumor_in:ro" \
        -v "$TRUTH_DIR:/truth_in:ro" \
        -w /work \
        "$image" \
        "$@"
}

# ── 1. Index reference (BWA + faidx) ────────────────────────────────────
# BWA produces five sibling files alongside the FASTA: .amb .ann .bwt .pac .sa.
# Skip if all are already present.
need_bwa_index=false
for ext in amb ann bwt pac sa; do
    [[ -f "${REFERENCE}.${ext}" ]] || need_bwa_index=true
done
if $need_bwa_index; then
    echo ">> Indexing reference for BWA (one-time, ~1 min per Gb)..."
    run_in "$BWA_IMG" bwa index "/refs/$REF_NAME"
else
    echo ">> BWA index present — skipping."
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo ">> Building .fai for reference..."
    run_in "$SAMTOOLS_IMG" samtools faidx "/refs/$REF_NAME"
else
    echo ">> .fai present — skipping."
fi

# GATK requires a sequence dictionary (.dict) next to the FASTA.
REF_DICT="${REFERENCE%.fa}.dict"
[[ -f "$REF_DICT" ]] || REF_DICT="${REFERENCE%.fasta}.dict"
if [[ ! -f "$REF_DICT" ]]; then
    echo ">> Building reference .dict for GATK..."
    run_in "$GATK_IMG" gatk CreateSequenceDictionary -R "/refs/$REF_NAME"
else
    echo ">> .dict present — skipping."
fi

# ── 2 & 3. Align FASTQ → sorted/indexed BAM ─────────────────────────────
# Single-end alignment (cancer_simulate.sh defaults to SE). The @RG SM tag
# becomes the sample name Mutect2 reads — must be distinct between normal
# and tumor for tumor/normal mode to work.
align_one() {
    local in_dir_alias="$1"  # /normal_in or /tumor_in
    local in_name="$2"        # FASTQ filename inside that dir
    local out_bam="$3"        # absolute host path under OUTPUT_DIR
    local sample_name="$4"    # NORMAL or TUMOR — used in @RG SM tag

    local out_name="$(basename "$out_bam")"

    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then
        echo ">> $out_name already present — skipping alignment."
        return 0
    fi

    echo ">> Aligning $sample_name FASTQ → $out_name..."
    local rg="@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${sample_name}\tPL:ILLUMINA"
    # Stream BWA → samtools sort: keeps the unsorted SAM out of memory
    # and avoids a multi-GB intermediate file. Two containers piped
    # together via a host fifo would be ideal but is awkward with
    # plain `docker run`; instead we do BWA → SAM in the BWA container,
    # then sort+index in the samtools container. The intermediate
    # unsorted BAM is ~2× the sorted BAM size but disposable.
    run_in "$BWA_IMG" sh -c \
        "bwa mem -t $THREADS -R '$rg' '/refs/$REF_NAME' '${in_dir_alias}/${in_name}' \
            > '/work/${out_name%.bam}.unsorted.sam'"
    run_in "$SAMTOOLS_IMG" sh -c \
        "samtools sort -@ $THREADS -o '/work/${out_name}' '/work/${out_name%.bam}.unsorted.sam' \
         && samtools index '/work/${out_name}' \
         && rm '/work/${out_name%.bam}.unsorted.sam'"
}

align_one "/normal_in" "$NORMAL_NAME" "$NORMAL_BAM" "NORMAL"
align_one "/tumor_in"  "$TUMOR_NAME"  "$TUMOR_BAM"  "TUMOR"

# ── 4. Mutect2 in tumor/normal mode ─────────────────────────────────────
if [[ -f "$RAW_VCF" ]]; then
    echo ">> Mutect2 unfiltered VCF already present — skipping."
else
    echo ">> Running Mutect2 (tumor/normal mode)..."
    run_in "$GATK_IMG" gatk Mutect2 \
        -R "/refs/$REF_NAME" \
        -I "/work/$(basename "$TUMOR_BAM")" \
        -I "/work/$(basename "$NORMAL_BAM")" \
        -tumor TUMOR \
        -normal NORMAL \
        --native-pair-hmm-threads "$THREADS" \
        -O "/work/$(basename "$RAW_VCF")"
fi

# ── 5. FilterMutectCalls ────────────────────────────────────────────────
if [[ -f "$FILTERED_VCF" ]]; then
    echo ">> FilterMutectCalls output already present — skipping."
else
    echo ">> Running FilterMutectCalls..."
    run_in "$GATK_IMG" gatk FilterMutectCalls \
        -R "/refs/$REF_NAME" \
        -V "/work/$(basename "$RAW_VCF")" \
        -O "/work/$(basename "$FILTERED_VCF")"
fi

# ── 6a. Filter truth VCF before scoring ─────────────────────────────────
# som.py treats every record in the truth VCF as an expected somatic
# call. The eidolon-emitted truth carries both germline-carry-through
# (EIDOLON_ORIGIN=shared) and somatic (EIDOLON_ORIGIN=somatic) records, so
# without filtering, Mutect2's (correct) refusal to call the germline
# carry-through variants shows up as ~99% false negatives. Filter the
# truth down to just the somatic ground-truth records before scoring.
SCORING_TRUTH="$OUTPUT_DIR/scoring_truth.vcf.gz"
SCORING_TRUTH_NAME="$(basename "$SCORING_TRUTH")"
if [[ -n "$TRUTH_FILTER" ]]; then
    echo ">> Filtering truth VCF with: $TRUTH_FILTER"
    run_in "$BCFTOOLS_IMG" sh -c \
        "bcftools view -i '$TRUTH_FILTER' -O z -o '/work/$SCORING_TRUTH_NAME' '/truth_in/$TRUTH_NAME' \
         && bcftools index -f -t '/work/$SCORING_TRUTH_NAME'"
    # If the filtered truth is empty the filter didn't match anything —
    # almost always means INFO/EIDOLON_ORIGIN isn't present on this truth
    # VCF (e.g. a non-eidolon truth set). Point the user at the escape
    # hatch rather than silently scoring against zero records.
    filtered_count=$(run_in "$BCFTOOLS_IMG" sh -c \
        "bcftools view -H '/work/$SCORING_TRUTH_NAME' | wc -l" | tr -d '\r\n')
    if [[ "$filtered_count" == "0" ]]; then
        echo "WARNING: truth filter '$TRUTH_FILTER' matched zero records." >&2
        echo "  If your truth VCF doesn't have INFO/EIDOLON_ORIGIN, pass --truth-filter ''" >&2
    else
        echo "    Filtered truth: $filtered_count records retained."
    fi
else
    echo ">> No truth filter applied — scoring against full truth VCF."
    cp "$TRUTH_VCF" "$SCORING_TRUTH"
    run_in "$BCFTOOLS_IMG" bcftools index -f -t "/work/$SCORING_TRUTH_NAME"
fi

# ── 6b. Score against the truth VCF with som.py ─────────────────────────
# som.py uses the reference for variant normalization (left-alignment,
# multi-allelic decomposition).
echo ">> Scoring filtered calls against truth VCF via som.py..."
"$DOCKER" run --rm \
    -v "$REF_DIR:/refs:ro" \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    "$HAPPY_IMG" \
    /opt/hap.py/bin/som.py \
        "/work/$SCORING_TRUTH_NAME" \
        "/work/$(basename "$FILTERED_VCF")" \
        -r "/refs/$REF_NAME" \
        -o "/work/$(basename "$SCORE_PREFIX")" \
        --feature-table generic

# ── Summary ─────────────────────────────────────────────────────────────
cat <<EOF

────────────────────────────────────────────────────────────────
Cancer benchmark complete

Inputs:
  Reference:    $REFERENCE
  Normal FASTQ: $NORMAL_FASTQ
  Tumor FASTQ:  $TUMOR_FASTQ
  Truth VCF:    $TRUTH_VCF

Outputs (in $OUTPUT_DIR):
  Normal BAM:        normal.bam   (+ .bai)
  Tumor BAM:         tumor.bam    (+ .bai)
  Raw Mutect2 calls: $(basename "$RAW_VCF")
  Filtered calls:    $(basename "$FILTERED_VCF")
  Score outputs:     $(basename "$SCORE_PREFIX").{stats.csv,metrics.json,vcf.gz}

The top-line numbers — recall, precision, F1 — are in:
  $SCORE_PREFIX.stats.csv

Look at the SNVs and INDELs rows in particular; the "Total" row is
the union.
────────────────────────────────────────────────────────────────
EOF
