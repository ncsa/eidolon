#!/usr/bin/env bash
# Score a somatic SV caller against an eidolon cancer-MVP truth VCF.
#
# This is the SV-caller companion to tools/cancer_benchmark.sh (which runs
# Mutect2 + som.py for SNVs/indels). It closes the stage-2 exit criterion
# in docs/cancer_simulator.md: "a current somatic SV caller (GRIDSS or
# Manta) calls translocation breakpoints from the simulated data with
# high enough recall that the simulator could plausibly be used as a
# benchmark fixture."
#
# Pipeline:
#   1. Index reference (bwa index + samtools faidx) if not already done.
#   2. Align the normal FASTQ → normal.bam, tumor FASTQ → tumor.bam
#      (sorted, indexed). Reuses BAMs in --output-dir if they already
#      exist — so this script chains cleanly with cancer_benchmark.sh,
#      which produces the same BAMs.
#   3. Run Illumina Manta in tumor/normal mode.
#   4. Filter the truth VCF to SV records (symbolic ALTs and long literal
#      indels with |ILEN| >= 50 — including v1.12.0's de novo INS records
#      that are emitted as literal-ALT) and optionally restrict to
#      INFO/EIDOLON_ORIGIN="somatic" so we measure somatic-SV-calling
#      performance, not germline-subtraction quality.
#   5. Score Manta's somatic SV calls against the filtered truth using
#      truvari bench.
#
# Every external tool runs in a Docker (or Podman-as-Docker) container
# with the image tags pinned, so the benchmark is reproducible across
# machines. Per-step outputs are idempotent (skipped if already present).
#
# Usage:
#   tools/cancer_sv_benchmark.sh \
#       --reference ~/code/data/chr22.fa \
#       --normal-fastq smoketest_normal_r1.fastq.gz \
#       --tumor-fastq  smoketest_merged_r1.fastq.gz \
#       --truth-vcf    smoketest_merged_truth.vcf.gz \
#       --output-dir   sv_benchmark_out
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Pinned image versions ───────────────────────────────────────────────
# Match cancer_benchmark.sh where they overlap so a single fetch covers both.
BWA_IMG="quay.io/biocontainers/bwa:0.7.18--he4a0461_0"
SAMTOOLS_IMG="quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
BCFTOOLS_IMG="quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
# Manta 1.6.0 is the last upstream release (project archived Sept 2024).
# Biocontainers tag includes the htslib + python2 deps Manta needs.
MANTA_IMG="quay.io/biocontainers/manta:1.6.0--h9ee0642_1"
# Truvari 4.3.1 — current stable. Used for SV comparison.
TRUVARI_IMG="quay.io/biocontainers/truvari:4.3.1--pyhdfd78af_0"

# ── Defaults ────────────────────────────────────────────────────────────
REFERENCE=""
NORMAL_FASTQ=""
NORMAL_FASTQ_R2=""
TUMOR_FASTQ=""
TUMOR_FASTQ_R2=""
TRUTH_VCF=""
OUTPUT_DIR="sv_benchmark_out"
THREADS="4"
DOCKER="docker"
# Default truth filter for SV scoring: keep only SYMBOLIC SV records
# (INFO/SVTYPE present — i.e. <DEL>/<DUP>/<INV>/<CNV>/<BND>/<INS> from the
# sv_model), AND restrict to somatic-tagged records so we measure somatic SV
# calling. We deliberately EXCLUDE literal large indels (the old filter's
# `abs(ILEN) >= 50` arm): those come from the SNP/indel side of the model, are
# scored by the SNV/indel benchmark (cancer_benchmark.sh), and a structural
# caller like Manta is not expected to emit them — leaving them in tanks the
# apparent recall (they all become false negatives). Set --truth-filter '' to
# score against the full truth (e.g. when the truth has no EIDOLON_ORIGIN/SVTYPE).
TRUTH_FILTER='INFO/SVTYPE!="." && INFO/EIDOLON_ORIGIN="somatic"'

usage() {
    cat <<'EOF'
Usage: cancer_sv_benchmark.sh [options]

Required:
  --reference        Decompressed reference FASTA (.fa, not .fa.gz). Must be
                     on a path the Docker daemon can mount. If a .fai is
                     missing it will be built; same for the BWA index.
  --normal-fastq     R1 of the pure-normal sample (typically <prefix>_normal_r1.fastq.gz)
  --normal-fastq-r2  R2 of the pure-normal sample (typically <prefix>_normal_r2.fastq.gz).
                     REQUIRED: Manta, DELLY, and GRIDSS all use pair-orientation
                     signal for SV detection and refuse to run on single-end data.
                     Re-run cancer_simulate.sh with `--paired-ended --fragment-mean
                     <N> --fragment-st-dev <M>` if your existing data is SE.
  --tumor-fastq      R1 of the purity-mixed tumor biopsy (typically
                     <prefix>_merged_r1.fastq.gz)
  --tumor-fastq-r2   R2 of the tumor sample (typically <prefix>_merged_r2.fastq.gz).
                     Also REQUIRED for the same reason.
  --truth-vcf        Origin-tagged truth VCF emitted by the merge step of
                     cancer_simulate.sh (typically <prefix>_merged_truth.vcf.gz)

Output:
  --output-dir     Directory for intermediate BAMs, Manta run dir, filtered
                   truth, and truvari scoring outputs (default: sv_benchmark_out).
                   Reuses existing normal.bam/tumor.bam if present — so this
                   script chains with cancer_benchmark.sh's output dir.

Resources:
  --threads        Threads passed to bwa / samtools / Manta (default: 4)

Truth filtering:
  --truth-filter   bcftools view -i expression applied to the truth VCF
                   before scoring. Default:
                     (TYPE="other" || abs(ILEN) >= 50) && INFO/EIDOLON_ORIGIN="somatic"
                   — keeps SV records that Manta in tumor/normal mode is
                   expected to call. Pass an empty string to skip filtering.

Reproducibility:
  --docker         Container runtime (default: docker; podman in compat mode also works)

  -h, --help       Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)     REFERENCE="$2"; shift 2 ;;
        --normal-fastq)     NORMAL_FASTQ="$2"; shift 2 ;;
        --normal-fastq-r2)  NORMAL_FASTQ_R2="$2"; shift 2 ;;
        --tumor-fastq)      TUMOR_FASTQ="$2"; shift 2 ;;
        --tumor-fastq-r2)   TUMOR_FASTQ_R2="$2"; shift 2 ;;
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
for var in REFERENCE NORMAL_FASTQ NORMAL_FASTQ_R2 TUMOR_FASTQ TUMOR_FASTQ_R2 TRUTH_VCF; do
    if [[ -z "${!var}" ]]; then
        flag="--$(echo "$var" | tr '[:upper:]_' '[:lower:]-')"
        echo "Missing required $flag" >&2
        if [[ "$var" == *"_R2" ]]; then
            echo "" >&2
            echo "Manta requires paired-end reads to determine fragment-size and" >&2
            echo "orientation distributions used for SV detection. If your existing" >&2
            echo "cancer_simulate.sh outputs are single-end, regenerate with:" >&2
            echo "  tools/cancer_simulate.sh ... \\" >&2
            echo "      --paired-ended --fragment-mean 350 --fragment-st-dev 50" >&2
        fi
        exit 2
    fi
done

for f in "$REFERENCE" "$NORMAL_FASTQ" "$NORMAL_FASTQ_R2" "$TUMOR_FASTQ" "$TUMOR_FASTQ_R2" "$TRUTH_VCF"; do
    [[ -f "$f" ]] || { echo "File not found: $f" >&2; exit 2; }
done

if [[ "$REFERENCE" == *.gz ]]; then
    echo "Reference must be uncompressed for BWA / Manta indexing." >&2
    echo "  zcat \"$REFERENCE\" > \"\${REFERENCE%.gz}\"" >&2
    exit 2
fi

if ! command -v "$DOCKER" >/dev/null 2>&1; then
    echo "$DOCKER not found on PATH; install it or pass --docker <runtime>" >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"

REFERENCE="$(realpath "$REFERENCE")"
NORMAL_FASTQ="$(realpath "$NORMAL_FASTQ")"
NORMAL_FASTQ_R2="$(realpath "$NORMAL_FASTQ_R2")"
TUMOR_FASTQ="$(realpath "$TUMOR_FASTQ")"
TUMOR_FASTQ_R2="$(realpath "$TUMOR_FASTQ_R2")"
TRUTH_VCF="$(realpath "$TRUTH_VCF")"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"

# Manta needs the R1 and R2 reads to live in directories the container
# can mount. Validate that R1 and R2 are siblings (typical cancer_simulate.sh
# layout) so we can mount one dir per sample instead of four.
if [[ "$(dirname "$NORMAL_FASTQ")" != "$(dirname "$NORMAL_FASTQ_R2")" ]]; then
    echo "Normal R1 and R2 must live in the same directory; got:" >&2
    echo "  R1: $NORMAL_FASTQ" >&2
    echo "  R2: $NORMAL_FASTQ_R2" >&2
    exit 2
fi
if [[ "$(dirname "$TUMOR_FASTQ")" != "$(dirname "$TUMOR_FASTQ_R2")" ]]; then
    echo "Tumor R1 and R2 must live in the same directory; got:" >&2
    echo "  R1: $TUMOR_FASTQ" >&2
    echo "  R2: $TUMOR_FASTQ_R2" >&2
    exit 2
fi

REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"
NORMAL_DIR="$(dirname "$NORMAL_FASTQ")"
NORMAL_NAME="$(basename "$NORMAL_FASTQ")"
NORMAL_NAME_R2="$(basename "$NORMAL_FASTQ_R2")"
TUMOR_DIR="$(dirname "$TUMOR_FASTQ")"
TUMOR_NAME="$(basename "$TUMOR_FASTQ")"
TUMOR_NAME_R2="$(basename "$TUMOR_FASTQ_R2")"
TRUTH_DIR="$(dirname "$TRUTH_VCF")"
TRUTH_NAME="$(basename "$TRUTH_VCF")"

NORMAL_BAM="$OUTPUT_DIR/normal.bam"
TUMOR_BAM="$OUTPUT_DIR/tumor.bam"
MANTA_DIR="$OUTPUT_DIR/manta_run"
SOMATIC_SV_VCF="$MANTA_DIR/results/variants/somaticSV.vcf.gz"
SCORING_TRUTH="$OUTPUT_DIR/scoring_truth_svs.vcf.gz"
TRUVARI_DIR="$OUTPUT_DIR/truvari_out"

# ── Helpers ─────────────────────────────────────────────────────────────
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
need_bwa_index=false
for ext in amb ann bwt pac sa; do
    [[ -f "${REFERENCE}.${ext}" ]] || need_bwa_index=true
done
if $need_bwa_index; then
    echo ">> Indexing reference for BWA..."
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

# ── 2 & 3. Align paired-end FASTQ → sorted/indexed BAM ──────────────────
# BWA-MEM with both R1 and R2 sets the proper PE flags (0x1 paired, 0x2
# proper-pair, 0x40 first-of-pair / 0x80 second-of-pair) that Manta needs
# to determine fragment-orientation statistics.
align_one() {
    local in_dir_alias="$1" in_name_r1="$2" in_name_r2="$3" out_bam="$4" sample_name="$5"
    local out_name
    out_name="$(basename "$out_bam")"

    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then
        echo ">> $out_name already present — skipping alignment."
        return 0
    fi

    echo ">> Aligning $sample_name FASTQ (PE) -> $out_name..."
    local rg="@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${sample_name}\tPL:ILLUMINA"
    run_in "$BWA_IMG" sh -c \
        "bwa mem -t $THREADS -R '$rg' '/refs/$REF_NAME' \
            '${in_dir_alias}/${in_name_r1}' '${in_dir_alias}/${in_name_r2}' \
            > '/work/${out_name%.bam}.unsorted.sam'"
    run_in "$SAMTOOLS_IMG" sh -c \
        "samtools sort -@ $THREADS -o '/work/${out_name}' '/work/${out_name%.bam}.unsorted.sam' \
         && samtools index '/work/${out_name}' \
         && rm '/work/${out_name%.bam}.unsorted.sam'"
}

align_one "/normal_in" "$NORMAL_NAME" "$NORMAL_NAME_R2" "$NORMAL_BAM" "NORMAL"
align_one "/tumor_in"  "$TUMOR_NAME"  "$TUMOR_NAME_R2"  "$TUMOR_BAM"  "TUMOR"

# ── 4. Manta tumor/normal SV calling ────────────────────────────────────
if [[ -f "$SOMATIC_SV_VCF" ]]; then
    echo ">> Manta somatic SV output already present — skipping."
else
    if [[ -d "$MANTA_DIR" ]]; then
        # Stale partial run from a previous attempt — clear it so configManta
        # doesn't refuse to overwrite.
        echo ">> Clearing stale Manta run dir $MANTA_DIR"
        rm -rf "$MANTA_DIR"
    fi
    echo ">> Configuring Manta..."
    run_in "$MANTA_IMG" configManta.py \
        --tumorBam "/work/$(basename "$TUMOR_BAM")" \
        --normalBam "/work/$(basename "$NORMAL_BAM")" \
        --referenceFasta "/refs/$REF_NAME" \
        --runDir "/work/$(basename "$MANTA_DIR")"
    echo ">> Running Manta workflow (~minutes for chr22-scale)..."
    run_in "$MANTA_IMG" "/work/$(basename "$MANTA_DIR")/runWorkflow.py" -j "$THREADS"
fi

# ── 5a. Filter truth VCF to SV records ──────────────────────────────────
SCORING_TRUTH_NAME="$(basename "$SCORING_TRUTH")"
if [[ -f "$SCORING_TRUTH" ]]; then
    echo ">> Filtered truth-SV VCF already present — skipping."
else
    # See cancer_benchmark.sh: bcftools rejects an undefined INFO tag in a -i
    # expression at parse time, so probe the header first — otherwise the advice
    # below is unreachable for a truth VCF that simply lacks the tag.
    if [[ -n "$TRUTH_FILTER" && "$TRUTH_FILTER" == *EIDOLON_ORIGIN* ]]; then
        truth_header="$(run_in "$BCFTOOLS_IMG" bcftools view -h "/truth_in/$TRUTH_NAME")"
        if ! grep -q '^##INFO=<ID=EIDOLON_ORIGIN,' <<<"$truth_header"; then
            if grep -q '^##INFO=<ID=NEAT_ORIGIN,' <<<"$truth_header"; then
                echo "ERROR: $TRUTH_NAME carries the pre-v3.0.0 NEAT_ORIGIN tag." >&2
                echo "  eidolon v3.0.0 renamed the emitted output tokens. Convert it once:" >&2
                echo "    tools/migrate_legacy_tokens.sh '$TRUTH_VCF' migrated_truth.vcf.gz" >&2
                echo "    $0 --truth-vcf migrated_truth.vcf.gz ..." >&2
            else
                echo "ERROR: $TRUTH_NAME has no INFO/EIDOLON_ORIGIN tag." >&2
                echo "  If this is a non-eidolon truth set, pass --truth-filter '' or a" >&2
                echo "  custom expression." >&2
            fi
            exit 1
        fi
    fi

    if [[ -n "$TRUTH_FILTER" ]]; then
        echo ">> Filtering truth VCF for SVs with: $TRUTH_FILTER"
        run_in "$BCFTOOLS_IMG" sh -c \
            "bcftools view -i '$TRUTH_FILTER' -O z -o '/work/$SCORING_TRUTH_NAME' '/truth_in/$TRUTH_NAME' \
             && bcftools index -f -t '/work/$SCORING_TRUTH_NAME'"
        filtered_count=$(run_in "$BCFTOOLS_IMG" sh -c \
            "bcftools view -H '/work/$SCORING_TRUTH_NAME' | wc -l" | tr -d '\r\n')
        if [[ "$filtered_count" == "0" ]]; then
            echo "WARNING: truth filter matched zero SV records." >&2
            echo "  If your truth VCF doesn't carry SVs (or doesn't have INFO/EIDOLON_ORIGIN)," >&2
            echo "  pass --truth-filter '' or a custom expression." >&2
        else
            echo "    Filtered truth: $filtered_count SV records retained."
        fi
    else
        echo ">> No truth filter applied — scoring against full truth VCF."
        cp "$TRUTH_VCF" "$SCORING_TRUTH"
        run_in "$BCFTOOLS_IMG" bcftools index -f -t "/work/$SCORING_TRUTH_NAME"
    fi
fi

# ── 5b. Score Manta calls against truth-SV VCF with truvari ─────────────
# Truvari requires both VCFs bgzipped + tabix-indexed and a reference for
# sequence-aware comparison. somaticSV.vcf.gz already meets that; our
# filtered truth was indexed in step 5a.
if [[ -d "$TRUVARI_DIR" ]]; then
    echo ">> Clearing previous truvari output dir $TRUVARI_DIR"
    rm -rf "$TRUVARI_DIR"
fi
echo ">> Scoring Manta somaticSV against filtered truth via truvari bench..."
# Flags tuned for REALISTIC symbolic SVs (PCAWG-derived sizes: DEL median ~115kb,
# events up to several Mb). truvari's defaults were calibrated for small,
# sequence-resolved SVs and silently drop ours:
#   --pctseq 0       symbolic <DEL>/<DUP>/... carry no ALT sequence to compare;
#                    the default 0.7 filters every unresolved SV (recall→0).
#   --sizemax 3e8    default 50,000 excludes our large SVs from BOTH truth and
#                    calls; raise above chr length so nothing is size-excluded.
#   --typeignore     Manta cross-classifies junctions (e.g. inversions emitted
#                    as BND pairs, BNDs as DEL/DUP); match on position+size, not
#                    label. Mirrors the v1.13.1 cross-type scoring rationale.
#   --refdist 500    breakpoint tolerance, matching v1.13.1's ±500bp.
# NOTE: truvari still does not score BND-type records — see the BND-specific
# positional pass (step 5c) below for translocation recall.
SOMATIC_SV_REL="${SOMATIC_SV_VCF#$OUTPUT_DIR/}"
run_in "$TRUVARI_IMG" truvari bench \
    -b "/work/$SCORING_TRUTH_NAME" \
    -c "/work/$SOMATIC_SV_REL" \
    -o "/work/$(basename "$TRUVARI_DIR")" \
    -f "/refs/$REF_NAME" \
    --passonly --pctseq 0 --typeignore --refdist 500 --pctsize 0.5 \
    --sizemax 300000000 --sizefilt 50

# ── 5c. Positional recall for BND + INV (truvari can't score either) ─────
# Translocations and inversions are detected by Manta but emitted as BND-form
# junctions, which truvari drops (and re-labels as DEL/DUP). Reproduce
# v1.13.1's positional matcher: a somatic truth record is "recovered" if Manta
# emits ANY PASS breakpoint within ±BP_WINDOW bp of one of its breakpoints
# (POS for BND; POS or END for INV, whose two junctions Manta calls
# separately). Pure awk on bcftools-query output.
BP_WINDOW=500
# Manta PASS breakpoints: POS for every call, plus END for span calls (DEL/DUP)
# — a re-labeled junction boundary still counts.
MANTA_BP_POS="$OUTPUT_DIR/manta_breakpoints.txt"
run_in "$BCFTOOLS_IMG" sh -c \
    "bcftools view -H -f PASS '/work/$SOMATIC_SV_REL' \
     | awk '{print \$1\"\t\"\$2; if(match(\$8,/END=[0-9]+/)){e=substr(\$8,RSTART+4,RLENGTH-4); print \$1\"\t\"e}}'" \
     > "$MANTA_BP_POS" || true
score_positional() {  # svtype  check_end(0|1)
    local svtype="$1" check_end="$2"
    local truth_pos="$OUTPUT_DIR/truth_${svtype}_pos.txt"
    run_in "$BCFTOOLS_IMG" sh -c \
        "bcftools view -H -i 'INFO/SVTYPE=\"$svtype\" && INFO/EIDOLON_ORIGIN=\"somatic\"' '/truth_in/$TRUTH_NAME' \
         | awk '{p=\$2; e=p; if(match(\$8,/END=[0-9]+/))e=substr(\$8,RSTART+4,RLENGTH-4); print \$1\"\t\"p\"\t\"e}'" \
        > "$truth_pos" || true
    awk -v W="$BP_WINDOW" -v CE="$check_end" '
        NR==FNR { split($0,a,"\t"); pos[$1] = pos[$1]" "a[2]; next }
        { tot++; split(pos[$1], cand, " "); found=0
          for (i in cand) if (cand[i]!="") {
              if ((($2-cand[i])<=W && (cand[i]-$2)<=W) ||
                  (CE=="1" && ($3-cand[i])<=W && (cand[i]-$3)<=W)) { found=1; break } }
          hit += found }
        END { if (tot>0) printf "%d/%d = %.0f%%", hit+0, tot, 100*hit/tot; else print "none in truth" }
    ' "$MANTA_BP_POS" "$truth_pos"
}
echo ">> Scoring somatic BND + INV recall positionally (±${BP_WINDOW}bp, cross-type)..."
BND_RECALL=$(score_positional BND 0)
INV_RECALL=$(score_positional INV 1)
echo "    Somatic BND positional recall: $BND_RECALL"
echo "    Somatic INV positional recall: $INV_RECALL"

# ── Summary ─────────────────────────────────────────────────────────────
cat <<EOF

────────────────────────────────────────────────────────────────
Cancer SV benchmark complete

Inputs:
  Reference:    $REFERENCE
  Normal FASTQ: $NORMAL_FASTQ
  Tumor FASTQ:  $TUMOR_FASTQ
  Truth VCF:    $TRUTH_VCF

Outputs (in $OUTPUT_DIR):
  Normal BAM:        normal.bam   (+ .bai)
  Tumor BAM:         tumor.bam    (+ .bai)
  Manta run dir:     manta_run/
  Manta somatic SVs: manta_run/results/variants/somaticSV.vcf.gz
  Filtered truth:    $SCORING_TRUTH_NAME
  Scoring outputs:   $(basename "$TRUVARI_DIR")/

Top-line numbers — precision, recall, F1 — are in:
  $TRUVARI_DIR/summary.json

Read the per-SV-type breakdown there (TP/FP/FN counts by SVTYPE);
truvari reports them separately so BND / INV / DEL / DUP / INS
performance can be compared.
────────────────────────────────────────────────────────────────
EOF
