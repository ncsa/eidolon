#!/usr/bin/env bash
# Simulate a tumor / normal mixture sample by running eidolon gen-reads twice
# over the same reference and concatenating the resulting FASTQs.
#
# Architecture (per docs/cancer_simulator.md):
#
#   Pass 1 (normal): coverage = (1 - purity) × total_coverage.
#     Produces a "germline" FASTQ + a golden VCF of randomly-generated
#     germline variants (the "patient's genome").
#
#   Pass 2 (tumor):  coverage = purity × total_coverage.
#     Re-uses pass 1's golden VCF as input_vcf — guaranteeing the tumor
#     cells carry the same germline variants as the normal cells. Adds
#     tumor-specific somatic mutations on top (via the tumor mutation
#     model and optional sv_rate_scale).
#
#   Merge: concatenate FASTQs (independent reads, just append). The two
#     per-pass golden VCFs survive for downstream truth-set work; merging
#     them with origin annotation is eidolon issue #185, not this script.
#
# Ported from NEAT 2.1's `genReadsTumorTutorial/simulate.sh` idiom.
#
# Usage:
#   tools/cancer_simulate.sh \
#       --reference /path/to/ref.fa \
#       --output-prefix /path/to/out/sample \
#       --total-coverage 30 \
#       --purity 0.7
#
# Run with --help for the full flag list.

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────
REFERENCE=""
OUTPUT_PREFIX=""
TOTAL_COVERAGE="30"
PURITY="0.5"
NORMAL_MODEL=""
TUMOR_MODEL=""
NORMAL_MUTATION_RATE=""
# Default the tumor pass to a realistic per-tumor SOMATIC rate rather than the
# model's fitted rate. The de-novo mutations added in the tumor pass are somatic
# (germline is carried through via input_vcf), and most tumor models — including
# the bundled COSMIC pan-cancer model (~5.5e-3) — carry a corpus-aggregated rate
# that overstates single-tumor burden by ~50-500x. 1e-5 ≈ a typical solid tumor.
# Pass `--tumor-mutation-rate model` to defer to the model's own fitted rate.
TUMOR_MUTATION_RATE="1e-5"
GERMLINE_VCF=""
READ_LEN="151"
PAIRED_END="false"
FRAGMENT_MEAN=""
FRAGMENT_ST_DEV=""
SV_RATE_SCALE="0.0"
RNG_SEED_ROOT="cancer-simulate"
EIDOLON_BIN="eidolon"

usage() {
    cat <<'EOF'
Usage: cancer_simulate.sh [options]

Required:
  --reference        FASTA reference (.fa or .fa.gz)
  --output-prefix    Path stem for outputs (directory must be writable; created if absent).
                     Produces:
                       <prefix>_normal_r1.fastq.gz  (and _r2 if paired)
                       <prefix>_normal.vcf.gz       (germline truth)
                       <prefix>_tumor_r1.fastq.gz   (and _r2 if paired)
                       <prefix>_tumor.vcf.gz        (germline + somatic truth)
                       <prefix>_merged_r1.fastq.gz  (and _r2 if paired)

Coverage / purity:
  --total-coverage   Combined coverage of the merged output (default: 30)
  --purity           Tumor cell fraction in [0.0, 1.0]                (default: 0.5)
                     Tumor pass gets purity × total, normal pass gets the rest.

Mutation models (optional; both default to eidolon's bundled MutationModel::default()
which is germline-derived. Supplying a cancer-trained tumor model is strongly
recommended for any non-toy run — see issue #186):
  --normal-model     Path to a .json.gz mutation model for the normal pass
  --tumor-model      Path to a .json.gz mutation model for the tumor pass
  --seq-error-model  Path to a .json.gz sequencing error model, used by BOTH passes.
                     Defaults to the model compiled into the binary.

Mutation-rate overrides:
  --normal-mutation-rate  Override the normal pass's per-base mutation rate.
                          Omitted → the normal model's fitted (germline) rate,
                          which is the right magnitude for germline burden.
  --tumor-mutation-rate   The tumor pass's per-base SOMATIC mutation rate.
                          DEFAULT: 1e-5 (typical solid tumor). The de-novo
                          mutations added in the tumor pass are somatic, and most
                          tumor models carry a corpus-aggregated rate that
                          overstates single-tumor burden by ~50–500× — e.g. the
                          bundled COSMIC pan-cancer model fits ~5.5e-3 per bp
                          ("fraction of bp mutated anywhere across the COSMIC
                          catalog"), which would emit ~hundreds of thousands of
                          SNVs per tumor. Realistic per-tumor rates: ~1e-5
                          (typical solid tumor) to ~1e-4 (high-burden / MSI).
                          Pass `--tumor-mutation-rate model` to instead defer to
                          the tumor model's own fitted rate.

Germline VCF (optional):
  --germline-vcf     A VCF to use as the shared germline. If omitted, pass 1
                     generates one de novo and pass 2 consumes it. Supply this
                     when you want a specific germline genotype (or when running
                     multiple tumor scenarios over the same germline).

Read parameters:
  --read-len         Read length in bp (default: 151)
  --paired-ended     Generate paired-end reads. Requires --fragment-mean and
                     --fragment-st-dev. (default: single-ended)
  --fragment-mean    Mean insert length, paired-end only
  --fragment-st-dev  Insert-length stddev, paired-end only

Tumor-pass SV generation:
  --sv-rate-scale    Multiplier for the tumor model's per_base_rate (default: 0,
                     SVs from the model are disabled). Set to 1.0 to use the
                     model's nominal rate; higher values stress-test SV callers.

                     As of v1.12.0, gen-reads also emits BND translocations,
                     INV inversions, and de novo INS (literal insertions with
                     novel sequence) when --sv-rate-scale > 0 AND the supplied
                     mutation model has those SV types in its pool.

                     As of v1.12.1, the bundled COSMIC tumor model
                     (tools/cosmic_v104_pancancer_model.json.gz) carries a
                     literature-derived (PCAWG, Li et al. 2020) sv_model
                     component covering all six SV types. The germline default
                     (default_sv_model) was already populated; before v1.12.1
                     the COSMIC model's sv_model was null, so the tumor pass
                     silently produced zero symbolic SVs. See #218 for the
                     follow-up plan to refit from a real cancer-SV corpus.

                     Known caveat: chimeric BND/INV junction reads are
                     generated on top of regular reference reads at the same
                     breakpoint positions, so a homozygous BND or INV ends up
                     with ~2x coverage at the junction (regular reads from the
                     unbroken reference + junction-spanning reads). Het lands
                     closer to correct. Downstream SV callers handle this fine
                     (they expect both reference and junction-supporting reads
                     at heterozygous loci), but coverage-depth plots will
                     show a bump.

Reproducibility / wiring:
  --rng-seed         Seed root (default: "cancer-simulate"); per-pass seeds
                     append "-normal" and "-tumor" so the two passes don't
                     share RNG state.
  --eidolon-bin        Path to the eidolon binary (default: "eidolon" on PATH)

  -h, --help         Print this message
EOF
}

# ── Arg parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)        REFERENCE="$2"; shift 2 ;;
        --output-prefix)    OUTPUT_PREFIX="$2"; shift 2 ;;
        --total-coverage)   TOTAL_COVERAGE="$2"; shift 2 ;;
        --purity)           PURITY="$2"; shift 2 ;;
        --normal-model)         NORMAL_MODEL="$2"; shift 2 ;;
        --tumor-model)          TUMOR_MODEL="$2"; shift 2 ;;
        --seq-error-model)      SEQ_ERROR_MODEL="$2"; shift 2 ;;
        --normal-mutation-rate) NORMAL_MUTATION_RATE="$2"; shift 2 ;;
        --tumor-mutation-rate)  TUMOR_MUTATION_RATE="$2"; shift 2 ;;
        --germline-vcf)     GERMLINE_VCF="$2"; shift 2 ;;
        --read-len)         READ_LEN="$2"; shift 2 ;;
        --paired-ended)     PAIRED_END="true"; shift ;;
        --fragment-mean)    FRAGMENT_MEAN="$2"; shift 2 ;;
        --fragment-st-dev)  FRAGMENT_ST_DEV="$2"; shift 2 ;;
        --sv-rate-scale)    SV_RATE_SCALE="$2"; shift 2 ;;
        --rng-seed)         RNG_SEED_ROOT="$2"; shift 2 ;;
        --eidolon-bin)        EIDOLON_BIN="$2"; shift 2 ;;
        -h|--help)          usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# ── Validation ──────────────────────────────────────────────────────────
[[ -z "$REFERENCE" ]] && { echo "Missing --reference" >&2; exit 2; }
[[ ! -f "$REFERENCE" ]] && { echo "Reference file not found: $REFERENCE" >&2; exit 2; }
[[ -z "$OUTPUT_PREFIX" ]] && { echo "Missing --output-prefix" >&2; exit 2; }

# Purity in (0, 1) — endpoints would zero one of the passes; refuse rather
# than silently skip a pass (use plain eidolon gen-reads for purity=0 / 1).
if ! awk -v p="$PURITY" 'BEGIN { exit !(p > 0 && p < 1) }'; then
    echo "Purity must be strictly between 0 and 1; got: $PURITY" >&2
    echo "  (For purity=0 or 1 use eidolon gen-reads directly — no mixing needed.)" >&2
    exit 2
fi

if [[ "$PAIRED_END" == "true" ]]; then
    [[ -z "$FRAGMENT_MEAN"   ]] && { echo "--paired-ended requires --fragment-mean"   >&2; exit 2; }
    [[ -z "$FRAGMENT_ST_DEV" ]] && { echo "--paired-ended requires --fragment-st-dev" >&2; exit 2; }
fi

if ! command -v "$EIDOLON_BIN" >/dev/null 2>&1; then
    echo "eidolon binary not found on PATH (override with --eidolon-bin): $EIDOLON_BIN" >&2
    exit 2
fi

# ── Per-pass coverage ────────────────────────────────────────────────────
# Round to nearest integer; reject if either pass would round to zero.
NORMAL_COVERAGE=$(awk -v t="$TOTAL_COVERAGE" -v p="$PURITY" 'BEGIN { printf "%d", t * (1 - p) + 0.5 }')
TUMOR_COVERAGE=$(awk -v t="$TOTAL_COVERAGE" -v p="$PURITY" 'BEGIN { printf "%d", t * p + 0.5 }')

if [[ "$NORMAL_COVERAGE" -lt 1 || "$TUMOR_COVERAGE" -lt 1 ]]; then
    echo "Per-pass coverage rounds to zero (normal=$NORMAL_COVERAGE, tumor=$TUMOR_COVERAGE)." >&2
    echo "  Bump --total-coverage or pick a less extreme --purity." >&2
    exit 2
fi

# ── Output layout ────────────────────────────────────────────────────────
OUTPUT_DIR="$(dirname "$OUTPUT_PREFIX")"
PREFIX_NAME="$(basename "$OUTPUT_PREFIX")"
mkdir -p "$OUTPUT_DIR"

NORMAL_PREFIX="${PREFIX_NAME}_normal"
TUMOR_PREFIX="${PREFIX_NAME}_tumor"
MERGED_PREFIX="${PREFIX_NAME}_merged"

# ── Pass-1 config writer ─────────────────────────────────────────────────
write_config() {
    local config_path="$1"
    local prefix_name="$2"
    local coverage="$3"
    local model_path="$4"
    local input_vcf="$5"
    local seed_suffix="$6"
    local sv_scale="$7"
    local mutation_rate="$8"

    cat > "$config_path" <<EOF
reference: $REFERENCE
read_len: $READ_LEN
coverage: $coverage
ploidy: 2
paired_ended: $PAIRED_END
EOF
    if [[ "$PAIRED_END" == "true" ]]; then
        cat >> "$config_path" <<EOF
fragment_mean: $FRAGMENT_MEAN
fragment_st_dev: $FRAGMENT_ST_DEV
EOF
    fi
    cat >> "$config_path" <<EOF
produce_fastq: true
produce_vcf: true
overwrite_output: true
output_dir: $OUTPUT_DIR
output_filename: $prefix_name
rng_seed: ${RNG_SEED_ROOT}-${seed_suffix}
sv_rate_scale: $sv_scale
EOF
    # NOTE: don't use `[[ -n "$x" ]] && echo ...` here. A false `[[ ]]` test
    # returns exit 1, and with `set -e` a function whose last command exits
    # non-zero aborts the whole script — even if the && short-circuit was the
    # intended behavior. Spell it out with `if`.
    if [[ -n "$model_path" ]]; then
        echo "mutation_model: $model_path" >> "$config_path"
    fi
    # Both passes get the same sequencing error model. It describes the instrument,
    # not the sample, so a tumor/normal pair that differed here would be comparing
    # two machines rather than two genomes.
    if [[ -n "${SEQ_ERROR_MODEL:-}" ]]; then
        echo "sequence_error_model: $SEQ_ERROR_MODEL" >> "$config_path"
    fi
    if [[ -n "$input_vcf" ]]; then
        echo "input_vcf: $input_vcf" >> "$config_path"
    fi
    if [[ -n "$mutation_rate" ]]; then
        echo "mutation_rate: $mutation_rate" >> "$config_path"
    fi
    return 0
}

# ── Pass 1: normal ───────────────────────────────────────────────────────
echo ">> Pass 1: normal at ${NORMAL_COVERAGE}× coverage"
NORMAL_CONFIG="${OUTPUT_DIR}/${NORMAL_PREFIX}.yml"
write_config "$NORMAL_CONFIG" "$NORMAL_PREFIX" "$NORMAL_COVERAGE" \
    "$NORMAL_MODEL" "$GERMLINE_VCF" "normal" "0.0" "$NORMAL_MUTATION_RATE"
"$EIDOLON_BIN" gen-reads -c "$NORMAL_CONFIG"

# If the user didn't supply a germline VCF, use pass 1's golden as the
# input for pass 2 — that's the whole point of running two passes.
if [[ -z "$GERMLINE_VCF" ]]; then
    GERMLINE_VCF="${OUTPUT_DIR}/${NORMAL_PREFIX}.vcf.gz"
    if [[ ! -f "$GERMLINE_VCF" ]]; then
        echo "Expected pass-1 golden VCF at $GERMLINE_VCF but file not found." >&2
        echo "  Check that produce_vcf is honored in your gen-reads build." >&2
        exit 3
    fi
fi

# ── Pass 2: tumor ────────────────────────────────────────────────────────
# `--tumor-mutation-rate model` opts out of the realistic-rate default and lets
# the tumor pass use whatever rate the model fitted (write_config omits the
# mutation_rate line when the value is empty).
if [[ "$TUMOR_MUTATION_RATE" == "model" ]]; then
    TUMOR_MUTATION_RATE=""
    echo ">> Tumor mutation rate: using the model's fitted rate (--tumor-mutation-rate model)"
else
    echo ">> Tumor mutation rate: ${TUMOR_MUTATION_RATE} per bp (somatic; override with --tumor-mutation-rate)"
fi
echo ">> Pass 2: tumor at ${TUMOR_COVERAGE}× coverage (germline from ${GERMLINE_VCF})"
TUMOR_CONFIG="${OUTPUT_DIR}/${TUMOR_PREFIX}.yml"
write_config "$TUMOR_CONFIG" "$TUMOR_PREFIX" "$TUMOR_COVERAGE" \
    "$TUMOR_MODEL" "$GERMLINE_VCF" "tumor" "$SV_RATE_SCALE" "$TUMOR_MUTATION_RATE"
"$EIDOLON_BIN" gen-reads -c "$TUMOR_CONFIG"

# ── Merge FASTQs ─────────────────────────────────────────────────────────
# gen-reads names reads by genomic position (`EIDOLON_generated_<contig>_
# <start>_<end>/<mate>`), so if the normal and tumor passes happen to
# sample reads at the same start coordinate they get IDENTICAL names. On
# chr22 at 30× combined coverage we observed ~456k collisions out of
# ~9.6M reads, which Mutect2 warns about and which Picard MarkDuplicates
# would silently drop (halving effective coverage at those loci).
#
# Fix: stream each per-pass FASTQ through awk and prefix the read name
# (line 1 of every 4-line record) with `N_` or `T_` before re-bgzipping.
# Quality lines (line 4) can start with `@` too — that's why we filter
# by line position within each record, not by leading character.
prefix_reads() {
    local in_fq="$1" out_fq="$2" tag="$3"
    zcat "$in_fq" | awk -v tag="$tag" '
        NR % 4 == 1 { sub(/^@/, "@" tag "_") }
        { print }
    ' | gzip -c > "$out_fq"
}

NORMAL_R1="${OUTPUT_DIR}/${NORMAL_PREFIX}_r1.fastq.gz"
TUMOR_R1="${OUTPUT_DIR}/${TUMOR_PREFIX}_r1.fastq.gz"
MERGED_R1="${OUTPUT_DIR}/${MERGED_PREFIX}_r1.fastq.gz"
NORMAL_R1_TAGGED="${OUTPUT_DIR}/${NORMAL_PREFIX}_r1.tagged.fastq.gz"
TUMOR_R1_TAGGED="${OUTPUT_DIR}/${TUMOR_PREFIX}_r1.tagged.fastq.gz"
prefix_reads "$NORMAL_R1" "$NORMAL_R1_TAGGED" "N"
prefix_reads "$TUMOR_R1"  "$TUMOR_R1_TAGGED"  "T"
cat "$NORMAL_R1_TAGGED" "$TUMOR_R1_TAGGED" > "$MERGED_R1"
rm -f "$NORMAL_R1_TAGGED" "$TUMOR_R1_TAGGED"

if [[ "$PAIRED_END" == "true" ]]; then
    NORMAL_R2="${OUTPUT_DIR}/${NORMAL_PREFIX}_r2.fastq.gz"
    TUMOR_R2="${OUTPUT_DIR}/${TUMOR_PREFIX}_r2.fastq.gz"
    MERGED_R2="${OUTPUT_DIR}/${MERGED_PREFIX}_r2.fastq.gz"
    NORMAL_R2_TAGGED="${OUTPUT_DIR}/${NORMAL_PREFIX}_r2.tagged.fastq.gz"
    TUMOR_R2_TAGGED="${OUTPUT_DIR}/${TUMOR_PREFIX}_r2.tagged.fastq.gz"
    prefix_reads "$NORMAL_R2" "$NORMAL_R2_TAGGED" "N"
    prefix_reads "$TUMOR_R2"  "$TUMOR_R2_TAGGED"  "T"
    cat "$NORMAL_R2_TAGGED" "$TUMOR_R2_TAGGED" > "$MERGED_R2"
    rm -f "$NORMAL_R2_TAGGED" "$TUMOR_R2_TAGGED"
fi

# ── Merge golden VCFs into a single origin-tagged truth set ─────────────
# Resolves eidolon's per-pass `INFO/EIDOLON_PROVENANCE` (which says "denovo or
# input") into the cancer-specific `INFO/EIDOLON_ORIGIN` (which says "germline,
# somatic, or shared"):
#
#   only in normal_golden                → EIDOLON_ORIGIN=germline
#   only in tumor_golden  (denovo there) → EIDOLON_ORIGIN=somatic
#   in both passes                       → EIDOLON_ORIGIN=shared
#
# bcftools is required for the set operations. If it isn't installed we
# leave the per-pass VCFs and warn — the merge is a post-processing step,
# not part of the simulation correctness.
NORMAL_VCF="${OUTPUT_DIR}/${NORMAL_PREFIX}.vcf.gz"
TUMOR_VCF="${OUTPUT_DIR}/${TUMOR_PREFIX}.vcf.gz"
MERGED_TRUTH_VCF="${OUTPUT_DIR}/${MERGED_PREFIX}_truth.vcf.gz"
MERGE_OK="false"

if ! command -v bcftools >/dev/null 2>&1; then
    echo ">> [skipped] bcftools not found — per-pass VCFs only."
    echo "   Install bcftools to enable the origin-tagged truth-VCF merge."
elif ! command -v bgzip >/dev/null 2>&1; then
    echo ">> [skipped] bgzip not found — install htslib for the truth-VCF merge."
else
    echo ">> Merging golden VCFs into origin-tagged truth set..."
    ISEC_DIR="${OUTPUT_DIR}/_isec_$$"
    mkdir -p "$ISEC_DIR"
    # gen-reads (v1.11.1+) emits bgzf-framed output directly via
    # noodles::bgzf, so `bcftools index -t` can tabix the per-pass VCFs
    # without a transcode step. Earlier (v1.11.0) builds wrote plain gzip
    # and required `bcftools view -O z` here.
    bcftools index -f -t "$NORMAL_VCF"
    bcftools index -f -t "$TUMOR_VCF"
    # Outputs:
    #   0000.vcf.gz — only in NORMAL_VCF (germline that didn't round-trip)
    #   0001.vcf.gz — only in TUMOR_VCF  (somatic de-novo)
    #   0002.vcf.gz — common, drawn from NORMAL_VCF
    #   0003.vcf.gz — common, drawn from TUMOR_VCF (preferred — preserves tumor's SV INFO)
    bcftools isec -p "$ISEC_DIR" -O z "$NORMAL_VCF" "$TUMOR_VCF" >/dev/null

    # Annotate each disjoint subset with EIDOLON_ORIGIN and re-bgzip.
    # awk owns the INFO-column rewrite because we want a stream-pure
    # operation that doesn't require building an annotations file.
    annotate_origin() {
        local in_vcf="$1" out_vcf="$2" origin="$3"
        zcat "$in_vcf" | awk -v origin="$origin" '
            BEGIN { FS = "\t"; OFS = "\t" }
            /^##/ { print; next }
            /^#CHROM/ {
                print "##INFO=<ID=EIDOLON_ORIGIN,Number=1,Type=String,Description=\"" \
                      "Origin of variant in tumor/normal mix: germline | somatic | shared\">"
                print; next
            }
            {
                tag = "EIDOLON_ORIGIN=" origin
                if ($8 == "." || $8 == "") $8 = tag
                else                       $8 = $8 ";" tag
                print
            }
        ' | bgzip -c > "$out_vcf"
    }

    annotate_origin "$ISEC_DIR/0000.vcf.gz" "$ISEC_DIR/normal_only.vcf.gz" germline
    annotate_origin "$ISEC_DIR/0001.vcf.gz" "$ISEC_DIR/tumor_only.vcf.gz"  somatic
    annotate_origin "$ISEC_DIR/0003.vcf.gz" "$ISEC_DIR/shared.vcf.gz"      shared

    # eidolon's gen-reads writer preserves INFO values from input_vcf
    # verbatim but doesn't add ##INFO header declarations for arbitrary
    # user-supplied fields (e.g. INFO/MATEID on BND records, INFO/SVTYPE
    # / END / SVLEN / CN on symbolic SVs). bcftools sort does internal
    # BCF translation that errors out on undeclared INFO fields. Inject
    # the standard SV-INFO declarations into each subset's header via
    # `bcftools annotate -h` before sorting so the merge succeeds even
    # when the per-pass VCFs carry BND or other SV records.
    SV_INFO_HDR="$ISEC_DIR/sv_info_decls.txt"
    cat > "$SV_INFO_HDR" <<'HDR'
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description="Read depth at local breakend (Manta-style)">
##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at remote breakend (Manta-style)">
HDR
    for f in "$ISEC_DIR/normal_only.vcf.gz" \
             "$ISEC_DIR/tumor_only.vcf.gz" \
             "$ISEC_DIR/shared.vcf.gz"; do
        bcftools annotate -h "$SV_INFO_HDR" -O z -o "$f.hdr" "$f"
        mv "$f.hdr" "$f"
        bcftools index -f -t "$f"
    done

    # bcftools concat -a (allow-overlaps) requires indexed inputs because
    # the three subset files all sit on overlapping contigs (different
    # positions per file, but the same chr names). Index, concat, sort.
    bcftools concat -a -O u \
        "$ISEC_DIR/normal_only.vcf.gz" \
        "$ISEC_DIR/tumor_only.vcf.gz" \
        "$ISEC_DIR/shared.vcf.gz" \
      | bcftools sort -O z -o "$MERGED_TRUTH_VCF"
    bcftools index -f -t "$MERGED_TRUTH_VCF"
    rm -rf "$ISEC_DIR"
    MERGE_OK="true"
fi

# ── Summary ──────────────────────────────────────────────────────────────
cat <<EOF

────────────────────────────────────────────────────────────────
Cancer simulation complete

Configuration:
  Reference:        $REFERENCE
  Total coverage:   ${TOTAL_COVERAGE}×
  Purity:           $PURITY  (tumor cell fraction)
    → normal pass:  ${NORMAL_COVERAGE}× coverage
    → tumor pass:   ${TUMOR_COVERAGE}× coverage
  Paired-ended:     $PAIRED_END
  Tumor SV scale:   $SV_RATE_SCALE

Outputs (in ${OUTPUT_DIR}):
  Normal FASTQ:     ${NORMAL_PREFIX}_r1.fastq.gz $( [[ "$PAIRED_END" == "true" ]] && echo ", ${NORMAL_PREFIX}_r2.fastq.gz" )
  Normal truth:     ${NORMAL_PREFIX}.vcf.gz             (germline only)
  Tumor FASTQ:      ${TUMOR_PREFIX}_r1.fastq.gz $( [[ "$PAIRED_END" == "true" ]] && echo ", ${TUMOR_PREFIX}_r2.fastq.gz" )
  Tumor truth:      ${TUMOR_PREFIX}.vcf.gz              (germline + somatic)
  Merged FASTQ:     ${MERGED_PREFIX}_r1.fastq.gz$( [[ "$PAIRED_END" == "true" ]] && echo ", ${MERGED_PREFIX}_r2.fastq.gz" )
$( [[ "$MERGE_OK" == "true" ]] && echo "  Merged truth:     ${MERGED_PREFIX}_truth.vcf.gz   (INFO/EIDOLON_ORIGIN = germline | somatic | shared)" )

The merged FASTQ is what a downstream somatic-variant pipeline should
consume.$( [[ "$MERGE_OK" == "true" ]] && echo " The merged truth VCF carries INFO/EIDOLON_ORIGIN tags that distinguish germline, somatic, and shared (germline-carried-through-tumor) variants — feed it to hap.py / som.py / vcfeval as the truth set for somatic-caller benchmarking." || echo " Per-pass golden VCFs are preserved; install bcftools + bgzip to enable the origin-tagged truth merge." )
────────────────────────────────────────────────────────────────
EOF
