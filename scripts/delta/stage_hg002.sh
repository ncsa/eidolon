#!/bin/bash
# SLURM job: stage REAL GIAB HG002 (human, GRCh38) data for the model-builder harness.
#
# Unlike soy, HG002 has a CURATED truth VCF (GIAB v4.2.1) that already matches GRCh38,
# so we don't call variants ourselves — we align a read subset for the BAM/FASTQ inputs
# and point gen-mut-model at the truth VCF. This is the human-scale counterpart to
# stage_soy.sh, and the run that tests the ~12-14 GB peak-RSS extrapolation recorded
# in docs/model_builder_baseline.md (soy ~1 Gb → ~4 GB; human ~3.1 Gb → ~12-14 GB?).
#
# Pipeline: download an Illumina 2x250 read chunk (GIAB FTP, no sra-tools) -> bwa-mem2
# align to GRCh38.fa -> coordinate-sorted BAM. FASTQ input = the downloaded reads;
# VCF input = the GIAB truth VCF (fetch_validation_data.sh DATA=hg002).
# Outputs (in $DATA_DIR): hg002.bam(+.bai), hg002_R1.fastq.gz, hg002_R2.fastq.gz.
#
# Prereqs:
#   * GRCh38.fa staged at $SCRATCH/neat_data/GRCh38.fa (same build already used for sim).
#   * truth VCF staged:  DATA=hg002 bash scripts/delta/fetch_validation_data.sh
#   * bwa-mem2 + samtools (bioinf conda env / modules). eidolon NOT needed here.
#
# HEAVY: ~12 GB read download per chunk; the bwa-mem2 index of GRCh38 needs ~64 GB RAM
# and ~30 GB disk (built once, reused); human alignment of one chunk is ~1-2 h.
#
# Usage:
#   sbatch scripts/delta/stage_hg002.sh
#   CHUNKS="L001:001 L001:002" MAX_PAIRS=20000000 sbatch scripts/delta/stage_hg002.sh
#
#   # reads only, for a quality-model fit -- no reference, no index, no alignment
#   ALIGN=0 sbatch scripts/delta/stage_hg002.sh
#
#   # a different GIAB library: give the directory and the two filename shapes
#   ALIGN=0 BASE=https://.../SOME_DATASET/reads \
#     R1_PATTERN='PREFIX_{lane}_R1_{id}.fastq.gz' \
#     R2_PATTERN='PREFIX_{lane}_R2_{id}.fastq.gz' \
#     DATA_DIR=$SCRATCH/neat_data/hg002_other sbatch scripts/delta/stage_hg002.sh
#
# The default library is NIST's 2x250, which is HiSeq 2500-era continuous quality scoring.
# There is no NovaSeq dataset under the GIAB AshkenazimTrio HG002 path -- its Illumina sets are
# NIST_Illumina_2x250bps and NIST_HiSeq_HG002_Homogeneity, both HiSeq. Staging a binned-quality
# library means pointing this at a dataset hosted elsewhere (#730).

#SBATCH --job-name=eidolon-stagehg002
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + setup_conda/conda_activate

D="${DATA_DIR:-$SCRATCH/neat_data/hg002}"
REF="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
VCF="${HG002_VCF:-$D/HG002_GRCh38_v4.2.1_benchmark.vcf.gz}"
# Space-separated LANE:CHUNK ids to download+concatenate. One chunk (~20M pairs) is
# plenty for the builders; add more for a bigger/heavier BAM.
CHUNKS="${CHUNKS:-L001:001}"
MAX_PAIRS="${MAX_PAIRS:-0}"                 # 0 = use all downloaded reads; else subsample
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# WHICH GIAB LIBRARY. Defaults to NIST's 2x250, which is what this script has always fetched.
# That set is HiSeq 2500-era CONTINUOUS quality scoring -- measured on the fitted model: 31
# options spanning Q2-Q40, near-contiguous, not binned (#677). Point BASE and the two patterns
# at another GIAB dataset to stage a different chemistry; the naming differs per dataset, so
# both the directory and the filename shape are parameters rather than one guessed template.
# `{lane}` and `{id}` are substituted from each CHUNKS entry.
BASE="${BASE:-https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads}"
# NOT `${R1_PATTERN:-...}`. A `${VAR:-default}` ends at the FIRST unescaped `}`, so a default
# CONTAINING braces is cut off there. `${R1_PATTERN:-D1_S1_{lane}_R1_{id}.fastq.gz}` yields
# `D1_S1_{lane` plus the literal tail `_R1_{id}.fastq.gz}` -- `{lane}` then never exists to be
# substituted, `{id}` survives and IS substituted, and a stray `}` lands on the end. Measured:
# the URL fetched was `D1_S1_{lane_R1_001.fastq.gz}`, a 404 that reads like a moved dataset.
# Assigning conditionally keeps every brace out of an expansion.
[[ -n "${R1_PATTERN:-}" ]] || R1_PATTERN='D1_S1_{lane}_R1_{id}.fastq.gz'
[[ -n "${R2_PATTERN:-}" ]] || R2_PATTERN='D1_S1_{lane}_R2_{id}.fastq.gz'

# Substitution as a function so a test can reach it. The defaults above are the part that
# broke, and a test that cannot see them is not a test of them -- PRINT_URLS below exposes
# both together.
expand_pattern() {  # <pattern> <lane> <id>
    local p="$1"
    p="${p//\{lane\}/$2}"
    p="${p//\{id\}/$3}"
    # Refuse a name still holding a placeholder rather than asking the server about it. A 404
    # from a mangled name is indistinguishable from a dataset that moved, which is how the
    # brace bug above cost a submit round trip to diagnose.
    case "$p" in
        *[{}]*) echo "FATAL: unsubstituted placeholder in '$p' (from '$1'). Patterns take {lane} and {id}; CHUNKS entries are LANE:ID." >&2; return 1 ;;
    esac
    printf '%s' "$p"
}

# ALIGN=0 stages the READS ONLY -- no reference index, no alignment, no BAM.
#
# `gen-seq-error-model` reads FASTQ and nothing else, so fitting a quality model needs none of
# section 2 or 3 below. Those exist for the model builders, which want a BAM. Skipping them
# turns a job that spends 1-2 h aligning (plus a one-time ~64 GB-RAM bwa-mem2 index) into a
# download, which is the difference between a plausible experiment and an expensive one.
ALIGN="${ALIGN:-1}"

# PRINT_URLS=1 prints what would be downloaded and exits, touching no filesystem and no
# network. Worth having on its own -- the failure this guards against was a 404 whose cause
# was invisible in the URL until it was percent-decoded.
if [[ -n "${PRINT_URLS:-}" ]]; then
    for spec in $CHUNKS; do
        lane="${spec%%:*}"; id="${spec##*:}"
        # Assigned first: a command substitution that fails inside `echo`'s arguments does not
        # abort under `set -e`, so the guard in expand_pattern would be advisory there.
        u1="$(expand_pattern "$R1_PATTERN" "$lane" "$id")"
        u2="$(expand_pattern "$R2_PATTERN" "$lane" "$id")"
        echo "$BASE/$u1"
        echo "$BASE/$u2"
    done
    exit 0
fi

mkdir -p "$D"
# Only the align path needs these. A reads-only run is wget and gzip, so loading them would
# make the job fail on an environment it never uses.
if [[ "$ALIGN" -eq 1 ]]; then
    module load samtools/1.22-cce19.0.0
    setup_conda
    conda_activate bioinf                   # bwa-mem2
fi

[[ "$ALIGN" -ne 1 || -f "$REF" ]] || { echo "GRCh38 reference not staged: $REF" >&2; exit 1; }
echo "=== stage HG002: chunks=[$CHUNKS]  max_pairs=$MAX_PAIRS  threads=$THREADS ==="
echo "    source: $BASE"
echo "    align:  $ALIGN  (0 = reads only, no BAM)"
[[ "$ALIGN" -eq 1 ]] && echo "    ref:    $REF"

# ── 1. download + concatenate read chunks (GIAB FTP, resumable) ─────────────────
R1="$D/hg002_R1.fastq.gz"; R2="$D/hg002_R2.fastq.gz"
if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    : > "$R1"; : > "$R2"
    for spec in $CHUNKS; do
        lane="${spec%%:*}"; id="${spec##*:}"
        f1="$(expand_pattern "$R1_PATTERN" "$lane" "$id")"
        f2="$(expand_pattern "$R2_PATTERN" "$lane" "$id")"
        echo "downloading $f1 / $f2 ..."
        wget -c -O "$D/$f1" "$BASE/$f1"
        wget -c -O "$D/$f2" "$BASE/$f2"
        cat "$D/$f1" >> "$R1"; cat "$D/$f2" >> "$R2"
    done
fi

# Optional subsample (head keeps R1/R2 in the same order → pairs stay matched).
if [[ "$MAX_PAIRS" -gt 0 && ! -s "$D/sub_1.fastq.gz" ]]; then
    echo "subsampling first $MAX_PAIRS pairs..."
    # head closes the pipe early → zcat gets SIGPIPE (141); expected, not an error.
    # pipefail would promote it to a failure and set -e would abort → disable it here.
    set +o pipefail
    zcat "$R1" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_1.fastq.gz"
    zcat "$R2" | head -n $(( MAX_PAIRS * 4 )) | gzip > "$D/sub_2.fastq.gz"
    set -o pipefail
fi
[[ "$MAX_PAIRS" -gt 0 ]] && { R1="$D/sub_1.fastq.gz"; R2="$D/sub_2.fastq.gz"; }

# Sections 2-4 build the BAM the model builders want. A quality-model fit needs none of it.
VCF_OK=0
if [[ "$ALIGN" -eq 1 ]]; then
    # ── 2. index reference (bwa-mem2 + faidx); -s guard against a 0-byte stub ────────
    [[ -f "$REF.fai" ]] || samtools faidx "$REF"
    if [[ ! -s "$REF.bwt.2bit.64" ]]; then
        echo "building bwa-mem2 index for GRCh38 (one-time, ~64 GB RAM, ~30 GB disk)..."
        bwa-mem2 index "$REF"
    fi

    # ── 3. align → coordinate-sorted BAM ─────────────────────────────────────────────
    if [[ ! -s "$D/hg002.bam" ]]; then
        echo "aligning HG002 reads to GRCh38..."
        bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" 2> "$D/bwa.log" \
            | samtools sort -@ "$THREADS" -o "$D/hg002.bam" -
        samtools index "$D/hg002.bam"
    fi

    # ── 4. sanity: does the truth VCF's contig naming match the reference? ──────────
    # gen-mut-model silently produces nothing if the VCF's contigs aren't in the FASTA
    # (the classic 'chr1' vs '1' trap). Warn loudly rather than let mut_model no-op.
    VCF_OK=1
    if [[ -s "$VCF" ]]; then
        vcf_ctg=$(zcat "$VCF" | grep -v '^#' | head -1 | cut -f1)
        if ! cut -f1 "$REF.fai" | grep -qxF "$vcf_ctg"; then
            echo "WARNING: truth VCF contig '$vcf_ctg' not found in $REF.fai — gen-mut-model" >&2
            echo "         would see zero variants. Check GRCh38 contig naming (chr-prefix?)." >&2
            VCF_OK=0
        fi
    else
        echo "WARNING: truth VCF not found at $VCF — run: DATA=hg002 bash scripts/delta/fetch_validation_data.sh" >&2
        VCF_OK=0
    fi
fi

nreads=$(( $(zcat "$R1" | wc -l) / 4 ))
echo
echo "════════════════════════════════════════════════════════════════"
if [[ "$ALIGN" -eq 1 ]]; then
    echo "HG002 staged: $nreads read pairs aligned to GRCh38"
    echo "Run the model-builders (human-scale stress; expect gc_bias/bam_models RSS in the GB range):"
    echo "  REFERENCE=$REF INPUT_BAM=$D/hg002.bam \\"
    echo "  INPUT_FASTQ=$R1 \\"
    [[ "$VCF_OK" -eq 1 ]] && echo "  INPUT_VCF=$VCF \\" || echo "  # (INPUT_VCF omitted — see contig-naming warning above; mut_model will be skipped)"
    echo "    sbatch scripts/delta/model_builders.sbatch"
else
    echo "HG002 reads staged: $nreads pairs (no BAM — ALIGN=0)"
    echo "  R1: $R1"
    echo "  R2: $R2"
    echo
    echo "Fit ONE model carrying both mates (#723). REBUILD THE BINARY FIRST: fastq_file_r2 is"
    echo "read by key lookup, so a binary predating it ignores the key and fits R1 twice."
    echo "  printf 'fastq_file: %s\\nfastq_file_r2: %s\\noutput_file: %s\\noverwrite_output: true\\nmax_reads: 0\\nqual_offset: 33\\n' \\"
    echo "    $R1 $R2 $D/hg002_qual.json.gz > $D/fit.yml"
    echo "  \$EIDOLON gen-seq-error-model -c $D/fit.yml"
    echo
    echo "Then read its shape against the shipped default:"
    echo "  bash scripts/delta/describe_quality_model.sh $D/hg002_qual.json.gz"
fi
echo "════════════════════════════════════════════════════════════════"
