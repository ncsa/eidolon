#!/bin/bash
# SLURM job: stage REAL cancer data — SEQC2 HCC1395 (tumor) / HCC1395BL (normal),
# a triple-negative breast cancer cell line + matched normal from the same donor —
# for a COMPOUND tumor model: germline layer from the normal + somatic layer from
# the tumor's real high-confidence truth VCF, composed by `eidolon gen-cancer-reads`.
#
# WHY / WHAT THE COMPOUND MODEL IS
#   gen-cancer-reads separates the germline process (normal_model + germline_vcf)
#   from the somatic process (tumor_model). This stages both from real data:
#     * somatic  → SEQC2 high-confidence somatic SNV+INDEL truth VCF (real TNBC
#                  signature; feeds gen-mut-model → tumor_model)
#     * germline → CALLED from the normal BAM here (this donor's germline; feeds
#                  gen-mut-model → normal_model + used as the shared germline_vcf)
#     * sequencing (seq-error/frag/gc/bam-models) → the real tumor/normal BAMs
#
# SEQC2 hosts pre-aligned + deduped NovaSeq BAMs (GRCh38), so NO alignment/index is
# needed. To respect scratch quota, we region-scope (default chr1): samtools pulls
# just that region from the remote BAM via range requests (MODE=remote), or set
# MODE=download to fetch the whole ~60-80 GB BAM first.
#
# Prereqs: GRCh38.fa staged at $SCRATCH/neat_data/GRCh38.fa (same build used for
#          HG002); samtools + bcftools (bioinf conda env / modules). eidolon NOT needed.
#
# Usage:
#   sbatch scripts/delta/stage_hcc1395.sh
#   REGION="chr1 chr2" sbatch scripts/delta/stage_hcc1395.sh
#   REGION=all MODE=download sbatch scripts/delta/stage_hcc1395.sh   # whole genome, full BAMs
#
# HEAVY: region chr1 ≈ ~5 GB pulled per BAM; germline calling on chr1 (~50x) is
# multi-hour. REGION=all downloads ~140 GB and calls genome-wide.

#SBATCH --job-name=eidolon-stagehcc1395
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + setup_conda/conda_activate

D="${DATA_DIR:-$SCRATCH/neat_data/hcc1395}"
REF="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
REGION="${REGION:-chr1}"                 # space-separated contigs, or "all" for whole genome
MODE="${MODE:-remote}"                    # remote = range-request subset; download = full BAM first
THREADS="${SLURM_CPUS_PER_TASK:-16}"
BAM_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS"
VCF_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest"
TUMOR_BAM="${TUMOR_BAM:-WGS_NS_T_1.bwa.dedup.bam}"    # NovaSeq replicate 1
NORMAL_BAM="${NORMAL_BAM:-WGS_NS_N_1.bwa.dedup.bam}"
# Which BAMs to pull. The realism panel needs only the normal, and each sample is a
# separate multi-GB range-request subset, so fetching both to use one doubles the download
# for nothing.
SAMPLES="${SAMPLES:-tumor normal}"
# Stage the BAMs and stop. Steps 3 and 4 exist to build model inputs: germline calling over
# the whole region at ~50x takes hours, and the tumor FASTQ gzips a multi-GB BAM. Neither
# is an input to anything that just needs real reads to measure against.
BAM_ONLY="${BAM_ONLY:-0}"
EIDOLON_HINT="${EIDOLON_BIN:-$SCRATCH/cargo-target/eidolon/release/eidolon}"   # for the printed next-step cmds

mkdir -p "$D"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf                     # bcftools (+ samtools with libcurl if module lacks it)

[[ -f "$REF" ]] || { echo "GRCh38 reference not staged: $REF" >&2; exit 1; }
[[ -f "$REF.fai" ]] || samtools faidx "$REF"
echo "=== stage HCC1395: region=[$REGION] mode=$MODE ref=$REF threads=$THREADS ==="

# Preflight (fail fast + loud): the SEQC2 BAMs and truth VCF are chr-prefixed (chr1).
# If REFERENCE is plain/Ensembl-named ("1"), mpileup emits N ref bases and germline
# calling silently yields ZERO variants — job 19901611 burned 55 min on exactly this.
if [[ "$REGION" != "all" ]]; then
    for ctg in $REGION; do
        cut -f1 "$REF.fai" | grep -qxF "$ctg" || {
            echo "FATAL: contig '$ctg' not found in $REF" >&2
            echo "       reference contigs look like: $(cut -f1 "$REF.fai" | head -3 | tr '\n' ' ')" >&2
            echo "       SEQC2 data is chr-prefixed — point REFERENCE at a chr-prefixed GRCh38" >&2
            echo "       (GDC GRCh38.d1.vd1, or GIAB no_alt analysis set)." >&2
            exit 1
        }
    done
fi

# ── 1. somatic truth VCF (SNV + INDEL, high-conf in high-conf regions; GRCh38) ──
if [[ "$BAM_ONLY" == "0" && ! -s "$D/somatic.vcf.gz" ]]; then
    for f in high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz \
             high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz; do
        [[ -s "$D/$f" ]] || wget -c -O "$D/$f" "$VCF_BASE/$f"
        bcftools index -f -t "$D/$f"
    done
    # combine SNV + INDEL into one somatic model input, sorted + bgzipped.
    bcftools concat -a "$D/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz" \
                       "$D/high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz" -Ou \
        | bcftools sort -Oz -o "$D/somatic.vcf.gz"
    bcftools index -t "$D/somatic.vcf.gz"
fi

# region arg forms: samtools/bcftools want space-separated contigs; "all" = no filter.
region_args=(); vcf_region=""
if [[ "$REGION" != "all" ]]; then
    region_args=($REGION)
    vcf_region=$(printf '%s,' $REGION); vcf_region="${vcf_region%,}"   # bcftools -r csv
fi

# ── 2. subset the tumor + normal BAMs to REGION (no alignment — pre-aligned) ─────
fetch_bam() {   # <remote_name> <local_prefix>
    local name="$1" out="$D/$2.bam"
    [[ -s "$out" ]] && return 0
    if [[ "$MODE" == "download" || "$REGION" == "all" ]]; then
        [[ -s "$D/$name" ]] || wget -c -O "$D/$name" "$BAM_BASE/$name"
        samtools view -@ "$THREADS" -b "$D/$name" "${region_args[@]}" -o "$out"
    else
        # Range-request subset straight from the remote BAM (needs htslib+libcurl).
        # samtools auto-discovers the remote .bai alongside the .bam.
        samtools view -@ "$THREADS" -b "$BAM_BASE/$name" "${region_args[@]}" -o "$out" \
            || { echo "remote subset failed (samtools may lack libcurl) — retry with MODE=download" >&2; exit 1; }
    fi
    samtools index "$out"
    [[ -s "$out" ]] && [[ "$(samtools view -c "$out")" -gt 0 ]] \
        || { echo "empty BAM after subset: $out (bad REGION contig name? check chr-prefix)" >&2; exit 1; }
}
for sample in $SAMPLES; do
    case "$sample" in
        tumor)  fetch_bam "$TUMOR_BAM"  tumor  ;;
        normal) fetch_bam "$NORMAL_BAM" normal ;;
        # Not a warning: a typo here would otherwise stage nothing and exit 0, and an
        # absent BAM is indistinguishable from one that was never asked for.
        *) echo "FATAL: unknown SAMPLES entry '$sample' (want: tumor, normal)" >&2; exit 1 ;;
    esac
done

if [[ "$BAM_ONLY" != "0" ]]; then
    echo
    echo "════════════════════════════════════════════════════════════════"
    echo "BAM_ONLY: staged reads only, over [$REGION]. No truth VCF, no germline"
    echo "calls, no FASTQ, no model configs — rerun without BAM_ONLY for those."
    for sample in $SAMPLES; do
        n="$(samtools view -c "$D/$sample.bam")"
        echo "  $D/$sample.bam — $n reads on [$REGION]"
        # A BAM that exists but holds nothing measures nothing while reading as staged.
        [[ "$n" -gt 0 ]] || { echo "FATAL: $sample.bam is empty" >&2; exit 1; }
    done
    echo "════════════════════════════════════════════════════════════════"
    exit 0
fi

# ── 3. call germline from the NORMAL bam (this donor's germline; per-contig || ) ─
[[ -s "$D/normal.bam" ]] || { echo "FATAL: step 3 needs normal.bam — add 'normal' to SAMPLES, or set BAM_ONLY=1" >&2; exit 1; }
if [[ ! -s "$D/germline.vcf.gz" ]]; then
    echo "calling germline from normal BAM over [$REGION]..."
    if [[ "$REGION" == "all" ]]; then
        contigs=$(cut -f1 "$REF.fai")
    else
        contigs="$REGION"
    fi
    parts="$D/germ_parts"; mkdir -p "$parts"
    printf '%s\n' $contigs | xargs -P "$THREADS" -I CTG bash -c '
        ctg="$1"; out="'"$parts"'/$ctg.vcf.gz"
        [[ -s "$out" ]] && exit 0
        # NOTE: mpileup stderr is NOT suppressed — a "sequence not found" ref/BAM
        # naming error must surface (the preflight above should already catch it).
        bcftools mpileup -r "$ctg" -f "'"$REF"'" "'"$D"'/normal.bam" \
            | bcftools call -mv -Oz -o "$out"
    ' _ CTG
    printf '%s\n' $contigs | sed "s#^#$parts/#; s#\$#.vcf.gz#" > "$D/germ_parts.list"
    bcftools concat -Oz -f "$D/germ_parts.list" -o "$D/germline.vcf.gz"
    bcftools index -t "$D/germline.vcf.gz"
    rm -rf "$parts" "$D/germ_parts.list"
fi

# ── 4. FASTQ for seq-error (tumor reads) + contig-naming sanity ─────────────────
[[ -s "$D/tumor.bam" ]] || { echo "FATAL: step 4 needs tumor.bam — add 'tumor' to SAMPLES, or set BAM_ONLY=1" >&2; exit 1; }
[[ -s "$D/tumor.fastq.gz" ]] || samtools fastq -n "$D/tumor.bam" 2>/dev/null | gzip > "$D/tumor.fastq.gz"

# Summary + sanity below use `head` and `grep -q`, which close their pipe early →
# the upstream bcftools/cut take SIGPIPE (exit 141), and set -o pipefail + set -e
# would turn that into a fatal false-failure AFTER all real staging finished (job
# 19901611 died here). Disable pipefail for this block. (Same footgun as stage_soy.)
set +o pipefail
som_ctg=$(bcftools view -H "$D/somatic.vcf.gz" 2>/dev/null | head -1 | cut -f1)
cut -f1 "$REF.fai" | grep -qxF "$som_ctg" \
    || echo "WARNING: somatic VCF contig '$som_ctg' not in $REF.fai — check GRCh38 contig naming" >&2

n_som=$(bcftools view -H "$D/somatic.vcf.gz" ${vcf_region:+-r "$vcf_region"} 2>/dev/null | wc -l)
n_germ=$(bcftools view -H "$D/germline.vcf.gz" 2>/dev/null | wc -l)
set -o pipefail

# Pre-write the two gen-mut-model configs to SHARED fs ($D) — NOT /tmp, which is
# node-local on Delta and invisible to an srun/sbatch task on a compute node.
for pair in "tumor_model:$D/somatic.vcf.gz" "normal_model:$D/germline.vcf.gz"; do
    name="${pair%%:*}"; vcf="${pair##*:}"
    printf 'reference: %s\nvcf_file: %s\noutput_file: %s\noverwrite_output: true\nbed_file: .\n' \
        "$REF" "$vcf" "$D/$name.json.gz" > "$D/$name.yml"
done

echo
echo "════════════════════════════════════════════════════════════════"
echo "HCC1395 staged (region [$REGION]): $n_som somatic, $n_germ germline variants"
echo "Build the COMPOUND model, then simulate (configs written to $D):"
echo "  # somatic tumor_model (real TNBC signature) + germline normal_model (this donor):"
echo "  srun --account=bhrd-delta-cpu -p cpu --mem=8G -t 00:30:00 $EIDOLON_HINT gen-mut-model -c $D/tumor_model.yml"
echo "  srun --account=bhrd-delta-cpu -p cpu --mem=8G -t 00:30:00 $EIDOLON_HINT gen-mut-model -c $D/normal_model.yml"
echo "  # sequencing models from the real tumor BAM/FASTQ:"
echo "  REFERENCE=$REF INPUT_BAM=$D/tumor.bam INPUT_FASTQ=$D/tumor.fastq.gz INPUT_VCF=$D/somatic.vcf.gz \\"
echo "    sbatch scripts/delta/model_builders.sbatch"
echo "  # compound tumor/normal simulation (gen-cancer-reads):"
echo "  #   reference: $REF ; normal_model: normal_model.json.gz ; tumor_model: tumor_model.json.gz"
echo "  #   germline_vcf: $D/germline.vcf.gz ; tumor_mutation_rate: model ; paired_ended + fragment_model"
echo "════════════════════════════════════════════════════════════════"
