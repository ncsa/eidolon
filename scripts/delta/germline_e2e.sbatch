#!/bin/bash
# SLURM job: germline recall comparison — eidolon vs NEAT 4 — on Delta.
#
# For EACH simulator (eidolon and NEAT 4), simulates SNP/indel-only germline reads
# and runs the SAME downstream pipeline (BWA-MEM2 → GATK HaplotypeCaller → hap.py),
# scoring the caller's recall/precision/F1 against THAT tool's own truth VCF. The
# two summaries are reported side-by-side: "how well does a standard germline
# caller recover each simulator's variants" — the fidelity half of goal 2
# (performance is benchmark.sbatch).
#
# Both tools use their DEFAULT mutation model (realistic Ts/Tv ~2.3 for each) and
# identical coverage/read-length/fragment params, so the comparison is fair.
#
# Stages per tool (all idempotent — existing outputs skipped on re-run):
#   1. simulate (SNP/indel only) -> reads + truth VCF
#   2. BWA-MEM2 align -> SAMtools sort/index
#   3. GATK HaplotypeCaller (gVCF -> genotyped VCF)
#   4. hap.py: called vs that tool's truth VCF
#
# Prerequisites: run scripts/delta/setup.sh once (builds eidolon, neat4 + bioinf
# conda envs, hap.py SIF). Submit from the repo root.
#
# Usage:
#   sbatch scripts/delta/germline_e2e.sbatch
#   REFERENCE=/path/ref.fa COVERAGE=30 TOOLS="eidolon neat" sbatch scripts/delta/germline_e2e.sbatch
#   TOOLS=eidolon sbatch ...        # eidolon only (skip the NEAT 4 arm)
#   SHUFFLE=1 TOOLS=eidolon sbatch ...   # randomize read order before alignment
#       (off by default; run once shuffled vs unshuffled to confirm the caller
#        is read-order-independent — same recall expected)
#   ADAPTERS=1 MEASURE_REALISM=1 TOOLS=eidolon sbatch ...   # 3' adapter readthrough
#       on, aligned raw, with adapter/insert/soft-clip QC captured (#125).
#   ADAPTERS=1 TRIM=1 MEASURE_REALISM=1 TOOLS=eidolon sbatch ...  # + fastp-trim
#       before alignment (the realistic Illumina pipeline).
#       Drive the off/on-untrimmed/on-trimmed A/B/C matrix with
#       scripts/delta/run_adapter_validation.sh (REPS for mean±sd).

#SBATCH --job-name=eidolon-germline
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

# ── Parameters ───────────────────────────────────────────────────────────
# Under sbatch, $0 is a spooled copy — derive the repo from the submit dir
# (run `sbatch` from the repo root) or override with EIDOLON_REPO.
REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # archive_run() + Delta path resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chr22.fa}"
OUTDIR="${OUTDIR:-$SCRATCH/germline_$SLURM_JOB_ID}"
TOOLS="${TOOLS:-eidolon neat}"          # which simulators to run + score
# Seeds. Under `sbatch --array=1-3` (triplicate), each replicate gets a distinct
# seed (independent draw) and SLURM spreads tasks across nodes; plain submit uses
# the base seeds.
SEED="${SEED:-germline fidelity delta}${SLURM_ARRAY_TASK_ID:+ rep$SLURM_ARRAY_TASK_ID}"
NEAT_SEED="${NEAT_SEED:-42}"
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] && NEAT_SEED=$((NEAT_SEED + SLURM_ARRAY_TASK_ID))

COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
# Read-order randomization before alignment (off by default). eidolon/NEAT emit
# reads in genome order; SHUFFLE=1 randomizes (seqkit, paired) to mimic real
# sequencer output and to demonstrate the caller is order-independent. Needs
# seqkit in the bioinf env (setup.sh installs it).
SHUFFLE="${SHUFFLE:-0}"
SHUFFLE_SEED="${SHUFFLE_SEED:-11}"

# ── 3' adapter readthrough (#125) + adapter-aware QC ──────────────────────
# All three default OFF → sim.yml carries no `adapters:` block and no extra
# stages run, so output is byte-identical to pre-adapter runs (the variety and
# regression jobs that call this script are unaffected).
#   ADAPTERS=1        enable eidolon 3' adapter readthrough (adapters.enabled=true)
#   ADAPTER_PRESET    truseq|nextera|custom (only used when ADAPTERS=1)
#   TRIM=1            fastp adapter-trim BEFORE alignment (the realistic Illumina
#                     pipeline). Off = align raw reads, exposing soft-clipped
#                     adapter tails so we can measure the untrimmed penalty.
#   MEASURE_REALISM=1 capture adapter/insert/soft-clip/mapping QC → realism.tsv
#                     (fastp report + samtools stats). Independent of TRIM.
ADAPTERS="${ADAPTERS:-0}"
ADAPTER_PRESET="${ADAPTER_PRESET:-truseq}"
TRIM="${TRIM:-0}"
MEASURE_REALISM="${MEASURE_REALISM:-0}"
# KEEP_SHORT=1 keeps short inserts as genomic reads WITHOUT adapter readthrough
# (eidolon keep_short_fragments) — the adapter-free short-insert control arm that
# isolates the coverage effect from the adapter effect (#125). Independent of
# ADAPTERS; leave both off for the default long-insert baseline.
KEEP_SHORT="${KEEP_SHORT:-0}"

# PRUNE=1 deletes the bulky FASTQ + BAM + gVCF intermediates after a SUCCESSFUL
# score (per-tool scored.summary.csv present), keeping the small scores/VCFs/QC.
# Matches cancer_pipeline.sbatch; default off so variety/regression runs are
# unchanged. The adapter matrix (3 conditions × reps, plus trimmed read copies)
# is data-heavy on the shared scratch PROJECT quota — run_adapter_validation.sh
# defaults it on.
PRUNE="${PRUNE:-0}"

HAPPY_SIF="${HAPPY_SIF:-$SCRATCH/sif/happy_v0.3.12.sif}"
NEAT_ENV="${NEAT_ENV:-neat4}"         # conda env with NEAT 4
THREADS=$SLURM_CPUS_PER_TASK
# fastp's reader/worker/writer pipeline is built for a small thread pool; at very
# high counts (e.g. THREADS=64 on a Delta cpu node) it livelocks — barely any CPU,
# zero output — so cap it. Observed 2026-07-01: fastp sat at 38% CPU with 0-byte
# trimmed output for 30+ min until capped. Override with FASTP_THREADS=.
FASTP_THREADS="${FASTP_THREADS:-$(( THREADS < 16 ? THREADS : 16 ))}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="$CARGO_TARGET_DIR/release/eidolon"

# ── Environment ─────────────────────────────────────────────────────────
source "$HOME/.cargo/env"
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1
setup_conda
conda_activate bioinf   # bwa-mem2, gatk4, bcftools (NEAT runs via `conda run -n $NEAT_ENV`)

mkdir -p "$OUTDIR"
REF_DIR="$(dirname "$REFERENCE")"
REF_NAME="$(basename "$REFERENCE")"

echo "=== germline fidelity comparison: $SLURM_JOB_ID ==="
echo "Reference: $REFERENCE   Coverage: ${COVERAGE}x   Read: ${READ_LEN}bp   Threads: $THREADS"
echo "Tools:     $TOOLS"
echo "Output:    $OUTDIR"

# ── Node-local staging (LOCAL_STAGE=1) ────────────────────────────────────
# Run all heavy per-tool I/O (FASTQ/BAM/gVCF) — plus the reference index reads,
# eidolon binary, and hap.py SIF — on the compute node's LOCAL disk, copying only
# the small durable artifacts (scores, VCFs, QC) back to $OUTDIR on Lustre at
# stage-out. A shared-filesystem OST dropping mid-run can then no longer wedge
# the job in uninterruptible I/O wait (the failure that stalled the 2026-07-01
# adapter matrix ~15h in the fastp step): Lustre is touched only briefly at
# stage-in/out. Default OFF (work directly in $OUTDIR, cross-submit resume
# preserved); run_adapter_validation.sh turns it on. Override the local root
# with LOCAL_SCRATCH=..., keep the local tree for debugging with KEEP_LOCAL=1.
LOCAL_STAGE="${LOCAL_STAGE:-0}"
if [[ "$LOCAL_STAGE" == "1" ]]; then
    LOCAL_ROOT="${LOCAL_SCRATCH:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}/eidolon_e2e_${SLURM_JOB_ID:-$$}"
    mkdir -p "$LOCAL_ROOT" || { echo "ERROR: cannot create local scratch $LOCAL_ROOT" >&2; exit 1; }
    [[ "${KEEP_LOCAL:-0}" == "1" ]] || trap 'rm -rf "$LOCAL_ROOT"' EXIT
    # Tool workdirs go under work/ so they can't collide with the staged assets
    # ($LOCAL_ROOT/ref, .../eidolon binary, .../happy.sif) — the "eidolon" tool dir
    # would otherwise clash with the staged "eidolon" binary of the same name.
    WORK_ROOT="$LOCAL_ROOT/work"
    echo "Local stage: $LOCAL_ROOT  (durable outputs copied back to $OUTDIR)"
else
    WORK_ROOT="$OUTDIR"
fi

# ── Shared one-time reference indexing ───────────────────────────────────
# Guard on -s (exists AND non-empty), not -f: an index build killed mid-write
# leaves a 0-byte .bwt.2bit.64 that -f treats as "present", so every later run
# skips indexing and bwa-mem2 dies with "Unexpected end of file". -s rebuilds it.
[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
if [[ ! -s "${REFERENCE}.bwt.2bit.64" && ! -s "${REFERENCE}.bwt" ]]; then
    echo "Indexing reference for BWA-MEM2 (one-time)..."
    bwa-mem2 index "$REFERENCE"
fi
REF_DICT="${REFERENCE%.fa}.dict"
[[ -f "$REF_DICT" ]] || gatk CreateSequenceDictionary -R "$REFERENCE" 2>&1 | tail -2

# ── Node-local stage-in of the reference + tools (LOCAL_STAGE=1) ───────────
# After the shared one-time Lustre index build above, copy the reference and all
# its sidecar indexes (fai/dict + bwa-mem2 .0123/.amb/.ann/.bwt.2bit.64/.pac) —
# plus the eidolon binary and hap.py SIF — onto local disk and repoint at them, so
# alignment/calling/scoring stream them from the node instead of Lustre.
if [[ "$LOCAL_STAGE" == "1" ]]; then
    echo "Staging reference + tools to $LOCAL_ROOT ..."
    mkdir -p "$LOCAL_ROOT/ref"
    cp -f "$REFERENCE" "$LOCAL_ROOT/ref/" 2>/dev/null || true
    for ext in .fai .0123 .amb .ann .bwt.2bit.64 .pac; do
        [[ -e "${REFERENCE}${ext}" ]] && cp -f "${REFERENCE}${ext}" "$LOCAL_ROOT/ref/" || true
    done
    [[ -e "$REF_DICT" ]] && cp -f "$REF_DICT" "$LOCAL_ROOT/ref/" || true
    REFERENCE="$LOCAL_ROOT/ref/$REF_NAME"; REF_DIR="$LOCAL_ROOT/ref"; REF_DICT="${REFERENCE%.fa}.dict"
    [[ -s "$REFERENCE" ]] || { echo "ERROR: reference stage-in failed ($REFERENCE)" >&2; exit 1; }
    [[ -x "$EIDOLON_BIN" ]] && cp -f "$EIDOLON_BIN" "$LOCAL_ROOT/eidolon"     && EIDOLON_BIN="$LOCAL_ROOT/eidolon"     || true
    [[ -f "$HAPPY_SIF" ]] && cp -f "$HAPPY_SIF" "$LOCAL_ROOT/happy.sif" && HAPPY_SIF="$LOCAL_ROOT/happy.sif" || true
    df -h "$LOCAL_ROOT" | tail -1
fi

# ── Simulators (each writes <wd>/sample_r1/r2.fastq.gz + a truth VCF) ─────
# eidolon: SNP/indel only via sv_rate_scale=0; truth = sample.vcf.gz
simulate_eidolon() {
    local wd="$1"
    cat > "$wd/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: 2
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $wd
output_prefix: sample
rng_seed: $SEED
sv_rate_scale: 0
num_threads: $THREADS
YAML
    # 3' adapter readthrough (#125). Appended only when enabled so the default
    # config is byte-identical to the pre-adapter sim.yml.
    if [[ "$ADAPTERS" == "1" ]]; then
        cat >> "$wd/sim.yml" <<YAML
adapters:
  enabled: true
  preset: $ADAPTER_PRESET
YAML
    elif [[ "$KEEP_SHORT" == "1" ]]; then
        # Adapter-free short-insert control: keep short fragments as genomic reads.
        # (Only when ADAPTERS!=1 — adapters.enabled already implies keeping them.)
        echo "keep_short_fragments: true" >> "$wd/sim.yml"
    fi
    "$EIDOLON_BIN" --log-level warn gen-reads -c "$wd/sim.yml"
}

# NEAT 4: default variant model is already SNP/indel only (no symbolic SVs);
# truth = sample_golden.vcf.gz. Runs in its own conda env.
simulate_neat() {
    local wd="$1"
    cat > "$wd/sim.yml" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: 2
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
rng_seed: $NEAT_SEED
overwrite_output: true
YAML
    conda run -n "$NEAT_ENV" neat read-simulator -c "$wd/sim.yml" -o "$wd" -p sample
}

# ── Shared align → call → score (scores against that tool's own truth) ────
# align_call_score <tool> <wd> <r1> <r2> <truth_vcf>
align_call_score() {
    local tool="$1" wd="$2" r1="$3" r2="$4" truth="$5"
    local bam="$wd/aligned.bam" gvcf="$wd/raw.g.vcf.gz" called="$wd/called.vcf.gz"

    if [[ -f "$bam" && -f "$bam.bai" ]]; then
        echo "  [$tool 2/3] BAM present — skipping alignment."
    else
        echo "  [$tool 2/3] BWA-MEM2 align..."
        local rg="@RG\tID:$tool\tSM:$tool\tLB:$tool\tPL:ILLUMINA"
        bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$r1" "$r2" 2>"$wd/bwa.log" \
            | samtools sort -@ "$THREADS" -o "$bam"
        samtools index "$bam"
    fi

    if [[ -f "$called" ]]; then
        echo "  [$tool 3/3] Called VCF present — skipping HaplotypeCaller."
    else
        echo "  [$tool 3/3] GATK HaplotypeCaller (gVCF -> genotyped)..."
        [[ -f "$gvcf" ]] || gatk HaplotypeCaller -R "$REFERENCE" -I "$bam" -O "$gvcf" \
            -ERC GVCF --native-pair-hmm-threads "$THREADS" 2>&1 | tail -5
        gatk GenotypeGVCFs -R "$REFERENCE" -V "$gvcf" -O "$called" 2>&1 | tail -3
    fi

    if [[ -f "$wd/scored.summary.csv" ]]; then
        echo "  [$tool score] hap.py output present — skipping."
    else
        echo "  [$tool score] hap.py (truth=$(basename "$truth"))..."
        [[ -f "$HAPPY_SIF" ]] || { echo "hap.py SIF not found: $HAPPY_SIF — run setup.sh" >&2; exit 1; }
        [[ -f "$truth.tbi" ]] || bcftools index -f -t "$truth"
        # Truth + called both live in $wd, bind it to /work.
        apptainer exec --bind "$REF_DIR:/refs:ro" --bind "$wd:/work" "$HAPPY_SIF" \
            /opt/hap.py/bin/hap.py \
                "/work/$(basename "$truth")" \
                "/work/$(basename "$called")" \
                -r "/refs/$REF_NAME" \
                -o "/work/scored" \
                --engine xcmp --no-roc
    fi
}

# ── Adapter-aware QC helpers (only used when MEASURE_REALISM/TRIM set) ─────
# fastp isn't a Delta module and may post-date the bioinf env; ensure it the
# same just-in-time way setup.sh ensures seqkit, so this runs even if setup.sh
# wasn't re-run after fastp was added.
ensure_fastp() {
    conda run -n bioinf which fastp >/dev/null 2>&1 && return 0
    echo "  [setup] installing fastp into bioinf env (one-time)..."
    conda install -y -n bioinf -c bioconda -c conda-forge fastp >/dev/null 2>&1
    conda run -n bioinf which fastp >/dev/null 2>&1
}

# run_fastp <wd> <r1> <r2>  — adapter detect/trim + QC report (idempotent).
# Always writes fastp.json/.log and trimmed_r{1,2}.fastq.gz; the caller decides
# whether to feed the trimmed reads to alignment (TRIM=1) or just keep the report.
run_fastp() {
    local wd="$1" r1="$2" r2="$3"
    local tr1="$wd/trimmed_r1.fastq.gz" tr2="$wd/trimmed_r2.fastq.gz"
    if [[ -f "$wd/fastp.json" && -s "$tr1" && -s "$tr2" ]]; then
        echo "  [fastp] report present — skipping."; return 0
    fi
    echo "  [fastp] adapter detect/trim + QC report..."
    conda run -n bioinf fastp -i "$r1" -I "$r2" -o "$tr1" -O "$tr2" \
        --detect_adapter_for_pe --thread "$FASTP_THREADS" \
        -j "$wd/fastp.json" -h "$wd/fastp.html" 2>"$wd/fastp.log"
}

# capture_realism <tool> <wd> — adapter/insert/soft-clip/mapping signals from
# samtools stats (BAM) + fastp.json into realism.tsv (key<TAB>value). These are
# the metrics adapter readthrough is expected to move; the collector diffs them
# off-vs-on. Soft-clip fraction = 1 - bases_mapped(cigar)/total_length, the clean
# signal for untrimmed adapter tails getting clipped at the read 3' end.
capture_realism() {
    local tool="$1" wd="$2"
    local bam="$wd/aligned.bam" st="$wd/samstats.txt" out="$wd/realism.tsv"
    [[ -f "$bam" ]] || { echo "  [$tool realism] no BAM — skipping." >&2; return 0; }
    [[ -s "$st" ]] || samtools stats -@ "$THREADS" "$bam" > "$st"

    # samtools-stats SN scalar by label (value is field 3, may carry a trailing comment).
    sn() { awk -F'\t' -v k="$1" '$1=="SN" && $2==k {v=$3; sub(/[ \t].*/,"",v); print v+0; exit}' "$st"; }
    local total mapped paired_pct errrate isize isd avlen avqual tlen bmapcig softclip
    total=$(sn "raw total sequences:")
    mapped=$(sn "reads mapped:")
    errrate=$(sn "error rate:")
    isize=$(sn "insert size average:")
    isd=$(sn "insert size standard deviation:")
    avlen=$(sn "average length:")
    avqual=$(sn "average quality:")
    tlen=$(sn "total length:")
    bmapcig=$(sn "bases mapped (cigar):")
    paired_pct=$(awk -v m="$(sn 'reads properly paired:')" -v t="$total" \
        'BEGIN{printf "%.4f", (t>0)? m/t : 0}')
    mapped=$(awk -v m="$mapped" -v t="$total" 'BEGIN{printf "%.4f", (t>0)? m/t : 0}')
    softclip=$(awk -v b="$bmapcig" -v t="$tlen" 'BEGIN{printf "%.4f", (t>0)? 1-(b/t) : 0}')

    # fastp.json scalars (robust JSON parse; python3 is always present on Delta).
    local fj="$wd/fastp.json"
    local adp_reads=0 adp_bases=0 dup=0 q30=0 gc=0 ipeak=0
    if [[ -f "$fj" ]]; then
        read -r adp_reads adp_bases dup q30 gc ipeak < <(python3 - "$fj" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
ac = d.get("adapter_cutting", {})
af = d.get("summary", {}).get("after_filtering", {})
print(ac.get("adapter_trimmed_reads", 0),
      ac.get("adapter_trimmed_bases", 0),
      d.get("duplication", {}).get("rate", 0),
      af.get("q30_rate", 0),
      af.get("gc_content", 0),
      d.get("insert_size", {}).get("peak", 0))
PY
)
    fi

    {
        printf '%s\t%s\n' \
            tool                 "$tool" \
            adapters             "$ADAPTERS" \
            trimmed              "$TRIM" \
            adapter_preset       "$ADAPTER_PRESET" \
            mapping_rate         "${mapped:-0}" \
            properly_paired_rate "${paired_pct:-0}" \
            error_rate           "${errrate:-0}" \
            softclip_frac        "${softclip:-0}" \
            insert_size_avg      "${isize:-0}" \
            insert_size_sd       "${isd:-0}" \
            avg_read_len         "${avlen:-0}" \
            avg_base_qual        "${avqual:-0}" \
            fastp_adapter_reads  "${adp_reads:-0}" \
            fastp_adapter_bases  "${adp_bases:-0}" \
            fastp_q30_rate       "${q30:-0}" \
            fastp_gc_content     "${gc:-0}" \
            fastp_dup_rate       "${dup:-0}" \
            fastp_insert_peak    "${ipeak:-0}"
    } > "$out"
    echo "  [$tool realism] wrote $(basename "$out")"
}

# ── Per-tool pipeline ─────────────────────────────────────────────────────
for tool in $TOOLS; do
    wd="$WORK_ROOT/$tool"; mkdir -p "$wd"
    r1="$wd/sample_r1.fastq.gz"; r2="$wd/sample_r2.fastq.gz"
    echo
    echo "──────── $tool ────────"

    case "$tool" in
        eidolon)
            truth="$wd/sample.vcf.gz"
            if [[ -f "$r1" && -f "$truth" ]]; then echo "  [eidolon 1/3] sim present — skipping."
            else echo "  [eidolon 1/3] simulating (SNP/indel only)..."; simulate_eidolon "$wd"; fi
            ;;
        neat)
            truth="$wd/sample_golden.vcf.gz"
            if ! conda run -n "$NEAT_ENV" which neat >/dev/null 2>&1; then
                echo "  [neat] NEAT 4 not found in env '$NEAT_ENV' — skipping NEAT arm." >&2
                continue
            fi
            if [[ -f "$r1" && -f "$truth" ]]; then echo "  [neat 1/3] sim present — skipping."
            else echo "  [neat 1/3] simulating (NEAT 4 default = SNP/indel)..."; simulate_neat "$wd"; fi
            ;;
        *) echo "  unknown tool '$tool' — skipping." >&2; continue ;;
    esac

    # Optional read-order randomization (SHUFFLE=1; off by default). eidolon (and
    # NEAT) write reads in genome-coordinate order, whereas a real sequencer's
    # FASTQ is in arbitrary order. The caller is order-independent (the BAM is
    # coordinate-sorted before calling, no order-based dup-marking), so this does
    # not change recall in principle — the flag exists to *demonstrate* that
    # (shuffled vs unshuffled diff) and to better mimic real input. seqkit's
    # per-index permutation with a fixed seed shuffles R1/R2 identically, keeping
    # mates in register.
    if [[ "$SHUFFLE" == "1" ]]; then
        sh1="$wd/shuf_r1.fastq.gz"; sh2="$wd/shuf_r2.fastq.gz"
        if [[ -f "$sh1" && -f "$sh2" ]]; then
            echo "  [$tool shuffle] shuffled FASTQ present — skipping."
        else
            echo "  [$tool shuffle] randomizing read order (seqkit, paired)..."
            seqkit shuffle -s "$SHUFFLE_SEED" -o "$sh1" "$r1"
            seqkit shuffle -s "$SHUFFLE_SEED" -o "$sh2" "$r2"
        fi
        r1="$sh1"; r2="$sh2"
    fi

    # Adapter-aware QC + optional trim. fastp runs once (report + trimmed reads);
    # TRIM=1 feeds the trimmed reads to alignment (realistic pipeline), TRIM=0
    # keeps the raw reads for alignment but still records the adapter report.
    if [[ "$MEASURE_REALISM" == "1" || "$TRIM" == "1" ]]; then
        ensure_fastp || { echo "  [$tool] fastp unavailable — skipping QC/trim." >&2; }
        if conda run -n bioinf which fastp >/dev/null 2>&1; then
            run_fastp "$wd" "$r1" "$r2"
            if [[ "$TRIM" == "1" && -s "$wd/trimmed_r1.fastq.gz" ]]; then
                echo "  [$tool trim] aligning fastp-trimmed reads."
                r1="$wd/trimmed_r1.fastq.gz"; r2="$wd/trimmed_r2.fastq.gz"
            fi
        fi
    fi

    align_call_score "$tool" "$wd" "$r1" "$r2" "$truth"
    # Stage the summary at the top level with a tool tag for the report.
    [[ -f "$wd/scored.summary.csv" ]] && cp -f "$wd/scored.summary.csv" "$OUTDIR/${tool}_scored.summary.csv"

    # Realism signals (depends on the aligned BAM, so after align_call_score).
    if [[ "$MEASURE_REALISM" == "1" ]]; then
        capture_realism "$tool" "$wd"
        [[ -f "$wd/realism.tsv" ]] && cp -f "$wd/realism.tsv" "$OUTDIR/${tool}_realism.tsv"
    fi

    # Stage durable per-tool artifacts back to Lustre when working locally: drop
    # the bulky FASTQ/BAM/gVCF (kept only on the node's local disk, discarded with
    # it), then copy the small remainder (scores, called/truth VCFs, QC, logs) to
    # $OUTDIR/$tool for provenance. The two files the collector needs were already
    # copied to $OUTDIR top-level above; this makes PRUNE=1 a no-op here.
    if [[ "$LOCAL_STAGE" == "1" ]]; then
        rm -f "$wd"/sample_r1.fastq.gz "$wd"/sample_r2.fastq.gz \
              "$wd"/trimmed_r1.fastq.gz "$wd"/trimmed_r2.fastq.gz \
              "$wd"/shuf_r1.fastq.gz "$wd"/shuf_r2.fastq.gz \
              "$wd"/aligned.bam "$wd"/aligned.bam.bai \
              "$wd"/raw.g.vcf.gz "$wd"/raw.g.vcf.gz.tbi
        mkdir -p "$OUTDIR/$tool"
        cp -rf "$wd/." "$OUTDIR/$tool/" 2>/dev/null || true
    fi
done

# ── Side-by-side summary ──────────────────────────────────────────────────
echo
echo "════════════════════════════════════════════════════════════════"
echo "Germline fidelity comparison — Job $SLURM_JOB_ID  (${COVERAGE}x, $REF_NAME)"
echo "Recall/precision/F1 of GATK recovering each simulator's truth set."
for tool in $TOOLS; do
    f="$OUTDIR/${tool}_scored.summary.csv"
    echo
    echo "### $tool"
    [[ -f "$f" ]] && column -t -s',' "$f" || echo "  (no summary — arm skipped or failed)"
done

# Persist report artifacts off scratch + record provenance/resources.
archive_run germline "$OUTDIR" \
    eidolon_scored.summary.csv neat_scored.summary.csv \
    eidolon_realism.tsv neat_realism.tsv

# Optional cleanup of bulky intermediates once scored + archived, so the
# adapter matrix / large sweeps don't exhaust the shared scratch project quota.
# Per tool, only prune when THAT tool scored successfully (scored.summary.csv
# present) so a failed/partial arm keeps its inputs for resume. Keeps the small
# durable bits: scores, called/truth VCFs, realism.tsv, samstats, fastp report.
if [[ "$PRUNE" == "1" && "$LOCAL_STAGE" != "1" ]]; then
    for tool in $TOOLS; do
        wd="$OUTDIR/$tool"
        [[ -f "$wd/scored.summary.csv" ]] || { echo "[prune] $tool not scored — keeping inputs."; continue; }
        echo "[prune] $tool — removing FASTQ/BAM/gVCF (keeping scores, VCFs, QC)..."
        rm -f "$wd"/sample_r1.fastq.gz "$wd"/sample_r2.fastq.gz \
              "$wd"/trimmed_r1.fastq.gz "$wd"/trimmed_r2.fastq.gz \
              "$wd"/shuf_r1.fastq.gz "$wd"/shuf_r2.fastq.gz \
              "$wd"/aligned.bam "$wd"/aligned.bam.bai \
              "$wd"/raw.g.vcf.gz "$wd"/raw.g.vcf.gz.tbi
    done
fi
