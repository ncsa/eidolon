#!/bin/bash
# SLURM job: exercise eidolon's MODEL-BUILDER subcommands on REAL, full-size,
# messy inputs — the half of eidolon the read-generation validation never touches.
#
# WHY
#   All Phase-1/2 Delta validation ran gen-reads / gen-cancer-reads against the
#   BUNDLED models. The builders that a real user runs on THEIR OWN data first —
#   gen-mut-model, gen-seq-error-model, gen-frag-length-model, gen-gc-bias-model,
#   gen-bam-models — are only covered by `eidolon/tests/model_parity.rs`, which runs
#   them on FIXED TINY fixtures (a 1.7 kb H1N1 reference, ~50 reads) and checks
#   byte-parity. That proves the algorithm doesn't drift; it proves NOTHING about
#   behavior on a 30 GB BAM, a real quality profile, a messy multi-thousand-contig
#   reference, or a VCF full of odd contigs — exactly where a first-contact user
#   failure (OOM, panic on an unexpected record, pathological runtime) hides.
#
# WHAT THIS DOES
#   1. Runs each builder for which an input is provided, under `/usr/bin/time -v`
#      (captures wall-clock + peak RSS — the resource behavior we have no data on).
#   2. Validates each output model: gzip-decodes, parses as JSON, is non-empty.
#   3. ROUND-TRIP: feeds the freshly-built models back into gen-reads and validates
#      the simulated FASTQ the way an aligner would (seq/qual length parity). This
#      is the real test — it proves a model built from real data is USABLE, not
#      merely parseable.
#
# INPUTS (provide whichever you have; each builder is skipped if its input is unset)
#   REFERENCE     reference FASTA                (gc-bias, mut-model, round-trip)
#   INPUT_BAM     a real aligned BAM             (frag-length, gc-bias, bam-models)
#   INPUT_FASTQ   a real FASTQ(.gz)              (seq-error-model)
#   INPUT_VCF     a real VCF(.gz) of variants    (mut-model)
#
# STRESS-TEST DATA — pick inputs that are big AND messy (see docs discussion):
#   * a full-size real Illumina BAM (tens of GB, real MAPQ/TLEN/dup/supplementary,
#     multiple read groups, unmapped reads) — not a toy alignment;
#   * a messy reference: many contigs / scaffolds, long N runs, IUPAC codes
#     (a less-vetted non-model organism assembly is ideal — soybean, or a
#     fragmented draft — precisely the reference that also stresses gen-reads);
#   * a population/clinical VCF with diverse variant types and contig naming.
#
# USAGE
#   REFERENCE=$SCRATCH/neat_data/genome.fa INPUT_BAM=$SCRATCH/data/real.bam \
#   INPUT_FASTQ=$SCRATCH/data/real_R1.fastq.gz INPUT_VCF=$SCRATCH/data/real.vcf.gz \
#     sbatch scripts/delta/model_builders.sbatch

#SBATCH --job-name=eidolon-modelbuild
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -uo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + archive_run

REFERENCE="${REFERENCE:-}"
INPUT_BAM="${INPUT_BAM:-}"
INPUT_FASTQ="${INPUT_FASTQ:-}"
INPUT_VCF="${INPUT_VCF:-}"
OUTDIR="${OUTDIR:-$SCRATCH/modelbuild_$SLURM_JOB_ID}"
MODELS="$OUTDIR/models"; mkdir -p "$MODELS"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="${EIDOLON_BIN:-$CARGO_TARGET_DIR/release/eidolon}"
source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0 2>/dev/null || true

echo "=== eidolon model-builder stress test: $SLURM_JOB_ID ==="
echo "eidolon: $("$EIDOLON_BIN" --version 2>/dev/null || echo '??')   git: $(git -C "$REPO_ROOT" describe --tags --always --dirty 2>/dev/null || echo '?')"
echo "ref=$REFERENCE  bam=$INPUT_BAM  fastq=$INPUT_FASTQ  vcf=$INPUT_VCF"
# An input that is SET but missing stops the job here. Each builder below runs only when its
# input exists, so a mistyped or purged path used to skip that builder without a word while
# the job carried on: job 22772429 had INPUT_BAM pointing at a deleted BAM, skipped three
# builders, and reported only the round-trip failure that followed from it.
missing=0
for var in REFERENCE INPUT_BAM INPUT_FASTQ INPUT_VCF; do
    path="${!var}"
    if [[ -n "$path" && ! -s "$path" ]]; then
        echo "ERROR: $var is set but missing or empty: $path" >&2
        missing=1
    fi
done
[[ "$missing" -eq 0 ]] || exit 1
[[ -x "$EIDOLON_BIN" ]] || { echo "eidolon binary missing: $EIDOLON_BIN" >&2; exit 1; }

SUMMARY="$OUTDIR/summary.tsv"; : > "$SUMMARY"
overall=0

# Run one builder: name, subcommand, config-body (heredoc content), output model path.
run_builder() {
    local name="$1" sub="$2" cfg_body="$3" out="$4"
    local cfg="$MODELS/$name.yml" tlog="$MODELS/$name.time.txt"
    printf '%s\n' "$cfg_body" > "$cfg"
    echo "── $name ($sub) ──"
    if /usr/bin/time -v -o "$tlog" "$EIDOLON_BIN" --log-level info "$sub" -c "$cfg" \
            > "$MODELS/$name.run.log" 2>&1; then
        local wall rss
        wall=$(awk -F': ' '/Elapsed \(wall/{print $2}' "$tlog")
        rss=$(awk -F': ' '/Maximum resident/{printf "%.0f", $2/1024}' "$tlog")
        # validate: gzip-decodes to non-empty valid JSON
        if [[ -f "$out" ]] && zcat "$out" 2>/dev/null | python3 -c 'import sys,json; json.load(sys.stdin)' 2>/dev/null; then
            echo "   PASS  wall=$wall  peakRSS=${rss}MB  -> $out"
            printf '%s\tPASS\t%s\t%s\n' "$name" "$wall" "${rss}MB" >> "$SUMMARY"
        else
            echo "   FAIL  builder exited 0 but output is missing/invalid JSON: $out" >&2
            printf '%s\tFAIL-OUTPUT\t%s\t%s\n' "$name" "$wall" "${rss}MB" >> "$SUMMARY"; overall=1
        fi
    else
        echo "   FAIL  builder exited non-zero (see $MODELS/$name.run.log):" >&2
        tail -8 "$MODELS/$name.run.log" >&2
        printf '%s\tFAIL-RUN\t-\t-\n' "$name" >> "$SUMMARY"; overall=1
    fi
}

# ── seq-error-model (FASTQ) ────────────────────────────────────────────────────
if [[ -n "$INPUT_FASTQ" && -f "$INPUT_FASTQ" ]]; then
    run_builder seq_error gen-seq-error-model \
"fastq_file: $INPUT_FASTQ
output_file: $MODELS/seq_error.json.gz
overwrite_output: true" "$MODELS/seq_error.json.gz"
fi

# ── frag-length-model (BAM) ─────────────────────────────────────────────────────
if [[ -n "$INPUT_BAM" && -f "$INPUT_BAM" ]]; then
    run_builder frag_length gen-frag-length-model \
"input_file: $INPUT_BAM
output_file: $MODELS/frag_length.json.gz
overwrite_output: true
min_reads: 100" "$MODELS/frag_length.json.gz"
fi

# ── gc-bias-model (reference + BAM) ─────────────────────────────────────────────
if [[ -n "$INPUT_BAM" && -f "$INPUT_BAM" && -n "$REFERENCE" && -f "$REFERENCE" ]]; then
    run_builder gc_bias gen-gc-bias-model \
"reference: $REFERENCE
bam_file: $INPUT_BAM
output_file: $MODELS/gc_bias.json.gz
overwrite_output: true
min_mapq: 20
bed_file: .
window_size: 1000
window_stride: 1000
min_windows_per_bin: 5" "$MODELS/gc_bias.json.gz"
fi

# ── mut-model (reference + VCF) ─────────────────────────────────────────────────
if [[ -n "$INPUT_VCF" && -f "$INPUT_VCF" && -n "$REFERENCE" && -f "$REFERENCE" ]]; then
    run_builder mut_model gen-mut-model \
"reference: $REFERENCE
vcf_file: $INPUT_VCF
output_file: $MODELS/mut_model.json.gz
overwrite_output: true
bed_file: ." "$MODELS/mut_model.json.gz"
fi

# ── bam-models (unified frag-length + GC in one pass over the BAM) ──────────────
if [[ -n "$INPUT_BAM" && -f "$INPUT_BAM" && -n "$REFERENCE" && -f "$REFERENCE" ]]; then
    run_builder bam_models gen-bam-models \
"bam_file: $INPUT_BAM
min_mapq: 20
frag_length:
  output_file: $MODELS/bam_frag.json.gz
  overwrite_output: true
  min_reads: 100
gc_bias:
  reference: $REFERENCE
  output_file: $MODELS/bam_gc.json.gz
  overwrite_output: true
  window_size: 1000
  window_stride: 1000
  min_windows_per_bin: 5" "$MODELS/bam_gc.json.gz"
fi

# ── ROUND-TRIP: simulate FROM the freshly-built models, validate the output ─────
# The real proof: a model built from real data is usable by gen-reads and yields
# valid reads. Wire in whichever models were built successfully.
if [[ -n "$REFERENCE" && -f "$REFERENCE" ]] && grep -q 'PASS' "$SUMMARY"; then
    echo "── round-trip: gen-reads using the built models ──"
    RT="$OUTDIR/roundtrip"; mkdir -p "$RT"
    [[ -f "$REFERENCE.fai" ]] || samtools faidx "$REFERENCE"
    read -r CTG LEN _ < "$REFERENCE.fai"
    END=$(( 2000000 < LEN ? 2000000 : LEN ))          # one small window is plenty
    printf '%s\t0\t%d\n' "$CTG" "$END" > "$RT/rt.bed"
    {
        echo "reference: $REFERENCE"
        echo "read_len: 151"
        echo "coverage: 30"
        echo "ploidy: 2"
        echo "paired_ended: true"
        # No fragment_mean/st_dev: the built fragment_model is the fragment-length
        # source (paired-end now accepts a fragment_model alone — #355). This makes
        # the round-trip exercise the realistic "user supplies only their model" path.
        echo "produce_fastq: true"
        echo "produce_vcf: true"
        echo "produce_bam: false"
        echo "overwrite_output: true"
        echo "output_dir: $RT"
        echo "output_filename: rt"
        echo "rng_seed: eidolon modelbuild roundtrip"
        echo "target_bed: $RT/rt.bed"
        echo "num_threads: 1"
        [[ -f "$MODELS/seq_error.json.gz" ]]  && echo "sequence_error_model: $MODELS/seq_error.json.gz"
        [[ -f "$MODELS/frag_length.json.gz" ]] && echo "fragment_model: $MODELS/frag_length.json.gz"
        [[ -f "$MODELS/gc_bias.json.gz" ]]    && echo "gc_bias_model: $MODELS/gc_bias.json.gz"
        [[ -f "$MODELS/mut_model.json.gz" ]]  && echo "mutation_model: $MODELS/mut_model.json.gz"
    } > "$RT/rt.yml"
    if "$EIDOLON_BIN" --log-level info gen-reads -c "$RT/rt.yml" > "$RT/run.log" 2>&1; then
        bad=$(zcat "$RT/rt_r1.fastq.gz" 2>/dev/null | awk 'NR%4==2{s=length($0)} NR%4==0{if(length($0)!=s)b++} END{print b+0}')
        reads=$(( $(zcat "$RT/rt_r1.fastq.gz" 2>/dev/null | wc -l) / 4 ))
        if [[ "$bad" == "0" && "$reads" -gt 0 ]]; then
            echo "   PASS  round-trip produced $reads valid read pairs from built models"
            printf 'roundtrip\tPASS\treads=%s\t-\n' "$reads" >> "$SUMMARY"
        else
            echo "   FAIL  round-trip output invalid (bad-length records=$bad reads=$reads)" >&2
            printf 'roundtrip\tFAIL\treads=%s\tbadlen=%s\n' "$reads" "$bad" >> "$SUMMARY"; overall=1
        fi
    else
        echo "   FAIL  gen-reads from built models exited non-zero (see $RT/run.log):" >&2
        tail -8 "$RT/run.log" >&2
        printf 'roundtrip\tFAIL-RUN\t-\t-\n' >> "$SUMMARY"; overall=1
    fi
fi

# No builder ran → every input was missing/unreadable. This is NOT a pass:
# a stress test that exercises nothing must fail loudly, or a mis-staged run
# (wrong paths, purged scratch) silently reports green. (see job 19887446)
if [[ ! -s "$SUMMARY" ]]; then
    echo "── no builder ran ──" >&2
    echo "   FAIL  all inputs missing/unreadable — nothing was exercised." >&2
    echo "         checked: ref=$REFERENCE bam=$INPUT_BAM fastq=$INPUT_FASTQ vcf=$INPUT_VCF" >&2
    printf 'NONE\tFAIL-NO-INPUT\t-\t-\n' >> "$SUMMARY"
    overall=1
fi

echo
echo "════════════════════════════════════════════════════════════════"
echo "Model-builder stress test summary (name / status / wall / peakRSS):"
column -t -s$'\t' "$SUMMARY" | sed 's/^/  /'
echo "  artifacts: $OUTDIR"
echo "  overall: $([ "$overall" -eq 0 ] && echo PASS || echo FAIL)"
echo "════════════════════════════════════════════════════════════════"
archive_run modelbuild "$OUTDIR" summary.tsv 2>/dev/null || true
exit "$overall"
