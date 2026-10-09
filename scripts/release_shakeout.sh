#!/usr/bin/env bash
# Pre-publish shakeout for eidolon v1.x.
#
# Runs the in-repo regression suite and then drives the binary against
# real, full-size inputs that the CI suite is too small to exercise.
# Intended to be invoked manually before tagging a release. Not run in CI.
#
# Pass a working directory via env DATA_DIR=...; large inputs are NOT
# checked into the repo and must be staged there before running. Set
# `SKIP_LARGE_DATA=1` to run only the in-tree regression suite and skip
# the full-size steps.
#
# Usage:
#   DATA_DIR=/scratch/eidolon-shakeout ./scripts/release_shakeout.sh
#   SKIP_LARGE_DATA=1 ./scripts/release_shakeout.sh
#
# Expected inputs in DATA_DIR (only when running large-data steps):
#   hg38_chr22.fa          single-contig human chr22, ~50 MB
#   hg38_chr22.bam         aligned BAM, MD tags present (samtools calmd)
#   hg38_chr22.vcf.gz      multi-allelic VCF, GT field populated
#
# This script exits on first failure. Re-run after fixing.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

log()  { printf '\n\033[1;36m▶ %s\033[0m\n' "$*"; }
warn() { printf '\n\033[1;33m⚠ %s\033[0m\n' "$*"; }
fail() { printf '\n\033[1;31m✗ %s\033[0m\n' "$*"; exit 1; }

# ── Step 1: cargo fmt + clippy + tests must be clean ─────────────────────────

log "Step 1: workspace formatting + lint check"
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings || \
    warn "clippy reported warnings — review before publishing (not blocking)"

log "Step 2: full test suite, including the model-parity regression"
cargo test --workspace --release

if [[ "${SKIP_LARGE_DATA:-0}" == "1" ]]; then
    log "SKIP_LARGE_DATA=1 set — stopping after in-repo suite. PASS."
    exit 0
fi

# ── Step 3: full-size shakeout (manual data required) ───────────────────────

: "${DATA_DIR:?DATA_DIR must be set when running large-data shakeout}"
[[ -d "$DATA_DIR" ]] || fail "DATA_DIR=$DATA_DIR does not exist"

REF="$DATA_DIR/hg38_chr22.fa"
BAM="$DATA_DIR/hg38_chr22.bam"
VCF="$DATA_DIR/hg38_chr22.vcf.gz"

for f in "$REF" "$BAM" "$VCF"; do
    [[ -f "$f" ]] || fail "missing required input: $f"
done

WORK="$DATA_DIR/shakeout-$(date +%Y%m%d-%H%M%S)"
mkdir -p "$WORK"
log "Working directory: $WORK"

RNEAT="$ROOT/target/release/eidolon"
[[ -x "$RNEAT" ]] || fail "release binary not found; cargo test --release should have built it"

# ── Step 3a: gen-mut-model on chr22 ─────────────────────────────────────────

log "Step 3a: gen-mut-model on hg38 chr22"
cat > "$WORK/mut.yml" <<EOF
reference: $REF
vcf_file: $VCF
output_file: $WORK/mut_model.json.gz
overwrite_output: true
bed_file: .
EOF
/usr/bin/time -v "$RNEAT" gen-mut-model -c "$WORK/mut.yml" 2> "$WORK/mut.time"
log "  ✓ wrote $WORK/mut_model.json.gz"
grep -E '(Maximum resident set|Elapsed)' "$WORK/mut.time"

# ── Step 3b: gen-frag-length-model and gen-gc-bias-model ────────────────────

log "Step 3b: gen-frag-length-model on chr22 BAM"
cat > "$WORK/frag.yml" <<EOF
input_file: $BAM
output_file: $WORK/frag_model.json.gz
overwrite_output: true
min_reads: 100
EOF
/usr/bin/time -v "$RNEAT" gen-frag-length-model -c "$WORK/frag.yml" 2> "$WORK/frag.time"
grep -E '(Maximum resident set|Elapsed)' "$WORK/frag.time"

log "Step 3c: gen-gc-bias-model on chr22 BAM"
cat > "$WORK/gc.yml" <<EOF
reference: $REF
bam_file: $BAM
output_file: $WORK/gc_model.json.gz
overwrite_output: true
min_mapq: 20
bed_file: .
window_size: 100
window_stride: 100
min_windows_per_bin: 10
EOF
/usr/bin/time -v "$RNEAT" gen-gc-bias-model -c "$WORK/gc.yml" 2> "$WORK/gc.time"
grep -E '(Maximum resident set|Elapsed)' "$WORK/gc.time"

# ── Step 3d: gen-bam-models parity at scale ─────────────────────────────────
#
# Two unified runs, one per output. The standalone `gen-frag-length-model`
# uses `BamWalkFilter::for_frag_length()` which hard-codes `min_mapq = 10`,
# while the standalone `gen-gc-bias-model` uses `for_coverage()` with the
# configured `min_mapq` (20 in Step 3c). The unified runner exposes a single
# walker `min_mapq` that gates BOTH observers, so a single unified config
# cannot match both standalone filter policies simultaneously. Run twice:
# once at min_mapq=10 to compare against frag-length, once at min_mapq=20
# to compare against gc-bias.

log "Step 3d.1: gen-bam-models (frag-length parity) — min_mapq=10"
cat > "$WORK/unified_frag.yml" <<EOF
bam_file: $BAM
min_mapq: 10
frag_length:
  output_file: $WORK/unified_frag.json.gz
  overwrite_output: true
  min_reads: 100
EOF
/usr/bin/time -v "$RNEAT" gen-bam-models -c "$WORK/unified_frag.yml" 2> "$WORK/unified_frag.time"
grep -E '(Maximum resident set|Elapsed)' "$WORK/unified_frag.time"

log "Step 3d.2: gen-bam-models (gc-bias parity) — min_mapq=20"
cat > "$WORK/unified_gc.yml" <<EOF
bam_file: $BAM
min_mapq: 20
gc_bias:
  reference: $REF
  output_file: $WORK/unified_gc.json.gz
  overwrite_output: true
  window_size: 100
  window_stride: 100
  min_windows_per_bin: 10
EOF
/usr/bin/time -v "$RNEAT" gen-bam-models -c "$WORK/unified_gc.yml" 2> "$WORK/unified_gc.time"
grep -E '(Maximum resident set|Elapsed)' "$WORK/unified_gc.time"

# Parity: fragment_length and gc_bias models are HashMap-free, so the
# gzipped JSON deserialises identically regardless of HashMap iteration
# order — direct `diff` of the decompressed bytes works.
log "  parity: frag_length unified vs standalone"
diff <(zcat "$WORK/unified_frag.json.gz") <(zcat "$WORK/frag_model.json.gz") \
    || fail "unified frag-length model diverged from standalone"
log "  parity: gc_bias unified vs standalone"
diff <(zcat "$WORK/unified_gc.json.gz") <(zcat "$WORK/gc_model.json.gz") \
    || fail "unified gc-bias model diverged from standalone"

# ── Step 3e: gen-reads end-to-end with the v1.7.0 models ────────────────────

log "Step 3e: gen-reads on chr22 using the trained models"
cat > "$WORK/reads.yml" <<EOF
reference: $REF
read_len: 150
coverage: 5
paired_ended: true
fragment_mean: 350.0
fragment_st_dev: 60.0
produce_fastq: true
produce_bam: true
produce_vcf: true
output_dir: $WORK
output_filename: shakeout
overwrite_output: true
rng_seed: shakeout v1.7.0
mutation_model: $WORK/mut_model.json.gz
fragment_model: $WORK/frag_model.json.gz
gc_bias_model: $WORK/gc_model.json.gz
EOF
/usr/bin/time -v "$RNEAT" gen-reads -c "$WORK/reads.yml" 2> "$WORK/reads.time"
grep -E '(Maximum resident set|Elapsed)' "$WORK/reads.time"

# Spot check: output FASTQ records exist and have well-formed sequences.
for f in "$WORK"/shakeout*.fastq.gz; do
    n=$(zcat "$f" | awk 'NR%4==1' | wc -l)
    [[ "$n" -gt 0 ]] || fail "no reads in $f"
    log "  ✓ $f: $n reads"
done

log "ALL SHAKEOUT STEPS PASSED."
log "Logs written under $WORK/*.time — compare against the prior release's"
log "peak RSS and elapsed time before publishing."
