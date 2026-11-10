#!/usr/bin/env bash
# FAIL-FAST smoke test: prove a build + config + reference actually produce VALID
# output on ONE small window before you commit a big exclusive-node array to it.
#
# Runs a single small gen-reads window (a few minutes, one core, no --exclusive)
# and then validates the output THE WAY A DOWNSTREAM ALIGNER WOULD — most
# importantly that every FASTQ record has len(seq) == len(qual). That exact defect
# (a malformed FASTQ that `zcat | wc -l` cannot see, but bwa-mem2 silently
# truncates on) cost a full at-scale run during adapter validation (#125). This
# catches it in minutes, for the price of one small window.
#
# Run it whenever you change the build, the config knobs, or the reference:
#   REFERENCE=$SCRATCH/neat_data/GRCh38.fa bash scripts/delta/wg_smoke.sh
#   REFERENCE=... COVERAGE=30 ADAPTERS=1 WINDOW_BP=2000000 bash scripts/delta/wg_smoke.sh
#
# Runnable on a small interactive allocation (salloc/srun --pty) or as the payload
# of a tiny sbatch. Exits non-zero (loudly) on any problem so you find out now.
set -uo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH resolution

REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
WINDOW_BP="${WINDOW_BP:-2000000}"        # size of the single test window
COVERAGE="${COVERAGE:-30}"
READ_LEN="${READ_LEN:-151}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-50}"
PLOIDY="${PLOIDY:-2}"
SV_RATE_SCALE="${SV_RATE_SCALE:-0}"
ADAPTERS="${ADAPTERS:-0}"                # 1 → exercise the 3' adapter readthrough path too
OUTDIR="${OUTDIR:-$SCRATCH/wg_smoke_$$}"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="${EIDOLON_BIN:-$CARGO_TARGET_DIR/release/eidolon}"
source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0 2>/dev/null || true

fail() { echo "SMOKE: FAIL — $*" >&2; exit 1; }

[[ -x "$EIDOLON_BIN" ]] || fail "eidolon binary not found/executable: $EIDOLON_BIN (run setup.sh / cargo build --release)"
[[ -f "$REFERENCE" ]] || fail "reference not found: $REFERENCE"
[[ -f "$REFERENCE.fai" ]] || samtools faidx "$REFERENCE" 2>/dev/null || fail "cannot index reference (need samtools)"

echo "=== eidolon smoke test ==="
echo "  binary:    $EIDOLON_BIN ($("$EIDOLON_BIN" --version 2>/dev/null || echo '??'))"
echo "  reference: $REFERENCE"
echo "  window:    first ${WINDOW_BP} bp of the first contig   cov=${COVERAGE}x  adapters=$ADAPTERS"

# One small window on the first contig — the cheapest possible representative run.
read -r CTG LEN _ < "$REFERENCE.fai"
END=$(( WINDOW_BP < LEN ? WINDOW_BP : LEN ))
mkdir -p "$OUTDIR"
BED="$OUTDIR/smoke.bed"
printf '%s\t0\t%d\n' "$CTG" "$END" > "$BED"

CFG="$OUTDIR/smoke.yml"
cat > "$CFG" <<YAML
reference: $REFERENCE
read_len: $READ_LEN
coverage: $COVERAGE
ploidy: $PLOIDY
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $OUTDIR
output_prefix: smoke
rng_seed: eidolon smoke $CTG
sv_rate_scale: $SV_RATE_SCALE
target_bed: $BED
num_threads: 1
YAML
# Exercise the adapter readthrough path (where the malformed-FASTQ bug lived) when
# asked. Short inserts (mean < read_len) are what actually trigger 3' readthrough.
if [[ "$ADAPTERS" == "1" ]]; then
    cat >> "$CFG" <<YAML
keep_short_fragments: true
adapters:
  enabled: true
  preset: truseq
YAML
    # re-point the fragment size short so readthrough fires (overrides above)
    sed -i 's/^fragment_mean:.*/fragment_mean: 120/; s/^fragment_st_dev:.*/fragment_st_dev: 60/' "$CFG"
fi

echo "  running gen-reads (log-level info)..."
if ! "$EIDOLON_BIN" --log-level info gen-reads -c "$CFG" > "$OUTDIR/run.log" 2>&1; then
    tail -20 "$OUTDIR/run.log" >&2
    fail "gen-reads exited non-zero (see $OUTDIR/run.log)"
fi

R1="$OUTDIR/smoke_r1.fastq.gz"; R2="$OUTDIR/smoke_r2.fastq.gz"; VCF="$OUTDIR/smoke.vcf.gz"
[[ -f "$R1" && -f "$R2" ]] || fail "expected FASTQ not produced ($R1 / $R2)"
[[ -f "$VCF" ]] || fail "expected truth VCF not produced ($VCF)"

# ── The critical check: seq/qual length parity on EVERY record ──────────────────
# (a too-long qual line is the malformed-FASTQ signature; 4-line count stays intact
#  so wc can't see it, but an aligner stops at the first offending record.)
check_fastq() {
    local f="$1" label="$2"
    zcat "$f" | awk -v L="$label" '
        NR%4==1 { if ($0 !~ /^@/) { print "  "L": record "int(NR/4)+1" header missing @" > "/dev/stderr"; bad++ } }
        NR%4==2 { s=length($0) }
        NR%4==0 { q=length($0); recs++; if (q!=s) { badlen++; if (badlen<=3) printf "  %s: record %d seq=%d qual=%d\n", L, recs, s, q > "/dev/stderr" } }
        END {
            printf "  %s: %d records, %d bad-length, %d bad-header\n", L, recs, badlen+0, bad+0
            if (recs==0)   { print "ZERO_RECORDS" }
            else if (badlen>0 || bad>0) { print "MALFORMED" }
            else           { print "OK" }
        }'
}
echo "  validating FASTQ record integrity..."
v1=$(check_fastq "$R1" R1); s1=$(echo "$v1" | tail -1); echo "$v1" | sed '$d'
v2=$(check_fastq "$R2" R2); s2=$(echo "$v2" | tail -1); echo "$v2" | sed '$d'
[[ "$s1" == "OK" && "$s2" == "OK" ]] || fail "FASTQ integrity check failed (R1=$s1 R2=$s2) — DO NOT run this config at scale"

# ── Sanity: reads and truth variants actually exist ─────────────────────────────
n_reads=$(( $(zcat "$R1" | wc -l) / 4 ))
n_vars=$(zcat "$VCF" | grep -vc '^#' || true)
(( n_reads > 0 )) || fail "zero read pairs produced"
(( n_vars  > 0 )) || echo "  NOTE: 0 truth variants in this small window (can be legitimate for a low-mutation region)"
echo "  reads=$n_reads pairs   truth variants=$n_vars"

echo
echo "════════════════════════════════════════════════════════════════"
echo "SMOKE: PASS — build+config+reference produce valid output."
echo "  Safe to scale up (make_shard_beds.sh + genome_array.sbatch)."
echo "  Artifacts: $OUTDIR   (rm -rf when done)"
echo "════════════════════════════════════════════════════════════════"
