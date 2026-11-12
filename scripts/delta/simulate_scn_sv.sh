#!/bin/bash
# SLURM job: simulate reads FROM the SV-inclusive SCN model and verify structural variants
# actually RENDER in the golden VCF. This closes the loop end-to-end:
#   real SCN SVs (Delly) → gen-mut-model sv_model → gen-reads renders SVs into the truth set.
# It is the ultimate proof that eidolon's SV machinery works on a real complex organism, not just
# that the model file has an sv_model field.
#
# gen-reads only emits SVs when `sv_rate_scale > 0` (default 0 = opt-in) AND the mutation model
# carries a populated sv_model; sv_rate_scale=1.0 reproduces the model's trained per-base SV rate.
#
# Prereqs: call_scn_sv.sh + model_builders already run with the SV-inclusive VCF (INPUT_VCF=
# …all.vcf.gz), producing a mut_model.json.gz whose sv_model is populated. Pass MODELS=<that
# model_builders output>/models. Reference staged by stage_scn.sh. eidolon built (setup.sh).
#
# Usage:
#   MODELS=$SCRATCH/modelbuild_<JOB>/models sbatch scripts/delta/simulate_scn_sv.sh
#   MODELS=… COV=5 sbatch scripts/delta/simulate_scn_sv.sh          # lighter/heavier coverage

#SBATCH --job-name=eidolon-scnsvsim
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

D="${DATA_DIR:-$SCRATCH/neat_data/scn}"
ACC="${ACC:-GCA_040805935.1}"
REF="${REFERENCE:-$D/${ACC}.fa}"
MODELS="${MODELS:?set MODELS=<model_builders output>/models (the SV-inclusive build)}"
COV="${COV:-8}"
SV_SCALE="${SV_RATE_SCALE:-1.0}"      # 1.0 = the model's trained SV rate
OUTDIR="${OUTDIR:-$SCRATCH/scn_svsim_${SLURM_JOB_ID:-manual}}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
EIDOLON_BIN="${EIDOLON_BIN:-$CARGO_TARGET_DIR/release/eidolon}"

source "$HOME/.cargo/env" 2>/dev/null || true
module load samtools/1.22-cce19.0.0
module load htslib/1.22-gcc13.3.1

MUT="$MODELS/mut_model.json.gz"
[[ -s "$REF" ]] || { echo "reference not staged: $REF (run stage_scn.sh)" >&2; exit 1; }
[[ -s "$MUT" ]] || { echo "mut_model not found: $MUT (run model_builders with the SV-inclusive VCF)" >&2; exit 1; }
[[ -x "$EIDOLON_BIN" ]] || { echo "eidolon binary not found: $EIDOLON_BIN (run setup.sh)" >&2; exit 1; }

# Preflight: the model must actually carry an sv_model, else this test is meaningless.
python3 - "$MUT" <<'PY' || { echo "ABORT: mut_model has no sv_model — build from the SV-inclusive VCF first." >&2; exit 1; }
import gzip, json, sys
d = json.load(gzip.open(sys.argv[1]))
sv = d.get("sv_model")
assert sv is not None, "sv_model is None"
print(f"sv_model present: per_base_rate={sv.get('per_base_rate')} types={list(sv.get('type_probabilities',{}))}")
PY

mkdir -p "$OUTDIR"
cat > "$OUTDIR/sim.yml" <<YML
reference: $REF
read_len: 151
coverage: $COV
ploidy: 2
paired_ended: true
sv_rate_scale: $SV_SCALE
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: $OUTDIR
output_basename: scn_svsim
rng_seed: scn sv render check
num_threads: $THREADS
mutation_model: $MUT
YML
# The other built models are optional here (we're testing SV rendering, not error/frag fidelity),
# but include them if present so the reads are realistic too.
[[ -s "$MODELS/seq_error.json.gz" ]]  && echo "sequence_error_model: $MODELS/seq_error.json.gz" >> "$OUTDIR/sim.yml"
[[ -s "$MODELS/frag_length.json.gz" ]] && echo "fragment_model: $MODELS/frag_length.json.gz"     >> "$OUTDIR/sim.yml"
[[ -s "$MODELS/gc_bias.json.gz" ]]    && echo "gc_bias_model: $MODELS/gc_bias.json.gz"           >> "$OUTDIR/sim.yml"

echo "=== gen-reads (sv_rate_scale=$SV_SCALE, cov=$COV) from $MUT ==="
"$EIDOLON_BIN" --log-level info gen-reads -c "$OUTDIR/sim.yml"

GOLDEN="$OUTDIR/scn_svsim.vcf.gz"
[[ -s "$GOLDEN" ]] || { echo "no golden VCF at $GOLDEN — did produce_vcf run?" >&2; exit 1; }

echo
echo "── SVs RENDERED in the simulated golden VCF (by SVTYPE) ──"
zcat "$GOLDEN" | grep -v '^#' | grep -oE 'SVTYPE=[A-Za-z]+' | sort | uniq -c | sort -rn || true
nsv=$(zcat "$GOLDEN" | grep -v '^#' | grep -c 'SVTYPE=' || true)
ntot=$(zcat "$GOLDEN" | grep -vc '^#' || true)
echo
echo "════════════════════════════════════════════════════════════════"
echo "Simulated golden VCF: $ntot variants, of which $nsv are SVs (symbolic ALT / SVTYPE)."
if [[ "$nsv" -gt 0 ]]; then
    echo "✅ eidolon RENDERED structural variants from the SCN-trained sv_model — loop closed."
else
    echo "⚠️ NO SVs rendered — check sv_rate_scale>0 and that sv_model.per_base_rate>0."
fi
echo "Golden VCF: $GOLDEN"
echo "════════════════════════════════════════════════════════════════"
