#!/usr/bin/env bash
# Mutational-signature fidelity check (ticket #317, realism epic #311).
#
# Does eidolon reconstruct the COSMIC signature its tumor model encodes? This runs on
# the *truth* somatic SNVs (what eidolon generated) — so it validates the simulator's
# trinucleotide/signature machinery directly, independent of any variant caller.
# Extracts somatic SNVs from an eidolon truth VCF and fits COSMIC SBS signatures with
# SigProfilerAssignment; the assigned profile should reflect the model, not a
# flat/random spectrum.
#
# Usage (from repo root; needs the `sigprofiler` env from setup.sh):
#   TRUTH_VCF=$SCRATCH/cancer_<id>/sample_merged_truth.vcf.gz \
#   OUTDIR=$SCRATCH/sigcheck_<id> bash scripts/delta/signature_check.sh
#
# Run via a compute node for the first invocation — it installs the SigProfiler
# reference genome (~3 GB, one-time, slow):
#   srun -A bhrd-delta-cpu -p cpu -N1 -n1 -c8 --mem=32G -t 02:00:00 \
#     bash -lc 'TRUTH_VCF=... OUTDIR=... bash scripts/delta/signature_check.sh'
#
# NOTE: this is a first-draft harness — SigProfiler's API/genome-naming may need a
# tweak on the first real run (like sv_pipeline did). SigProfiler GRCh38 uses
# Ensembl-style contig names (1/2/3); a chr-prefixed truth VCF (chr22.fa runs) may
# need `bcftools annotate --rename-chrs` first.
set -euo pipefail

REPO_ROOT="${EIDOLON_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}}"
source "$REPO_ROOT/scripts/delta/lib_report.sh"

TRUTH_VCF="${TRUTH_VCF:?set TRUTH_VCF to an eidolon truth VCF (INFO/EIDOLON_ORIGIN-tagged)}"
OUTDIR="${OUTDIR:-$SCRATCH/sigcheck_$(basename "$(dirname "$TRUTH_VCF")")}"
GENOME="${GENOME:-GRCh38}"
SIGPROFILER_ENV="${SIGPROFILER_ENV:-sigprofiler}"
# Minimum somatic SNVs before a signature fit is meaningful. Keeps the historical 50
# so this change enforces the existing bar rather than silently moving it — but note
# SBS-96 has 96 contexts, so 50 is a hard floor and not a sufficiency threshold.
MIN_SNVS="${MIN_SNVS:-50}"

setup_conda
mkdir -p "$OUTDIR/vcf"

# 1. Extract somatic SNVs only (SigProfiler SBS works on single-base substitutions).
conda_activate bioinf   # bcftools
bcftools view -v snps -i 'INFO/EIDOLON_ORIGIN="somatic"' \
    -O v -o "$OUTDIR/vcf/eidolon_somatic_snvs.vcf" "$TRUTH_VCF"
n=$(grep -vc '^#' "$OUTDIR/vcf/eidolon_somatic_snvs.vcf" || true)
echo "Somatic SNVs for signature fitting: $n"
# Fatal, not a warning. This harness produces exactly one thing — a signature fit —
# so too few SNVs is not a degraded result, it is not a result: SBS-96 has 96
# trinucleotide contexts, and a spectrum built from fewer mutations than contexts
# cannot populate them. SigProfiler will still return an assignment, and that
# assignment would be reported as though it measured something. The report's own
# signature analysis (§3.9) used 6,225 SNVs; MIN_SNVS is a floor, not a target.
if [[ "$n" -lt "$MIN_SNVS" ]]; then
    echo "ERROR: only $n somatic SNV(s) — below MIN_SNVS=$MIN_SNVS." >&2
    echo "  A COSMIC SBS fit needs to populate 96 trinucleotide contexts; below that" >&2
    echo "  floor the assignment is noise that would be reported as a measurement." >&2
    if [[ "$n" -eq 0 ]]; then
        echo "  0 SNVs usually means the truth VCF has no somatic SNVs at all — raise" >&2
        echo "  TUMOR_MUTATION_RATE or coverage. (A pre-v3.0.0 truth VCF uses" >&2
        echo "  NEAT_ORIGIN and would have failed the filter above instead; convert it" >&2
        echo "  with tools/migrate_legacy_tokens.sh.)" >&2
    else
        echo "  Use a larger / higher-coverage run, or set MIN_SNVS to accept the" >&2
        echo "  noise deliberately." >&2
    fi
    exit 1
fi

# 2. SigProfiler: install the reference genome once (no-op if present), then fit.
conda_activate "$SIGPROFILER_ENV"
python3 - "$GENOME" <<'PY'
import sys
genome = sys.argv[1]
try:
    from SigProfilerMatrixGenerator import install as genInstall
    genInstall.install(genome, rsync=False, bash=True)   # no-op if already installed
except Exception as e:
    print(f"(genome install note: {e})")
PY

python3 - "$OUTDIR" "$GENOME" <<'PY'
import sys
outdir, genome = sys.argv[1], sys.argv[2]
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples=f"{outdir}/vcf",
                   output=f"{outdir}/sigprofiler",
                   input_type="vcf",
                   genome_build=genome,
                   collapse_to_SBS96=True)
print("SigProfilerAssignment cosmic_fit complete.")
PY

# 3. Report assigned signatures.
echo
echo "════════════════════════════════════════════════════════════════"
echo "Signature fidelity check — $OUTDIR  ($n somatic SNVs)"
echo "Assigned COSMIC signature activities:"
act=$(find "$OUTDIR/sigprofiler" -name '*Activities*.txt' 2>/dev/null | head -1)
[[ -f "$act" ]] && column -t "$act" || echo "  (activities table not found — see $OUTDIR/sigprofiler)"
cat <<'EOF'

Interpretation: the dominant assigned signature(s) should reflect the tumor model's
COSMIC composition, not a flat/uniform spectrum. A sensible COSMIC fit = eidolon
reproduces the signature it was given (trinucleotide machinery faithful). Per-tissue
SNV signature is pan-cancer in the current models (only the SV side is tissue-
specific), so expect a pan-cancer-like profile rather than tissue-defining SBS.
EOF
