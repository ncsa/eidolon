#!/usr/bin/env bash
# Print the SHAPE of one or more sequencing error models: option set, read length, binning.
#
# WHY. #677 asks whether the shipped default quality model still describes a current
# instrument. The shipped default is 101 bp with 42 continuous scores (0-41) and
# `binned_scores: false`; modern instruments commonly emit BINNED scores at 151 bp or more.
# #677's suggested first step is "build a quality model from a modern library and compare its
# shape against the shipped default" -- and #720 already did exactly that, twice, on HG002 R1
# and R2. This reads those models rather than fitting anything new.
#
# It answers the structural half only. The per-cycle quality PROFILE of the real library is
# already measured (measure_quality_degradation.sh section 4); what is missing is what the
# fitted MODEL's option set looks like, which is what decides whether a binned default is the
# representative one.
#
# Read-only. Reads the model files and writes nothing.
#
# USAGE
#   bash scripts/delta/describe_quality_model.sh $SCRATCH/qualdeg_r1/model_degraded.json.gz ...
#
# Compare against the shipped default, which is in the repo:
#   bash scripts/delta/describe_quality_model.sh eidolon-core/src/models/model_data/default_sequencing_error_model.json.gz

#SBATCH --job-name=eidolon-modelshape
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=modelshape-%j.log

set -euo pipefail

[[ $# -ge 1 ]] || { echo "usage: describe_quality_model.sh <model.json.gz> [...]" >&2; exit 1; }
for f in "$@"; do
    [[ -f "$f" ]] || { echo "FATAL: not a file: $f" >&2; exit 1; }
done

python3 - "$@" <<'PY'
import gzip, json, sys

# Phred+33 encodes 31 as '@'. eidolon never emits it as the FIRST character of a quality line
# (README, "Q31 and the first quality character"), substituting the nearest other option. That
# is a deliberate, bounded departure from the fitted model, and this reports its size.
AT_SYMBOL = 31

def qsm_of(m):
    # Accepts either a full SequencingErrorModel or a bare QualityScoreModel, because the
    # shipped default is stored as the latter and a fitted model as the former. Guessing
    # wrong would silently describe the wrong object, so key on a field only the outer
    # struct has.
    if "quality_score_model" in m:
        return m["quality_score_model"], m
    return m, None

for path in sys.argv[1:]:
    with gzip.open(path, "rt") as fh:
        raw = json.load(fh)
    q, outer = qsm_of(raw)
    opts = q["quality_score_options"]
    n = len(opts)
    contiguous = opts == list(range(opts[0], opts[0] + n))
    gaps = [b - a for a, b in zip(opts, opts[1:])]

    print(f"=== {path}")
    print(f"  kind                 {'SequencingErrorModel' if outer else 'QualityScoreModel (bare)'}")
    print(f"  assumed_read_length  {q['assumed_read_length']}")
    print(f"  binned_scores flag   {q.get('binned_scores')}")
    print(f"  options              {n} values, Q{opts[0]}-Q{opts[-1]}")
    print(f"  contiguous?          {contiguous}" + ("" if contiguous else f"   gaps: {sorted(set(gaps))}"))
    if not contiguous:
        print(f"  option values        {opts}")
    print(f"  transition positions {len(q['distros_from_one'])}")
    # Position 1, which is where the Q31 -> '@' suppression acts. `weights` is a stored CDF
    # (starts at 0.0, ends at 1.0), so the mass on a value is its step, and the first value's
    # mass is its own entry. Getting that off by one would misreport the size of the deviation,
    # which is the whole point of printing it.
    sd = q["seed_dist"]
    vals, cdf = sd["values"], sd["weights"]
    mass = [cdf[0]] + [b - a for a, b in zip(cdf, cdf[1:])]
    top = sorted(zip(mass, vals), reverse=True)[:4]
    print("  seed (position 1)    " + ", ".join(f"Q{v} {100*m:.2f}%" for m, v in top))
    if AT_SYMBOL in vals:
        at = mass[vals.index(AT_SYMBOL)]
        sub = min((v for v in vals if v != AT_SYMBOL),
                  key=lambda v: (abs(v - AT_SYMBOL), v), default=None)
        print(f"  Q{AT_SYMBOL} at position 1     {100*at:.3f}% of reads, rewritten to Q{sub}"
              f"  <- the deviation from the fit")
    else:
        print(f"  Q{AT_SYMBOL} at position 1     not an option; the rewrite never fires")

    deg = q.get("degradation")
    if deg is not None:
        print(f"  degradation          present, read_fraction {deg['read_fraction']:.4f}, "
              f"{len(deg['degraded_distros'])} positions")
    else:
        print(f"  degradation          none")
    if outer is not None:
        print(f"  error_rate           {outer.get('error_rate')}")
        print(f"  indel_probability    {outer.get('indel_probability')}")
    print()
PY
