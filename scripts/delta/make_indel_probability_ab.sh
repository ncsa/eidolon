#!/usr/bin/env bash
# Build two sequencing error models that differ in ONE field: `indel_probability`.
#
# WHY. `somatic_snv_recall` sits 4.5 sigma below its baseline on every run since #660
# (issue #680), replicated on two builds. #660 corrected `indel_probability` from 0.4 to
# 0.01, and because an error is either an indel or a substitution, that raised the
# SUBSTITUTION share of sequencing errors from 0.60 to 0.99 -- a 1.65x higher SNV noise
# floor. Somatic calling works close to that floor; germline calling does not, which fits
# the observation that germline SNV recall is unchanged while somatic dropped.
#
# That is a hypothesis. This script makes it testable without real data: run the cancer
# pipeline twice, identical in every respect except this one field.
#
# ALREADY RUN, AND THE HYPOTHESIS FAILED (jobs 21823970 / 21823971). Reverting
# indel_probability to 0.4 moved SNV recall by +0.0136 -- five variants of 368 -- against
# a 0.081 gap, and both arms landed near 0.85 rather than the 0.921 baseline. #660 is not
# the cause. The script is kept because the A/B is the right shape for the next candidate;
# change the field it patches and the guidance below still applies.
#
# BOTH arms share one quality-score model -- the SHIPPED one, swapped into the built model --
# so the only difference between them is the field under test, and both sit at the operating
# point the tier baseline was measured at. The arms must not differ in the quality model, or
# the comparison says nothing about indel_probability.
#
# This previously pinned `error_rate` to the shipped value instead, which did nothing at all:
# that field is a fitted summary and generation never reads it (#724). Both arms ran at the
# built model's much quieter distribution, so the guard against an underpowered null was not
# in effect for the run recorded above.
#
# Usage:
#   bash scripts/delta/make_indel_probability_ab.sh <training.fastq.gz> <outdir>
#
# Then, for each arm:
#   SEQ_ERROR_MODEL=<outdir>/seq_error_p0.01.json.gz OUTDIR=... \
#     REFERENCE=$SCRATCH/neat_data/chr22.fa TOTAL_COVERAGE=30 PURITY=0.6 PRUNE=1 \
#     sbatch scripts/delta/cancer_pipeline.sbatch
set -euo pipefail

FASTQ="${1:?usage: make_indel_probability_ab.sh <training.fastq.gz> <outdir>}"
OUTDIR="${2:?usage: make_indel_probability_ab.sh <training.fastq.gz> <outdir>}"
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EIDOLON_BIN="${EIDOLON_BIN:-$CARGO_TARGET_DIR/release/eidolon}"

[[ -f "$FASTQ" ]]       || { echo "ERROR: training FASTQ not found: $FASTQ" >&2; exit 1; }
[[ -x "$EIDOLON_BIN" ]] || { echo "ERROR: eidolon not built at $EIDOLON_BIN" >&2; exit 1; }
mkdir -p "$OUTDIR"

BASE="$OUTDIR/seq_error_base.json.gz"
CFG="$OUTDIR/build.yml"
cat > "$CFG" <<EOF
fastq_file: $FASTQ
output_file: $BASE
overwrite_output: true
EOF

echo "[1/3] building the shared base model from $FASTQ ..."
"$EIDOLON_BIN" gen-seq-error-model -c "$CFG"
[[ -s "$BASE" ]] || { echo "ERROR: gen-seq-error-model produced nothing" >&2; exit 1; }

echo "[2/3] writing the two arms ..."
# The shipped default is one whole sequencing-error model as of v3.4.0; its quality half is a
# field inside it rather than its own file.
SHIPPED_QSM="$REPO_ROOT/eidolon-core/src/models/model_data/default_sequencing_error_model.json.gz"
[[ -f "$SHIPPED_QSM" ]] || { echo "ERROR: shipped model not found: $SHIPPED_QSM" >&2; exit 1; }

python3 - "$BASE" "$OUTDIR" "$SHIPPED_QSM" <<'PY'
import gzip, json, sys
base, outdir, shipped_qsm = sys.argv[1], sys.argv[2], sys.argv[3]
with gzip.open(base, "rt") as fh:
    model = json.load(fh)

# Fail loudly rather than writing two identical files: a silent no-op patch and a real one
# produce the same exit status, and the A/B would then compare a model against itself.
if "indel_probability" not in model:
    sys.exit("ERROR: built model has no indel_probability field; cannot patch it")
built = model["indel_probability"]
print(f"       built model carries indel_probability = {built}")

# Both arms are put at the shipped OPERATING POINT by giving them the shipped quality model.
# A model built from a modern library lands far quieter than the inherited default --
# 0.000392 against 0.006638 on HCC1395 normal, ~17x -- so arms built straight from real reads
# would run well below the point where the effect appears, and a null result could not be told
# apart from an underpowered one.
#
# THIS USED TO PIN `error_rate` INSTEAD, WHICH DID NOTHING (#724). That field is a fitted
# summary of the quality histogram and is never read during generation: errors are injected per
# base from the base's own quality score. Both arms therefore ran at the BUILT model's quiet
# distribution, which is exactly what the pin was meant to correct for -- so the guard against
# an underpowered null was not in effect for jobs 21823970 / 21823971 (#680).
#
# The quality model is what sets the rate, so swap that. `error_rate` is recomputed to stay a
# true description of the model it now describes rather than of the one it came from.
built_rate = model.get("error_rate")
with gzip.open(shipped_qsm, "rt") as fh:
    shipped = json.load(fh)
# v3.4.0 made the shipped default a whole SequencingErrorModel, so take its quality half.
if "quality_score_model" not in shipped:
    sys.exit("ERROR: shipped model has no quality_score_model field")
model["quality_score_model"] = shipped["quality_score_model"]
qsm = model["quality_score_model"]
# The field now has to describe the model that just replaced it. That value is known exactly --
# it is the shipped default's own summary, recorded in model_data/README.md -- so it is copied
# rather than re-derived. Re-deriving the whole-histogram rate means marginalizing the
# transition tensor over positions, and a position-1 approximation written into a field
# documented as the histogram's summary would be a wrong number rather than a missing one.
SHIPPED_ERROR_RATE = shipped.get("error_rate", 0.003774138199848705)
print(f"       built model error_rate = {built_rate}")
print(f"       swapped in the shipped quality model "
      f"({qsm['assumed_read_length']} bp, {len(qsm['quality_score_options'])} options)")
print(f"       error_rate set to the shipped model's own {SHIPPED_ERROR_RATE}")
model["error_rate"] = SHIPPED_ERROR_RATE

for value, tag in ((0.01, "p0.01"), (0.40, "p0.40")):
    model["indel_probability"] = value
    path = f"{outdir}/seq_error_{tag}.json.gz"
    with gzip.open(path, "wt") as fh:
        json.dump(model, fh)
    print(f"       wrote {path}  (indel_probability = {value})")
PY

echo "[3/3] verifying the two arms differ in exactly one field ..."
python3 - "$OUTDIR" <<'PY'
import gzip, json, sys
outdir = sys.argv[1]
a = json.load(gzip.open(f"{outdir}/seq_error_p0.01.json.gz", "rt"))
b = json.load(gzip.open(f"{outdir}/seq_error_p0.40.json.gz", "rt"))
if sorted(a) != sorted(b):
    sys.exit("ERROR: the two arms have different fields")
differing = [k for k in a if json.dumps(a[k], sort_keys=True) != json.dumps(b[k], sort_keys=True)]
if differing != ["indel_probability"]:
    sys.exit(f"ERROR: arms differ in {differing}, expected only ['indel_probability']")
print(f"       OK: differ only in indel_probability ({a['indel_probability']} vs {b['indel_probability']})")
# The substitution share is what the hypothesis is about; state it so the run can be read.
for m, tag in ((a, "p0.01"), (b, "p0.40")):
    p = m["indel_probability"]
    print(f"       {tag}: {1 - p:.2%} of sequencing errors are substitutions")
PY

cat <<EOF

Two arms are in $OUTDIR. Run the cancer pipeline once per arm, changing nothing else:

  for arm in p0.01 p0.40; do
    SEQ_ERROR_MODEL=$OUTDIR/seq_error_\$arm.json.gz \\
    OUTDIR=\$SCRATCH/somatic_ab_\$arm \\
    REFERENCE=\$SCRATCH/neat_data/chr22.fa TOTAL_COVERAGE=30 PURITY=0.6 PRUNE=1 \\
      sbatch $REPO_ROOT/scripts/delta/cancer_pipeline.sbatch
  done

Read somatic SNV recall from each run's scored.stats.csv:

  python3 -c "
import csv,sys
r={(x.get('type') or '').strip().lower():x for x in csv.DictReader(open(sys.argv[1]))}
print(r['snvs']['recall'], r['snvs']['precision'])" <outdir>/scored.stats.csv

COMPARE THE TWO ARMS TO EACH OTHER. Both arms now carry the SHIPPED quality model, so
they sit at the operating point the tier baseline was measured at -- which is the point of
the swap in step 2, and is what the old `error_rate` pin was trying and failing to do.
Everything else in the model still comes from your FASTQ, so treat a comparison against
historical runs as indicative and the arm-to-arm difference as the result.

  arms differ materially  -> indel_probability affects somatic calling at this operating
                             point, and the question becomes which value is right.
  arms within noise       -> it does not, and the hypothesis under test is eliminated.

For scale: 368 truth SNVs means one variant is ~0.003 of recall. A difference of a few
variants is not a result. Run replicates before concluding anything from a small gap.
EOF
