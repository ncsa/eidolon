#!/usr/bin/env bash
# Pool truvari results across replicate sv_pipeline.sbatch runs.
#
# WHY REPLICATES. Recall/precision at the model's NOMINAL SV rate needs enough events to
# be meaningful, and the honest way to get them is more runs, not a higher rate.
# `SV_RATE_SCALE=40` was used earlier to make a geometry distribution readable; that
# inflates SV density forty-fold, so every recall figure from such a run characterises
# the mechanism rather than caller performance on realistic input (ACCESS §3.7 says so).
# At `SV_RATE_SCALE=1.0`, GRCh38 yields ~25 translocations per run — so ~8 replicates.
#
# Pooling is TP/FN/FP summed across replicates, NOT a mean of per-rep recalls: a mean of
# ratios weights a rep with 3 events the same as one with 40.
#
# Usage:
#   scripts/delta/aggregate_sv_reps.sh --expect 8 $SCRATCH/sv_12345_*
#   scripts/delta/aggregate_sv_reps.sh --expect 8 --dirs-from manifest.txt
set -euo pipefail

EXPECT=0
DIRS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --expect)    EXPECT="$2"; shift 2 ;;
        --dirs-from) mapfile -t DIRS < "$2"; shift 2 ;;
        -h|--help)   sed -n '2,16p' "$0"; exit 0 ;;
        *)           DIRS+=("$1"); shift ;;
    esac
done

if [[ ${#DIRS[@]} -eq 0 ]]; then
    echo "ERROR: no replicate directories given" >&2; exit 2
fi

# ── Coverage of the aggregator's own inputs ─────────────────────────────────────
# A pooling script that silently averages 3 of 8 replicates is precisely the failure
# this repo keeps finding: a metric reported over an unstated denominator. Every
# replicate is accounted for by name before any number is printed.
present=(); missing=()
for d in "${DIRS[@]}"; do
    if compgen -G "$d/truvari_*/summary.json" > /dev/null; then present+=("$d"); else missing+=("$d"); fi
done

printf 'Replicates: %d given, %d with truvari output, %d without\n' \
    "${#DIRS[@]}" "${#present[@]}" "${#missing[@]}"
for d in "${missing[@]}"; do echo "  NO OUTPUT: $d" >&2; done

if [[ "$EXPECT" -gt 0 && "${#present[@]}" -ne "$EXPECT" ]]; then
    echo "ERROR: expected $EXPECT replicate(s) with output, found ${#present[@]}." >&2
    echo "  Pooling a partial set would report a real-looking number over the wrong" >&2
    echo "  denominator. Fix or re-run the missing replicates, or lower --expect" >&2
    echo "  deliberately if some are known-dead." >&2
    exit 3
fi
[[ "${#present[@]}" -gt 0 ]] || { echo "ERROR: no replicate produced truvari output" >&2; exit 3; }

# ── Pool ────────────────────────────────────────────────────────────────────────
python3 - "${present[@]}" <<'PY'
import glob, json, os, sys
from collections import defaultdict

dirs = sys.argv[1:]
acc = defaultdict(lambda: {"TP": 0, "FN": 0, "FP": 0, "reps": 0, "seen_in": []})

for d in dirs:
    for path in sorted(glob.glob(os.path.join(d, "truvari_*", "summary.json"))):
        label = os.path.basename(os.path.dirname(path))[len("truvari_"):]
        try:
            s = json.load(open(path))
        except Exception as e:                      # a corrupt summary is a hard stop,
            print(f"ERROR: unreadable {path}: {e}", file=sys.stderr)  # not a skipped rep
            sys.exit(4)
        a = acc[label]
        a["TP"] += s.get("TP-base") or 0
        a["FN"] += s.get("FN") or 0
        a["FP"] += s.get("FP") or 0
        a["reps"] += 1
        a["seen_in"].append(os.path.basename(d))

n = len(dirs)
print(f"\n{'scorer':<22}{'reps':>5}{'TP':>8}{'FN':>8}{'FP':>8}{'recall':>9}{'prec':>8}{'f1':>8}")
print("-" * 76)
uneven = []
for label in sorted(acc):
    a = acc[label]
    tp, fn, fp = a["TP"], a["FN"], a["FP"]
    rec = tp / (tp + fn) if tp + fn else float("nan")
    pre = tp / (tp + fp) if tp + fp else float("nan")
    f1 = 2 * rec * pre / (rec + pre) if rec + pre else float("nan")
    print(f"{label:<22}{a['reps']:>5}{tp:>8}{fn:>8}{fp:>8}{rec:>9.3f}{pre:>8.3f}{f1:>8.3f}")
    if a["reps"] != n:
        uneven.append(f"{label}: present in {a['reps']} of {n} replicate(s)")

# A scorer that ran in only some replicates is pooled over a different denominator
# than its neighbours in the same table. Say so rather than let the rows be compared.
if uneven:
    print("\nWARNING — these scorers did not run in every replicate, so their rows are")
    print("not directly comparable with the others:")
    for u in uneven:
        print(f"  {u}")

print(f"\nPooled over {n} replicate(s). Recall/precision are computed from SUMMED")
print("TP/FN/FP, not averaged per-rep ratios — a mean of ratios would weight a")
print("replicate with 3 events the same as one with 40.")
PY