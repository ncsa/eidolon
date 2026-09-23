#!/usr/bin/env bash
# Run a regression tier on the CURRENT build and gate it against the baseline.
# Thin orchestrator over the existing harnesses — no new simulation logic.
#
#   MODE=run     TIER=0|1   — ONE COMMAND: submit the tier's jobs + a dependency
#                            "gate" job that auto-collects and writes a PASS/FAIL
#                            verdict when they finish. Fire-and-forget; no babysitting.
#   MODE=submit  TIER=0|1   — submit only; write the manifest (manual collect later).
#   MODE=collect MANIFEST=… — extract metrics from finished jobs → stdout (candidate TSV).
#   MODE=gate    MANIFEST=… — collect + run regression_gate.sh; write the verdict file
#                            (this is what the dependency job runs).
#
# One-command usage on a feature branch (from repo root):
#   MODE=run TIER=1 LABEL=adapters bash scripts/delta/regression_suite.sh
#   # ... jobs run; when the gate job finishes, read the verdict file it prints ...
#
# Tier 2 (variety, SV per-type reps, sweeps) is heavier and seam-specific — drive it
# with run_input_variety.sh / run_cancer_sweeps.sh (REPS=3) + collect_* and append rows.
# See docs/regression_protocol.md.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO_ROOT/scripts/delta/lib_report.sh"   # $SCRATCH + RESULTS_DIR

MODE="${MODE:-run}"
TIER="${TIER:-1}"
LABEL="${LABEL:-cand}"
DATA="${DATA:-$SCRATCH/neat_data}"
ACCT="${ACCESS_ACCOUNT:-bhrd-delta-cpu}"
RD="${RESULTS_DIR:-$HOME}"
MANIFEST="${MANIFEST:-$RD/regression_${LABEL}.manifest}"
CAND="$RD/regression_${LABEL}.candidate.tsv"
VERDICT="$RD/regression_${LABEL}.verdict.txt"
BASELINE="$REPO_ROOT/scripts/delta/baseline_metrics.tsv"
D="$REPO_ROOT/scripts/delta"

# ── submit the tier's harness jobs (writes manifest: kind prefix jobid outdir) ──
do_submit() {
    mkdir -p "$(dirname "$MANIFEST")"
    echo "# kind prefix jobid outdir" > "$MANIFEST"
    sub() { local kind="$1" prefix="$2" out="$3" script="$4"; shift 4
        local jid; jid=$(env "$@" OUTDIR="$out" sbatch --parsable "$D/$script")
        echo "$kind $prefix $jid $out" | tee -a "$MANIFEST"; }
    sub order    -      "$SCRATCH/reg_${LABEL}_order"      run_order_independence.sbatch REFERENCE="$DATA/yeast.fa"
    sub germline ecoli  "$SCRATCH/reg_${LABEL}_germ_ecoli" germline_e2e.sbatch           REFERENCE="$DATA/ecoli.fa" TOOLS=eidolon
    # Perf arm trimmed to EXACTLY what do_collect reads back: ecoli and chr22, eidolon,
    # one thread. The full benchmark also sweeps yeast, GRCh38 and a 5-point thread curve
    # against NEAT 4, none of which reaches the gate — ~4 h of discarded compute that timed
    # out the 6 h job (21766280) and produced no perf rows at all. NEAT_MAX_GENOME_MB=0
    # drops the NEAT arm; SCALING_THREAD_MODES=1 collapses Phase B. Use benchmark.sbatch
    # directly, unoverridden, for the head-to-head publication numbers.
    sub perf     -      "$SCRATCH/reg_${LABEL}_perf"       benchmark.sbatch \
        GENOMES="ecoli:$DATA/ecoli.fa chr22:$DATA/chr22.fa" \
        NEAT_MAX_GENOME_MB=0 SCALING_THREAD_MODES=1
    if [[ "$TIER" -ge 1 ]]; then
        sub germline chr22 "$SCRATCH/reg_${LABEL}_germ_chr22" germline_e2e.sbatch   REFERENCE="$DATA/chr22.fa" TOOLS=eidolon
        sub cancer  somatic "$SCRATCH/reg_${LABEL}_cancer"    cancer_pipeline.sbatch REFERENCE="$DATA/chr22.fa" TOTAL_COVERAGE=30 PURITY=0.6 PRUNE=1
    fi
}

# ── extract metrics from finished jobs → stdout (candidate TSV) ──
do_collect() {
    [[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; return 1; }
    while read -r kind prefix jobid outdir; do
        [[ "$kind" == \#* || -z "${kind:-}" ]] && continue
        case "$kind" in
          order)
            local f="$REPO_ROOT/eidolon-orderindep_${jobid}.out"
            [[ -f "$f" ]] || f="eidolon-orderindep_${jobid}.out"
            # Emits MISSING, never an empty string, when the label is absent. An empty
            # value string-compares to a FAIL on an `exact` gate and is then
            # indistinguishable from a measured failure -- which is how the two label
            # mismatches below went unnoticed while the gate reported them as regressions.
            get(){
                local v; v="$(grep -m1 "$1:" "$f" 2>/dev/null | grep -oE 'PASS|FAIL' | head -1)"
                if [[ -z "$v" ]]; then
                    echo "WARNING: '$1:' not found in $f — collector and harness labels have drifted" >&2
                    v=MISSING
                fi
                printf '%s' "$v"
            }
            # thread_invariant is the harness's A1-vs-B1 comparison (1 thread vs N on the
            # SAME build). It was grepped as `multithread_vs_1thread:.*same`, a label the
            # harness has not printed since d0cd060 renamed it -- so the local `ti=FAIL`
            # default was never overwritten and this metric reported FAIL on every run
            # without reading a measurement. Baseline says PASS, so every run since that
            # redesign showed a false regression here.
            local ti; ti="$(get thread_invariant)"
            if [[ "$ti" == PASS && "$(get determinism_1thread)" != PASS ]]; then
                # Thread invariance is meaningless if the build is not deterministic at a
                # fixed thread count: two identical hashes prove nothing when neither run
                # is reproducible.
                ti=FAIL
            fi
            printf '%s\t%s\n' determinism_thread_invariant "$ti"
            printf '%s\t%s\n' shard_order_independent "$(get shard_order_independent)"
            printf '%s\t%s\n' shard_disjoint          "$(get shard_disjoint)"
            # Metric name kept for baseline pairing; the harness label is
            # contig_name_invariant, which is what it actually measures.
            printf '%s\t%s\n' contig_order_independent "$(get contig_name_invariant)"
            ;;
          germline)
            python3 - "$outdir/eidolon_scored.summary.csv" "$prefix" <<'PY'
import csv, sys
f, pre = sys.argv[1], sys.argv[2]
try: rows = {r["Type"]: r for r in csv.DictReader(open(f)) if r.get("Filter") == "PASS"}
except Exception: rows = {}
def g(t, k):
    try: return f"{float(rows[t][k]):.4f}"
    except Exception: return "NA"
print(f"{pre}_snp_recall\t{g('SNP','METRIC.Recall')}")
print(f"{pre}_snp_precision\t{g('SNP','METRIC.Precision')}")
print(f"{pre}_indel_recall\t{g('INDEL','METRIC.Recall')}")
print(f"{pre}_indel_precision\t{g('INDEL','METRIC.Precision')}")
try: print(f"{pre}_tstv\t{float(rows['SNP']['QUERY.TOTAL.TiTv_ratio']):.3f}")
except Exception: print(f"{pre}_tstv\tNA")
PY
            ;;
          perf)
            local m="$outdir/median.tsv"
            for gg in ecoli chr22; do
              awk -F'\t' -v g="$gg" '$1==g && $4=="eidolon" && $3==1 {print g"_wall_s\t"$5"\n"g"_rss_mb\t"$6}' "$m" 2>/dev/null
            done
            ;;
          cancer)
            python3 - "$outdir/scored.stats.csv" <<'PY'
import csv, sys
f = sys.argv[1]
try: rows = {(r.get("type") or "").strip().lower(): r for r in csv.DictReader(open(f))}
except Exception: rows = {}
def g(types, k):
    for t in types:
        r = rows.get(t)
        if not r: continue
        for c in r:
            if c.strip().lower() == k:
                try: return f"{float(r[c]):.4f}"
                except Exception: return "NA"
    return "NA"
print(f"somatic_snv_recall\t{g(('snvs','snp','snps'),'recall')}")
print(f"somatic_indel_recall\t{g(('indels','indel'),'recall')}")
PY
            ;;
        esac
    done < "$MANIFEST"
}

case "$MODE" in
  submit)
    do_submit
    echo; echo "Submitted TIER=$TIER. Manifest: $MANIFEST"
    echo "When done: MODE=collect MANIFEST=$MANIFEST bash $0 > candidate.tsv ; then regression_gate.sh"
    ;;
  collect)
    do_collect
    ;;
  gate)
    do_collect > "$CAND"
    echo "candidate metrics → $CAND"
    set +e
    bash "$D/regression_gate.sh" "$BASELINE" "$CAND" "$TIER" | tee "$VERDICT"
    rc=${PIPESTATUS[0]}
    set -e 2>/dev/null || true
    echo "verdict ($([ "$rc" -eq 0 ] && echo PASS || echo FAIL)) → $VERDICT"
    exit "$rc"
    ;;
  run)
    do_submit
    deps=$(awk 'NR>1 && $3 ~ /^[0-9]+$/ {print $3}' "$MANIFEST" | paste -sd:)
    [[ -n "$deps" ]] || { echo "no jobs submitted — aborting gate" >&2; exit 1; }
    gjid=$(sbatch --parsable --job-name="reg-gate-${LABEL}" --account="$ACCT" --partition=cpu \
        --time=00:20:00 --mem=4G --output="$RD/reg-gate_${LABEL}_%j.out" \
        --dependency=afterany:"$deps" \
        --wrap="cd '$REPO_ROOT' && MODE=gate TIER=$TIER LABEL='$LABEL' MANIFEST='$MANIFEST' bash scripts/delta/regression_suite.sh")
    echo
    echo "TIER=$TIER submitted; gate job $gjid runs after them (afterany)."
    echo "Verdict will be written to: $VERDICT   (gate log: $RD/reg-gate_${LABEL}_${gjid}.out)"
    echo "Watch: squeue --me"
    ;;
  *) echo "unknown MODE=$MODE (run|submit|collect|gate)" >&2; exit 1 ;;
esac
