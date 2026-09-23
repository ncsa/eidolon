# eidolon re-validation / regression protocol

When a feature changes eidolon (e.g. adapter modeling, #125), this protocol confirms
the change does not:

- **(M) blow up memory** — peak RSS regression,
- **(S) slow down** — wall-clock regression,
- **(D) emit errant data** — fidelity or determinism regression.

It reuses the Delta validation harnesses already built, and gates a candidate
build against a frozen baseline with **statistically-grounded tolerances**.

## Why this is now possible: the replication study gives the noise floor

The Phase-2 replication runs measured the *run-to-run* standard deviation of every
metric (see [`replication_audit.md`](https://github.com/ncsa/eidolon/blob/develop/docs/archive/replication_audit.md)). That sd is exactly
what separates a **real regression** from normal scatter:

| metric | measured run-to-run sd | gate uses |
|---|---|---|
| SNV recall | ~0.01 | tight band |
| indel recall | ~0.05 | wider band |
| SV per-type recall | ~0.09–0.11 | widest band |
| Ts/Tv | ~0.14 | plausibility band 2.0–2.6 |
| determinism (md5) | 0 (exact) | zero tolerance |
| wall-clock | HPC-noisy → use n=3 median | % band |
| peak RSS | stable | % band + ceiling |

A gate fires when a candidate moves **more than k·sd below baseline** (default
k = 2 ≈ 95 % confidence). Without the reps you could not set these — you'd be
guessing whether a 0.97→0.94 dip is a bug or noise.

## The suite: failure mode → metric → harness → gate

| Mode | Metric | Harness | Gate |
|---|---|---|---|
| **M** | peak RSS | `benchmark.sbatch` (`/usr/bin/time -v`) | ≤ baseline × 1.15 **and** below an absolute ceiling |
| **S** | wall-clock (n=3 median) | `benchmark.sbatch` | median ≤ baseline × 1.20 (absorbs HPC noise) |
| **D** | determinism + thread-invariance | `run_order_independence.sbatch` | **byte-identical (exact)** — any drift fails |
| **D** | germline SNP/indel recall, precision, Ts/Tv | `germline_e2e.sbatch` | recall ≥ base − 2 sd; precision ≥ base − 2 sd; Ts/Tv ∈ [2.0, 2.6] |
| **D** | somatic SNV/indel recall, precision | `cancer_pipeline.sbatch` | ≥ base − 2 sd |
| **D** | SV per-type recall | `sv_pipeline.sbatch` (REPS) | ≥ base − 2 sd (only for SV-touching changes) |

All harnesses already emit the needed numbers via the `collect_*` scripts and
`archive_run` (which also records git SHA + core-hours for provenance).

## Tiers — run the cheapest tier that covers the change

| Tier | Scope | Cost | When |
|---|---|---|---|
| **0 — Smoke** | determinism (exact) + **E. coli** germline fidelity + wall/RSS vs baseline | minutes, 1 small node | **every commit / PR** (cheap enough to be near-CI) |
| **1 — Standard** | chr22 germline + cancer fidelity + perf (n=3) + order-independence | ~1 h | **every feature PR** — the default gate |
| **2 — Deep** | input variety (multi-genome) + SV per-type (REPS=3) + cancer sweeps | hours | **risky/perf/SV features, or pre-release** |

Tier 0 alone catches the big three failure modes grossly (crash, memory blowup,
non-determinism, fidelity collapse). Tiers 1–2 add resolution and breadth.

## Feature-touch matrix — which extra checks a change demands

| If the change touches… | Require |
|---|---|
| read generation / error / quality / **adapters** | Tier 0 + 1, plus read-level checks (below) |
| RNG / threading / sharding / seeding | **order-independence is mandatory** (Tier 1) |
| SV / CNV / cancer mutation model | **Tier 2 SV per-type (REPS=3)** |
| performance / memory / allocation | **benchmark Tier 1 + RSS ceiling mandatory** |
| VCF / truth-set output | fidelity + **truth-VCF determinism (md5)** |

## Worked example — adapter modeling (#125)

Adapters add read-level sequence (3′ adapter readthrough, etc.). Two question sets:

**Regression (must NOT change):**
- Determinism still byte-identical (exact).
- Wall-clock / peak RSS within tolerance — adapter logic must not inflate either.
- Truth VCF unchanged — adapters are a read artifact, not variants (md5 of truth
  VCF vs no-adapter run must match).
- **Post-trim fidelity == baseline:** simulate with adapters → trim (cutadapt /
  Trim Galore) → align → call; SNP/indel recall/precision must match the
  no-adapter baseline within 2 sd. Adapters that survive trimming and corrupt
  calls would show up here.

**New-feature checks (adapter-specific) — built as a Tier-2 seam:**
`run_adapter_validation.sh` (CONTROL=1) runs a four-arm matrix that isolates the
adapter effect from the short-insert coverage effect (`off` rejects short inserts,
so a naïve off-vs-on comparison is confounded):

| arm | inserts | adapter |
|---|---|---|
| `off` | long (short rejected) | none — baseline |
| `short_ctrl` | short (kept) | none — genomic reads (`keep_short_fragments`) |
| `on_raw` | short | readthrough, aligned raw |
| `on_trim` | short | readthrough, fastp-trimmed |

The callability check rests on the **adapter-only contrasts** (`on_trim − short_ctrl`,
`on_raw − on_trim`), not the confounded off-vs-on. `collect_adapter_validation.sh`
with `CANDIDATE_TSV=…` emits `adapter_*` metrics that `regression_gate.sh` checks
against the `adapter_*` rows in `baseline_metrics.tsv` — gating callability (recall
stays high → catches a malformed-FASTQ/collapse regression), readthrough presence
(`on_raw` soft-clip elevated), trim recovery, and the null adapter effect.

```bash
# after the 12 jobs finish:
CANDIDATE_TSV=$RESULTS_DIR/adapter_candidate.tsv \
  bash scripts/delta/collect_adapter_validation.sh $RESULTS_DIR/adapter_validation.manifest
bash scripts/delta/regression_gate.sh scripts/delta/baseline_metrics.tsv \
  $RESULTS_DIR/adapter_candidate.tsv 2
```

**Result (#125):** adapter readthrough is detectable and fully callable — the pure
adapter effect is within replication noise; the only fidelity difference vs the
long-insert baseline is the short-insert coverage effect, present with adapters
entirely absent (`short_ctrl`). **Lesson:** the collapse bug that motivated this
seam (a malformed FASTQ record) was invisible to `zcat`/`wc` read counts and only
surfaced through a real aligner — so the gate keys on *aligned/called* metrics, not
read counts.

**Tier plan for adapters:** Tier 0/1 (default path, adapters off) on every commit /
PR; run the Tier-2 adapter seam above when the change touches the adapter or
fragment-generation code, or pre-release.

## Baseline management

- Freeze baseline metrics on `develop` into a committed `baseline_metrics.tsv`
  (one row per metric: value, sd, n, git SHA, date). Re-freeze on each release.
- A candidate run produces `candidate_metrics.tsv` in the same shape.
- `regression_gate.sh` diffs the two, applies the gates + tolerances above, prints
  a PASS/FAIL table, and **exits non-zero on any FAIL** (CI / `--dependency`-friendly).

### Re-freezing PERF rows without disturbing fidelity (the v1.19.1 logging fix)

A change can be **output-preserving but performance-changing** — the v1.19.1
logging fix (#340) is the canonical case. Every Phase-1/2 run used the old default
log level `trace`, whose per-base debug events fire ~`coverage × read_len ×
ref_bp` times and burned most of the wall-clock on log I/O; so the frozen
`*_wall_s` / `*_rss_mb` rows were captured **under that handicap and are now
stale** (the fixed build is faster). The fix does **not** change output bytes, so
every fidelity / determinism / recall row stays valid — only the perf rows must be
re-measured and re-frozen. Two steps:

1. **Prove the change is output-preserving** (before trusting the split above). Run
   `baseline_capture.sbatch` on the fixed build and compare its `baseline.json`
   `fastq_md5` / `vcf_records_md5` against a pre-fix run's:
   ```bash
   bash scripts/delta/rebaseline_perf.sh --check-identity OLD.json NEW.json
   ```
   If this is not `PASS`, the change altered output and the fidelity rows must be
   re-validated too — it is not a pure perf/logging change.

2. **Re-freeze only the perf rows** from a fresh benchmark on the fixed build:
   ```bash
   MODE=submit TIER=0 LABEL=logfix bash scripts/delta/regression_suite.sh   # TIER=1 also re-freezes chr22_*
   # …when the perf job finishes…
   MODE=collect MANIFEST=$HOME/regression_logfix.manifest \
     bash scripts/delta/regression_suite.sh > $HOME/logfix.candidate.tsv
   bash scripts/delta/rebaseline_perf.sh $HOME/logfix.candidate.tsv          # rewrites *_wall_s / *_rss_mb only
   ```
   `rebaseline_perf.sh` touches **only** rows ending in `_wall_s` / `_rss_mb`
   (preserving their gate + tier), writes a `.bak`, and leaves every other row
   untouched — so a stray fidelity value in the candidate can never overwrite a
   fidelity baseline. Commit the result with the fixed build's git SHA.

## Automation (built)

Three files implement this, thin wrappers over the existing harnesses — no new
simulation logic:

- **`scripts/delta/baseline_metrics.tsv`** — the frozen baseline, *seeded from the
  Phase-1/2 validation + replication numbers* (value, sd, gate, tier). Re-freeze on
  each release by running the collect step on `develop`.
- **`scripts/delta/regression_suite.sh`** — `MODE=submit TIER=n` submits the tier's
  jobs (order-independence + germline + perf, plus chr22 germline/cancer at Tier 1)
  and writes a manifest; `MODE=collect MANIFEST=…` extracts metrics into a candidate
  TSV. (Tier 2 — variety, SV per-type reps, sweeps — is driven by
  `run_input_variety.sh` / `run_cancer_sweeps.sh REPS=3` and appended.)
- **`scripts/delta/regression_gate.sh baseline candidate [max_tier]`** — applies the
  gates, prints a PASS/FAIL table, and exits non-zero on any FAIL (CI-friendly).

**One command (fire-and-forget)** on a feature branch, from repo root:

```bash
MODE=run TIER=1 LABEL=adapters bash scripts/delta/regression_suite.sh
```

This submits the tier's jobs **plus a dependency gate job** (`--dependency=afterany`)
that auto-collects and writes a PASS/FAIL verdict to
`$RESULTS_DIR/regression_<label>.verdict.txt` when they finish — no babysitting.
A crashed harness surfaces as a `FAIL` (non-numeric), not a silent pass.

Manual three-step (if you want to inspect between stages):

```bash
MODE=submit  TIER=1 LABEL=cand bash scripts/delta/regression_suite.sh
MODE=collect MANIFEST=<path>   bash scripts/delta/regression_suite.sh > candidate.tsv
bash scripts/delta/regression_gate.sh scripts/delta/baseline_metrics.tsv candidate.tsv 1
```

There is no automatic trigger (git hook / CI) — GitHub runners can't reach Delta,
so a run is kicked off deliberately (the convention: `MODE=run TIER=1` on the
feature branch before merging; `TIER=2` for SV/perf-touching changes).

The gate logic and all collect parsers are unit-tested; the `sbatch` submit path
runs on Delta. Capturing a fresh baseline is just the same collect step pointed at
a `develop` run.
