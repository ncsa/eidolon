# eidolon — working practices for agents

Hand-written guidance below; the GitNexus block that follows is auto-generated
(regenerated between its `gitnexus` markers — keep edits **above** the markers).

## Before you touch anything
- **Sync first.** `git fetch --all` and confirm the branch is current before working.
  `develop` is the integration branch (PRs target `develop`, not `main`) and has
  silently diverged from `main` before — don't build on a stale checkout.
- **After a PR merges, verify your commits actually landed.** Late pushes to an open
  PR have missed the merge more than once. Check with
  `git merge-base --is-ancestor <sha> origin/develop`; recover a missed commit with
  `git cherry-pick`.

## Don't trust a green result — read the artifact
- Harness "overall: PASS" / summary lines can be **false passes** — e.g. a run that
  found no inputs still printed PASS over an empty `summary.tsv` (job 19887446). Open
  the actual output (per-item rows, counts, sizes) before believing a success.

### The recurring failure mode: harnesses that don't check their own coverage
Every quiet failure found so far has the same shape — the harness reports a *metric*
without asserting it actually **measured everything it planted**. Verified instances:

| what it reported | what was true |
|---|---|
| `VERDICT: PASS`, bias/MAE in range (#450) | 160 of 567 planted sites silently excluded — the whole lowest-VAF cluster, i.e. the case the harness exists to test |
| `BND recall=0.000` (#451) | truth was emitted unpaired/MATEID-less, so BND was unmatchable by construction; the reads were fine |
| `nsom > 0` abort guard passed | the VAF values were the malformed string `AF=AF=0.3000`; the guard counted **records**, not content |

**Rules when writing or reviewing a harness:**
- Assert on **coverage of its own inputs**, not just on the metric: `n_scored` vs
  `n_planted`, and report the shortfall per stratum. A metric over an unknown
  denominator is not a result.
- A zero or unexpectedly-small denominator is a **hard failure**, never a `WARNING`.
  `cancer_pipeline.sbatch` only warned on an empty somatic truth, so it would have
  scored every caller against nothing and printed results.
- Guards must check **content**, not counts. Record counts pass while values are garbage.
- If a step drops data deliberately (filters, LoD, min-depth), log how much it dropped.

**This is not purely a bash problem.** The two deepest defects (#450, #451) were
`bcftools` genotype semantics and Rust record emission — a Rust harness that likewise
never asked "did I score all my planted sites?" would have failed just as silently. Bash
contributed the *fragility* (see the footguns below), not the blind spots.

### Bash footguns actually hit in this repo
- `set -euo pipefail` + `zcat … | head` makes `zcat` take SIGPIPE (exit 141) and `set -e`
  aborts a step that succeeded. Wrap in `set +o pipefail` … `set -o pipefail`.
- **`$(cmd; echo $?)` does not capture a failing status** — the inherited `set -e` aborts
  the subshell before `echo` runs, yielding an empty string. Use `cmd || rc=$?`.
- **A process substitution's failure is invisible to `set -e`** — `while read … done <
  <(cmd)` just runs zero iterations if `cmd` dies, so downstream logic silently sees "no
  results". Assign to a variable first, or derive the data from something already captured.
- Don't `case`/`if` on only the statuses you expect: handle `*)` explicitly, or an
  unexpected rc gets treated as the happy path (`cancer_pipeline.sbatch` read every rc
  except 10 as "already current").
- Hardcoded offsets over token names (`substr($8,RSTART+9,…)`, `line[17..]`) break the
  moment a token is renamed. Anchor on the token (`sub(/^EIDOLON_VAF=/,"",v)`) or use
  `PREFIX.len()`.

## Delta / HPC (`scripts/delta/`)
- Real-data validation runs on **NCSA Delta** (SLURM, account `bhrd-delta-cpu`). The
  cluster filesystem is **not reachable from this workstation** — the user runs jobs
  and pastes results. Artifacts archive to **two** roots — the project was renamed
  partway through, so Phase 1 / pre-v2.0.0 runs are under
  `/projects/bhrd/jallen17/rneat-access-results/` (germline, cancer, sv, benchmark,
  baseline, tune, threadscale, ppn, modelbuild — the corpus behind the ACCESS report's
  §3.1–3.11) and v2.0.0-onward runs under
  `/projects/bhrd/jallen17/eidolon-access-results/`. `lib_report.sh` writes only to the
  latter (`RESULTS_DIR`, overridable), and `collect_report.sh` scans one root at a time —
  so pass `RESULTS_DIR=` explicitly when looking for historical runs, or you will
  conclude they are missing.
- Staging: `fetch_validation_data.sh` (references), `stage_soy.sh` (align + call a
  self-consistent ref/BAM/VCF; `FULL_GENOME=1` for the whole-genome stress vs the
  fast single-chromosome default). `model_builders.sbatch` exercises the builders.
- **Shell footgun:** `set -euo pipefail` + `zcat … | head` makes `zcat` take SIGPIPE
  (exit 141) when `head` closes the pipe early, and `set -e` then aborts a step that
  actually succeeded. Wrap head-truncated pipelines in `set +o pipefail` … `set -o pipefail`.

## Testing (Rust)
- `cargo test` builds a **fresh** binary (`assert_cmd::cargo_bin`) — no staleness worry.
  Integration tests live in `eidolon/tests/` (`mod common;` for shared helpers).
- `model_parity.rs` pins builder output byte-for-byte on the H1N1 fixture
  (`BLESS_BASELINES=1` to regenerate after an intentional change).
- **Model fidelity** (built model file actually shapes gen-reads output) is covered by
  `model_fragment_fidelity.rs` and `model_output_fidelity.rs`; `docs/model_builder_baseline.md`
  records the Delta resource envelope and the fidelity status.

## Git / GitHub mechanics
- **`gh pr edit` is broken on this repo** (deprecated Projects-classic GraphQL). Patch a
  PR with `gh api -X PATCH repos/ncsa/eidolon/pulls/<N> -f body=…`. `gh issue edit` works.
- End commit messages with `Co-Authored-By: Claude <noreply@anthropic.com>`.

## GitNexus block (below)
- Auto-generated code-intelligence MCP guidance. Genuinely useful for impact analysis
  when editing Rust **symbols**; treat its "MUST" language as scoped to symbol edits —
  it does not apply to shell scripts, docs, or test-only changes.

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **eidolon** (2845 symbols, 7093 relationships, 241 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

> If any GitNexus tool warns the index is stale, run `npx gitnexus analyze` in terminal first.

## Always Do

- **MUST run impact analysis before editing any symbol.** Before modifying a function, class, or method, run `gitnexus_impact({target: "symbolName", direction: "upstream"})` and report the blast radius (direct callers, affected processes, risk level) to the user.
- **MUST run `gitnexus_detect_changes()` before committing** to verify your changes only affect expected symbols and execution flows.
- **MUST warn the user** if impact analysis returns HIGH or CRITICAL risk before proceeding with edits.
- When exploring unfamiliar code, use `gitnexus_query({query: "concept"})` to find execution flows instead of grepping. It returns process-grouped results ranked by relevance.
- When you need full context on a specific symbol — callers, callees, which execution flows it participates in — use `gitnexus_context({name: "symbolName"})`.

## Never Do

- NEVER edit a function, class, or method without first running `gitnexus_impact` on it.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis.
- NEVER rename symbols with find-and-replace — use `gitnexus_rename` which understands the call graph.
- NEVER commit changes without running `gitnexus_detect_changes()` to check affected scope.

## Resources

| Resource | Use for |
|----------|---------|
| `gitnexus://repo/eidolon/context` | Codebase overview, check index freshness |
| `gitnexus://repo/eidolon/clusters` | All functional areas |
| `gitnexus://repo/eidolon/processes` | All execution flows |
| `gitnexus://repo/eidolon/process/{name}` | Step-by-step execution trace |

## CLI

| Task | Read this skill file |
|------|---------------------|
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus/gitnexus-cli/SKILL.md` |

<!-- gitnexus:end -->