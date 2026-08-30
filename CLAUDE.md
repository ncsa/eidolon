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

## Vetting standard (standing requirement)
"It ran and produced output" is not evidence of correctness, and it is the bar this repo
keeps accidentally settling for. **Nothing is "done" until there is evidence it works as
intended; until then it is in progress** — and that governs how work is *reported*, not
just how it is tested. A merged PR is not evidence. A tagged release is not evidence.
Passing CI is evidence about the tests, which is only as good as the tests. Say what
level of evidence exists and what is missing: *"vetted on a known-answer fixture; not yet
run on real data"* beats a checkmark.

Each rule below was earned by a defect that shipped green. Case histories live in
`docs/claude_engineering_audit.md` — **add new ones there, not here.**

1. **Assert content, not existence.** A file existing, a non-zero count, or exit 0 is
   necessary and not sufficient. The `nsom -gt 0` guard passed while every value was the
   malformed string `AF=AF=0.3000`. If a wrong value would still pass, the test is
   decoration.
2. **Prove non-vacuity by mutation.** Break the code a test covers and watch it fail; if
   it still passes it is not a test. Free for a fix (revert it), deliberate for a
   feature. A coverage claim with no mutation experiment behind it is an opinion.
   **Verify the mutation was applied.** A surviving mutant and an unapplied edit produce
   identical output, and the failure is silent: a `sed` or `str.replace` whose pattern does not
   match changes nothing and says nothing. Assert the file changed, or diff it. The tell is a
   "survivor" whose numbers match the baseline *exactly* — thirteen decimal places of agreement
   is not tolerance, it is the same code running twice.
3. **Vet the premise, not just the implementation.** `bnd_proximity.py` was correct code
   that should not have existed — built on an inherited "truvari cannot benchmark
   breakends", true of v4 and reversed in v5.0.0. Check a third-party tool's actual
   version before working around its supposed limits, and **grep for every instance of a
   claim**: a stale copy of that same belief survived one retraction by a day.
4. **Check coverage of the inputs, not just the metric.** Report `n_scored` vs
   `n_planted` per stratum; a metric over an unknown denominator is not a result. A zero
   or unexpectedly-small denominator is a **hard failure, never a `WARNING`**. If a step
   drops data deliberately (filters, LoD, min-depth), log how much.
5. **Chase the evidence past the first plausible story.** `BND recall=0.000` drew three
   confident explanations before the real cause. Each was plausible enough to stop at.
6. **Say what was NOT verified** — checked by hand rather than CI, on a fixture rather
   than real data. Current example: `scripts/delta/tests/` covers 2 of
   `sv_pipeline.sbatch`'s 14 functions; `score_caller` and `check_denominator` produce
   every recall figure in ACCESS §3.5–3.7 and are untested, as is `sbs96_compare.py` (#466).

**The recurring shape**, every quiet failure so far — a harness reporting a metric
without asserting it measured everything it planted:

| reported | true |
|---|---|
| `VERDICT: PASS`, bias in range (#450) | 160 of 567 planted sites silently excluded — the whole lowest-VAF cluster |
| `BND recall=0.000` (#451) | truth emitted unpaired/MATEID-less, unmatchable by construction; the reads were fine |
| `nsom > 0` guard passed | values were `AF=AF=0.3000`; the guard counted **records**, not content |

**Not a bash problem.** #450 and #451 were `bcftools` genotype semantics and Rust record
emission. A Rust harness that never asks "did I score all my planted sites?" fails just
as silently. Bash contributes fragility (footguns below), not blind spots.

### The bar is HIGHER for features than fixes
A fix has a known-bad baseline; a feature has none, so the negative case must be built
deliberately.

- **State the correctness criterion before implementing**, and how it will be falsified.
  If that cannot be written down, the feature is not ready to build.
- **Include a known-answer fixture** — correct output computable independently of the
  code under test. A 7-alt-in-100-reads BAM has VAF 0.070 whatever the implementation thinks.
- **Include a case where it must NOT fire.** Most defects were things matching or
  counting when they should not have.
- **Verify the whole chain.** #405's subclonal VAF passed unit tests and shipped while
  the harness validating it excluded the very sites it existed to test.

### Test adequacy (not test presence)
**New code ships with tests that would catch it being wrong. No exemption for "small",
"obvious", or "hard to test".** Hard to test usually means the seam is wrong — narrow the
inputs until it is testable (`get_bnd_pieces` took a whole `ContigContext` but used only
contig lengths).

- **Test the path that can break, not the path that works.** `bnd_fastq.rs` "covered" BND
  generation via the *input-VCF* path while the *de novo* path shipped a truth VCF
  contradicting its own reads from v1.13.1 to v3.1.0.
- **A function that makes a decision needs a test of that decision.** `get_bnd_pieces`
  chooses which piece is reverse-complemented — the entire semantics of a breakend — and
  had zero tests.
- **Invariants spanning two components need their own test.** Neither side was wrong:
  `sv_model.rs` emitted a correct ALT, `runner.rs` honoured its flags. They disagreed and
  nothing asserted they must agree. Use the *same* helper both sides use so they cannot drift.

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

## Languages: Rust, bash, and (reluctantly) Python
- **The shipped artifact is pure Rust** and must stay that way: the binary invokes no
  interpreter, `Cargo.lock` has no `pyo3`/`cpython`, and `conda-recipe/meta.yaml`
  declares **no runtime requirements** (build-only). Verify before adding anything that
  would change that.
- **Do not introduce new Python.** Product logic is Rust; harness orchestration is bash.
  Python is acceptable only where an **external tool forces it** — We will keep existing parsers for truvari and
  SigProfiler are Python packages, but avoid adding more.
- **Vet what exists** (all validation/prep only, none shipped):
  `scripts/delta/sbs96_compare.py` is a per-validation measurement helper. It parses
  SigProfiler's output in SigProfiler's own env, which is the one justification for
  Python here, and it is **not covered by CI** (#466). `scn_af_compare.py` is gone —
  ported to `eidolon compare-af`, which needed no external tool and so should never have
  been Python; the port is pinned by golden fixtures reproducing the Python's output byte
  for byte plus known-answer tests. `tools/{build_pcawg_sv_vcf,normalize_pcawg_sv_model,graft_sv_model}.py` are
  offline corpus prep. `tools/inject_cancer_sv_model.py` is **dead** — deprecated in
  v1.14.0, no live caller, safe to delete.

## Delta / HPC (`scripts/delta/`)
- Real-data validation runs on **NCSA Delta** (SLURM, account `bhrd-delta-cpu`). The
  cluster filesystem is **not reachable from this workstation** — the user runs jobs
  and pastes results. Artifacts archive to **one** root:
  `/projects/bhrd/jallen17/eidolon-access-results/`. `rneat-access-results/` no longer
  exists — the pre-v2.0.0 runs were merged in when the project rename settled (confirmed
  gone 2026-08-02), so **do not check for it**; a search there wastes a round trip and
  its absence is not evidence a run is missing.

- **Building on Delta: link with `gcc`, never the Cray `cc` wrapper.** Delta's site-wide
  `default` module set loads `PrgEnv-gnu` alongside **`craype-accel-nvidia80`**, which points
  the Cray driver at NVIDIA A100s. `cc` then builds its link line from pkg-config's
  `virtual:world`, dragging in a CUDA runtime plus `cray-mpich`, `cray-libsci`, `cray-dsmml`
  and `libfabric`. **eidolon uses none of them** — it is pure Rust plus zlib-ng/libdeflate via
  cmake. After the RHEL/PE upgrade that world named `cray-sdk-cudatoolkit-25.3_11.8`, whose
  `.pc` file is gone, and `cc` refused to link anything at all — a hello-world included:
  ```
  Package 'cray-sdk-cudatoolkit-25.3_11.8', required by 'virtual:world', not found
  ```
  Nobody loaded those modules; they come from `default`, which is why this appeared overnight
  with no local change (2026-08-12). Build with:
  ```bash
  CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc \
    CARGO_TARGET_DIR=$SCRATCH/cargo-target/eidolon cargo build --release
  ```
  `setup.sh` now exports this. **Confirmed working 2026-08-12.** Note rustc already links with
  its own bundled `lld`; `cc` only ever supplied a driver front-end.
  What does NOT work: **`module unload cudatoolkit/25.3_11.8` alone — measured, still fails.**
  It removes the package while leaving `craype-accel-nvidia80`, the module that demands it.
  Unloading `craype-accel-nvidia80` instead is untested; it should work by the same reasoning,
  but it would not persist anyway since `default` reloads at every login.
  **Symptom to recognise:** every dependency *build script* fails to link on a fresh
  `CARGO_TARGET_DIR`, while a cached target dir gets all the way to the final binary and fails
  there. Same cause; the cache only changes where it surfaces.
  **`CARGO_TARGET_DIR` must be `$SCRATCH/cargo-target/eidolon`** — the pipeline reads
  `$CARGO_TARGET_DIR/release/eidolon` and nothing else. Building into any other path (e.g. the
  pre-rename `cargo-target/rneat-*`) leaves the old binary in place; #514's provenance guard
  now catches that at job start rather than letting it produce stale numbers.

- **The Delta checkout is a separate working copy and goes stale silently.** A queued job
  runs the script as it was when the job *started*, which on a busy queue can be days after
  it was submitted — so a result can predate a fix you believe is in it. Establish the SHA
  **before** interpreting any number:
  `cd /projects/bhrd/jallen17/eidolon && git log -1 --format='%h %ci' -- <script>`, then
  `git merge-base --is-ancestor <sha> origin/develop`. Job 20853511 was diagnosed for two
  findings before the checkout turned out to be three days behind (`6e4c288`, predating
  #508), which invalidated the run regardless. Ask for the SHA in the same round trip as
  the results, not after.

- **Disk is the binding constraint on Delta, and quota views lie about it.** `df` inside a
  project-quota'd directory reports the *quota*, not the mount: `df -h /work/nvme/bhrd`
  showed `500G/500G/0` while `df -h /work/nvme` showed 8.5P at 56%. The `quota` table and
  `df` also disagreed for `/scratch` (1000G vs a 500G mount). What actually stops a job is
  the **project quota** — check it with
  `lfs quota -h -p <projid> /work/nvme`, and read `Disk quota exceeded` literally.
  Measured footprint for `sv_pipeline.sbatch` on GRCh38 at 30x: **~214 GB peak per
  replicate** — the pipeline's disk gate scales this by genome size and coverage, so chr22 at
  30x asks for the 5 GB floor rather than 214 (a constant would have refused every smoke run
  once `/scratch` passed 336 GB used). Itemised from job 20884022 — FASTQ 116 GB (merged 29+29, tumor 17+17,
  normal 12+12) **plus** BAMs 98 GB (the figure `prune_bams` itself reported across
  campaign 20925151), which coexist because pruning only runs at the end. An earlier
  "~113 GB" figure here counted the FASTQ only and was wrong; every capacity estimate built on it was ~half of reality. With a ~171 GB fixed
  baseline (`neat_data` + bwa-mem2 index) **one replicate at a time is all that fits** in a
  500 GB quota (171 + 214 = 385 GB), so serialize with `%1` — `%2` cannot work at any array
  size.
  **A replicate that FAILS keeps all ~214 GB**: its FASTQ prune never runs. Job 20884022 died
  incomplete, held 203 GB indefinitely, and starved array 20904141 — all five tasks failed,
  task 4 spending 8 h 45 m producing 8 KB against a full filesystem. Clear a failed
  replicate's directory before the next campaign, and clear finished ones down to
  `truvari_*/summary.json`, which is all `aggregate_sv_reps.sh` reads.
- **Smoke first, and never learn a mechanism from a 16-hour run.** `sv_pipeline.sbatch`
  defaults `REFERENCE` to `$SCRATCH/neat_data/chr22.fa` — use that default with
  `SV_RATE_SCALE=30` for a fast run that exercises every stage. Job 21025737 did the whole
  pipeline in **4 minutes and 25 core-hours**, against 1215–1445 core-hours for GRCh38.
  **Caveat, measured on that run: chr22 alone cannot plant BND at all.** De novo BND is
  inter-chromosomal only (P3a), so a single-contig reference yields `truth BND: 0` and a
  COVERAGE HOLE for BND by construction. Use a ≥2-contig reference to smoke the BND path. Only then spend GRCh38 at `SV_RATE_SCALE=1.0`, which is for *numbers*,
  not for finding out whether the machinery works. Faster still: most read-level questions are
  answerable on the H1N1 fixture in **seconds** via `eidolon/tests/sv_support_matrix.rs` — that
  is where #516 was finally caught, after being chased through three multi-hour campaigns.
  Ask "can a local test answer this?" before submitting anything long.

- **But H1N1 answers "does the machinery work", never "what is the number".** It is 13.5 kb
  across EIGHT contigs, longest 2280 bp — a shape no real genome has, and four separate
  defects in one day (2026-08-26/27) came from measuring on it and believing the result:
  - **A guard that cannot fire at scale fired constantly.** `generate_fragments` refuses a
    region smaller than `read_length + max_del_len * 2`. With a 500 bp deletion that is
    1150 bp, so splitting a 2280 bp contig left both halves under it and **every read on the
    contig vanished** — silently, via `debug!`. On a real chromosome that guard never fires.
  - **A 2280 bp contig saturates.** At `sv_rate_scale=40` it packs ~33 SVs of up to 570 bp,
    so overlap rejection dominates and #603's budget change was invisible (33.12 → 29.62,
    inside noise). It only became measurable at an unsaturated rate.
  - **Small windows are almost all noise.** A 201 bp window at 30x holds ~60 independent
    reads: sigma ~13%. #499's "+8% over-delivery" was measured on a 1200 bp event, squarely
    in that regime, and does not reproduce — a 10 kb event on a 1 Mb reference gives 0.98
    at every multiplier with sigma ~1.9%.
  - **Eight contigs in 13.5 kb make cross-contig contamination trivial.** Painting every
    contig's reads into one depth array read 415x against a configured 60x, and reported a
    working deletion as 0.98 of control. Hit twice in one day, in two different tests.

  **So H1N1 is a smoke fixture, and for SV work it is ONLY that.** Its longest contig is
  shorter than the guards, comparable to the fragments, and saturates at any interesting SV
  rate, so an SV number measured there is not a weak measurement — it is not a measurement.

  **Use `eidolon/test_data/references/ecoli.fa` for SV tests.** It is already in the repo,
  4.6 Mb in one contig, and no SV test currently uses it. Measured on the #590 deletion
  fixture, same assertion, 17 s per run at 30x:

  | | H1N1 (2280 bp contig) | ecoli (4.6 Mb contig) |
  |---|---|---|
  | deleted-span ratio | 0.00 | 0.000 |
  | largest `D` op | 500 | 500 |
  | must-not-fire tolerance needed | **±20%** | **1.0017** |

  That last row is the whole argument: on H1N1 the "coverage outside the event is unchanged"
  guard needed a 20% tolerance to survive window noise, which is loose enough to miss a real
  regression. On ecoli the same guard reads 0.17%. The single contig means BND still needs
  H1N1 (>= 2 contigs) or Delta — that is the one thing H1N1 is genuinely better at.

  Anything that will be **quoted** still needs Delta (#607). A local pass says "the mechanism
  works"; it is not the number.

- **Evaluate the BAMs before pruning them.** `prune_bams` reclaims ~98 GB per replicate, and
  for three campaigns the pipeline produced BAMs, scored VCF-against-VCF, and deleted them
  without ever asking whether the reads contained what the truth VCF declared. That is exactly
  how #516 survived. `verify_planted_ins` now probes the reads for a 30-mer from the **middle**
  of each planted insertion before anything is deleted (the middle, not the head — a partially
  realized insertion has a perfectly good head). Any new SV type added to the pipeline wants
  the same treatment: **a caller recall of 0 is uninterpretable until you know the evidence was
  there to find.** `PRUNE_BAM=0` keeps the BAMs when a run is specifically diagnostic.

- Staging: `fetch_validation_data.sh` (references), `stage_soy.sh` (align + call a
  self-consistent ref/BAM/VCF; `FULL_GENOME=1` for the whole-genome stress vs the
  fast single-chromosome default). `model_builders.sbatch` exercises the builders.

## Where the evidence lives (so this file can stay a rulebook)
- `docs/claude_engineering_audit.md` — defect taxonomy, timeline, and the case studies
  behind the rules above. **Add new war stories there, not here.**
- `docs/sv_support_matrix.md` + `eidolon/tests/sv_support_matrix.rs` — what eidolon can
  reproduce from `input_vcf` vs generate de novo, measured, with the broken cells pinned.
- `docs/pcawg_sv_measurement.md` — the PCAWG corpus measurements the SV model rests on.

## Testing (Rust)
- `cargo test` builds a **fresh** binary (`assert_cmd::cargo_bin`) — no staleness worry.
  Integration tests live in `eidolon/tests/` (`mod common;` for shared helpers).
- `model_parity.rs` pins builder output byte-for-byte on the H1N1 fixture
  (`BLESS_BASELINES=1` to regenerate after an intentional change).
- **Model fidelity** (built model file actually shapes gen-reads output) is covered by
  `model_fragment_fidelity.rs` and `model_output_fidelity.rs`; `docs/model_builder_baseline.md`
  records the Delta resource envelope and the fidelity status.
- **The toolchain is pinned** in `rust-toolchain.toml`, so `cargo clippy -D warnings` here
  gives the same verdict CI does. It was not always so: three PRs passed clippy locally and
  failed on the runner (#618 `useless_borrows_in_formatting`, 23 sites; #621; #626
  `byte_char_slices`) purely because the workstation was four minor versions behind. If you
  see a lint you cannot reproduce, check `rustup show` before believing it is flaky.
  Bump the pin deliberately and in its own PR — new lints are a feature, not an ambush.
  **`dtolnay/rust-toolchain@stable` sets `RUSTUP_TOOLCHAIN` and silently overrides the pin.**
  That broke the v3.2.0 release build: the action installed `targets:` into *stable* while
  `rust-toolchain.toml` pointed cargo at the pinned version, which had no cross target — so
  `x86_64-apple-darwin` died with "can't find crate for `core`" and the release shipped four
  binaries instead of five. Only genuine cross-compiles are affected; native targets do not
  need `target add` at all, which is why the other four passed. `rust_binaries.yml` now reads
  the channel out of `rust-toolchain.toml` and installs that toolchain with its target.
  It is **tag-only**, so it gets no PR-time signal — use its `workflow_dispatch` to test any
  change to it before believing it.

## Git / GitHub mechanics
- **`gh pr edit` failing is a STALE CLIENT, not this repo.** Fixed by upgrading gh; verified
  working on 2.98.0 (2026-08-30). It has nothing to do with any project board.
  gh <= 2.45 requested `projectCards` unconditionally when fetching a PR, GitHub retired
  Projects classic for everyone, and gh treated that GraphQL deprecation error as fatal —
  so it aborted before applying the edit. Current gh gates the field behind feature
  detection (`ProjectsV1()` returns `Unsupported` for any non-Enterprise host).
  **Ubuntu/Pop 24.04 ships 2.45.0 and will never move**: LTS freezes upstream versions and
  backports only security fixes, which is what `2.45.0-1ubuntu0.3` means. That model suits
  libraries, not clients of a remote API that changes underneath them. Install gh from
  `cli.github.com/packages`, not the distro archive.
  Fallback if stuck on an old gh: `gh api -X PATCH repos/ncsa/eidolon/pulls/<N> -f body=…`.
- **`gh pr view --json commits` has served a STALE commit list** — it showed only the first
  commit of a branch while `gh api repos/ncsa/eidolon/pulls/<N>/commits` was correct. When
  checking whether a push made it into a PR, trust the API, not `gh pr view`. (This compounds
  the missed-merge hazard above: both the "did it land" checks can lie in the same direction.)
- End commit messages with `Co-Authored-By: Claude <noreply@anthropic.com>`.

## GitNexus block (below)
- Auto-generated code-intelligence MCP guidance. Genuinely useful for impact analysis
  when editing Rust **symbols**; treat its "MUST" language as scoped to symbol edits —
  it does not apply to shell scripts, docs, or test-only changes.

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **eidolon** (3759 symbols, 9214 relationships, 300 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

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
