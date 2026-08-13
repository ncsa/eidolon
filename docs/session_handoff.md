# Session handoff — 2026-08-12

Written for a fresh agent session picking this up. Read `CLAUDE.md` first; this file is
the delta on top of it.

## Do this before anything else

```bash
git fetch --all --tags
git log --oneline -1 origin/develop
grep -m1 '^version' Cargo.toml     # expect 3.1.0 or later
```

**A stale checkout is the live hazard in this repo, not a hypothetical.** This session
opened on a `develop` that was **231 commits behind** `origin/develop` — v1.20.1 against
v3.1.0 — and the project had been **renamed `rusty-neat` / `rneat` → `eidolon`** in v2.0.0
in the meantime. Every file path, crate name, and binary name in older notes is wrong:
`common/` is now `eidolon-core/`, the binary is `eidolon`, and output tokens are
`EIDOLON_*` rather than `NEAT_*` / `RNEAT_*`. The git remote is still
`github.com:ncsa/rusty-neat.git` and resolves by redirect to `ncsa/eidolon`, so remote URLs
are not a reliable signal of which name is current.

## Open PRs from this session

All target `develop`. None merged as of writing.

| PR | Branch | Contents |
|---|---|---|
| [#523](https://github.com/ncsa/eidolon/pull/523) | `docs/neat-comparison-refresh` | README NEAT-comparison table: corrected, refreshed, expanded. Docs only |
| [#524](https://github.com/ncsa/eidolon/pull/524) | `chore/ignore-junie` | One `.gitignore` line for `/.junie/`. Trivial |
| [#525](https://github.com/ncsa/eidolon/pull/525) | `test/transition-matrix-content` | Real assertions for the SNP transition matrix + a Delta plan. One accessor added, no behavior change |

`#522` (`origin/docs/delta-cray-linker`) predates this session and was not touched.

**Verify merges landed.** Late pushes have missed merges here more than once:
`git merge-base --is-ancestor <sha> origin/develop`. Also note `gh pr view` served a
**stale commit list** during this session — the API (`gh api repos/ncsa/rusty-neat/pulls/N/commits`)
was correct while `gh pr view --json commits` showed only the first commit. Trust the API.

## Where the work stands

### Done and pushed

- **The README's NEAT comparison table is now accurate** (#523). Two claims in it were
  false and are corrected with the evidence recorded inline: NEAT 4.6.1 does **not**
  generate structural variants (the `Inversion` / `Duplication` / `Translocation` /
  `Transposition` / `CopyNumberVariant` classes exist under `neat/variants/` but are
  exported from nothing, referenced nowhere outside their own files, and absent from
  `VariantTypes.types`), and it cannot render one from an input VCF either (symbolic ALTs
  fall to `UnknownVariant`, whose `get_alt()` and `get_ref_len()` both raise). NEAT 2.1 has
  no SV code at all.
- **The transition matrix is now genuinely tested** (#525), with non-vacuity established by
  four mutations rather than asserted. See `docs/transition_matrix_validation_plan.md`.

### Next, in rough priority order

1. **Run the Delta tier for the transition matrix.** Fully specified in
   `docs/transition_matrix_validation_plan.md` — inputs, steps, pass criteria, failure
   modes, and the guards that matter. Nothing has been run; the script
   (`scripts/delta/run_transition_matrix_validation.sh`) does not exist yet. **The user runs
   Delta jobs and pastes results — the cluster filesystem is not reachable from the
   workstation.**
2. **[#526](https://github.com/ncsa/eidolon/issues/526) — evaluate SVE / GATK-SV, secondarily
   LUMPY / CNVnator, as additional SV-caller verification.** From Joao; see the issue for
   why it matters and what to check first. Step one is whether the tools still run at all.
3. **Upstream courtesy:** NEAT 4.6.1's README recommends this tool as `rneat` and links
   `ncsa/rusty-neat` — pre-rename, working only by redirect. Worth a PR to `ncsa/NEAT`.
   Separately, the two `UnknownVariant` defects found above (`__init__` never sets `alt`;
   the call site passes `kwargs=data`, nesting metadata a level deeper than the accessors
   read) are real upstream bugs and were not reported.

## Findings worth not rediscovering

- **A symmetric test fixture cannot catch a transposed matrix.** `ref=ACGT` / `read=CATG`
  produces A→C and C→A in equal number — a count matrix equal to its own transpose. My
  first version of that test passed the transpose mutation. Single-entry rows fail for a
  related reason: normalization puts one nonzero cell at probability 1.0 wherever it sits.
  Asymmetric, multi-target fixtures only.
- **A forced transition matrix does not drive output to 100%.** `indel_probability` is 0.4
  and `insertion_bias` is uniform over ACGT, so ~40% of sequencing errors are indels whose
  inserted bases never consult the matrix. An absolute `>98%` assertion **fails on correct
  code**; the fidelity test is differential against a control run instead (default matrix
  T share 0.198, forced A→T 0.867). Insertions show up as CIGAR `I` rather than mismatches,
  so this floor should *not* contaminate a real-data comparison — confirm rather than assume.
- **`cargo test --package eidolon --lib` fails** — there is no library target in the
  `eidolon` package. Unit tests live in the binary: `--bin eidolon`.
- **eidolon is behind NEAT 4 in two places**, both worth keeping visible rather than
  quietly dropping: default-model Ts/Tv is 2.21 against NEAT's 2.33 (#410, default model
  only), and NEAT ships **named** quality-bin presets (`--quality-preset novaseq`) where
  eidolon requires the bins spelled out in `binned_quality_bins`.
- NEAT source is readable without cloning:
  `gh api -H "Accept: application/vnd.github.raw" repos/ncsa/NEAT/contents/<path>?ref=4.6.1`.
  A shallow clone at a tag is easier for anything more than one file.

## Artifact

A published comparison page, **The NEAT Lineage** —
https://claude.ai/code/artifact/970b8fe1-d7ad-482e-82ca-84f9f332eab4 — covering the
capability table, a subcommand map across all three generations, the v1.21.0 → v3.1.0
release ledger, and an explicit "how to read these numbers" section. To update it, pass
that URL as `url` to the Artifact tool; publishing without it creates a *separate*
artifact rather than updating this one.