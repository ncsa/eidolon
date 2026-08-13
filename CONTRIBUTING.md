# Contributing to eidolon

Thanks for your interest. Please review the [Code of Conduct](CODE_OF_CONDUCT.md) first.

eidolon is a Rust port of [NEAT](https://github.com/ncsa/NEAT), extended with a native
tumour/normal cancer workflow. See the [README](README.md) for what it does, or
[`ONBOARDING.md`](ONBOARDING.md) if you are new to the project and want the tour rather than
the reference.

## Getting set up

```bash
cargo build --release          # binary at target/release/eidolon
cargo test --workspace         # the full suite, ~30s
cargo fmt --all                # before every commit
```

Unit tests live in the binary crate, so `cargo test --package eidolon --lib` **fails** — there
is no library target. Use `--bin eidolon`.

CI gates three things: `cargo fmt --all -- --check`, `cargo build --workspace --all-targets`
with `RUSTFLAGS="-D warnings"` (a disabled test shows up as a rustc warning, and that has
silently switched off a regression guard before), and `cargo test --workspace`.

## Where changes go

**PRs target `develop`, not `main`.** `main` carries release merges and tags; `develop` is the
integration branch. Fetch before you start — `develop` moves quickly.

After a PR merges, **confirm your commits actually landed**:

```bash
git merge-base --is-ancestor <sha> origin/develop
```

Late pushes have missed a merge more than once here. Recover with `git cherry-pick`.

## The testing bar

This is the part worth reading properly. These rules were each earned by a defect that shipped
green, and the case histories are in
[`docs/claude_engineering_audit.md`](docs/claude_engineering_audit.md).

**Assert content, not existence.** A file existing, a non-zero count, or exit 0 is necessary and
not sufficient. A guard once checked `count > 0` while every value was the malformed string
`AF=AF=0.3000`. If a wrong value would still pass, the test is decoration.

**Prove non-vacuity by mutation.** Break the code a test covers and watch it fail. If it still
passes, it is not a test. This is free for a bug fix — you already have the broken state, so
just revert your fix and check. For a new feature it has to be deliberate. A coverage claim with
no mutation experiment behind it is an opinion.

**Test the path that can break.** BND generation was "covered" through the input-VCF path while
the *de novo* path shipped a truth file contradicting its own reads for eight releases.

**Report the denominator.** A metric over an unknown denominator is not a result. If a step
drops data — filters, thresholds, minimum depth — say how much. A zero or unexpectedly small
denominator is a failure, not a warning.

**Say what you did not verify.** "Vetted on a known-answer fixture; not yet run on real data"
beats a checkmark. A merged PR is not evidence, a passing CI run is evidence about the tests,
and *nothing is done until there is evidence it works as intended*.

The bar is **higher for a new feature than for a fix**, because a fix has a known-bad baseline
and a feature has none. A feature wants a known-answer fixture (an output computable
independently of the code under test) and a case where it must *not* fire — most defects found
here were things matching or counting when they should not have.

`CLAUDE.md` has the full version, including the bash footguns this repo has actually hit. It is
written for AI agents but the engineering content applies to everyone.

## Languages

**Product logic is Rust. Harness and orchestration are bash.** The shipped artifact invokes no
interpreter, and the conda recipe declares no runtime requirements.

Please do not add Python. The few Python files that exist are there because an external tool
forces it — SigProfiler and truvari are Python packages — and they are all offline analysis, none
shipped.

## Commits and PRs

Explain **why**, not what — the diff shows what. If a change fixes a defect, say what the defect
produced, because that is what makes the test reviewable. Note explicitly what you did *not*
verify.

Two mechanics worth knowing: `gh pr edit` is broken on this repo (deprecated Projects-classic
GraphQL), so patch a PR body with
`gh api -X PATCH repos/ncsa/eidolon/pulls/<N> -f body=...`. `gh issue edit` works fine.

## Reporting problems

Open a GitHub issue. Bugs: what you ran, what you expected, what happened, and the version
(`eidolon --version`). There is also a **Feedback** issue type for things that are not quite
bugs — a rough edge, a confusing message, or a good experience worth hearing about.

If you would rather talk first, message Joshua Allen (links on the repo page).

## License

BSD 3-Clause — see [LICENSE.md](LICENSE.md). Because some code derives from Biopython, the
Biopython License Agreement is included there as well.
