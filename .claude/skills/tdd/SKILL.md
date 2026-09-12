---
name: tdd
description: "Use when writing new eidolon code or fixing a defect and you want the test to come first — the ordered procedure for red/green, proving a test is non-vacuous by mutation, and building the known-answer + must-not-fire cases a feature needs. Examples: \"add X, test-first\", \"write a test for this bug before fixing it\", \"is this test actually testing anything?\", \"TDD the new SV type\"."
---

# Test-driven development in eidolon

This skill turns the vetting standard in `CLAUDE.md` into an ordered procedure. That
standard says what "done" means; this says what to type, in what order, to get there.

**The one rule underneath all of it:** a test that would still pass if the code were wrong
is decoration. Every step below exists to make that failure mode visible *before* the work
is reported as finished.

Run the whole sequence without being asked to enumerate it. The user should say "add X,
test-first" and get steps 0–5.

---

## Step 0 — Write the correctness criterion down, before any code

One or two sentences, in the test file's `//!` header, stating **what would be true if this
worked and what observation would falsify it**. Not "tests the deletion path" — that is a
topic, not a criterion.

> A large deletion must not remove coverage anywhere except where it is.
> — `eidolon/tests/deletion_region_padding.rs`

If you cannot write that sentence, the feature is not ready to build. Say so and stop there
rather than producing a test that ratifies whatever the code happens to do.

For a **fix**, the criterion is free: the known-bad baseline is the defect. Write down the
observed wrong value.

For a **feature**, there is no baseline, so the negative case has to be built deliberately —
see step 4.

## Step 1 — Red: write the failing test first

Put it in the right place:

| kind | location |
|---|---|
| unit, on a pure function | `#[cfg(test)] mod tests` beside the code |
| anything that runs the binary | `eidolon/tests/<topic>.rs`, with `mod common;` |

Shared helpers already exist — do not re-roll them. `eidolon/tests/common/mod.rs` has
`eidolon()` (an `assert_cmd` `Command`), `GenReadsConfig`, `fresh_workdir()`,
`h1n1_reference()`, `read_gzip_fastq_lines()`, `revcomp()`, `synthetic_insert()`.
`eidolon/tests/common/gate2.rs` has the align-and-analyze path: `generate_reads()`,
`align()`, `analyse()`, `run_gate()`.

**Assert content, not existence.** A file existing, a count above zero, or exit 0 is
necessary and never sufficient. `nsom -gt 0` passed while every value was the malformed
string `AF=AF=0.3000` — the guard counted records, not content. Assert the *value*.

Name the test as the claim it makes, so a failure reads as a sentence:
`a_large_deletion_does_not_silence_the_rest_of_the_contig`, not `test_deletion_2`.

Then **run it and watch it fail, for the reason you expect**:

```bash
cargo test --workspace --no-fail-fast <test_name> -- --nocapture
```

A test that fails to compile, or fails on a missing fixture, has not gone red — it has
gone nowhere. Read the assertion message.

## Step 2 — Green: the smallest change that makes it pass

Then the gates, which are exactly what CI runs:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --verbose --workspace --no-fail-fast
```

`--no-fail-fast` matters: plain `cargo test` stops at the first failing test *binary*, so a
later break stays hidden. The toolchain is pinned in `rust-toolchain.toml` (1.98.0), so a
clippy verdict here is the verdict CI gives. If you see a lint you cannot reproduce, check
`rustup show` before believing it is flaky.

## Step 3 — Prove the test is non-vacuous, by mutation

**This step is not optional and it is not a formality.** Break the code the test covers and
watch the test fail. If it still passes, it is not a test.

For a fix this is free — revert the fix and re-run. For a feature it is deliberate: pick the
single line carrying the decision and invert it (flip a comparison, drop a `+ 1`, return the
other branch, replace a computed value with a constant).

`cargo-mutants` is **not installed** on this workstation, so this is done by hand. Use the
helper, which closes the trap that makes this step lie:

```bash
.claude/skills/tdd/mutate.sh start eidolon-core/src/thing.rs
#   ... now make the mutation with the Edit tool ...
.claude/skills/tdd/mutate.sh run eidolon-core/src/thing.rs <test_name>
```

`run` refuses to proceed unless the file actually differs from the snapshot, runs the test,
reports SURVIVOR or KILLED, and restores the file either way.

**The trap it closes:** a `sed` or `str.replace` whose pattern does not match changes nothing
and says nothing, so an unapplied edit and a surviving mutant produce identical output. The
tell is a "survivor" whose numbers match the baseline *exactly* — thirteen decimal places of
agreement is not tolerance, it is the same code running twice.

Record the result in the PR body: which line was mutated, and that the test failed. A
coverage claim with no mutation experiment behind it is an opinion.

## Step 4 — For a feature, three cases are required, not one

A fix inherits its negative case from the defect. A feature has none, so build all three:

1. **A known-answer fixture** — correct output computable independently of the code under
   test. A 7-alt-in-100-reads BAM has VAF 0.070 whatever the implementation thinks. If the
   expected value was produced *by* the code, it is a snapshot, not a known answer; say so.
2. **A case where it must NOT fire.** Most defects in this repo were things matching or
   counting when they should not have. `deletion_region_padding.rs` measures the untouched
   part of the contig, not just the deleted part.
3. **The denominator.** Report `n_scored` vs `n_planted`, not just the metric. A metric over
   an unknown denominator is not a result, and a zero or unexpectedly-small denominator is a
   **hard failure, never a warning**. #450 reported `VERDICT: PASS` while 160 of 567 planted
   sites — the whole lowest-VAF cluster — were silently excluded.

Also: **a function that makes a decision needs a test of that decision.** `get_bnd_pieces`
chooses which piece is reverse-complemented, which is the entire semantics of a breakend,
and had zero tests. And **an invariant spanning two components needs its own test** — neither
`sv_model.rs` nor `runner.rs` was wrong on its own; they disagreed and nothing asserted they
must agree. Have both sides call the *same* helper so they cannot drift.

If a thing is hard to test, the seam is usually wrong. Narrow the inputs until it is testable
— `get_bnd_pieces` took a whole `ContigContext` and used only contig lengths.

## Step 5 — Choose the fixture deliberately, and say which one you used

| fixture | what it can answer | what it cannot |
|---|---|---|
| `eidolon/test_data/references/ecoli.fa` — 4.6 Mb, one contig | **Default for SV work.** Real window statistics: the "coverage outside the event is unchanged" guard reads 0.17% here | anything inter-chromosomal — single contig |
| `H1N1.fa` — 13.5 kb, 8 contigs, longest 2280 bp | "does the machinery work"; BND, because it has ≥2 contigs | **any number.** Its longest contig is shorter than the guards, saturates at any interesting SV rate, and a 201 bp window at 30x has sigma ~13% |
| Delta (`scripts/delta/`) | anything that will be **quoted** | fast iteration — GRCh38 is 1215–1445 core-hours |

Four separate defects in one day came from measuring on H1N1 and believing the result. The
same must-not-fire guard needed a ±20% tolerance on H1N1 and 1.0017 on ecoli — loose enough,
on H1N1, to miss a real regression. Use ecoli unless you need a second contig.

Ask **"can a local test answer this?"** before submitting anything long. #516 was chased
through three multi-hour Delta campaigns and finally caught in seconds by
`eidolon/tests/sv_support_matrix.rs`.

---

## Reporting: what to say when you are done

Report the level of evidence, not a checkmark. "Vetted on a known-answer fixture; not yet run
on real data" beats "✅ done". A merged PR is not evidence; a tagged release is not evidence;
passing CI is evidence about the tests, which is only as good as the tests.

State explicitly what was **not** verified — checked by hand rather than in CI, on a fixture
rather than real data, one path covered and the sibling path not.

Until there is evidence it works as intended, the work is **in progress**, and that governs
how it is *reported*, not just how it is tested.

## Anti-patterns, each earned by a defect that shipped green

- **Testing the path that works instead of the path that can break.** `bnd_fastq.rs` covered
  BND generation via the *input-VCF* path while the *de novo* path shipped a truth VCF
  contradicting its own reads, from v1.13.1 to v3.1.0.
- **Writing the test after the code, from the code.** It will encode current behavior,
  including the bug. This is the whole reason step 0 comes before step 1.
- **Stopping at the first plausible explanation.** `BND recall=0.000` drew three confident
  explanations before the real cause; each was plausible enough to stop at.
- **Working around a third-party limitation without checking its version.** `bnd_proximity.py`
  was correct code that should not have existed — built on "truvari cannot benchmark
  breakends", true of v4 and reversed in v5.0.0. And grep for *every* instance of the claim:
  a stale copy of that belief survived one retraction by a day.
- **Skipping tests for "small" or "obvious" changes.** There is no exemption.
- **Leaving a mutation in the working tree.** `mutate.sh run` restores automatically; if you
  mutated by hand, `git diff` before committing.

New war stories go in `docs/claude_engineering_audit.md`, not here.
