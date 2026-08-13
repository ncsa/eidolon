# Welcome to eidolon

A tour of the project for new teammates. **Everything marked *Deeper* is safely skippable** —
the top of each section carries the whole idea, and the deep dives are there for when you need
them.

If you are here to *contribute code*, read [`CONTRIBUTING.md`](CONTRIBUTING.md) after this.

---

## The whole idea, in one paragraph

**eidolon makes realistic fake DNA sequencing data where we already know the right answer.**

You hand it a reference genome and a description of what genetic changes to plant. It writes out
sequencing files that look like they came off a real instrument — and, alongside them, a "truth"
file listing exactly what it planted and where. You then run your normal analysis software on the
fake data and check whether it found what was actually there.

That last step is the point, and it is something you **cannot** do with real patient data: for a
real genome, nobody knows the complete correct answer. Simulated data is the only place you can
measure "how much did my pipeline miss?" rather than guess at it.

## Why anyone needs that

- **Benchmarking.** Two variant callers disagree on a real sample. Which is right? On simulated
  data you can just check.
- **Pipeline validation.** Before running a study on real samples, confirm the pipeline actually
  recovers what it should — and quantify what it drops.
- **Cancer.** A tumour biopsy is a *mixture*: some cancer cells, some normal tissue, often several
  competing cancer subpopulations. eidolon can build that mixture at a known composition, so you
  can measure how well a tool copes with it. This is the part of the project with the most
  original work in it.
- **Development.** Small, fast, known-answer datasets to test against.

## The mental model

```
   reference genome                          your analysis pipeline
   + trained models      ┌──────────┐        (aligner, variant caller)
   + what to plant  ───► │ eidolon  │ ───►  FASTQ reads  ───►  called variants
                         └──────────┘                                │
                                     └───►  TRUTH file  ────────────┐│
                                            (what was planted)      ▼▼
                                                              compare  ───►  "it found 86%
                                                                               of deletions"
```

Two outputs matter. The **reads** are what your pipeline consumes; they are meant to be
indistinguishable from real data in their statistical behaviour. The **truth file** is the answer
key, and your pipeline never sees it.

> ### Deeper: what "trained models" means *(skip freely)*
>
> Real sequencing data has quirks — some machines make more errors at the end of a read, some
> genome regions get read more often than others, fragment lengths follow a particular
> distribution. If simulated data ignores that, it is too clean and every tool looks better than
> it is.
>
> So eidolon *learns* those quirks from your own real data first, then reproduces them. Each
> `gen-*-model` subcommand builds one reusable model file: mutation patterns from a VCF,
> sequencing-error behaviour from FASTQ, fragment lengths and GC bias from a BAM. You build them
> once and reuse them across runs.

## What you can run

| Command | In plain terms |
|---|---|
| `gen-reads` | The main event: make reads + truth from a reference |
| `gen-cancer-reads` | Make a tumour/normal mixture at a chosen purity |
| `gen-mut-model` | Learn mutation patterns from real variant data |
| `gen-seq-error-model` | Learn a sequencer's error behaviour from real reads |
| `gen-frag-length-model` / `gen-gc-bias-model` / `gen-bam-models` | Learn physical/coverage properties from a real alignment |
| `compare-vcfs` | Score a caller's output against the truth file |
| `compare-af` | Check that per-variant allele fractions came out as requested |
| `validate` | Check an emitted file against what downstream tools will accept |
| `filter-reads` | Post-filter generated reads |

Everything is driven by a YAML config file:
`eidolon gen-reads -c my_config.yml`.

## Try it in five minutes

```bash
cargo build --release                  # or: conda install -c bioconda eidolon
./target/release/eidolon --help
./target/release/eidolon gen-reads --help
```

[`docs/cancer_howto.md`](docs/cancer_howto.md) is the best starting point for a real worked
example — copy-paste configs with explanations.

## Where things live

| Path | What is in it |
|---|---|
| `eidolon/` | The command-line program: argument parsing, per-subcommand runners |
| `eidolon-core/` | The engine: models, variant types, file readers/writers |
| `eidolon/tests/` | Integration tests — the ones that run the real binary end to end |
| `docs/` | Design documents, how-tos, and measurement records |
| `scripts/delta/` | The HPC validation harness (see below) |
| `CLAUDE.md` | Working practices, written for AI coding agents but useful to anyone |

## The unusual thing about this project

Most simulators are validated by "we ran it and it produced output." This one is not, and that is
deliberate — because a simulator that quietly generates *wrong* data is worse than one that
crashes. A crash you notice. Wrong truth data silently invalidates every benchmark built on it,
and you may not find out for months.

So the project holds an explicit standard: **nothing is "done" until there is evidence it works as
intended, and reports say which kind of evidence exists.** "Checked on a small known-answer
fixture, not yet on real data" is a normal and respectable thing to write here.

> ### Deeper: how that plays out *(skip freely)*
>
> Two habits do most of the work.
>
> **Assert content, not existence.** A test that checks "a file was produced" passes just as
> happily when the file is garbage. Real example from this repo: a check confirmed a count was
> above zero while every value in it was the malformed text `AF=AF=0.3000`.
>
> **Prove a test can fail.** Deliberately break the code a test covers and confirm the test
> notices. If it still passes, it was never testing anything. This has repeatedly caught tests
> that looked thorough and asserted nothing — in one case, eleven separate ways of breaking the
> code that the whole suite ignored.
>
> `docs/claude_engineering_audit.md` collects the case histories, including the ones where the
> tooling fooled us for weeks. It is unusually candid for a project document and it is the best
> single read for understanding how this codebase thinks.

> ### Deeper: HPC validation *(skip freely)*
>
> Small tests prove the machinery runs. They cannot prove the output looks like real data at
> genome scale. So there is a second tier on NCSA Delta: simulate a whole tumour genome, run it
> through production variant callers (Manta, Delly, GATK), and score their output against the
> truth with an independent tool (truvari).
>
> That harness lives in `scripts/delta/` and is written in bash. It carries its own positive and
> negative controls — it scores the truth against *itself* (must be perfect) and against a
> deliberately shifted copy (must find nothing), because a scoring configuration loose enough to
> match anything would otherwise report excellent results. Findings from those runs feed back into
> `docs/sv_support_matrix.md`, which records exactly which capabilities are measured, which are
> broken, and which are unverified.

## How the work is organised

Changes go through pull requests targeting the **`develop`** branch. `main` holds releases.
Issues on GitHub track bugs and planned work — including a **Feedback** issue type for things
that are not quite bugs, like a confusing error message or a rough edge.

Versions follow [Semantic Versioning](https://semver.org) as of v3.0.0. The project was called
`rusty-neat` / `rneat` before v2.0.0, so older notes and links use those names; the emitted output
tokens changed at v3.0.0 as well, which matters if you have scripts that parse eidolon's output.

## Questions

Open an issue, or message **Joshua Allen** (links on the repo page). Questions that turn out to be
gaps in the documentation are useful — say so and it gets fixed.

Worth knowing if you are wondering whether to ask: nobody here expects you to have read the whole
codebase, and "I could not tell from the docs" is treated as a docs bug rather than a
you-problem.
