# The eidolon project

[![eidolon-tests](https://github.com/ncsa/eidolon/actions/workflows/eidolon-tests.yml/badge.svg?branch=main)](https://github.com/ncsa/eidolon/actions/workflows/eidolon-tests.yml)
[![DOI](https://zenodo.org/badge/765847780.svg)](https://doi.org/10.5281/zenodo.20100558)

> **Formerly `rusty-neat` / `rneat`** — renamed to `eidolon` in v2.0.0. Same tool, same
> NEAT lineage; the `rneat` command still works as a deprecated alias for one transition
> release. See `CHANGELOG.md`.

> **Upgrading from 2.0.0 → 3.0.0? The names of emitted output tokens changed.**
> See [Upgrading from 2.0.0](docs-site/src/upgrading/from-2-0-0.md) — **read it if you have
> any script that parses eidolon VCFs or FASTQ/BAM read names.**

`eidolon` is a Rust port of [NEAT](https://github.com/ncsa/neat): it simulates FASTQ that
looks like it came off a sequencer and carries your data's statistical properties, alongside
a golden BAM with ideal alignments and a truth VCF saying exactly what was planted. It adds
"noise" in the form of sequencing errors as it writes. Train models on your own data and
`eidolon` will reproduce that dataset's statistics — which is what makes it useful for tuning
alignment and variant-calling software.

Recent work targets cancer genetics: structural variants (CNV, BND, INV, INS), a native
tumor/normal workflow at configurable purity with an origin-tagged truth VCF, per-tissue
somatic models, and trinucleotide-context-aware SNP placement so mutational signatures
reproduce. Memory stays low and flat, and output is byte-identical for a given seed
regardless of thread count.

Tell us about your real-world experience by opening a Feedback issue — bugs, or things that
are not quite bugs. See `CHANGELOG.md` for the full release history.

## Install

```bash
conda install -c bioconda eidolon
```

A prebuilt binary with dependencies handled, no Rust toolchain required. Release binaries and
build-from-source instructions are in the docs.

## A minimal run

```yaml
# my_config.yml
reference: /path/to/reference.fa
read_len: 151
coverage: 10
ploidy: 2
output_dir: /path/to/output
output_filename: my_run
produce_fastq: true
produce_bam: true
produce_vcf: true
```

```bash
eidolon gen-reads -c my_config.yml
```

That writes FASTQ, a coordinate-sorted golden BAM, and a truth VCF. `eidolon --help` lists
every subcommand; `eidolon <subcommand> --help` covers one.

## Documentation

The full guide — every subcommand's config keys, the model builders, cancer simulation,
targeting, parallelism, HPC, and the versioning policy — is an mdBook site under
[`docs-site/`](docs-site/), with a sidebar and search.

Build and read it locally:

```bash
cargo install mdbook
mdbook serve docs-site --open
```

The pages are plain Markdown and readable directly on GitHub. Start at
[`docs-site/src/SUMMARY.md`](docs-site/src/SUMMARY.md) for the table of contents, or jump to:

| | |
|---|---|
| [Installing eidolon](docs-site/src/getting-started/installation.md) | install, build from source, the CLI tour |
| [Scope: germline is general, somatic is human](docs-site/src/overview/scope.md) | what is and is not claimed |
| [How `eidolon` compares to NEAT](docs-site/src/overview/comparison-with-neat.md) | feature, speed and memory comparison |
| [Cancer simulation how-to](docs/cancer_howto.md) | copy-paste tumor/normal guide |
| [Model builders](docs-site/src/models/mutation-model.md) | mutation, sequencing error, GC bias, fragment length |
| [Versioning and the public API](docs-site/src/reference/versioning.md) | what a MAJOR bump protects |
| [Upgrading from 2.0.0](docs-site/src/upgrading/from-2-0-0.md) | the v3.0.0 token rename |

## Citing

NEAT: Stephens et al. (2016), *PLOS ONE* 11(11):e0167047,
[doi:10.1371/journal.pone.0167047](https://doi.org/10.1371/journal.pone.0167047); and
Allen et al. (2026), *Journal of Open Source Software* 11(121):9056,
[doi:10.21105/joss.09056](https://doi.org/10.21105/joss.09056). `eidolon`:
[doi:10.5281/zenodo.20100558](https://doi.org/10.5281/zenodo.20100558).
