# The eidolon project

[![eidolon-tests](https://github.com/ncsa/eidolon/actions/workflows/eidolon-tests.yml/badge.svg?branch=main)](https://github.com/ncsa/eidolon/actions/workflows/eidolon-tests.yml)

> **Formerly `rusty-neat` / `rneat`** — renamed to `eidolon` in v2.0.0. Same tool, same
> NEAT lineage; the `rneat` command still works as a deprecated alias for one transition
> release. See `CHANGELOG.md`.

> **Upgrading from 2.0.0 → 3.0.0? The names of emitted output tokens changed.**
> See [Upgrading from 2.0.0](upgrading/from-2-0-0.md) below — **read it if you have any
> script that parses eidolon VCFs or FASTQ/BAM read names.**

Welcome to `eidolon`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates fastq that appear to be from sequencers, carry the same statistical properties as your data, and generate a golden bam and fastq that gives you ideal alignments and what variants were inserted. In addition, `eidolon` generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment and variant calling software to your data. Training models on your data will allow `eidolon` to faithfully reproduce the statistical properties of your dataset.

We have spent some dedicated time toward gearing the current version of `eidolon` to simulate cancer genetics, including adding structural variant simulations (CNV, BND, SVs, and others), and creating a wrapper that simulates purity levels and then stitches the results back together. We've geared this software with an aim at keeping memory usage as low and CPU time as short as possible. Let us know your real world experience by creating a Feedback issue, if you have something that's not quite a bug, or have a positive experience to share. As always, let us know if you find a bug and give as many details as you can to help us troubleshoot.

`eidolon` trains reusable models from your own data (mutation, sequencing-error, fragment-length, GC-bias, and BAM/alignment models), outputs a golden BAM with ideal alignments alongside the FASTQ and truth VCF, accepts a BED to target read creation to regions, and can read custom variants from an input VCF — including per-variant allele frequencies, to reproduce a continuous AF spectrum for pooled or somatic data. More recent releases add native tumor/normal simulation (`eidolon gen-cancer-reads`), structural variants and copy number, per-tissue somatic models, and trinucleotide-context-aware SNP placement so context-specific mutational signatures reproduce. A file-streaming writer keeps disk I/O and footprint low. See `CHANGELOG.md` for the full release history, and please open a Feedback issue with your real-world experience.

Find us on Zenodo:
[![DOI](https://zenodo.org/badge/765847780.svg)](https://doi.org/10.5281/zenodo.20100558)
