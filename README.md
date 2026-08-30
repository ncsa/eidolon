# The eidolon project

> **Formerly `rusty-neat` / `rneat`** — renamed to `eidolon` in v2.0.0. Same tool, same
> NEAT lineage; the `rneat` command still works as a deprecated alias for one transition
> release. See `CHANGELOG.md`.

> **Upgrading from 2.0.0 → 3.0.0? The names of emitted output tokens changed.**
> See [Upgrading from 2.0.0](#upgrading-from-200) below — **read it if you have any
> script that parses eidolon VCFs or FASTQ/BAM read names.**

Welcome to `eidolon`, a Rust port of NEAT (https://github.com/ncsa/neat), a genetic simulation program that creates fastq that appear to be from sequencers, carry the same statistical properties as your data, and generate a golden bam and fastq that gives you ideal alignments and what variants were inserted. In addition, `eidolon` generates "noise" in the form of sequencing errors as it is writing out files. These features can help you hone in your alignment and variant calling software to your data. Training models on your data will allow `eidolon` to faithfully reproduce the statistical properties of your dataset.

We have spent some dedicated time toward gearing the current version of `eidolon` to simulate cancer genetics, including adding structural variant simulations (CNV, BND, SVs, and others), and creating a wrapper that simulates purity levels and then stitches the results back together. We've geared this software with an aim at keeping memory usage as low and CPU time as short as possible. Let us know your real world experience by creating a Feedback issue, if you have something that's not quite a bug, or have a positive experience to share. As always, let us know if you find a bug and give as many details as you can to help us troubleshoot.

`eidolon` trains reusable models from your own data (mutation, sequencing-error, fragment-length, GC-bias, and BAM/alignment models), outputs a golden BAM with ideal alignments alongside the FASTQ and truth VCF, accepts a BED to target read creation to regions, and can read custom variants from an input VCF — including per-variant allele frequencies, to reproduce a continuous AF spectrum for pooled or somatic data. More recent releases add native tumor/normal simulation (`eidolon gen-cancer-reads`), structural variants and copy number, per-tissue somatic models, and trinucleotide-context-aware SNP placement so context-specific mutational signatures reproduce. A file-streaming writer keeps disk I/O and footprint low. See `CHANGELOG.md` for the full release history, and please open a Feedback issue with your real-world experience.

Find us on Zenodo:
[![DOI](https://zenodo.org/badge/765847780.svg)](https://doi.org/10.5281/zenodo.20100558)

## Scope: germline is general, somatic is human

These two halves of `eidolon` have deliberately different reach, and it is worth knowing which
one you are using.

**Germline simulation is species-agnostic.** `gen-reads` needs a reference and, optionally,
models trained from your own data — nothing in that path assumes a human genome. It has been
exercised on H1N1 (13.5 kb, eight contigs), *E. coli*, yeast, soy, and human references from a
single chromosome up to whole-genome GRCh38. If you are simulating germline variation for any
organism, that is a supported use.

**Somatic / cancer simulation is specialized to human.** The bundled mutational models are
derived from human cancer cohorts — COSMIC signatures and PCAWG structural-variant
distributions — and the per-tissue models (BRCA, skin, lung) are human tissue types. SV size
and type distributions, mutation rates, and signature weights all carry that provenance.

`gen-cancer-reads` will *run* against a non-human reference and produce output, but the model
driving it has no meaning there: a yeast genome given a human pan-cancer SV spectrum produces
numbers, not biology. **We do not claim cancer simulation for non-human references**, and the
validation campaigns behind those models are human-only. Use another organism's reference for
cancer work only if you have trained your own models against matched data, and treat the
bundled ones as human-specific.

This split is also why the two halves are validated differently: germline correctness is
checked across deliberately varied contig shapes and sizes, while somatic realism is measured
against real human sequencing data.

## Cancer simulation

`eidolon` simulates tumor / normal sequencing data end-to-end. The
`eidolon gen-cancer-reads -c <config.yml>` subcommand runs two `gen-reads` passes —
one normal-genotype, one tumor — over the same reference and merges them at a
configurable purity into a single "tumor biopsy" FASTQ that downstream somatic
callers (Mutect2, Strelka, Manta, …) consume directly, plus an origin-tagged truth
VCF (`INFO/EIDOLON_ORIGIN ∈ {germline, somatic, shared}`) for scoring. It also
generates the foundational structural-variant types cancer SVs depend on — `<BND>`
translocations, `<INV>` inversions, and de novo `<INS>` — and ships bundled
pan-cancer and per-tissue (BRCA / skin / lung) models plus Docker-based benchmark
pipelines.

**See [`docs/cancer_howto.md`](docs/cancer_howto.md) for a copy-paste guide** with
worked examples, output reference, benchmarking, and model training. Design
rationale and calibration caveats live in
[`docs/cancer_simulator.md`](docs/cancer_simulator.md).

## How `eidolon` compares to NEAT

`eidolon` is a Rust port of NEAT that tracks the NEAT feature set while adding a
native cancer workflow, a low and flat memory footprint, and reproducible
output. The table below compares the original Python 2 NEAT (the 2.x "genReads"
line), the current Python 3 NEAT 4.x, and `eidolon`.

|                                            | **NEAT 2.x** (genReads)        | **NEAT 4.x**                              | **`eidolon`**                                              |
| ------------------------------------------ | ------------------------------ | ----------------------------------------- | -------------------------------------------------------- |
| Latest version                             | 2.1                            | 4.6.1                                      | 3.2.0                                                    |
| Language                                   | Python 2                       | Python 3                                  | Rust                                                     |
| FASTQ reads (single / paired)              | ✅                             | ✅                                        | ✅                                                       |
| Golden BAM + VCF truth set                 | ✅                             | ✅                                        | ✅                                                       |
| SNPs + indels                              | ✅                             | ✅                                        | ✅                                                       |
| Structural variants                        | ❌ none                        | ❌ placeholder classes only, not wired in (see note) | ✅ Native: BND, INV, INS, **CNV**              |
| Native tumor / normal cancer workflow      | ❌ (tutorial script only)      | ❌ (stated roadmap)                       | ✅ `gen-cancer-reads` (purity mix + origin-tagged truth) |
| Cancer / per-tissue models                 | ❌                             | ❌                                        | ✅ pan-cancer + BRCA / skin / lung (COSMIC / PCAWG)      |
| Intra-tumor heterogeneity (subclones)      | ❌                             | ❌                                        | ✅ CCF mixtures; `purity × dosage × CCF`; PyClone-VI / DPClust ingest; `EIDOLON_CCF` / `EIDOLON_VAF` truth |
| Empirical mutation model (trinucleotide)   | ✅                             | ✅                                        | ✅ (incl. positional trinucleotide bias — SBS-96 cosine 0.99) |
| Sequencing error model                     | ✅                             | ✅                                        | ✅                                                       |
| SNP transition matrix source               | fit from VCF                   | fit from VCF                              | ✅ VCF, **or inferred from a BAM's MD tags**, or a custom 4×4 TSV |
| Quality-score binning                      | ❌                             | ✅ + named instrument presets (`--quality-preset novaseq`) | ✅ explicit bin list (`binned_quality_bins`); no named presets |
| Default het/hom ratio (bundled model)      | ~99 (`homozygous_frequency` 0.010) | ~999 (0.001)                          | ✅ ~2.0 (0.333) — human-realistic, locked by a test      |
| Fragment-length + GC-bias models           | ✅                             | ✅                                        | ✅ (incl. one-pass `gen-bam-models`)                     |
| BED targeting / VCF variant insertion      | ✅                             | ✅                                        | ✅                                                       |
| Continuous per-variant allele fraction     | ❌ (genotype `{0.5, 1.0}`)     | ❌ (genotype `{0.5, 1.0}`)                | ✅ input-VCF `AF`/`AD` → matched AF spectrum (pooled / somatic) |
| Long reads (ONT / PacBio)                  | ❌                             | ⚠️ PacBio-like single-end "given a model" | ⚠️ `long_reads:` fragment mode; no long-read error model yet (#319) |
| Parallelism                                | Manual job sharding (`--job`)  | Multiprocessing (`--threads`): genome split into ~8 chunks/thread, then stitched | Multithreading (rayon) |
| Deterministic output                       | Not guaranteed                 | Not guaranteed across thread counts       | ✅ byte-identical for a given seed, across thread counts (`determinism.rs` pins single-threaded and rayon-default) |
| Speed (single thread, vs NEAT 4.6.1)       | —                              | 1× (baseline)                             | **~10–13× faster** (E. coli 10.4× · yeast 11.2× · chr22 13.3×) |
| Peak memory (same runs)                    | —                              | 1× (baseline)                             | **3–7× less**, and flatter (chr22 227 MB vs 1525 MB)      |
| VCF comparison tooling                      | Bundled scripts                | ✅ `compare-vcfs`                          | ✅ `compare-vcfs`                                        |
| I/O / memory                               | Temp files                     | Temp files                                | Streaming writes, low-memory focus                       |
| Versioning policy                          | ad hoc                         | ad hoc                                    | SemVer 2.0.0 with a declared public API (as of v3.0.0)   |
| Distribution                               | GitHub source                  | GitHub / PyPI                             | GitHub + binaries; Bioconda (`conda install -c bioconda eidolon`) |

Speed and memory are single-thread medians (n=3, 10× coverage) on one exclusive
Delta node against NEAT 4.6.1; see `docs/access_report_draft.md` §3.2. Germline
fidelity is statistically equivalent between the two tools through an identical
GATK pipeline (§3.3), which is the intended result for a port.

**On NEAT 4 and structural variants.** NEAT 4.6.1 ships `Inversion`, `Duplication`,
`Translocation`, `Transposition` and `CopyNumberVariant` classes under
`neat/variants/`, but none is reachable: none is exported from that package's
`__init__.py`, none is referenced anywhere outside its own file, and
`VariantTypes.types` lists only `Insertion`, `Deletion`, `SingleNucleotideVariant`
and `UnknownVariant`. The mutation model can emit only SNVs, insertions and
deletions. The input-VCF path classifies records by REF/ALT length arithmetic, so a
symbolic `<INV>`/`<DUP>` or a breakend ALT becomes an `UnknownVariant` — whose
`get_alt()` and `get_ref_len()` both raise, since its `__init__` sets no `alt` and
its metadata is nested one level deeper than those accessors expect. NEAT 4 does not
generate structural variants and does not render them from an input VCF either.
NEAT 2.1 has no structural-variant code at all. Verified against tags `2.1` and
`4.6.1` on 12 Aug 2026.

**Citations.** NEAT: Stephens et al. (2016), *PLOS ONE* 11(11):e0167047,
[doi:10.1371/journal.pone.0167047](https://doi.org/10.1371/journal.pone.0167047);
and Allen et al. (2026), *Journal of Open Source Software* 11(121):9056,
[doi:10.21105/joss.09056](https://doi.org/10.21105/joss.09056). `eidolon`:
[doi:10.5281/zenodo.20100558](https://doi.org/10.5281/zenodo.20100558).

**Where `eidolon` fits.** NEAT 4.x is a capable, actively developed simulator, and
the two tools share most of their core feature set. `eidolon` is the right choice
when you want:

- **Cancer simulation** — a native tumor/normal workflow (`gen-cancer-reads`)
  with configurable purity and an origin-tagged truth VCF, plus CNVs and
  per-tissue cancer models. This is `eidolon`-only.
- **Intra-tumor heterogeneity** — a tumor as a mixture of subclones at distinct
  cancer-cell fractions, either specified inline, fitted from PyClone-VI /
  DPClust output, or replayed from a real tumor's observed VAF. Validated on
  SEQC2 HCC1395 at r = 0.99 against the intended per-site VAF. Also
  `eidolon`-only — the axis the NEAT lineage cannot represent at all.
- **Low, flat memory** — streaming FASTQ writes keep peak RSS small and roughly
  constant across genome size and thread count, which matters on shared HPC
  allocations.
- **Reproducibility** — the same seed yields byte-identical FASTQ regardless of
  the thread count.
- **Easy deployment** — a single self-contained binary with no Python
  environment to manage, installable via Bioconda
  (`conda install -c bioconda eidolon`).

## Versioning and the public API

As of **v3.0.0**, eidolon follows [Semantic Versioning 2.0.0](https://semver.org): a
MAJOR bump means something you may depend on changed incompatibly, a MINOR bump adds
functionality compatibly, and a PATCH bump is fixes only. Releases before v3.0.0 were
versioned less strictly — notably v1.11.0 and v1.12.0 changed the FASTQ read-name format
in minor bumps.

SemVer requires a project to say what its public API is. For eidolon:

**Public API — a change here means a MAJOR bump:**
- Names and semantics of emitted VCF INFO tags, the VCF sample-column name, and the
  `FILTER` / `FORMAT` fields eidolon writes
- FASTQ/BAM **read-name (QNAME) format** — the `EIDOLON_generated_` /
  `EIDOLON_chimeric_` prefixes and the positional fields encoded after them
- CLI subcommand names, flags, and configuration-YAML keys
- Model-file compatibility (a model built by one version staying readable by the next)

**Not public API — may change in a MINOR or PATCH release:**
- The `eidolon-core` Rust library surface; it exists to serve the binary
- Log output, progress reporting, and human-readable messages
- The exact simulated *content* for a given seed — reads are a random draw, so sampler
  and model changes legitimately alter output while preserving the format

If you parse eidolon's output, the first list is what you are relying on, and a major
version bump is your signal to check this file before upgrading.

## Upgrading from 2.0.0

v3.0.0 renamed the tokens eidolon **emits**. There is no behavior change — the same
reads and variants are produced — but anything that *parses* eidolon output needs
attention. In-tree consumers were all migrated; your own scripts were not.

| Surface | v2.0.0 and earlier | v3.0.0+ | Migration |
|---|---|---|---|
| VCF INFO tags | `NEAT_ORIGIN`, `NEAT_PROVENANCE`, `NEAT_REASON`, `NEAT_CCF`, `NEAT_VAF` | `EIDOLON_*` | **converter provided** |
| VCF sample column | `NEAT_simulated_sample` | `EIDOLON_simulated_sample` | **converter provided** |
| FASTQ read names | `@RNEAT_generated_*`, `RNEAT_chimeric_*` | `@EIDOLON_generated_*`, `EIDOLON_chimeric_*` | **not converted** — see below |
| BAM read names (QNAME) | `RNEAT_generated_*` | `EIDOLON_generated_*` | **not converted** — see below |

### VCFs — convert them

```bash
tools/migrate_legacy_tokens.sh old_truth.vcf.gz new_truth.vcf.gz
tools/migrate_legacy_tokens.sh --check old_truth.vcf.gz   # report only; exit 10 = legacy
```

Idempotent: safe to run on a file that's already current, or on a non-eidolon VCF, so
you can call it unconditionally in a pipeline. It rewrites the header declarations and
the record fields together, so the result is a valid VCF.

Note you **cannot** work around this with a dual-name bcftools filter — bcftools rejects
an undefined tag in a `-i` expression at parse time, so even
`INFO/EIDOLON_ORIGIN="somatic" || INFO/NEAT_ORIGIN="somatic"` fails outright on a file
that carries only one of the two. Convert the file instead.

### Read names — NOT converted, and this is the part that can bite you

Read names are **not** rewritten, in either FASTQ or BAM. Rewriting files that are
routinely hundreds of GB to change ~19 bytes per record isn't a reasonable ask, so
instead `eidolon filter-reads` accepts **both** prefixes natively — legacy FASTQs keep
working with no conversion step, and it warns once when it sees an old prefix.

**Nothing protects your own scripts.** A `grep '^@RNEAT_generated_'`, `awk` field split,
or regex over read names will match **zero records and exit 0** against v3.0.0 output —
a silent empty result, not an error. Audit for the old prefix before upgrading:

```bash
grep -rn 'RNEAT_generated_\|RNEAT_chimeric_' your_scripts/
```

#### Best fix: make your parser prefix-agnostic

**Only the prefix changed.** Everything after `_generated_` —
`<contig>_<start>_<end>_<uniq>/<mate>` — is byte-for-byte identical between versions
(verified across every read of a golden BAM re-prefixed both ways). So rather than
rewriting files, strip the prefix and your script works against *either* version, and
against any future rename:

```bash
# version-proof: drop whatever prefix is present, keep the encoded fields
samtools view x.bam | cut -f1 | sed -E 's/^[A-Z]+_generated_//'
zcat x_r1.fastq.gz  | awk 'NR%4==1' | sed -E 's/^@[A-Z]+_generated_//'
```

Same trick for chimeric reads: `sed -E 's/^@?[A-Z]+_chimeric_//'`.

If you are joining an old artifact against a new one on read name, normalize both sides
this way. One unrelated gotcha for QNAME joins: eidolon's golden BAM **keeps** the `/1`
`/2` mate suffix, while bwa strips it from aligned BAMs — so strip that too
(`sed 's|/[12]$||'`) when joining a golden BAM to an aligner's output. That mismatch
predates this release.

#### If you must rewrite the files instead

```bash
# FASTQ (expensive — a full re-compress; prefer the prefix-agnostic parse above)
zcat old_r1.fastq.gz | sed 's/^@RNEAT_generated_/@EIDOLON_generated_/' | gzip > new_r1.fastq.gz

# BAM QNAMEs (QNAME is field 1, and no header line contains the prefix, so ^ is safe)
samtools view -h old.bam | sed 's/^RNEAT_generated_/EIDOLON_generated_/' | samtools view -b -o new.bam
```

BAM QNAMEs are covered by no in-tree tooling — `filter-reads` reads FASTQ, and the
converter is bcftools-based so it cannot touch a BAM. Nothing eidolon ships parses BAM
read names either, so this only affects your own scripts.

# How to use `eidolon`

## Prerequisites

The easiest way to install `eidolon` is via [Bioconda](https://bioconda.github.io/):

```
conda install -c bioconda eidolon
```

This pulls a prebuilt binary with all dependencies handled — no Rust toolchain
required. If you prefer to build from source or grab a release binary, read on.

You will need to install the rust toolchain to compile `eidolon`, including `cargo`. Check the cargo documentation for instructions (https://doc.rust-lang.org/cargo/getting-started/installation.html). Alternatively, you can try one of the binaries on the release page. Select the one that matches your system and let us know if you run into errors. During compilation, you may run into errors, such as cmake not found. Some of the packages `eidolon` uses have these dependencies. For Debian/Ubuntu this should be a simple `sudo apt install cmake` and for RHEL/Rocky type distros this should be `sudo dnf install cmake`. There may be some other requirements. Drop a comment if you need specific help.

Download the executable in the release (current version 3.2.0).

```bash
$ eidolon --help
Usage: eidolon [OPTIONS] [SUB-COMMAND]

SUB-COMMANDS:
  gen-reads              Generates reads for an input dataset
  filter-reads           Filters the output of gen-reads
  gen-mut-model          Generates a mutation model from real VCF data
  gen-seq-error-model    Generates a sequencing error model from real FASTQ data
  gen-frag-length-model  Generates a fragment length model from a BAM or SAM file
  gen-gc-bias-model      Generates a GC bias model from a reference FASTA and aligned BAM
  gen-bam-models         Builds multiple models (frag-length, GC bias) from one BAM in a single pass
  compare-vcfs           Compares a NEAT-simulated golden VCF against a downstream variant-caller VCF
  gen-cancer-reads       Simulates a tumor/normal mixture (two gen-reads passes + merge)
  validate               Checks an emitted FASTQ/VCF/BAM against what its consumers accept
  compare-af             Per-allele AF correlation between a truth VCF and a simulated one
  help                   Print this message or the help of the given subcommand(s)

Options:
      --log-level [<log_level>]  Verbosity of the written .neat.log. The on-screen log is always
                                 info. [default: info] [possible values: trace, debug, info,
                                 warn, error, off]
      --log-dest <log_dest>      Sets the log destination (full path with full filename) for the written log
  -h, --help                     Print help
```

`validate` and `compare-af` take file paths and flags directly rather than a
configuration YAML — both are diagnostic/measurement runs invoked ad hoc and from
harness scripts, where writing a YAML per check is friction.

To check options for a subcommand:

```bash
$ eidolon gen-reads --help
Generates reads for an input dataset

Usage: eidolon gen-reads [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

To run filter reads, check the help menu.

```bash
$ eidolon filter-reads --help
Filters the output of gen-reads

Usage: eidolon filter-reads --configuration-yaml <configuration_yaml>

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

To run gen-mut-model:

```bash
$ eidolon gen-mut-model --help
Generates a mutation model from input VCF data

Usage: eidolon gen-mut-model [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

Use the help menu to see the available options and leave an issue if you find something bad happening.

To compile and run `eidolon` yourself, besides the rust toolchain, you will need `git` installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

```bash
~/$ git clone git@github.com:ncsa/eidolon.git
~/$ cd eidolon
~/eidolon/$
```

For Windows and Mac users, try the binary packaged with the latest version of `eidolon`.

Once in the repo, you can build the program either in debug (default) or release mode. The main difference is how much info it gives you if there is an error. Release mode also has some optimizations to run it faster.

```bash
~/eidolon/$ cargo build --release
```

If you prefer to run the package directly without using the binary, you can also use
```bash
~/eidolon/$ cargo run -- gen-reads -c my_config.yml
```

Rust will download any required packages. Compiling Rust code is the slowest part of the process. The final binary will be built and the program run immediately after in the second case. To run the program manually, from the repo main dir, run

```bash
~/eidolon/$ ./target/release/eidolon -h
```

`eidolon` uses a configuration file to read values it needs for the run. A command line execution might look like this:

```bash
~/eidolon/$ ./target/release/eidolon -c /path/to/filled/in/config.yml
```
If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results.

Fastq Output
============
The fastq output will have a key name that identifies the block where the read was drawn from, for quick comparisons in alignments. The output BAM file will contain the original sequence and cigar string. The name will have the format `EIDOLON_generated_<contig short name>_<fragment_start>_<fragment_end>_<uniq>/1` (or `/2` for the second read in a pair), where start and end are zero-padded to 10 digits and `<uniq>` is a 16-digit hex per-fragment tag that keeps same-position fragments from colliding (see #210).

```bash
@EIDOLON_generated_Chromosome_0000000000_0000000353_000000000000002a/1
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTTAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGC
+
>AC7<GDEGGGGEFGA<GFCGG;GGGGGF>GEGGEGGGFGFEFGCEGGGGGGCG:AEFGFFGG>FG;GDGA9$GGAGF=GFG=EFFCGGGFGGGGGC$BFGEFFAGGG9F7E@>?GFGGGG>EBFGFDG)DGC6DEDFA2EG:EGG%FFB
```
The above example comes from the chromosome in the reference (which was simply named "Chromosome") between 0 and 353, which is the size
of a fragment in the simulated DNA. If there is a pair with this read, it will have the same coordinates, though it started at index 352 instead.

BAM Output
==========
`eidolon` can write a golden BAM file alongside the FASTQ output. The BAM contains the same reads as the FASTQ — same sequences, same quality scores, same variants and sequencing errors applied — with alignment information included. To enable it, set `produce_bam: true` in your `gen-reads` config. The output path is derived automatically from `output_filename` (e.g. `output_filename: my_run` → `my_run.bam`). Please note that the BAM has a higher overhead than the fastq, and may take longer to produce.

```yaml
produce_bam: true
```

The CIGAR strings in the BAM reflect the full ground-truth alignment:
- Genomic variants (SNPs, insertions, deletions) are encoded as `M`, `I`, and `D` ops.
- Sequencing error indels are also encoded: deletion errors add `D` ops for skipped reference bases; insertion errors add `I` ops for inserted bases.

The BAM is written in coordinate-sorted order. No post-processing sort is required before indexing:

```bash
samtools index my_run.bam
```

Note: heterozygous variants are applied probabilistically, identical to the FASTQ. A read drawn to the reference allele will not show the variant in either the FASTQ or the BAM.

Targeted Generation with a BED File
=====================================
For exome, targeted-panel, or any run where you only want reads over specific regions, you can supply a BED file at generation time:

```yaml
target_bed: /path/to/targets.bed
```

When `target_bed` is set, `gen-reads` skips contigs absent from the BED entirely — no blocks are read, no variants are placed, and no reads are generated for those contigs. Within covered contigs, reads and variants are generated only over the intersecting non-N regions. This is the recommended approach for targeted runs on large genomes; it is far more efficient than generating genome-wide reads and post-filtering with `filter-reads`.

The BED contig names must match the short names derived from the reference FASTA (the text after `>` up to the first whitespace character). Both `.bed` and `.bed.gz` inputs are accepted.

### `coverage` on a targeted run is a fragment budget, not on-target depth

**On a targeted run, expect on-target depth to come in *below* the `coverage` you asked for, and scale your request up to compensate.** This is a real property of targeted sequencing, not a defect, but the arithmetic is not obvious and `coverage` does not currently correct for it.

A fragment is *owned* by the target its start falls in, but its far end runs into the flanking sequence beyond the target — as it must: a 250 bp fragment cannot fit inside a 178 bp exon. Those spilled bases are real reads (they are exactly the off-target reads a real capture kit produces), but they are not on-target depth. So:

```
delivered on-target depth  ≈  coverage × on-target base fraction
```

and that fraction is governed by target width against fragment length. Measured on a 436-target BED (mean width 178 bp) over real GRCh38 chr22 sequence, `read_len: 100`, `fragment_mean: 250`, `coverage: 60`:

| | value |
|---|---|
| fragments placed | 23,462 (101% of the requested budget) |
| on-target base fraction | 40.5% |
| **delivered on-target depth** | **24.5×** (against a requested 60×) |

To get ~60× on target with those parameters, ask for roughly `60 / 0.405 ≈ 148×`. The fraction rises quickly as targets get wider relative to fragments, so a panel with 1 kb targets loses far less than a 178 bp exome, and a whole-genome run (no BED) loses nothing.

Two related notes:

- **Do not read a 100%-on-target simulation as "better".** eidolon before v3.1.x confined every fragment inside its target, which yields 100% on-target and is physically impossible for targets narrower than the fragment length. It also under-delivered the fragment budget itself (60% of requested, in the run above) and collapsed the realized fragment-length distribution near target edges.
- Making `coverage` mean on-target depth directly — i.e. having eidolon scale the budget for you — is tracked in [#578](https://github.com/ncsa/eidolon/issues/578). Until then, scale it yourself.

Custom Mutation Rate Regions with a BED File
============================================
In addition to targeting certain regions, you can use a BED file with a custom field entered in any column after the third. 

```yaml
mutation_regions: /path/to/mutation_regions.bed
```

The text pattern is "mut_rate=0.0001" followed by a delimiter (end of line, space, semicolon, comma, bar) where the number after the equal sign is any float, and will be the mutation rate for the region defined in columns 1-3. 
```text
chr1    39930   39957   Bed_info    mut_rate=0.002
chr2    0   111199  Other_bed_info  mut_rate=0.02
```
Any region outside of the regions defined in the mutation regions BED will be assigned the default mutation rate, set by the mutation model, which can be overridden a custom default defined in the config file. For example, if you only want to have mutations in exomes, but you want reads from the full genome, you could set the `mutation_rate` to 0.0 and use a bed file defining all exomes, with a column appende `mut_rate=0.0011` or whatever mutation rate you desire. This BED works in conjunction with the `target_bed`, and any region outside the target bed will be excluded from reads anyway, and thus will have no variants. 

Input Variants VCF
==================
You can supply a VCF of variants to force into the simulation:

```yaml
input_vcf: /path/to/variants.vcf.gz
```

`eidolon` will place every variant from the VCF into the corresponding position in the simulated reads and the output VCF. Random variants are still generated at `mutation_rate` for all positions not covered by the input VCF; set `mutation_rate: 0.0` to disable random variants entirely and output only the provided set of variants.

**Requirements**

- The VCF must be single-sample (one sample column).
- Every record must include `GT` in the FORMAT field. `eidolon` uses the genotype to determine whether to apply the variant to all reads covering the position (homozygous, e.g. `1/1`) or only a probabilistic subset (heterozygous, e.g. `0/1`). Records without `GT` are rejected.
- Contig names must match the short names derived from the reference FASTA (text after `>` up to the first whitespace character). Variants on unrecognised contigs are skipped with a warning.
- Both `.vcf` and `.vcf.gz` files are accepted.

**Supported variant types**

| Type | Condition | Handled |
|------|-----------|---------|
| SNP | REF and ALT both single base | Yes |
| Insertion | single-base REF, multi-base ALT | Yes |
| Deletion | multi-base REF, single-base ALT | Yes |
| Symbolic SV | ALT is `<DEL>` / `<DUP>` / `<CNV>` / `<INS>` / `<INV>` / breakend / other `<TAG>` | Yes — see "Symbolic / structural variants" below |
| Literal complex | multi-base REF **and** multi-base ALT (literal bases) | **No** — skipped with warning |

**Symbolic / structural variants**

Symbolic ALTs (VCF 4.2 §1.4) are accepted and round-tripped to the output VCF verbatim, with `INFO/END`, `INFO/SVLEN`, and `INFO/CN` preserved. As of v1.10, `eidolon` can also generate symbolic SVs *de novo* from a learned model — opt in by setting `sv_rate_scale: 1.0` (or higher) in your gen-reads YAML; see "De novo SV generation" below.

| SV | Effect on read depth | Effect on read sequence |
|----|----------------------|-------------------------|
| `<DEL>` | hom → ×0, het → ×(ploidy−1)/ploidy; or ×CN/ploidy if `INFO/CN` is set | none — bases in the deleted span just stop producing reads |
| `<DUP>` | hom → ×2, het → ×(ploidy+1)/ploidy; or ×CN/ploidy if `INFO/CN` is set | none — extra reads come from the forward-strand reference |
| `<CNV>` | ×CN/ploidy when `INFO/CN` is set; otherwise warned and passed through with no depth change | none |
| `<INS>` | none (insertion at a single anchor base — no span to modulate) | none — the inserted sequence is not synthesized |
| `<INV>` | none in this release | **not modeled** — reads come from the forward-strand reference, not the inverted sequence |
| Breakends, unknown `<TAG>` | none | none — round-tripped only |

DEL anchor convention: for `<DEL>`, POS is the unaffected base immediately before the deletion (per VCF 4.2), so the modulated span is `[POS+1, END]` in 1-based coords. For `<DUP>` / `<CNV>` / `<INV>`, POS is the first base of the affected region, so the span is `[POS, END]`.

When an SV zeroes out coverage (hom `<DEL>` or `INFO/CN=0`), the mutation rate over the same span is also zeroed so de-novo SNPs don't pollute the output VCF with variants that never appear in reads.

**Current caveats**

- *Multi-allelic records*: only the first ALT allele is used; additional alleles are silently ignored. Split multi-allelic records with `bcftools norm -m -` before passing to `eidolon`. (This is true for both literal and symbolic ALTs — `<DEL>,<DUP>` on the same line is treated as the first ALT only.)
- *REF allele verification*: `eidolon` does not check that the REF field matches the reference sequence at that position. Mismatches will produce biologically incorrect output without any warning.
- *Literal complex variants*: records whose REF and ALT are both multi-base strings (and the ALT is literal bases, not a `<TAG>`) are skipped with a logged warning and do not appear in the output.
- *`<INV>` interior is forward-strand*: inversion **junction** reads are emitted (the chimeric signal callers use to detect the inversion), but reads from the *interior* of an inverted span are still transcribed from the forward-strand reference. The `<INV>` record round-trips to the output VCF.
- *Breakends*: `<BND>` junction reads are generated (chimeric reads spanning the mate locus; validated against Manta — see the ACCESS report), and the record round-trips to the output VCF.

**De novo SV generation**

`eidolon` can sample SVs directly from a learned `SvModel` rather than relying on user-supplied records: `<DEL>`, `<DUP>`, `<CNV>`, `<INV>` and literal insertions are drawn per contig, and `<BND>` translocations are drawn separately as mate pairs across two contigs. **`<BND>` therefore needs a reference with at least two contigs** — on a single-contig reference the correct yield is zero rather than a same-contig junction, so a chr22-only run reports `truth BND: 0` by construction. Generation is **off by default** (opt-in via `sv_rate_scale`) so v1.9 pipelines remain unchanged.

To enable, add a single line to your gen-reads YAML:

```yaml
sv_rate_scale: 1.0
```

- `0.0` (the default) disables de novo SV generation entirely — only `input_vcf` records flow through.
- `1.0` reproduces the rate from the trained model.
- Larger values scale the rate proportionally for stress testing.

When enabled, `eidolon` consults `MutationModel.sv_model`:

1. **If you trained your own model** with `eidolon gen-mut-model -c <yaml>` against an SV-rich VCF (e.g. a gnomAD-SV slice), the model file includes a fitted SV component covering type / length / copy-number / homozygous frequencies. Pass it via `mutation_model: /path/to/your.json.gz` in the gen-reads YAML.
2. **If you don't supply a model**, gen-reads loads the bundled default. The default carries **approximate** gnomAD-SV v2.1 parameters — useful for kicking the tires, but not a substitute for retraining on data that matches your downstream use case. See `eidolon-core/src/models/sv_model_defaults.rs` for the parameter sources.

Per-contig sampling: Poisson count from `per_base_rate × contig_len × sv_rate_scale` → weighted type pick → log-normal length (rejected outside `[50bp, contig_len / 4]`) → uniform anchor with overlap and N-gap rejection → `INFO/CN` draw for `<CNV>` → Bernoulli for genotype.

De novo records merge with `input_vcf` SVs and flow through the same depth-modulation path, so simulated coverage and the golden VCF round-trip behave identically for both sources. `compare-vcfs` skips both into the `skipped.symbolic` bucket — symbolic ALTs aren't byte-comparable.

Caveats:
- The bundled default is a literature-derived approximation, not a refit on the actual gnomAD VCFs. Don't use it for distribution-faithful benchmarking — retrain.
- `<INV>` and breakends are not generated (the read content for either isn't modeled yet).
- A trained model's `sv_model` field is `None` if the training VCF lacked sufficient SV observations (< 2 per type after filtering). Loading that model with `sv_rate_scale > 0` is harmless — generation just no-ops.

FASTQ Shuffling
===============
`eidolon` writes reads in contig order. To shuffle the output, use `seqkit shuffle` as a post-processing step:

```bash
# single-ended
seqkit shuffle sample.fastq.gz -o sample_shuffled.fastq.gz

# paired-ended (keeps mates in sync)
seqkit shuffle -2 sample_R1.fastq.gz sample_R2.fastq.gz \
    -o sample_R1_shuffled.fastq.gz -o sample_R2_shuffled.fastq.gz
```

`seqkit` uses reservoir sampling and streams from disk, keeping memory use bounded regardless of file size (`seqkit` is an open source toolkit for FASTQ files: https://github.com/shenwei356/seqkit).

3′ Adapter Readthrough
======================
Real Illumina libraries whose insert is shorter than the read length "read through" the
fragment into the 3′ sequencing adapter — the read's tail is adapter sequence, not genome.
`eidolon` can simulate this so short-insert data looks realistic and exercises downstream
adapter-trimming QC. **Disabled by default** — output is byte-identical to prior versions
when the `adapters` block is omitted.

```yaml
adapters:
  enabled: true
  preset: truseq        # truseq | nextera | custom
  # for preset: custom, give explicit 5′→3′ uppercase-ACGT sequences:
  # r1: AGATCGGAAGAGC...
  # r2: AGATCGGAAGAGC...
```

When enabled, any fragment whose insert is shorter than `read_len` is kept (instead of
being resampled away) and the read is padded to `read_len` at its 3′ end with adapter
sequence — R1 gets `r1`, R2 gets `r2` (appended after R2 is reverse-complemented, matching
real read orientation). Adapter bases carry the same quality/error models as genomic bases,
are marked soft-clipped (`S`) in the golden BAM, and never introduce variants. **How much
readthrough you get is controlled by how many inserts fall below `read_len`** — i.e. by
`fragment_mean` / `fragment_st_dev` relative to `read_len`. With a long-insert setting
(`fragment_mean` ≫ `read_len`) essentially no reads reach the adapter and enabling this
changes nothing.

Presets:
- `truseq` — Illumina TruSeq (`AGATCGGAAGAGC…`), the most common.
- `nextera` — Nextera / transposase (`CTGTCTCTTATACACATCT`).
- `custom` — supply your own `r1` / `r2`.

Not compatible with `long_reads: true` (adapter readthrough is a short-read phenomenon).

**When to enable it:**
- **Validating adapter-trimming pipelines** — generate reads with known adapter content to
  confirm fastp / cutadapt / Trim Galore settings actually detect and remove it.
- **Short-insert libraries** — small-RNA / cfDNA / degraded-DNA style data where inserts sit
  near or below the read length and readthrough is expected; the simulated reads then match
  what the sequencer really produces.
- **Aligner soft-clip benchmarking** — untrimmed adapter tails should be soft-clipped by the
  aligner; use this to check that behaviour end-to-end.
- **Realistic FASTQ→VCF tests** — exercise the whole trim → align → call path on data that
  carries adapters instead of idealized adapter-free reads.

To keep short inserts as **genomic** reads *without* modeling adapters (e.g. to study the
coverage behaviour of a short-insert library on its own), set `keep_short_fragments: true`
instead of enabling `adapters` — short fragments are then emitted as insert-length reads
with no adapter padding.

**Validated behaviour** (paired germline benchmark: chr22, 30×, TruSeq, 3 reps/arm,
BWA-MEM2 → GATK HaplotypeCaller vs the golden VCF): adapter readthrough is **detectable**
(fastp adapter content and aligner soft-clip fraction rise only when it is enabled) and
**fully callable** — SNP and indel recall/precision/F1 are unchanged within replication
noise whether the adapters are fastp-trimmed or left for BWA-MEM2 to soft-clip. An
adapter-free short-insert control showed the only fidelity difference versus a long-insert
baseline is the reduced effective coverage of short inserts (mate overlap), not the adapters.

Parallel Processing
===================
`eidolon gen-reads` processes contigs in parallel by default using rayon's work-stealing thread pool. Each contig is an independent unit of work — variant generation, fragment sampling, and FASTQ/BAM writing all happen concurrently across contigs — so references with many contigs scale well across cores.

Output is **byte-identical regardless of `num_threads`** (the same seed always produces the same reads in the same order), so you can change the thread count freely without affecting results.

A note on scaling: read generation is largely **memory-bandwidth bound**, so on references dominated by one or a few large chromosomes the wall-time gain from extra cores is limited — the bottleneck is memory throughput, not CPU. An experimental sub-contig **chunk size** knob (below) can split a single chromosome across cores, but it is **disabled by default** because in our benchmarks it did not improve wall time and occasionally regressed it.

**Thread count:**

By default `eidolon` uses all available logical cores. You can cap the thread count with the `num_threads` config key:

```yaml
# use 4 threads instead of all available cores
num_threads: 4
```
or disable parallelism entirely (useful for debugging or reproducibility testing)
```yaml
# disable parallelism entirely (useful for debugging or reproducibility testing)
num_threads: 1
```

Omit `num_threads` (or set it to `.`) to restore the default all-cores behaviour.

**Thread count and hardware — fewer threads can be faster:**

Read generation moves a lot of data relative to the arithmetic it does, so it is
largely **memory-bandwidth bound**. On a typical desktop or laptop (which has
only a couple of memory channels), a single `eidolon` thread can already saturate
much of the available memory bandwidth — so adding cores yields little speedup
and, past a point, can even run *slower* as threads contend for the memory bus.
In our desktop benchmarks, large references ran about as fast (sometimes faster)
on **1–4 threads** as on all 8.

Practical guidance:

- **Desktop / laptop:** try `num_threads: 1` (or a small number) and compare — it
  is often as fast or faster than all-cores for a single large run, and leaves
  cores free for other work. If you are simulating **many samples**, running
  several single-threaded `eidolon` jobs in parallel typically beats one
  many-threaded job.
- **HPC nodes** with many memory channels (and multiple sockets) have far more
  aggregate bandwidth, so higher `num_threads` scales better there.
- Output is byte-identical regardless of `num_threads`, so it is safe to tune the
  thread count purely for speed on your hardware.

**Chunk size (experimental, opt-in):**

By default each contig is one unit of work. Setting `chunk_size` splits contigs
into sub-contig chunks so a single large chromosome can be worked by several
cores at once. It is **off by default**: read generation is memory-bandwidth
bound, so chunking did not improve wall time in our benchmarks (and added a
little overhead). It is kept for CPU-bound or very-many-core scenarios where it
may help. The size is in base pairs and independent of the thread count, so
output stays byte-identical regardless of `num_threads`.

```yaml
# default — disabled (one chunk per whole contig); omit or set to . or 0
chunk_size: 0

# opt in: fixed chunk size in base pairs
chunk_size: 5000000
```

*When might it help?* Only when read generation is **CPU-bound** rather than
memory-bandwidth bound — the opposite of what we measured on a typical desktop.
Consider trying it when **both** are true:

1. Your reference is dominated by **one or a few large contigs**, so the default
   per-contig parallelism leaves most cores idle (e.g. a single-chromosome
   assembly, or a genome where one chromosome dwarfs the rest).
2. You have **memory bandwidth to spare relative to cores** — a multi-socket or
   many-memory-channel HPC node rather than a commodity desktop — and/or you run
   compute-heavier settings (GC-bias-weighted coverage, long-read mode, very high
   coverage) where per-read CPU work dominates memory traffic.

It is **not** worth enabling on a typical workstation, or on references with many
contigs (those already parallelize per contig). Always benchmark on your own
hardware: start with `chunk_size ≈ longest_contig_bp / (4 × cores)` (a few chunks
per core), compare wall time against the default, and keep it only if it is
actually faster.

**BAM output is fully parallel:**

`eidolon` uses a per-contig temp-file strategy: each contig worker writes its alignment records to a private temporary BAM body file, then a single concatenation pass assembles them in reference order into the final coordinate-sorted BAM.

**Reproducibility:**

Each contig's random number generator is derived deterministically from the parent seed and the contig's position in the reference, so output is identical across runs with the same seed even when the number of threads changes.

Reading Bed Data
================
`eidolon` can now read bed data and filter the reads based on regions specified. This comes with some caveats. First is that rust can't handle any header lines or non-header rows in the bed file, because of potential range of possibilities. Second, the names must match what is in the fasta. Third, `eidolon` can only filter, and does not use the rest of the bed information. It treats each record as a region of interest only.

Bed data will have some challenges. For example, if the contig names in the bed file don't match the assumed contig name from the fasta file, how will rust be able to know which is what? To make things work, we are going to assume, for now, that the bed file contig names match the names as derived by `eidolon`. To determine this:

First, in a linux terminal, use your input fasta to find the names of the contigs like this:

```bash
$ grep "^>" input.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
>etc
```
This will give you the full fasta name. `eidolon` determines the short name of the input fasta by skipping the initial '>' character, then taking the string up to the first whitespace delimiter. In the above example, the short names are "NC_001133.9", "NC_001134.8", etc. 

Next, we check the bed. This example awk command will look for unique values in the first column of your bed file and print out what it finds. 

```bash
$ awk '!seen[$1]++ {print $1}' file.bed
1
2
3
etc
```
In this case, the bed file uses a simple numbering scheme to number the chromosomes. We can see a mapping with the roman numeral chromosomes in the name, but it's tricky to handle consistently. So, the following command will allow you, the user, to map these chromosomes and thus instruct `eidolon`. Yours will vary based on the inputs.

```bash
awk -F'\t' '{$1=($1=="1"?"NC_001133.9":$1); $1=($1=="2"?"NC_001134.8":$1); print}' OFS='\t' file.bed > file_renamed.bed
```
You will have to tailor this command to your dataset, but this is one way to map the names from the bed to the fasta in a way `eidolon` will understand.

Filtering Your Data
====================

To run filter-reads, you must first copy the filter_reads_template.yml file from eidolon/template_config/ to a directory of your choosing. Then you can edit the file in your favorite editor. The configuration has four fields:

```bash
bed_file: /path/to/my.bed # required
```
A filename (full path is best) for a bed file with regions to target for filtering. It should be tab-separated in standard bed format.

In configuration:
```bash
files_to_filter: [
  # required, in list form
]
```
Filenames (full path is best) for files to filter by the bed file. For example, if you ran paired-ended fastqs with the vcf, your list would look like:

```bash
files_to_filter: [
   /my/home/output/bacteria_r1.fastq.gz,
   /my/home/output/bacteria_r2.fastq.gz,
   /my/home/output/bacteria.vcf.gz,
]
```
Note the full extension is important, but filter-reads should be able to handle both gzipped and unzipped files, if you decide to unzip them.

```bash
filter_key: .
```
This is the key appended to the filename, before extensions. For example, if your filename is "miscanthus_r1.fastq.gz", and you set the key as "_only_genes", the output file would be "miscanthus_r1_only_genes.fastq.gz". The default is "_filter".

```bash
overwrite_output: false # optional, default false
```
Set to `true` to allow filter-reads to overwrite existing output files. When `false`, the run will error if the output file already exists.

Once your files are entered and your config is saved, you can run eidolon:

```bash
$ eidolon filter-reads -c my_config.yml
```
Your output will give you the filenames and success status.

IMPORTANT: This feature only works on fastq files generated by `eidolon`. It is untested on generic VCF files, but should work. The reason it only works with `eidolon`-generated fastq files is that `eidolon` uses a naming scheme that the parser can use to identify where the read is located (enabling filtering), whereas an average fastq will not have that info.

Generating a Mutation Model
====================
`eidolon` now includes the ability to read in real data and learn the parameters to reproduce it in a simulation run. So far, this is limited to the mutation model, using the `eidolon gen-mut-model` subcommand.

```bash
$ eidolon gen-mut-model -c gen_mut_model_config.yml
```

The inputs are a single-sample VCF and the reference FASTA the VCF was called against. `eidolon` computes statistics for indels and SNPs and builds the trinucleotide model of SNP generation, as in the original `eidolon`. The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `mutation_model` config key.

Copy `template_config/gen_mut_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to the reference FASTA used to call the VCF
reference: /path/to/reference.fa

# required: path to a single-sample VCF file
# The VCF must include FORMAT and SAMPLE columns, and GT must be present in FORMAT.
# QUAL=. is accepted and treated as 0.
vcf_file: /path/to/variants.vcf

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/output_model.json.gz

# optional: restrict model learning to regions in this BED file
# set to . (dot) to use the entire reference (default)
bed_file: .

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# optional: custom 4x4 SNP transition matrix TSV (rows/columns: A C G T).
# A single header line is allowed. Diagonal values are zeroed automatically.
# Overrides the transition matrix inferred from VCF SNP data.
transition_matrix_file: /path/to/matrix.tsv
```

**VCF requirements:**
- Single sample only — multi-sample VCFs are not yet supported (tracked in #412)
- Each variant record must have `GT` in the `FORMAT` column; `eidolon` hard-errors if GT is missing
- `QUAL=.` is accepted and treated as quality score 0

Caveats: Only one sample can be read at this point (#412). Currently, high-mutation regions and common variants features from Python NEAT are not yet implemented (#413).

Try it out and let us know if you run into any issues!

We will continue to improve this process in the future. If you have suggestions for parsing names to facilitate this, let us know in the Issues section!

Generating a Sequencing Error Model
====================
`eidolon` can also learn a sequencing error model from real FASTQ data using the `eidolon gen-seq-error-model` subcommand. This reads per-base quality scores from the FASTQ to build a Markov quality score model and computes an average base error rate. Optionally, you can supply a BAM file or a custom TSV to set the SNP transition matrix (which base errors are most likely to become which other bases).

```bash
$ eidolon gen-seq-error-model -c gen_seq_error_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `seq_error_model` config key.

Copy `template_config/gen_seq_error_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to input FASTQ file (.fastq or .fastq.gz)
fastq_file: /path/to/reads.fastq.gz

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/seq_error_model.json.gz

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# optional: maximum number of reads to use for model learning; 0 = unlimited (default: 0)
max_reads: 0

# optional: quality score ASCII offset; 33 for Illumina 1.8+/Sanger, 64 for older Illumina (default: 33)
qual_offset: 33

# optional: list of Q-score bins to quantize the learned model to (e.g. NovaSeq 6000 uses
# [2, 12, 23, 37]). When set, each observed Q-score is snapped to its nearest bin before
# the model is built, and gen-reads will emit only these values. Omit for a continuous model.
# Bins must be < 94 and cannot include 31 (encodes to '@' under Phred+33).
# binned_quality_bins: [2, 12, 23, 37]

# optional: aligned BAM file to infer the SNP transition matrix from read-vs-reference mismatches.
# The BAM must have MD tags. If your aligner did not add them, run:
#   samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
# Ignored if transition_matrix_file is also set.
# If the BAM yields no mismatches at all, eidolon errors out rather than quietly using the
# default matrix. To use the default deliberately, omit this key.
bam_file: /path/to/aligned.bam

# optional: custom 4x4 SNP transition matrix TSV (rows/columns: A C G T).
# A single header line is allowed. Diagonal values are ignored.
# Takes precedence over bam_file.
transition_matrix_file: /path/to/matrix.tsv
```

**SNP transition matrix priority:**
1. `transition_matrix_file` (explicit TSV) — highest priority
2. `bam_file` (inferred from MD-tagged BAM mismatches)
3. Built-in default matrix (inherited from Python NEAT) — used when neither is provided

Supplying `bam_file` is a request to *fit* the matrix from data, so if the BAM yields no
read-vs-reference mismatches, eidolon **errors out** instead of falling back to the default. A
model that looks trained and is actually the default is indistinguishable from a trained one
downstream. **Omitting `bam_file` is how you ask for the default matrix.**

**BAM MD tag requirement:**
The BAM path requires MD tags to identify reference bases at mismatch positions. MD is an
*optional* SAM tag, so check your aligner: BWA-MEM writes it, `minimap2` needs `--MD`, and STAR
needs `--outSAMattributes MD`. If your BAM lacks them, generate them with:
```bash
samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
samtools index aligned_with_md.bam
```
Without MD tags the inference finds nothing, which is an error rather than a silent fallback.

**TSV format:**
The transition matrix TSV has 4 data rows (one per reference base: A, C, G, T) and 4 whitespace-separated columns (one per read base: A, C, G, T). The first row is skipped if it is non-numeric (treated as a header). Diagonal values are zeroed automatically. Rows are re-normalized to sum to 1.

```
A    C    G    T
0.0  0.5  0.3  0.2
0.5  0.0  0.3  0.2
0.4  0.3  0.0  0.3
0.3  0.3  0.4  0.0
```

**Binned quality scores:**
Modern Illumina platforms (NovaSeq 6000, NextSeq 1000/2000, the simplified MiSeq reporting modes) no longer emit the full continuous Q0–Q40 range. Instead, they quantize each per-base quality into a small set of discrete bins — for example, NovaSeq 6000 emits only `{2, 12, 23, 37}`. Set `binned_quality_bins` in the config to mimic this behaviour:

```yaml
binned_quality_bins: [2, 12, 23, 37]
```

When this field is set, `gen-seq-error-model` snaps each observed Q-score from the input FASTQ to its nearest bin (ties round down) before learning seed, transition, and global-frequency counts. The resulting model file is flagged `binned_scores: true` and lists the bin values as its `quality_score_options`. `gen-reads` then samples only those bin values when emitting reads — no extra flag needed on the read-generation side.

Validation rules:
- Bins must be non-empty integers in `[0, 94)` (entries are sorted and deduped automatically).
- Value `31` is rejected — under Phred+33 it encodes to `@`, which would corrupt FASTQ output. If you need a bin near Q31, pick 30 or 32.
- Bins with no observed counts in the training FASTQ are still kept in `quality_score_options`; their transition rows fall back to a uniform distribution, and a warning is logged listing the empty bins.

Some common platform bin sets (consult your sequencer's documentation for the authoritative list):

| Platform                       | Suggested bins     |
|--------------------------------|--------------------|
| NovaSeq 6000 (4-bin)           | `[2, 12, 23, 37]`  |
| NextSeq 1000/2000 (3-bin)      | `[2, 15, 35]`      |
| NextSeq 1000/2000 (4-bin)      | `[2, 12, 23, 37]`  |

Generating a GC Bias Model
====================
`eidolon` can learn a GC bias model from a reference FASTA and an aligned BAM file using the `eidolon gen-gc-bias-model` subcommand. It walks the BAM once to accumulate per-base reference coverage, tiles the reference in fixed-size windows, computes the GC% of each window, looks up the mean coverage over that window, and accumulates the data into a 101-bin weight table (one bin per integer GC percentage, 0–100%). The resulting model can be passed to `gen-reads` to make fragment start positions favour regions whose GC content matches the coverage bias observed in your data.

```bash
$ eidolon gen-gc-bias-model -c gen_gc_bias_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `gc_bias_model` config key.

Copy `template_config/gen_gc_bias_model.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: reference FASTA used to compute GC content
# Both plain and gzipped (.fa.gz) references are accepted.
reference: /path/to/reference.fa

# required: aligned BAM. Coverage is accumulated in process per reference base.
bam_file: /path/to/aligned.bam

# optional: minimum mapping quality applied while accumulating coverage.
# Default 0 matches `samtools depth` defaults.
min_mapq: 0

# optional: BED file restricting model inference to target regions
# Use "." (dot) to use the entire reference (default)
bed_file: .

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/gc_bias_model.json.gz

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# Window size in base pairs for GC% and coverage averaging.
# Should match the typical read length (or fragment length for long reads).
# Default: 100
window_size: 100

# Step between successive windows. Must not exceed window_size.
# Use the same value as window_size for non-overlapping windows (recommended).
# Values smaller than window_size produce overlapping windows.
# Default: window_size
window_stride: 100

# Bins with fewer windows than this receive a neutral weight of 1.0 rather than
# a learned weight. Increase this to require more evidence before trusting a bin.
# Default: 10
min_windows_per_bin: 10
```

**BAM filtering:** unmapped, secondary, and supplementary records are skipped automatically. Reads with mapping quality `<= min_mapq` are dropped before contributing to coverage. The default `min_mapq: 0` reproduces `samtools depth` defaults; set it to 20 if your downstream `gen-reads` runs assume mapq-filtered coverage.

**Large genomes:** The tool walks the BAM once and allocates a per-contig depth array (`u32` per reference base) only for contigs that receive at least one record. Peak memory for hg38 at 30× is on the order of all-contigs-summed depth arrays — see the HPC section below for sizing.

**Long reads:** Set `window_size` to approximately your typical read length (e.g. 5000–50000 for ONT or PacBio). Larger windows mean fewer windows per contig and faster model building. Non-overlapping windows (`window_stride: window_size`) are recommended so each observation is independent.

**If the model has no effect in gen-reads:** If every GC% bin in your reference has fewer observations than `min_windows_per_bin`, all weights will be neutral (1.0) and no bias will be applied. This is logged as a warning. To diagnose, lower `min_windows_per_bin` or increase the region covered (use a larger reference or remove the BED restriction).

**Using the model in gen-reads:**

```yaml
# in gen_reads_config.yml
gc_bias_model: /path/to/gc_bias_model.json.gz

# optional: inflate fragment count to compensate for low-weight regions
# true (default): total coverage stays close to the requested depth
# false: coverage in low-GC regions will be lower than requested
gc_bias_normalize_coverage: true
```

Generating a Fragment Length Model
====================
`eidolon` can learn a fragment length model from real paired-end alignment data using the `eidolon gen-frag-length-model` subcommand. It reads a BAM or SAM file, collects the template lengths (TLEN) of confidently-mapped concordant read pairs, filters out rare and extreme-outlier lengths, then fits a normal distribution to the result and writes a `FragmentLengthModel` that `gen-reads` can use.

```bash
$ eidolon gen-frag-length-model -c gen_frag_length_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `fragment_model` config key.

Copy `template_config/gen_frag_length_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to input BAM or SAM file
# The file must contain paired-end aligned reads. BAM is strongly recommended.
input_file: /path/to/aligned.bam

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/fragment_length_model.json.gz

# optional: minimum number of reads for a fragment length to be included in the model
# Default is 2 to handle smaller datasets. Set to 0 to disable filtering entirely.
# For larger datasets (e.g. whole-genome), try 100 and adjust from there.
min_reads: 2

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false
```

**Read filtering rules:**  
Only reads that satisfy all of the following are used:
- Paired-end and first in pair
- Mapped (not unmapped, secondary, or supplementary)
- Mate mapped to the same reference sequence
- Mapping quality > 10

**Outlier filtering:**  
Fragment lengths that exceed `median + 10 × MAD` (median absolute deviation) are removed before fitting. This mirrors the filtering in Python NEAT's fragment length modeler. If no lengths survive the filter, lower `min_reads` or set it to `0`.

Building multiple BAM-derived models in one pass
====================
If you need both a fragment length model and a GC bias model from the same BAM, `eidolon gen-bam-models` walks the BAM once and feeds both observers from the single pass — half the I/O of running the two per-tool commands back-to-back.

```bash
$ eidolon gen-bam-models -c gen_bam_models_config.yml
```

Copy `template_config/gen_bam_models.yml` and fill in the sections for whichever models you want:

```yaml
bam_file: /path/to/aligned.bam
min_mapq: 0   # walker-level filter; FragLengthObserver enforces MAPQ > 10 internally

frag_length:
  output_file: /path/to/frag_length_model.json.gz
  overwrite_output: false
  min_reads: 2

gc_bias:
  reference: /path/to/reference.fa
  output_file: /path/to/gc_bias_model.json.gz
  overwrite_output: false
  bed_file: .
  window_size: 100
  window_stride: 100
  min_windows_per_bin: 10
```

Both sections are optional but at least one must be present. The output models are interchangeable with what `gen-frag-length-model` and `gen-gc-bias-model` would produce on their own, so the resulting files plug into `gen-reads` exactly the same way.

**One filter for the shared walk.** The top-level `min_mapq` is the only MAPQ knob exposed to the unified walker, and it gates *both* observers. The standalone tools differ: `gen-frag-length-model` hard-codes its own MAPQ > 10 cutoff (via `BamWalkFilter::for_frag_length()`), while `gen-gc-bias-model` honors the `min_mapq` you pass it. If you want the unified runner's frag-length output to match the standalone command byte-for-byte, set `min_mapq: 10`; if you want the GC bias output to match a standalone run with `min_mapq: 20`, set `min_mapq: 20`. When the two policies need to differ, run `gen-bam-models` twice — one config per output — and the per-observer BAM iteration savings still beat running the standalone commands.

Comparing a downstream caller's VCF against the golden VCF
==========================================================
After running `gen-reads`, you can validate a downstream variant caller against the simulated truth using `eidolon compare-vcfs`. The subcommand classifies every variant into TP / FN / FP, catches denotation-different alternates that exact matching would mis-classify as both an FN and an FP (e.g., a left-aligned indel in the called VCF versus the same indel right-aligned in the golden), and attributes the surviving FNs to specific reasons drawn from the simulator's configuration.

```bash
$ eidolon compare-vcfs -c compare_vcfs_config.yml
```

Copy `template_config/compare_vcfs_template.yml` and fill in:

```yaml
golden_vcf: /path/to/golden.vcf.gz
called_vcf: /path/to/called.vcf.gz
reference:  /path/to/reference.fa
output_dir: /path/to/report_dir
overwrite_output: false

# Optional: restrict comparison to these regions.
target_bed: .
# Optional: comma-separated list. "." treats every called contig as simulated.
contigs_simulated: .
# Optional: BED used by the simulator. Drives FN attribution.
mutation_bed: .
# Optional: TSV remapping BED chrom names (e.g., 1<TAB>chr1).
chrom_aliases: .

# Filters
include_homs: false       # ALT == REF rows
include_filtered: false   # FILTER != PASS / "."

# Equivalence sweep (NEAT 2.1 algorithm, ported)
equivalence_window: 50    # ±N bp around each FP
fast: false               # skip the sweep entirely

# Outputs
write_fp_vcf: false       # also produce FP.vcf
```

The run produces four files under `output_dir`:

- `comparison_summary.json` — schema-versioned (currently `1.3.0`), machine-readable. Includes TP/FN/FP totals + per-contig breakdown, precision / recall / F1, the FN attribution roll-up (`outside_simulated_contigs`, `outside_mutation_bed`, `outside_target_bed`, `unknown`), per-VCF skip counters (`multiallelic`, `homozygous_ref`, `filtered`, `outside_target_bed`, `outside_simulated_contigs`, `symbolic`), and any chrom-naming-mismatch warnings. Symbolic / structural ALTs (`<DEL>`, `<DUP>`, `<CNV>`, ...) are byte-incomparable, so they're counted into `skipped.symbolic` and excluded from TP/FN/FP classification.
- `comparison_summary.txt` — same content, human-readable.
- `FN_with_reasons.vcf` — every surviving FN, annotated with a `EIDOLON_REASON` INFO tag listing the attribution reasons.
- `FP.vcf` (optional) — every surviving FP, as-is. Off by default; enable with `write_fp_vcf: true`.

**Equivalence detection.** For each false-positive variant, `compare-vcfs` takes a ±`equivalence_window` bp window of the reference and applies both the FP set and the FN set within that window. If the resulting byte sequences are identical, the two sets are alternative spellings of the same edit and every consumed FN is promoted to TP. Set `fast: true` to skip this pass; the report's `totals.equivalents_promoted` counts how many TPs were rescued by it.

**FN attribution.** Each surviving FN is tagged with at least one reason. If the FN's contig wasn't in `contigs_simulated`, that reason is reported alone (the BED checks are skipped because they presuppose the contig was simulated). Otherwise the configured BEDs are checked and `unknown` is the fallback. Aliases configured via `chrom_aliases` are applied to BED chrom names at load time, so you can compare against a reference that uses one convention (e.g., `chr1`) when your BED uses another (e.g., `1`). When a BED has any chrom that doesn't appear in the reference's FASTA contig list (post-alias), `compare-vcfs` surfaces a warning in the report — useful for catching single-typo BEDs that would otherwise misattribute every variant on that contig.

Running on HPC
==============

`eidolon` runs as a single process and fits naturally onto a single HPC compute node. No MPI, distributed computing, or special environment setup is required. The notes below cover resource budgets for whole-genome human-scale runs.

### gen-reads

`gen-reads` is already multi-threaded via rayon (see [Parallel Processing](#parallel-processing) above). Set `num_threads` in your config to match your CPU allocation:

```yaml
num_threads: 16
```

Memory scales with `num_threads × largest_contig_size`, not total genome size — each thread processes one contig at a time. For hg38, chr1 is ~249 Mbp, so 16 simultaneous threads need roughly 16 GB for sequence data plus overhead. Budget additional scratch space for per-contig temp files before the final assembly: roughly 3× the expected output size.

| Genome | Threads | Recommended RAM | Scratch space | Typical wall time |
|--------|---------|-----------------|---------------|-------------------|
| Bacterial (~5 Mbp) | 4 | 4 GB | 5 GB | < 1 min |
| Human hg38, 10× SE | 16 | 32 GB | 100 GB | ~20 min |
| Human hg38, 30× PE | 16 | 48 GB | 300 GB | ~60 min |

Example SLURM header:

```bash
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --tmp=300G
```

### gen-gc-bias-model

Single-threaded; walks the BAM once and accumulates per-base reference coverage as a `Vec<u32>` per observed contig. Peak memory is roughly `4 bytes × total reference bases that received at least one record`. For a full human genome (~3.2 Gbp) this is ~13 GB. For a single chromosome or a targeted/exome BAM it is much smaller.

If memory is tight on whole-genome runs, split the BAM by chromosome and run `gen-gc-bias-model` once per chromosome on a representative subset that spans your GC range of interest, then average or combine the resulting models.

Recommended SLURM header (whole-genome BAM):

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
```

### gen-mut-model

Single-threaded; processes the VCF one chromosome at a time. A full human genome WGS VCF (~10 GB) typically completes in 30–60 minutes with peak RSS under 4 GB.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:30:00
```

### gen-seq-error-model

Single-threaded; streams through the FASTQ. For a full genome FASTQ (~600 M reads), expect 30–60 minutes. For most purposes, training on a representative subset produces a statistically equivalent model in a fraction of the time:

```yaml
max_reads: 5000000   # 5 M reads is sufficient; set to 0 for unlimited
```

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00   # with max_reads=5M; scale up for unlimited
```

### gen-frag-length-model

Single-threaded; streams TLEN fields from the BAM. A full human genome BAM at 30× typically completes in under 15 minutes with peak RSS under 2 GB.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
```

### Environment notes

- The `eidolon` binary has no runtime dependencies beyond a standard C library (glibc), which is present on all Linux HPC systems.
- No module loads or conda environments are required — copy the release binary to your scratch or project directory and run it directly.
- If you compile from source on the cluster, ensure `cmake` is available (`module load cmake` on most systems) before running `cargo build --release`.

**Current Benchmarks**
We ran some benchmarks on `eidolon` on my home desktop, a very basic Linux desktop, using data pulled from public resources on the internet (mostly ncbi). The "yeast" is brewer's yeast. These are the results:

```text
[16:26:39] Build complete.
eidolon Benchmark Report
Run:        20260510_162639
Machine:    pop-os  |  8 logical CPUs
Coverage:   10x    |  Read length: 151 bp
Binary:     /home/joshfactorial/code/eidolon/target/release/eidolon
──────────────────────────────────────────────────────────────────────
SECTION 1: Single-ended read generation
Genome          Size(MB)   Wall time  Peak RSS(MB)      CPU%
──────────────────────────────────────────────────────────────────────
[16:26:40] Single-ended: ecoli  (4.6 MB)
ecoli                4.6     0:14.83         416.9       99%
[16:26:54]   Done: wall=0:14.83  peak_rss=416.9 MB  cpu=99%
[16:26:54] Single-ended: pneumonia  (2.1 MB)
pneumonia            2.1     0:05.39         197.0       99%
[16:27:00]   Done: wall=0:05.39  peak_rss=197.0 MB  cpu=99%
[16:27:00] Single-ended: yeast  (12 MB)
yeast                 12     0:10.94         981.7      352%
[16:27:11]   Done: wall=0:10.94  peak_rss=981.7 MB  cpu=352%
[16:27:11] Single-ended: c_elegans  (97 MB)
c_elegans             97     7:35.23        8591.2      387%
[16:34:46]   Done: wall=7:35.23  peak_rss=8591.2 MB  cpu=387%
──────────────────────────────────────────────────────────────────────
Section 1 — Passed: 4 / 4  |  Failed: 0
──────────────────────────────────────────────────────────────────────
SECTION 2: Paired vs single comparison  (genome: yeast, 10x)
  Paired-end uses fragment_mean=400, fragment_st_dev=50
  Same reference, same coverage, same seed — only the read-generation mode differs.
Mode               Wall time  Peak RSS(MB)      CPU%
──────────────────────────────────────────────────────────────────────
[16:34:46] Paired comparison — single-ended pass
single-ended         0:10.48         970.4      365%
[16:34:57]   Single done: wall=0:10.48  peak_rss=970.4 MB  cpu=365%
[16:34:57] Paired comparison — paired-ended pass
paired-ended         0:13.30        1126.8      354%
[16:35:10]   Paired done: wall=0:13.30  peak_rss=1126.8 MB  cpu=354%
──────────────────────────────────────────────────────────────────────
  Paired/single wall-time ratio: 1.27x  (10.48s → 13.30s)
  Paired/single peak-RSS ratio:  1.16x  (970.4 MB → 1126.8 MB)
──────────────────────────────────────────────────────────────────────
```
