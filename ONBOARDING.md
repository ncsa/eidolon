# eidolon for NEAT users

Written for people who already know NEAT — what `gen_reads` does, why you train models from your
own data, what a golden BAM is for. This is the **delta**, not an introduction.

Short version: eidolon is NEAT reimplemented in Rust, tracking the NEAT feature set, plus a
native tumour/normal cancer workflow that has no NEAT counterpart. The conceptual model is
unchanged — reference in, models trained from real data, FASTQ + golden BAM + truth VCF out.

Sections marked *Deeper* are skippable.

---

## What is the same

If you've tried the latest versions of NEAT, you'll find the user interface very familiar. 
* Subcommand shape (`gen-reads`, `gen-mut-model`, `gen-seq-error-model`, `gen-frag-length-model`, 
`gen-gc-bias-model`)
* YAML config
* models trained once, stored, and reused as input
* BED targeting
* variant insertion from an input VCF
* golden BAM + truth VCF alongside the FASTQ, single or paired end.

The outputs are very similar, with some added features in the VCF output, and potential for much greater complexity.

## What is different in practice

| Change | Practical consequence |
|---|---|
| **~10–13× faster, 3–7× less memory** | chr22 at 10×: 61 s / 227 MB against NEAT 4.6.1's 810 s / 1525 MB, single thread. Shared-allocation runs stop being painful |
| **Deterministic for a given seed at any thread count** | Byte-identical output. NEAT's varies with `--threads`, so a rerun is not a rerun |
| **Output token names changed at v3.0.0** | `EIDOLON_*`, not `NEAT_*` / `RNEAT_*`. **Breaks any script that parses eidolon VCFs or read names** — see "Upgrading from 2.0.0" in the README |
| **Quality binning is explicit** | You spell out `binned_quality_bins`. NEAT 4 has named presets (`--quality-preset novaseq`); eidolon does not, which is a real ergonomic loss |
| **Default model Ts/Tv is 2.21** | NEAT 4's is 2.33. Default model only — a trained model uses your data ([#410](https://github.com/ncsa/eidolon/issues/410)) |
| **Structural variants are actually generated** | NEAT 4.6.1 ships `Inversion`/`Duplication`/`Translocation`/`CopyNumberVariant` classes that are exported from nothing and referenced nowhere; `VariantTypes.types` lists only SNV/insertion/deletion/unknown. It also cannot render an SV from an input VCF — symbolic ALTs become `UnknownVariant`, whose accessors raise. Verified against tag `4.6.1` |
| **SemVer since v3.0.0** | With a declared public API, so a MAJOR bump is a real signal |

## What is new

Cancer work is where the original effort went.

- **`gen-cancer-reads`** — two `gen-reads` passes (normal-genotype and tumour) over one reference,
  merged at a configurable purity into a single biopsy FASTQ that Mutect2 / Strelka / Manta
  consume directly, plus a truth VCF tagged `INFO/EIDOLON_ORIGIN ∈ {germline, somatic, shared}`
  so you can score germline and somatic separately.
- **Intra-tumour heterogeneity** — a tumour as subclones at distinct cancer-cell fractions,
  specified inline, fitted from PyClone-VI / DPClust output, or replayed from a real tumour's
  observed VAF. Per-site VAF is `purity × dosage × CCF`, and the truth carries `EIDOLON_CCF` and
  `EIDOLON_VAF`. Validated on SEQC2 HCC1395 at r = 0.99 against intended VAF.
- **Per-tissue somatic models** — bundled pan-cancer plus BRCA / skin / lung, fitted from
  COSMIC and PCAWG.
- **Continuous per-variant allele fraction** from an input VCF's `AF`/`AD`, rather than genotype
  `{0.5, 1.0}`. This is what makes pooled and somatic AF spectra reproducible.
- **Trinucleotide-context-aware SNP placement**, so SBS-96 signatures reproduce (cosine 0.99). Originally 
  trimmed from `rneat` for simplicity and because it's effects were difficult to detect, we discovered this 
  came into play in `eidolon` while looking for cancer signals, as cancer detection looks for high-SNV regions, 
  which were missing in the fully random placement. Signals began to show up once the placement became 
  context-aware again.
- **Sequencing-error nucleotide substitution matrix trainable from data** — from a BAM's MD tags, or a custom 
  4×4 TSV. NEAT 4.6.1 hardcodes this matrix (model_sequencing_error/runner.py, with its own TODO incorporate 
  these into the calculations); its trainable Markov model covers quality-score transitions from FASTQ, 
  a different axis.
- **`compare-af`** (per-allele AF correlation, truth vs simulated) and **`validate`** (checks an
  emitted file against what downstream tools actually accept, before you spend a pipeline run
  finding out).

## What to trust, and what not to

The part worth reading before you design an experiment.

**[`docs/sv_support_matrix.md`](docs/sv_support_matrix.md) is the authority on structural
variants** — every row is a measurement against read-level evidence, and the broken cells are
pinned by tests rather than quietly omitted. Current state worth knowing:

- **Insertions above ~a read length are not fully realised** and the engine now refuses to plant
  them de novo rather than write a truth record the reads cannot support
  ([#516](https://github.com/ncsa/eidolon/issues/516)). Consequence: the realised INS rate sits
  below the model's, and any insertion benchmark is limited to short insertions for now.
- **DEL / DUP / CNV coverage runs ~8% high** whenever the copy-number multiplier is not 0 or 1
  ([#499](https://github.com/ncsa/eidolon/issues/499)) — fine for detection, not for dosage
  estimation.
- **Symbolic `<INS>` from an input VCF is a silent no-op**, and a single unpaired breakend
  destroys local coverage ([#500](https://github.com/ncsa/eidolon/issues/500)).
- BND, INV and literal indels are measured good.

Caller recall figures from the HPC tier live in `docs/access_report_draft.md`. Treat its section
caveats as load-bearing — several sections say explicitly that they are not yet submittable and
why.

> ### Deeper: how the project decides something is verified *(skip freely)*
>
> A simulator that generates *wrong* truth data is worse than one that crashes: a crash you
> notice, wrong truth silently invalidates every benchmark built on it. So the standard is that
> nothing is "done" until there is evidence, and reports state which kind of evidence exists —
> "known-answer fixture, not yet real data" is a normal thing to read here.
>
> In practice two habits do the work. **Assert content, not existence**: a check that a file
> appeared passes just as happily when the file is garbage, and one here confirmed a non-zero
> count while every value was the malformed string `AF=AF=0.3000`. **Prove a test can fail**:
> deliberately break the code it covers and confirm it notices. That has repeatedly exposed
> thorough-looking tests asserting nothing.
>
> `docs/claude_engineering_audit.md` collects the case histories, including the ones where the
> tooling misled us for weeks. Unusually candid for a project document, and the best single read
> for how this codebase thinks.

> ### Deeper: the HPC validation tier *(skip freely)*
>
> Unit tests prove the machinery runs; they cannot prove genome-scale output resembles real data.
> So `scripts/delta/` simulates a whole tumour genome on NCSA Delta, runs production callers
> (Manta, Delly, GATK CNV) and scores them against the truth with truvari.
>
> It carries its own controls, which is the interesting part: the truth is scored against
> **itself** (must be perfect) and against a **deliberately displaced copy** (must find nothing),
> because a matching configuration loose enough to accept anything would otherwise report
> excellent recall. Findings feed back into the SV support matrix.

## Reading the outputs

Truth-VCF INFO tags you will actually encounter:

| Tag | Meaning |
|---|---|
| `EIDOLON_PROVENANCE` | `denovo` (sampled from a model) vs supplied via `input_vcf` |
| `EIDOLON_ORIGIN` | `germline` \| `somatic` \| `shared` — cancer runs only |
| `EIDOLON_CCF` | Cancer-cell fraction of the subclone carrying this variant |
| `EIDOLON_VAF` | Intended observed VAF after purity mixing |
| `EIDOLON_REASON` | Written by `compare-vcfs` on FN/FP records, explaining the classification |

Read names carry `EIDOLON_chimeric` for reads spanning a structural-variant junction, which is
how you find them without an aligner.

## Where things live

| Path | What is in it |
|---|---|
| `docs/cancer_howto.md` | **Start here** — copy-paste configs with worked examples |
| `docs/sv_support_matrix.md` | What is measured, broken, and unverified |
| `docs/cancer_simulator.md` | Cancer design rationale and calibration caveats |
| `docs/hpc_guide.md` | Running on a cluster |
| `eidolon/` | CLI: argument parsing, per-subcommand runners |
| `eidolon-core/` | Engine: models, variant types, readers/writers |
| `scripts/delta/` | The HPC validation harness (bash) |
| `CLAUDE.md` | Working practices — written for AI agents, useful to anyone |

## Getting it

```bash
conda install -c bioconda eidolon        # or cargo build --release
eidolon --help
```

## Working on it

PRs target **`develop`**; `main` holds releases. GitHub issues track bugs and planned work, and
there is a **Feedback** issue type for things that are not quite bugs — a confusing message, a
rough edge, or something that worked well. [`CONTRIBUTING.md`](CONTRIBUTING.md) has the testing
bar if you are writing code.

The project was `rusty-neat` / `rneat` before v2.0.0, so older notes, links and issue titles use
those names.

## Questions

Open an issue or message **Joshua Allen** (links on the repo page). If something was unclear here,
that is a docs bug — say so and it gets fixed.
