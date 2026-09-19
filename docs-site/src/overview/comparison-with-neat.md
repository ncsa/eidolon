# How `eidolon` compares to NEAT

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
