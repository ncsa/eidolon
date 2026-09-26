# Shipped default models

These are the models `gen-reads` falls back on when a config supplies none. They are
compiled into the binary with `include_bytes!`, so they ship with every release.

After analysis of real sequencing data, we have updated the default error model for
`eidolon`. This file records the source of each value, for future reference and for
repeatability.

## Status

| model / parameter | source | measured |
|---|---|---|
| `default_fragment_length_model.json.gz` | GIAB HG002 2x250, fitted | yes — see below |
| `error_rate` | GIAB HG002 2x250, fitted | yes — 0.003774, measured |
| `indel_probability` | NEAT2 static default | no |
| `insertion_fraction` | NEAT2 static default | yes — confirmed at 0.387 |
| indel-error lengths | HCC1395 normal | yes |
| homopolymer context curve | HCC1395 normal | yes |
| `default_sequencing_error_model.json.gz` | GIAB HG002 2x250, fitted | yes — see below |
| `default_mutation_model.json.gz` | GIAB HG002 v4.2.1 truth VCF, fitted | yes — see below |
| `default_indel_model.json.gz` | the mutation default's own indel model | yes, with it |
| `default_trinuc_model.json.gz` | the mutation default's own trinucleotide model | yes, with it |
| `default_mutation_model_bkup.json.gz` | NEAT2 `MutModel_NA12878.p.gz`, first conversion | not loaded; provenance only |

## `default_fragment_length_model.json.gz`

| | |
|---|---|
| **Source** | GIAB HG002, `NIST_Illumina_2x250bps`, aligned to GRCh38 (21x) — the library the quality model comes from |
| **Pairs used** | 61,588,264; 0.19% trimmed as outliers (discordant and chimeric pairs) |
| **Built with** | `eidolon gen-frag-length-model`, `min_reads: 100`, `distribution: discrete` |
| **Fitted by** | `scripts/delta/fit_hg002_defaults.sbatch`, jobs 22443316 / 22444161 (#752) |

Shape: 982 bins over 3–984 bp, no gaps. Mean 408.5, sd 93.2, skew +0.207,
p05/p50/p95/p99 = 262/404/567/649.

**Checked against the BAM with samtools**, not eidolon (`scripts/delta/validate_frag_model.sh`,
same filter as the builder): over the model's support the BAM reads mean 408.55, sd 93.19,
skew +0.207, and the model agrees to 0.00% on mean, sd and p99. That checks the builder. It is
not a held-out check: the model is fitted from the whole genome, so there is no unseen
chromosome to test it on.

At 250 bp reads, 4.71% of the mass falls below the sampler's `read_len + 10` floor.

**Chemistry is older than the model it replaced.** The previous default came from HCC1395's
NovaSeq normal (mean 431.8, sd 112.3, skew +0.528). This one is HiSeq 2500. It was chosen so
the shipped defaults describe one library rather than the individually newest source per
component (#752).

## Sequencing error model

Built inline in `eidolon-core/src/models/sequencing_error_model.rs` rather than as a data
file. It ships with every run that supplies no model of its own, so it gets the same
accounting as the files above.

### `error_rate` (0.006638164688495656)

Fitted, unlike the constants below. It is the `avgError` of NEAT2's bundled
`errorModel_toy.p`, computed from the sequencing data that model was built on. The
originating sample is not recorded upstream.

**It is a summary, not a setting.** `gen-reads` never reads it. Sequencing errors are
injected per base from that base's own quality score (`convert_score`, `10^(-q/10)`), so
the quality model is what determines the rate and this field describes it. Measured: three
models differing only in this field, including `0.0` and `0.5`, generate byte-identical
FASTQ.

It is recorded because it makes models comparable — the shipped 0.006638 against HG002 R1
at 0.002261 and R2 at 0.004454 (#695) is a statement about three quality models. **To
change how many errors a run produces, change the quality model**; there is no runtime
scale factor, and #725 scopes what one would have to be.

### Inherited constants (#660)

We initially tried to match NEAT2 as closely as possible. NEAT2 shipped the following
defaults hardcoded in its sequencing error model.

| parameter | value | NEAT2 name |
|---|---|---|
| `indel_probability` | 0.01 | `SIE_RATE` — odds a sequencing error is an indel |
| `insertion_fraction` | 0.4 | `SIE_INS_FREQ` — odds such an indel is an insertion |
| insertion base composition | uniform over ACGT | `SIE_INS_NUCL` |
| substitution transitions | 0.4918 / 0.3377 / 0.1705 … | `SSE_PROB` — now only the fallback for a fit without `bam_file`; the shipped matrix is fitted (below) |

**Two of these were initially mistranslated in the Rust port.** The insertion fraction was
used as the indel rate, the real indel rate was dropped, and the insertion split was
replaced by a hardcoded `0.5` — about **40x too many** sequencing errors made into indels.
Both were restored to their source values in #660.

`insertion_fraction` is confirmed by measurement: 668 of 1,726 low-support indels in
HCC1395 normal are insertions, a fraction of **0.387**.

`indel_probability` has not been measured. On Illumina data the indel error rate is around
1e-5/base; at Q35 this constant gives ~3.2e-6. Changing it needs its own measurement.

### Indel-error lengths

| | |
|---|---|
| **Source** | HCC1395 matched **normal**, SEQC2 Somatic Mutation WG reference sample |
| **Reference** | GRCh38, chr20 + chr21 + chr22 |
| **Region** | the ten 400 kb loci a realism-panel run placed (`realism_21795898/regions.bed`) |
| **Events** | 1,726 low-support indels — 1,058 deletions and 668 insertions |
| **Measured by** | `scripts/delta/indel_context.sbatch`, Delta job 21801707 |

| \|len\| | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | tail |
|---|---|---|---|---|---|---|---|---|---|---|---|
| deletions (n) | 762 | 135 | 42 | 47 | 11 | 8 | 5 | 9 | 5 | 6 | 11,12,13,15,16,17,19,20,21,22,23,25,27,34,38,45 |
| insertions (n) | 426 | 120 | 22 | 55 | 13 | 8 | 1 | 6 | 1 | 5 | 11,12,13,15,18,19,22,27,30 |

Indels were split from variants by support fraction: below 10% of local depth is slippage,
at or above 25% is a variant. Only the **low-support** side feeds this model; variant indel
size is a different population and belongs to placement.

Deletions: n = 1,058 over 26 bins. Insertions: n = 668 over 19 bins. The two arms did not
observe the same lengths, so each carries its own value list.

The previous distribution was a first-order approximation — `[0.999, 0.001]` over lengths
`[1, 2]` — which produced indel errors of a single base 99.9% of the time. The data showed
a significant (**16.4%**) number of slippage events at 3 bp or more, and a small but
measurable number (**0.70%**) at 20 bp or more. We built the new baseline from this data
set, as its error properties were more robust.

Deletions and insertions are not the same shape — **72.0%** of deletions are 1 bp against
**63.8%** of insertions — so we separated the distributions, similar to how the mutation
model treats insertions and deletions as independent events.

### Homopolymer context curve (#661)

We added a further refinement based on our findings of concentrated indels in homopolymer
regions. `indel_context_curve` scales `indel_probability` by the length of the homopolymer
run the base sits in, rather than spreading it uniformly.

| | |
|---|---|
| **Source** | HCC1395 matched **normal**, SEQC2 Somatic Mutation WG reference sample |
| **Reference** | GRCh38, chr20 + chr21 + chr22 at 46x |
| **Background** | 3,999,990 reference bases, exact (N runs excluded) |
| **Events** | 1,726 slippage errors |
| **Measured by** | `scripts/delta/indel_context.sbatch`, Delta job 21674484 |
| **Raw data** | `/projects/bhrd/jallen17/eidolon-access-results/indelctx/job_21674484/` |

| run | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | ≥10 |
|---|---|---|---|---|---|---|---|---|---|---|
| propensity | 0.64 | 0.76 | 0.82 | 1.11 | 1.58 | 1.84 | 5.64 | 12.16 | 24.24 | 39.20 |

**Each entry is a normalized enrichment** — the share of indel errors at that run length
divided by the share of reference bases at that run length. That makes the curve
1.0-centered by construction over its human background, so applying it *redistributes*
`indel_probability` rather than raising it: the genome-wide total on human is unchanged.
On a reference with different homopolymer composition the realized total moves with that
composition, which is the intended behavior — measured at **0.745x** on the 4.6 Mb E. coli
fixture and **0.734x** for an idealized 50% GC random sequence. A genome with fewer
homopolymers shows less slippage overall.

This is the sequencing-error curve. Variants carry their own, steeper propensity — 60.44x
at runs ≥10 against 39.20x here — which belongs to variant placement (#378).

## `default_sequencing_error_model.json.gz`

The shipped default, as of v3.4.0. `SequencingErrorModel::default()` deserializes this file
whole; `QualityScoreModel::default()` takes its R1 half, so the two cannot drift.

| | |
|---|---|
| **Source** | GIAB HG002, `NIST_Illumina_2x250bps`, chunk `L001:001` |
| **URL** | `ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/` |
| **Instrument** | HiSeq 2500. **Continuous quality scoring, not binned.** |
| **Reads fitted** | 3,391,610 per mate, taken 1-in-10 across the file |
| **Degraded cut** | Q<25 averaged over the last 50 bases |
| **Fitted by** | `scripts/delta/fit_hg002_pair.sbatch`, job 22233888 |
| **Substitution matrix** | the HG002 BAM with its GIAB v4.2.1 variants masked: `scripts/delta/fit_hg002_defaults.sbatch`, job 22447983 |

### Substitution matrix

Fitted from the same library's aligned reads by `gen-seq-error-model` with `bam_file:` and
`known_variants_vcf:` (#752). The refit reproduces job 22233888's `error_rate` and degraded
fractions exactly, so the matrix is the only field that differs from that fit.

| from → to | A | C | G | T |
|---|---|---|---|---|
| **A** | — | 0.392 | 0.310 | 0.299 |
| **C** | 0.437 | — | 0.227 | 0.337 |
| **G** | 0.344 | 0.216 | — | 0.440 |
| **T** | 0.295 | 0.311 | 0.395 | — |

**Masking.** 50,240,073 of 418,177,394 mismatches (12.01%) sat at the sample's own variant
positions and were left out. They were HG002's variants, not sequencing errors: masking them
lowered the transition share of every row by about 0.04, the direction a Ti/Tv 2.1 germline
contribution predicts.

**Strand symmetry, which the fit is not told about.** Reads align to both strands, so an
error profile should equal its own reverse complement. Every pair agrees to within 0.007
(A→C 0.392 / T→G 0.395, C→A 0.437 / G→T 0.440), and a test holds it to 0.015. The largest
substitutions are C→A and G→T.

**Not masked:** HG002 variants outside GIAB's high-confidence regions, which the truth VCF
does not list. Their share is unmeasured; measuring it needs the BAM mismatch count
restricted to a BED, which the fitter does not yet support.

The matrix it replaces is NEAT2's `SSE_PROB` (see above). That one is not strand symmetric,
for example C→T 0.210 against G→A 0.375, and it remains the fallback for a fit without a BAM.

### Shape

250 bp reads; 31 observed scores spanning Q2–Q40, non-contiguous (gaps of 1, 2 and 8); a
249-position transition tensor. Both mates are present, each with its own degraded
population: R1 `read_fraction` 0.1102, R2 0.2598. Fitted `error_rate` 0.003774.

R2 is measurably worse than R1 — 2.07x on fitted error rate, 2.36x on degraded fraction —
which is why the model carries the two separately (#723).

### What has been measured

The two-population fit reproduces the library it was fitted from. Generated against the same
subsample, with a one-population fit of the same reads by the same binary as the control:

| | real | this model | one-population control |
|---|---|---|---|
| R1 tail Q<25 | 11.12% | 11.95% | 0.54% |
| R1 tail Q<20 | 4.92% | 5.91% | **0.00%** |
| R2 tail Q<25 | 26.22% | 26.36% | 11.74% |
| R2 tail Q<20 | 12.48% | 14.00% | 0.32% |

### What it is NOT

**Not current chemistry.** HiSeq 2500 is 2020-era, and its scores are continuous. Modern
instruments commonly emit binned scores. GIAB's HG002 path has no NovaSeq library at all —
its Illumina sets are 2x250 HiSeq, HiSeq homogeneity, a mate-pair set and an exome — so this
is the best provenanced option available from that source, not the most modern one. Sourcing
a binned/current-chemistry library is #730.

**Not a claim about your data.** It is a starting point. Fit your own model with
`gen-seq-error-model` whenever you can; that is what the tooling is for.

**Fitted at 250 bp.** Generating at another read length rescales the curve to fit, which is
what NEAT2 did and is an approximation. `gen-reads` logs a warning naming both lengths when
it does (#742). `read_len` defaults to 250 so that taking both defaults needs no rescaling.

Measured cost of rescaling, R1 at 151 bp against this model's native 250 bp: the per-cycle
shape is preserved to a tenth of a Q at the 25%, 50% and 75% marks, and only the end of the
read moves — the last cycle reads Q30.6 instead of Q23.2, and reads whose last 50 bases
average below Q20 fall from 5.70% to 1.28%. The direction is conservative: a rescaled read is
cleaner than a native one, never dirtier. For comparison, the pre-v3.4.0 default produced
0.00% of those reads at any length. A 151 bp model is #744.

## The variant models: GIAB HG002

`default_mutation_model.json.gz` is fitted by `gen-mut-model` from the GIAB HG002 v4.2.1
GRCh38 truth VCF, restricted to its `noinconsistent` high-confidence BED, against an unmasked
GRCh38 (`scripts/delta/fit_hg002_defaults.sbatch`, job 22443316, #752).
`default_indel_model.json.gz` and `default_trinuc_model.json.gz` are its own indel and
trinucleotide components, extracted unchanged, so every mutation default describes the same
sample. A test enforces that.

| | |
|---|---|
| `mutation_rate` | 0.0015161 |
| `homozygous_frequency` | 0.3884 |
| SNP / insertion / deletion | 0.8727 / 0.0617 / 0.0655 |
| CpG context weight | 4.48x the mean context |

**Checked against the VCF without eidolon.** `bedtools intersect` of the biallelic records
(by POS) with the BED gives 3,364,039 SNPs, 237,674 insertions and 252,894 deletions over
2,542,242,843 bp. The model counts 353 fewer, all at BED edges (#770). Each of the 64
per-context SNP rows, rebuilt from those SNPs and the reference, agrees with the model to
3.4e-4. The insertion and deletion length distributions agree to 3e-5. Ti/Tv of the source
SNPs is 2.102.

**What it is not:** a callset. Rates are per base of GIAB's high-confidence regions, which
exclude the hardest parts of the genome. Multi-allelic sites (47,781) are not used.

### Previous default

NEAT2's `MutModel_NA12878.p.gz`, converted: `mutation_rate` 0.0010987, `homozygous_frequency`
1/3 (NEAT2 used 0.01), SNP / insertion / deletion 0.95 / 0.03 / 0.02. That file is kept as
`eidolon-core/src/models/test_fixtures/neat2_mutation_model.json.gz`, the fixture for loading
pre-stamp and pre-#763 models. `default_mutation_model_bkup.json.gz` is NEAT2's first
conversion and is kept for provenance only.

## Building a custom model

The default models are based on public data and may not match the properties of the data
that you want to simulate.

Fragment length distribution can vary across different labs and preparation methods. To
match your own data as closely as possible, `eidolon` provides `gen-frag-length-model`, a
tool that can generate a model based on your data, with the library prep you want to
simulate:

```bash
eidolon gen-frag-length-model -c your_config.yml   # or gen-bam-models for frag + GC together
```

Indel distributions can vary between datasets, and so our default model may not fit your
use case. To build one from your own reads, use `gen-seq-error-model`:

```bash
eidolon gen-seq-error-model -c your_config.yml
```

That fits the quality-score model from your FASTQ. The indel parameters above are static
defaults that this tool does not fit; #662 tracks making the context curve fittable from
a BAM.
