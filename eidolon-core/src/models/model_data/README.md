# Shipped default models

These are the models `gen-reads` falls back on when a config supplies none. They are
compiled into the binary with `include_bytes!`, so they ship with every release.

After analysis of real sequencing data, we have updated the default error model for
`eidolon`. This file records the source of each value, for future reference and for
repeatability.

## Status

| model / parameter | source | measured |
|---|---|---|
| `default_fragment_length_model.json.gz` | HCC1395 normal | yes |
| `error_rate` | GIAB HG002 2x250, fitted | yes — 0.003774, measured |
| `indel_probability` | NEAT2 static default | no |
| `insertion_fraction` | NEAT2 static default | yes — confirmed at 0.387 |
| indel-error lengths | HCC1395 normal | yes |
| homopolymer context curve | HCC1395 normal | yes |
| `default_sequencing_error_model.json.gz` | GIAB HG002 2x250, fitted | yes — see below |
| `default_mutation_model.json.gz` (+ `_bkup`) | NEAT2 `MutModel_NA12878.p.gz` | no — source verified, values unmeasured |
| `default_indel_model.json.gz` | NEAT2 `MutModel_NA12878.p.gz` (`INDEL_FREQ`) | no — source verified, values unmeasured |
| `default_trinuc_model.json.gz` | NEAT2 `MutModel_NA12878.p.gz` (`TRINUC_MUT_PROB`) | no — source verified, values unmeasured |

## `default_fragment_length_model.json.gz`

| | |
|---|---|
| **Source** | HCC1395 matched **normal**, SEQC2 Somatic Mutation WG reference sample |
| **Read group** | `WGS_NS_N_1` (NovaSeq replicate 1, `WGS_NS_N_1.bwa.dedup.bam`) |
| **Origin** | `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS` |
| **Reference** | GRCh38, chr-prefixed |
| **Region** | chr20 + chr21 + chr22 |
| **Pairs used** | 32,627,236 of 32,669,084 collected (0.13% trimmed as outliers) |
| **Built with** | `eidolon gen-frag-length-model`, `min_reads: 100`, default `distribution: discrete` |
| **Built at** | eidolon `3.2.1+2f98bb6`, 2026-08-30 |

Shape: 1087 bins over 8–1094 bp, no gaps. Mean 431.8, sd 112.3, **skew +0.528**,
p05/p50/p95/p99 = 258/424/623/746.

**Cross-validated against a different chromosome.** A model built from chr20/21/22 was
checked against chr1's fragments from the same library (27.6M independent pairs) with
`scripts/delta/validate_frag_model.sh`: mean within **0.34%**, sd **0.12%**, skew **0.011**,
p99 **0.13%** — against tolerances of 2% / 5% / 0.15 / 5%. A model built from chr1 itself
did only marginally better (0.01% / 0.03% / 0.001 / 0.00%).

Updating this model was motivated by careful analysis of public data. The previous default
produced left-skewed (−0.434) fragments where the real data we analyzed was consistently
right-skewed.

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
| substitution transitions | 0.4918 / 0.3377 / 0.1705 … | `SSE_PROB` |

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
what NEAT2 did and is an approximation — NEAT2 warned about it and eidolon does not yet
(#742). `read_len` defaults to 250 so that taking both defaults needs no rescaling.

Measured cost of rescaling, R1 at 151 bp against this model's native 250 bp: the per-cycle
shape is preserved to a tenth of a Q at the 25%, 50% and 75% marks, and only the end of the
read moves — the last cycle reads Q30.6 instead of Q23.2, and reads whose last 50 bases
average below Q20 fall from 5.70% to 1.28%. The direction is conservative: a rescaled read is
cleaner than a native one, never dirtier. For comparison, the pre-v3.4.0 default produced
0.00% of those reads at any length. A 151 bp model is #744.

## The variant models: NEAT2's NA12878

`default_mutation_model.json.gz` (+ `_bkup`), `default_indel_model.json.gz` and
`default_trinuc_model.json.gz` all derive from **one** upstream file:
`~/code/neat2/models/MutModel_NA12878.p.gz`.

That was recorded as "unrecorded provenance" until it was checked. It is not a cancer model
and not a toy: NA12878 is the CEPH/GIAB germline reference sample, which makes it a defensible
source for germline defaults. What is genuinely unrecorded is how NEAT2 built it — reference
build, aligner and caller are not stated upstream, and given NEAT2's age it is likely GRCh37
era.

Verified by value, not by reading:

| eidolon | upstream field | agreement |
|---|---|---|
| `mutation_rate` 0.0010987132390211135 | `AVG_MUT_RATE` | exact, 16 digits |
| `_bkup`'s `variant_dist` first weight 0.886404192662459 | `SNP_FREQ` | exact, 16 digits |
| `default_indel_model` `ins_dist` (70 lengths) | positive `INDEL_FREQ` keys, renormalized | max diff 5.4e-17 |
| `default_indel_model` `del_dist` (72 lengths) | negative `INDEL_FREQ` keys, renormalized | max diff 1.1e-16 |
| `default_trinuc_model` `snp_distro` (64) | `TRINUC_MUT_PROB`, normalized, alphabetical (AAA, AAC, AAG, …) indexed 0–63 | max diff 5.7e-17 |

### Two deliberate departures, and one unexplained number

**`homozygous_frequency` 0.01 → 0.3333.** `_bkup` carries NEAT2's 0.01, which implies a
het/hom ratio near 99. The live model uses 1/3, a ratio near 2.0, which is what human data
shows. This is the difference the NEAT comparison table in the README refers to, and it is
locked by a test.

**`variant_dist` is keyed by name** (`SNP`/`Insertion`/`Deletion`) where `_bkup` keys by
integer. Presentation only; the weights are unchanged.

**`insertion_probability` 0.4538979885714955 does not derive from the pickle by any obvious
route.** The upstream insertion share of total indel mass is 0.4768593189964158 and the ratio
of distinct insertion to deletion lengths is 0.4929577464788732. Neither is the shipped value,
while the two length distributions beside it match to floating-point epsilon. Recorded as
unexplained rather than reverse-engineered into a plausible story.

### What is still missing

A measurement. Knowing the source is not knowing whether the values describe the data eidolon
is asked to simulate, and none of the three has been checked against a modern human callset.
That is the germline-defaults review (#752).

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
