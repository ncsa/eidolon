# Shipped default models

These are the models `gen-reads` falls back on when a config supplies none. They are
compiled into the binary with `include_bytes!`, so they ship with every release.

**Updated Defaults** After analysis of real sequencing data, we have updated the default error 
model for `eidolon`. These are the parameters and the value that they are assigned in the absence
of a custom fragment model. This file records the source of each value, for future reference and for 
repeatability.

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

Shape: 1087 bins over 8–1094 bp, no gaps. mean 431.8, sd 112.3, **skew +0.528**,
p05/p50/p95/p99 = 258/424/623/746.

**Cross-validated against a different chromosome.** A model built from chr20/21/22 was
checked against chr1's fragments from the same library (27.6M independent pairs) with
`scripts/delta/validate_frag_model.sh`: mean within **0.34%**, sd **0.12%**, skew **0.011**,
p99 **0.13%** — against tolerances of 2% / 5% / 0.15 / 5%. A model built from chr1 itself
did only marginally better (0.01% / 0.03% / 0.001 / 0.00%).

### Indel-error lengths

`ins_length_distribution` / `del_length_distribution` were updated based on HCC1395.

| | |
|---|---|
| **Source** | HCC1395 matched **normal**, SEQC2 Somatic Mutation WG reference sample |
| **Reference** | GRCh38, chr20 + chr21 + chr22 |
| **Region** | the ten 400 kb loci a realism-panel run placed (`realism_21795898/regions.bed`) |
| **Events** | 1,726 low-support indels — 1,058 deletions and 668 insertions |
| **Measured by** | `scripts/delta/indel_context.sbatch`, Delta job 21801707 |

Indels were split from variants by support fraction: below 10% of local depth is slippage,
at or above 25% is a variant. Only the **low-support** side feeds this model; variant indel
size is a different population and belongs to placement.

| \|len\| | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | tail |
|---|---|---|---|---|---|---|---|---|---|---|---|
| deletions (n) | 762 | 135 | 42 | 47 | 11 | 8 | 5 | 9 | 5 | 6 | 11,12,13,15,16,17,19,20,21,22,23,25,27,34,38,45 |
| insertions (n) | 426 | 120 | 22 | 55 | 13 | 8 | 1 | 6 | 1 | 5 | 11,12,13,15,18,19,22,27,30 |

Deletions: n = 1,058 over 26 bins. Insertions: n = 668 over 19 bins. The two arms did not
observe the same lengths, so each carries its own value list.

### The indel context curve (#661)

We added a further refinement to the model, based on our findings of concentrated indels
in homopolymer regions. `indel_context_curve` scales `indel_probability` by the length 
of the homopolymer run the base sits in, rather than spreading it uniformly.

| | |
|---|---|
| **Source** | HCC1395 matched **normal**, SEQC2 Somatic Mutation WG reference sample |
| **Reference** | GRCh38, chr20 + chr21 + chr22 at 46x |
| **Background** | 3,999,990 reference bases, exact (N runs excluded) |
| **Events** | 1,726 slippage errors |
| **Measured by** | `scripts/delta/indel_context.sbatch`, Delta job 21674484 |
| **Raw data** | `/projects/bhrd/jallen17/eidolon-access-results/indelctx/job_21674484/` |

| run | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | >=10 |
|---|---|---|---|---|---|---|---|---|---|---|
| propensity | 0.64 | 0.76 | 0.82 | 1.11 | 1.58 | 1.84 | 5.64 | 12.16 | 24.24 | 39.20 |

### Building a Custom Model

The default models are based on public data and may not match the properties of the data
that you want to simulate.

Fragment length distribution can vary across different labs and preparation methods. To 
match your own data as closely  as possible, `eidolon` provides `gen-frag-length-model`, 
a tool that can generate a model based on your data, with the library prep you want to 
simulate:

```bash
eidolon gen-frag-length-model -c your_config.yml   # or gen-bam-models for frag + GC together
```

Indel distributions can vary between datasets, and so our default model may not fit your
use case. To simulate your data, use `gen-seq-err

### Analysis findings

Updating the models was motivated by careful analysis of public data. We found the current 
model produced left-skewed (−0.434) data where the real data we analyzed was consistently 
right-skewed.

## The sequencing error model

The model was built inline in `models/sequencing_error_model.rs`. It ships with every run 
that does not supply its own, so it gets the same accounting as the files above.

### The NEAT2 constants (#660)

We initially tried to match NEAT2 as closely as possible. NEAT2 shipped with the following
defaults hardcoded in its sequencing error model.

| parameter | value | NEAT2 name |
|---|---|---|
| `indel_probability` | 0.01 | `SIE_RATE` — odds a sequencing error is an indel |
| `insertion_fraction` | 0.4 | `SIE_INS_FREQ` — odds such an indel is an insertion |
| length distribution | `[0.999, 0.001]` over lengths `[1, 2]` | — |
| insertion base composition | uniform over ACGT | — |

**These were initially mistranslated in the Rust port.** The insertion fraction was used 
as the indel rate, the real indel rate was dropped, and the insertion split was replaced
by a hardcoded `0.5` — about **40x too many** sequencing errors made into indels.
#660 restored the NEAT2 defaults. 

**Analysis of the inherited defaults** We found on Illumina data, an indel error rate around
1e-5/base. The original error model was a first-order approximation using a simple length 
distribution that produced indel errors of a single base 99% of the time. The data showed
a significant (16.4%) number of slippage events at 3bp or more, and a small but measurable
(at scale) number of events (0.70%) 20 base pairs or more in length. We decided to use 
this data set to build our new baseline error model, as it's error properties were more 
robust.

We found in the data that deletions and insertions are not the same shape — 72.0% of 
deletions are 1 bp against 63.8% of insertions — so we separated the distributions, similar
to how the mutation model treats insertions and deletions as independent events.

The data showed an `insertion_fraction` 0.387, comparable to our existing default of 0.4.

**Each entry is a normalized enrichment** — the share of indel errors at that run length
divided by the share of reference bases at that run length. That makes the curve
1.0-centered by construction over its human background, so applying it *redistributes*
`indel_probability` rather than raising it: the genome-wide total on human is unchanged.
On a reference with different homopolymer composition the realized total moves with that
composition, which is the intended behavior — measured at **0.745x** on the 4.6 Mb E. coli
fixture and **0.734x** for an idealized 50% GC random sequence. A genome with fewer
homopolymers shows less slippage overall.

**It is deliberately not the variant curve.** #378 measures a *separate* homopolymer
propensity for germline and somatic variants, which reaches 60.44x at runs >= 10 where
errors reach 39.20x. The two are measurably different.
