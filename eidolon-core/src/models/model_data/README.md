# Shipped default models

These are the models `gen-reads` falls back on when a config supplies none. They are
compiled into the binary with `include_bytes!`, so they ship with every release.

**A default is a claim about what real data looks like.** Anyone who runs eidolon without
building their own models gets these, and every number they measure is downstream of them.
This file records where each came from, because for most of them nobody currently knows —
and that is how a fragment model that no real library resembles stayed in place for years.

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
did only marginally better (0.01% / 0.03% / 0.001 / 0.00%). Fragment length is a property
of the library preparation, not of the region, so measuring it on three chromosomes is not
a limitation of this file.

### What this default is NOT

It is one library, one instrument, one prep. **Fragment length is set by chemistry** —
size selection, bead ratios, kit version — so a newer library may well sit elsewhere, and
recent chemistries in particular may carry more of the long tail than this one does. If
your data differs, build your own; that is a two-minute job and the whole reason
`gen-frag-length-model` exists:

```bash
eidolon gen-frag-length-model -c your_config.yml   # or gen-bam-models for frag + GC together
```

The choice here was deliberate: **data we can account for beats data that looks wrong.**

### What it replaced, and why

The previous default was left-skewed (**−0.434**) where every real size-selected library is
right-skewed; centered at p50 554 against this library's 424; truncated at 799; carried 33
integer lengths inside its own range with **no bin at all**; and held an isolated spike at
fragment length **1** with a 30-wide hole above it — one stray read that survived a filter.
Its provenance is unknown; it predates the Rust port.

Nothing tested it. `test_discrete_default` pinned all 766 of its values as a literal array,
which asserted the bytes had not changed and said nothing about whether they were usable.
`the_shipped_default_is_a_usable_fragment_distribution` in `fragment_length.rs` now checks
the properties that actually matter, and it rejects the old file on its first gap.

## The sequencing error model

Not a file — it is built inline in `models/sequencing_error_model.rs`, small enough that it
was never worth a `.json.gz`. It ships with every run that does not supply its own, so it
gets the same accounting as the files above.

### Inherited parameters (#660)

Four parameters carry over from NEAT2's `genSeqErrorModel.py`, where they are static
defaults rather than fitted values — that tool fits the quality-score model from FASTQ and
leaves the indel parameters fixed.

| parameter | value | source name |
|---|---|---|
| `indel_probability` | 0.01 | `SIE_RATE` — odds a sequencing error is an indel |
| `insertion_fraction` | 0.4 | `SIE_INS_FREQ` — odds such an indel is an insertion |
| insertion base composition | uniform over ACGT | `SIE_INS_NUCL` |
| substitution transition matrix | 0.4918 / 0.3377 / 0.1705 … | `SSE_PROB` |

`indel_probability` and `insertion_fraction` were transposed during the Rust port: 0.4 was
applied as the indel rate and the insertion split was fixed at 0.5, giving roughly 40x the
intended indel-error rate. #660 corrected both.

**Status.** 0.01 matches its source but has not been measured against real data; at Q35 it
gives ~3.2e-6 indel errors per base against a real Illumina rate near 1e-5. Changing it
needs its own measurement. `insertion_fraction` **has** been measured and holds — see the
indel-error length section below.

### `error_rate` (0.006638164688495656)

Fitted, not a static default. It is the `avgError` of NEAT2's bundled `errorModel_toy.p`,
computed from the sequencing data that model was built on. The originating sample is not
recorded upstream. Listed separately here because it has different standing from the four
above.

### Indel-error lengths

`ins_length_distribution` / `del_length_distribution` are **measured**, not NEAT2's.

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

The source ships the **raw counts** as weights — `DiscreteDistribution` normalizes them — so
the constants can be checked against the job output directly.

**Not truncated.** Every observed length is carried, to 45 bp for deletions and 30 bp for
insertions. The tail bins hold single observations, so the *shape within* the tail is
uncertain — but its total mass (12 of 1,726 events, **0.70%** at >= 20 bp) is an ordinary
estimate, and dropping it would remove the only part of this distribution able to produce a
candidate breakpoint, deciding #672 by construction rather than measuring it.

The model can therefore emit a 45 bp deletion as a sequencing error, at p = 0.00095. If that
turns out to be a mapping artifact rather than slippage, the fix belongs in the classifier
upstream, not in a cutoff chosen here.

Deletions and insertions are **not the same shape** — 74.0% of deletions are 1 bp against
64.8% of insertions — and are now separate distributions. They were previously one
distribution cloned twice.

**What it replaced.** NEAT2's `[0.999, 0.001]` over lengths `[1, 2]`: 99.9% of indel errors
at a single base, and nothing above 2 bp emittable at all, against a measured **16.4%** of
slippage events at 3 bp or more. Same provenance as the constants above — the
`if PILEUP == None` branch, never fitted.

**Independently, this run confirms `insertion_fraction`.** 668 of 1,726 low-support indels
are insertions, a fraction of **0.387** against NEAT2's `SIE_INS_FREQ = 0.4`. That is the
first of these NEAT2 placeholders to be checked against real data and hold.

### The indel context curve (#661)

`indel_context_curve` scales `indel_probability` by the length of the homopolymer run the
base sits in — the mechanism that concentrates a realistic total where real slippage
actually happens, rather than spreading it uniformly.

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

Monotone, crossing 1.0 at run 4. Reproduced on a second sample with a different aligner:
candidate clip boundaries sit in homopolymers at 2.26x on HCC1395/bwa and 2.56x on
NA12878/novoalign, with controls at chance in both.

**Each entry is a normalized enrichment** — the share of indel errors at that run length
divided by the share of reference bases at that run length. That makes the curve
1.0-centered **by construction** over its human background, so applying it *redistributes*
`indel_probability` rather than raising it: the genome-wide total on human is unchanged.
On a reference with different homopolymer composition the realized total moves with that
composition, which is the intended behavior — measured at **0.745x** on the 4.6 Mb E. coli
fixture and **0.734x** for an idealized 50% GC random sequence. A genome with fewer
homopolymers really does slip less.

**This is a default, not a measurement of your data.** Same status as the fragment-length
model above: one sample, one instrument, one prep. Slippage depends on chemistry and on
the aligner's gap placement, so if yours differs this curve will not describe it. #662
makes it fittable from a BAM.

**This is the sequencing-error curve.** Variants have their own, measurably steeper,
propensity — 60.44x at runs >= 10 against 39.20x here. That one belongs to variant
placement; see #378.

## The others

### `default_quality_score_model.json.gz`

Converted from NEAT2's bundled `errorModel_toy.p` — its `initQ1` seed vector and `probQ1`
transition tensor, verified to agree to floating-point epsilon (max abs diff 5.6e-16 across
sampled cells). Same standing as `error_rate` above: fitted from the sequencing data that
model was built on, originating sample not recorded upstream.

Shape: 101-base reads, 42 continuous scores (0–41), a 100 x 42 x 42 position-by-previous-score
transition tensor. **That describes an older chemistry** — current instruments commonly emit
binned quality scores at 151 bp. The model supports binned scores and other read lengths;
this default does not exercise either. See #677.

### The others

**Provenance unknown.** All predate the Rust port and none is validated by anything beyond
round-trip serialization:

- `default_mutation_model.json.gz` (+ `_bkup`)
- `default_indel_model.json.gz`
- `default_trinuc_model.json.gz`

Each deserves the same treatment: a known source, a measurement against real data, and a
test asserting it is usable rather than unchanged. Until then, treat any result that leans
on them as resting on an unaudited assumption.
