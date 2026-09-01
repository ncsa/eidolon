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
right-skewed; centred at p50 554 against this library's 424; truncated at 799; carried 33
integer lengths inside its own range with **no bin at all**; and held an isolated spike at
fragment length **1** with a 30-wide hole above it — one stray read that survived a filter.
Its provenance is unknown; it predates the Rust port.

Nothing tested it. `test_discrete_default` pinned all 766 of its values as a literal array,
which asserted the bytes had not changed and said nothing about whether they were usable.
`the_shipped_default_is_a_usable_fragment_distribution` in `fragment_length.rs` now checks
the properties that actually matter, and it rejects the old file on its first gap.

## The others

**Provenance unknown.** All predate the Rust port and none is validated by anything beyond
round-trip serialization:

- `default_mutation_model.json.gz` (+ `_bkup`)
- `default_indel_model.json.gz`
- `default_quality_score_model.json.gz`
- `default_trinuc_model.json.gz`

Each deserves the same treatment: a known source, a measurement against real data, and a
test asserting it is usable rather than unchanged. Until then, treat any result that leans
on them as resting on an unaudited assumption.
