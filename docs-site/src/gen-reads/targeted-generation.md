# Targeted Generation with a BED File
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
