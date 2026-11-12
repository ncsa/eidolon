# Targeted Generation with a BED File
For exome, targeted-panel, or any run where you only want reads over specific regions, you can supply a BED file at generation time:

```yaml
target_bed: /path/to/targets.bed
```

When `target_bed` is set, `gen-reads` skips contigs absent from the BED entirely — no blocks are read, no variants are placed, and no reads are generated for those contigs. Within covered contigs, reads and variants are generated only over the intersecting non-N regions. This is the recommended approach for targeted runs on large genomes; it is far more efficient than generating genome-wide reads and post-filtering with `filter-reads`.

The BED contig names must match the short names derived from the reference FASTA (the text after `>` up to the first whitespace character). Both `.bed` and `.bed.gz` inputs are accepted.

### `coverage` on a targeted run is the mean on-target depth

From v5.0.0, `coverage: 30` with a `target_bed` gives a mean depth of about 30× over the target bases. Earlier versions treated it as a fragment budget and delivered less (see below).

A fragment is *owned* by the target its start falls in, but its far end runs into the flanking sequence beyond the target, as it must: a 250 bp fragment cannot fit inside a 178 bp exon. Those spilled read bases are real off-target reads, but they are not on-target depth. So `gen-reads` estimates, once per run, what fraction of its read bases will land on target, from the BED geometry and the fragment-length model, and divides the fragment budget by it. It logs the result:

```
target_bed: 51.4% of read bases are expected on target, so the fragment budget is scaled 1.946x ...
```

A base counts as on-target if *any* target covers it, so targets that touch or sit close together lose less to spill and are scaled less. A window much wider than a fragment (the shards in a sharded whole-genome run) is scaled by less than 0.01% (pinned by a unit test). The estimate is a sample of a million fragments, with a standard error I estimate at under 0.1% of the depth; that was not measured separately. The scale changes how many fragments are placed, not where, so the truth VCF and the BAM record layout are as before.

**Measured** with `samtools depth -b` over the output BAM, on `coverage: 30`:

| BED | mean on-target depth |
|---|---|
| 200 isolated 180 bp targets, synthetic 1 Mb reference (before: 11.4×) | 30.2× |
| 250 isolated 80 bp targets, 100 bp reads | 29.9× |
| 400 touching 1 kb intervals | 30.2× (no inflation) |
| 4,294 chr22 targets of the IDT exome design from the [DeepVariant exome case study](https://storage.googleapis.com/deepvariant/exome-case-study-testdata/idt_capture_novogene.grch38.bed) (GRCh38; 203,161 targets, mean width 191 bp), real chr22 sequence, 150 bp reads, 300 bp fragments | 30.2× |

**The mean is corrected, not each target.** A narrow target captures fewer on-target bases per fragment than a wide one, so depth still varies with target width: in the chr22 run, targets of 150 bp or less (63% of the targets, 32% of the target bases) averaged 13.8×, and wider ones 37.8×.

**Targets narrower than a read.** Before v5.0.0 a target no longer than `read_len` owned no fragments and was covered only by spill from its neighbors (1.1× against a requested 30× for those targets on the chr22 run). They now own fragments like any other target.

- **Do not read a 100%-on-target simulation as "better".** eidolon before v3.1.x confined every fragment inside its target, which yields 100% on-target and is physically impossible for targets narrower than the fragment length. It also under-delivered the fragment budget itself and collapsed the realized fragment-length distribution near target edges.
- To reproduce a pre-v5.0.0 run's fragment count, multiply its `coverage` by the on-target fraction from the log line.

**Not verified:** one real exome design on one chromosome; N-gap-adjacent targets (the estimate treats only contig ends as termini); other fragment-length models than a normal.
