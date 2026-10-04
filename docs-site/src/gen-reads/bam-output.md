# BAM Output
`eidolon` can write a golden BAM file alongside the FASTQ output. The BAM contains the same reads as the FASTQ — same sequences, same quality scores, same variants and sequencing errors applied — with alignment information included. To enable it, set `produce_bam: true` in your `gen-reads` config. The output path is derived automatically from `output_filename` (e.g. `output_filename: my_run` → `my_run.bam`). Please note that the BAM has a higher overhead than the fastq, and may take longer to produce.

```yaml
produce_bam: true
```

The CIGAR strings in the BAM reflect the full ground-truth alignment:
- Genomic variants (SNPs, insertions, deletions) are encoded as `M`, `I`, and `D` ops.
- Sequencing error indels are also encoded: deletion errors add `D` ops for skipped reference bases; insertion errors add `I` ops for inserted bases.

Every mapped record carries an `NM` tag, its edit distance to the reference: mismatched aligned bases plus inserted plus deleted bases. It is counted the way `samtools calmd` counts it: soft-clipped adapter bases are not edits, and an `N` on either side counts as a mismatch. `samtools stats` and MultiQC take their mismatch count and error rate from `NM`, so they report the simulated error rate. Unmapped records carry no `NM`, and no `MD` tag is written; `samtools calmd` can add one.

Reads spanning an SV junction (BND, INV, and the junctions of DEL and DUP) are written as one alignment at the position of their first piece, although part of each read comes from elsewhere. Their `NM` is correct for that alignment, so it is high, and when SVs are simulated these reads raise the error rate `samtools stats` reports above the sequencing error rate alone.

The BAM is written in coordinate-sorted order. No post-processing sort is required before indexing:

```bash
samtools index my_run.bam
```

Note: heterozygous variants are applied probabilistically, identical to the FASTQ. A read drawn to the reference allele will not show the variant in either the FASTQ or the BAM.
