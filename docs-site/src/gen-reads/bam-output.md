# BAM Output
`eidolon` can write a golden BAM file alongside the FASTQ output. The BAM contains the same reads as the FASTQ — same sequences, same quality scores, same variants and sequencing errors applied — with alignment information included. To enable it, set `produce_bam: true` in your `gen-reads` config. The output path is derived automatically from `output_filename` (e.g. `output_filename: my_run` → `my_run.bam`). Please note that the BAM has a higher overhead than the fastq, and may take longer to produce.

```yaml
produce_bam: true
```

The CIGAR strings in the BAM reflect the full ground-truth alignment:
- Genomic variants (SNPs, insertions, deletions) are encoded as `M`, `I`, and `D` ops.
- Sequencing error indels are also encoded: deletion errors add `D` ops for skipped reference bases; insertion errors add `I` ops for inserted bases.

The BAM is written in coordinate-sorted order. No post-processing sort is required before indexing:

```bash
samtools index my_run.bam
```

Note: heterozygous variants are applied probabilistically, identical to the FASTQ. A read drawn to the reference allele will not show the variant in either the FASTQ or the BAM.
