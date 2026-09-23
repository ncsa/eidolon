# FASTQ Shuffling
`eidolon` writes reads in contig order. To shuffle the output, use `seqkit shuffle` as a post-processing step:

```bash
# single-ended
seqkit shuffle sample.fastq.gz -o sample_shuffled.fastq.gz

# paired-ended (keeps mates in sync)
seqkit shuffle -2 sample_R1.fastq.gz sample_R2.fastq.gz \
    -o sample_R1_shuffled.fastq.gz -o sample_R2_shuffled.fastq.gz
```

`seqkit` uses reservoir sampling and streams from disk, keeping memory use bounded regardless of file size (`seqkit` is an open source toolkit for FASTQ files: https://github.com/shenwei356/seqkit).
