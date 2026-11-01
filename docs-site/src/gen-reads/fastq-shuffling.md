# FASTQ Shuffling
`eidolon` writes reads in genome order: contig by contig, by position within each contig. A
real sequencer's FASTQ is in effectively random order. Variant calling does not depend on this,
since the BAM is coordinate-sorted before calling. Tools that read the FASTQ in batches,
aligners in particular, do see the difference: a batch can be filled entirely with reads from
one region, such as a long satellite array. Shuffle before aligning if that matters to you.
#870 tracks shuffling in `eidolon` itself.

With `seqkit`, shuffle each mate separately **with the same seed**. seqkit has no paired mode,
but its permutation depends only on the seed and the number of records, and R1 and R2 have the
same count. So the two files come out in the same order and mates stay in register:

```bash
# single-ended
seqkit shuffle -s 11 sample.fastq.gz -o sample_shuffled.fastq.gz

# paired-ended: same -s on both mates
seqkit shuffle -s 11 sample_r1.fastq.gz -o sample_r1_shuffled.fastq.gz
seqkit shuffle -s 11 sample_r2.fastq.gz -o sample_r2_shuffled.fastq.gz
```

**Memory:** `seqkit shuffle` reads every record into memory, so expect memory on the order of
the uncompressed FASTQ, per mate. Its low-memory `-2/--two-pass` mode supports FASTA only.
Do not pass both mates to one `seqkit shuffle` call: they would be read as one input and shuffled
together into a single output.

(`seqkit` is an open-source FASTQ toolkit: https://github.com/shenwei356/seqkit.)
