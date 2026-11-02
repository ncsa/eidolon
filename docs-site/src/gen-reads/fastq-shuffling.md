# FASTQ Shuffling
`eidolon` writes its FASTQs in **random order** by default (`shuffle_fastq: true`), as a
sequencer does. Mates stay in register: line *n* of R1 and line *n* of R2 are the two reads of
one pair.

Reads are generated chunk by chunk along the genome, so without shuffling the FASTQ comes out in
genome order: contig by contig, by position within each contig. Variant calling does not depend
on this, because the BAM is coordinate-sorted before calling. Tools that read the FASTQ in
batches do see the difference, aligners in particular: a batch can be filled entirely with reads
from one region, such as a long satellite array. Measured on a 21x simulation of T2T-CHM13
(#870), bwa-mem2 peaked at 166 GB on reads in genome order and at 90 GB on the same reads
shuffled, with the same alignments either way.

```yaml
shuffle_fastq: true    # default; false writes reads in genome order
```

**What it does and does not change.**
- Only the order of FASTQ records changes. The reads themselves, the truth VCF and the golden BAM
  are identical with the shuffle on or off. The shuffle draws from its own branch of the random
  number generator, so it does not disturb the streams that generate reads.
- Output is deterministic: the same seed gives byte-identical FASTQs at any `num_threads`.
- Turning the shuffle off restores genome order exactly as earlier versions wrote it.

**Cost.** The shuffle is a two-pass external sort. Records are scattered at random into bucket
files in the run's temp directory, then each bucket is shuffled in memory and written out. Buckets
are shuffled in parallel only as far as they fit in about 4 GB, so memory stays bounded at any
genome size. Temp disk stays close to one copy of the output, because each input is deleted as it
is read.

Measured on T2T-CHM13 chr21 at 20x (about 3M pairs, 8 threads), against `shuffle_fastq: false`:
- wall time 112 s against 99 s;
- peak memory 2.3 GB against 0.5 GB;
- **gzipped FASTQs about 60% larger.** Reads in genome order overlap their neighbors, which gzip
  exploits; shuffled reads, like a real sequencer's, compress less well.

**Shuffling output from older versions, or from `shuffle_fastq: false`.** With `seqkit`, shuffle
each mate separately **with the same seed**. seqkit has no paired mode, but its permutation
depends only on the seed and the number of records, and R1 and R2 have the same count:

```bash
seqkit shuffle -s 11 sample_r1.fastq.gz -o sample_r1_shuffled.fastq.gz
seqkit shuffle -s 11 sample_r2.fastq.gz -o sample_r2_shuffled.fastq.gz
```

`seqkit shuffle` reads every record into memory; its low-memory `-2/--two-pass` mode supports
FASTA only. Do not pass both mates to one call: they would be shuffled together into a single
output. (`seqkit`: https://github.com/shenwei356/seqkit.)
