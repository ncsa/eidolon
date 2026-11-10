# Upgrading from 4.x

v5.0.0 removes duplicate config-key names and the pre-v3.0.0 read-name prefix, and changes
the order of `gen-reads` FASTQ records by default. A removed name stops the run before any
output is written, and the error names the replacement.

| What you see | Cause | What to do |
|---|---|---|
| `gen-reads`: `` `output_filename` (removed in v5.0.0; use `output_basename`) `` | Each key now has one name, shared across subcommands. | Rename `output_filename:` to `output_basename:`. Same value, same output names. |
| `gen-cancer-reads`: the same message for `output_filename`, `coverage` or `rng_seed_root` | These were second names for `output_basename`, `total_coverage` and `rng_seed`. | Use the name the error gives. |
| `filter-reads`: "read name uses the pre-v3.0.0 prefix `RNEAT_generated_`" | Output from v2.0.0 or earlier. 4.x still read it, with a warning. | Rename the reads first (below). |
| `gen-reads` or `gen-cancer-reads`: `` `output_prefix` (renamed `output_basename` in v5.0.0) `` | `output_prefix` held a file name with no directory, which "prefix" suggests otherwise. | Rename `output_prefix:` to `output_basename:`. A value containing `/` is now refused; put the directory in `output_dir`. |
| `gen-reads` FASTQs in random order, and about 60% larger | Records are shuffled by default (#870), as a sequencer's are. | Nothing, usually. Set `shuffle_fastq: false` for genome order. |

To update every config under a directory:

```bash
grep -rlE '^output_(filename|prefix):' configs/ \
  | xargs sed -i -E 's/^output_(filename|prefix):/output_basename:/'
```

For `gen-cancer-reads` configs, also rename `coverage:` to `total_coverage:` and
`rng_seed_root:` to `rng_seed:`. Do not run those two over `gen-reads` configs, where
`coverage` is correct.

## Old read names in `filter-reads`

Rename the prefix, then filter the result:

```bash
zcat old_r1.fastq.gz | sed 's/^@RNEAT_generated_/@EIDOLON_generated_/' | gzip > new_r1.fastq.gz
```

[Upgrading from 2.0.0](from-2-0-0.md) covers the other v3.0.0 renames.

## Shuffled FASTQs

Genome-ordered FASTQs made bwa-mem2 need 166 GB instead of 90 GB on a 21x T2T-CHM13
simulation, with identical alignments. [FASTQ shuffling](../gen-reads/fastq-shuffling.md)
has the measurements, the run-time and disk cost, and the gzip size increase. The truth VCF
and golden BAM are the same with the shuffle on or off.
