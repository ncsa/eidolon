# Generating a Sequencing Error Model
`eidolon` can also learn a sequencing error model from real FASTQ data using the `eidolon gen-seq-error-model` subcommand. This reads per-base quality scores from the FASTQ to build a Markov quality score model and computes an average base error rate. Optionally, you can supply a BAM file or a custom TSV to set the SNP transition matrix (which base errors are most likely to become which other bases).

```bash
$ eidolon gen-seq-error-model -c gen_seq_error_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `seq_error_model` config key.

Copy `template_config/gen_seq_error_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to input FASTQ file (.fastq or .fastq.gz)
fastq_file: /path/to/reads.fastq.gz

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/seq_error_model.json.gz

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# optional: maximum number of reads to use for model learning; 0 = unlimited (default: 0)
max_reads: 0

# optional: quality score ASCII offset; 33 for Illumina 1.8+/Sanger, 64 for older Illumina (default: 33)
qual_offset: 33

# optional: largest read length to model in bp; 0 = no limit (default: 1000)
# One transition matrix is stored per read position (~69 KiB each), so this bounds memory:
# 1000 bp is about 67 MiB. A longer read is an error, not a silent truncation. Raise it for
# genuinely longer reads; note that long-read error models are not supported yet (#319).
max_model_read_length: 1000

# optional: list of Q-score bins to quantize the learned model to (e.g. NovaSeq 6000 uses
# [2, 12, 23, 37]). When set, each observed Q-score is snapped to its nearest bin before
# the model is built, and gen-reads will emit only these values. Omit for a continuous model.
# Bins must be < 94 and cannot include 31 (encodes to '@' under Phred+33).
# binned_quality_bins: [2, 12, 23, 37]

# optional: aligned BAM file to infer the SNP transition matrix from read-vs-reference mismatches.
# The BAM must have MD tags. If your aligner did not add them, run:
#   samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
# Ignored if transition_matrix_file is also set.
# If the BAM yields no mismatches at all, eidolon errors out rather than quietly using the
# default matrix. To use the default deliberately, omit this key.
bam_file: /path/to/aligned.bam

# optional: custom 4x4 SNP transition matrix TSV (rows/columns: A C G T).
# A single header line is allowed. Diagonal values are ignored. Each row needs four
# finite, non-negative values with some weight off the diagonal; a malformed row stops
# the run and names its line.
# Takes precedence over bam_file.
transition_matrix_file: /path/to/matrix.tsv
```

## SNP transition matrix priority
1. `transition_matrix_file` (explicit TSV) — highest priority
2. `bam_file` (inferred from MD-tagged BAM mismatches)
3. Built-in default matrix (inherited from Python NEAT) — used when neither is provided

Supplying `bam_file` is a request to *fit* the matrix from data, so if the BAM yields no
read-vs-reference mismatches, eidolon **errors out** instead of falling back to the default. A
model that looks trained and is actually the default is indistinguishable from a trained one
downstream. **Omitting `bam_file` is how you ask for the default matrix.**

## BAM MD tag requirement
The BAM path requires MD tags to identify reference bases at mismatch positions. MD is an
*optional* SAM tag, so check your aligner: BWA-MEM writes it, `minimap2` needs `--MD`, and STAR
needs `--outSAMattributes MD`. If your BAM lacks them, generate them with:
```bash
samtools calmd -b aligned.bam reference.fa > aligned_with_md.bam
samtools index aligned_with_md.bam
```
Without MD tags the inference finds nothing, which is an error rather than a silent fallback.

## TSV format
The transition matrix TSV has 4 data rows (one per reference base: A, C, G, T) and 4 whitespace-separated columns (one per read base: A, C, G, T). The first row is skipped if it is non-numeric (treated as a header). Diagonal values are zeroed automatically. Rows are re-normalized to sum to 1.

```
A    C    G    T
0.0  0.5  0.3  0.2
0.5  0.0  0.3  0.2
0.4  0.3  0.0  0.3
0.3  0.3  0.4  0.0
```

## Binned quality scores
Modern Illumina platforms (NovaSeq 6000, NextSeq 1000/2000, the simplified MiSeq reporting modes) no longer emit the full continuous Q0–Q40 range. Instead, they quantize each per-base quality into a small set of discrete bins — for example, NovaSeq 6000 emits only `{2, 12, 23, 37}`. Set `binned_quality_bins` in the config to mimic this behaviour:

```yaml
binned_quality_bins: [2, 12, 23, 37]
```

When this field is set, `gen-seq-error-model` snaps each observed Q-score from the input FASTQ to its nearest bin (ties round down) before learning seed, transition, and global-frequency counts. The resulting model file is flagged `binned_scores: true` and lists the bin values as its `quality_score_options`. `gen-reads` then samples only those bin values when emitting reads — no extra flag needed on the read-generation side.

Validation rules:
- Bins must be non-empty integers in `[0, 94)` (entries are sorted and deduped automatically).
- Value `31` is rejected — under Phred+33 it encodes to `@`, which would corrupt FASTQ output. If you need a bin near Q31, pick 30 or 32.
- Bins with no observed counts in the training FASTQ are still kept in `quality_score_options`; their transition rows fall back to a uniform distribution, and a warning is logged listing the empty bins.

## Q31 (`@`) and the first quality character
Phred+33 encodes Q31 as `@` — the same byte that begins a FASTQ record header. eidolon never writes it as the **first** character of a quality line, because a reader that scans for `@` to find record boundaries would mis-frame the file. (A reader that consumes fixed four-line records, including eidolon's own, is unaffected.)

What that means depends on the model:

- **Binned models** reject `31` outright, at config time — see the validation rules above.
- **Continuous models** may legitimately learn Q31, since real data contains it, and will emit it anywhere in a read *except* position 1. There, eidolon substitutes the nearest other quality score the model actually contains, with ties going to the lower score. A model learned over the usual Q0–Q40 range therefore substitutes Q30. The rest of the read is untouched.

This is a deliberate, bounded departure from the fitted model: the distribution at position 1 is shifted off Q31 by design, and nowhere else is affected. A model whose *only* learned score is Q31 has nothing to substitute and is reported as an error rather than emitting a quality line that cannot be read.

Some common platform bin sets (consult your sequencer's documentation for the authoritative list):

| Platform                       | Suggested bins     |
|--------------------------------|--------------------|
| NovaSeq 6000 (4-bin)           | `[2, 12, 23, 37]`  |
| NextSeq 1000/2000 (3-bin)      | `[2, 15, 35]`      |
| NextSeq 1000/2000 (4-bin)      | `[2, 12, 23, 37]`  |
