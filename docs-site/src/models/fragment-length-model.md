# Generating a Fragment Length Model
`eidolon` can learn a fragment length model from real paired-end alignment data using the `eidolon gen-frag-length-model` subcommand. It reads a BAM or SAM file, collects the template lengths (TLEN) of confidently-mapped concordant read pairs, trims extreme outliers, and writes a `FragmentLengthModel` that `gen-reads` can use.

**Since v3.3.0 the model keeps the observed shape by default.** Real size-selected libraries are right-skewed — a steep left edge from size selection, and a long tail of larger fragments that escaped it. A normal distribution is symmetric by construction and cannot hold that tail, and the tail is exactly what a paired-end SV caller thresholds against when it decides an insert is "larger than expected". The builder therefore fits a *discrete* (empirical) distribution over the observed lengths. Sparse histograms are smoothed so the support has no holes in it, and a model too sparse to have a shape is refused rather than written. Set `distribution: normal` for the pre-v3.3.0 two-parameter fit, which is still the right choice for small or targeted BAMs (exome, amplicon) where there is not enough data to estimate a shape.

**The shipped default** — used when a config supplies neither `fragment_model` nor `fragment_mean` — is built from HCC1395 matched normal (SEQC2 `WGS_NS_N_1`, NovaSeq, GRCh38 chr20/21/22, 32.6M read pairs): 1087 bins over 8–1094 bp, mean 431.8, sd 112.3, skew +0.528. It is cross-validated against chr1 of the same library, agreeing within 0.34% on the mean and 0.011 on the skew, which is why three chromosomes suffice. Fragment length is set by library chemistry, so treat it as one real library rather than a universal shape: if yours differs, build your own with the command above. Full provenance is in `eidolon-core/src/models/model_data/README.md`.

```bash
$ eidolon gen-frag-length-model -c gen_frag_length_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `fragment_model` config key.

Copy `template_config/gen_frag_length_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to input BAM or SAM file
# The file must contain paired-end aligned reads. BAM is strongly recommended.
input_file: /path/to/aligned.bam

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/fragment_length_model.json.gz

# optional: `discrete` (default) keeps the measured shape, including the right tail.
# `normal` collapses it to mean + standard deviation, which is the pre-v3.3.0 behaviour
# and is the right choice for sparse input (exome, amplicon, a small targeted BAM).
distribution: discrete

# optional: floor on the TOTAL number of observations, below which the model is refused.
# Default is 2 to handle smaller datasets. Set to 0 to disable filtering entirely.
# NOTE: before v3.3.0 this also DELETED any individual length seen fewer than min_reads
# times, which punched holes in the distribution — hardest in the tail, where counts are
# lowest. Sparse lengths are now handled by smoothing instead.
min_reads: 2

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false
```

## Read filtering rules
Only reads that satisfy all of the following are used:
- Paired-end and first in pair
- Mapped (not unmapped, secondary, or supplementary)
- Mate mapped to the same reference sequence
- Mapping quality > 10

## Outlier filtering
Fragment lengths that exceed `median + 10 × MAD` (median absolute deviation) are removed before fitting. This mirrors the filtering in Python NEAT's fragment length modeler. If no lengths survive the filter, lower `min_reads` or set it to `0`.
