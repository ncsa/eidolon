# Generating a GC Bias Model
`eidolon` can learn a GC bias model from a reference FASTA and an aligned BAM file using the `eidolon gen-gc-bias-model` subcommand. It walks the BAM once to accumulate per-base reference coverage, tiles the reference in fixed-size windows, computes the GC% of each window, looks up the mean coverage over that window, and builds a 101-bin weight table (one bin per integer GC percentage, 0–100%). Each bin's weight is its **median** window coverage divided by the median over all fitted windows. The median is used because a few pileup windows (collapsed repeats) can move a bin's mean arbitrarily far, and a simulator working from a single-copy reference cannot reproduce them. Every fit also writes `<output_file>.bins.tsv`, listing each bin's window count, median, mean, weight, and whether the weight was fitted or interpolated. The resulting model can be passed to `gen-reads` to make fragment start positions favour regions whose GC content matches the coverage bias observed in your data.

```bash
$ eidolon gen-gc-bias-model -c gen_gc_bias_model_config.yml
```

The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `gc_bias_model` config key.

Copy `template_config/gen_gc_bias_model.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: reference FASTA used to compute GC content
# Both plain and gzipped (.fa.gz) references are accepted.
reference: /path/to/reference.fa

# required: aligned BAM. Coverage is accumulated in process per reference base.
bam_file: /path/to/aligned.bam

# optional: minimum mapping quality applied while accumulating coverage.
# Default 0 matches `samtools depth` defaults.
min_mapq: 0

# optional: BED file restricting model inference to target regions
# Use "." (dot) to use the entire reference (default)
bed_file: .

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/gc_bias_model.json.gz

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false

# Window size in base pairs for GC% and coverage averaging.
# Should match the typical read length (or fragment length for long reads).
# Default: 100
window_size: 100

# Step between successive windows. Must not exceed window_size.
# Use the same value as window_size for non-overlapping windows (recommended).
# Values smaller than window_size produce overlapping windows.
# Default: window_size
window_stride: 100

# Bins with fewer windows than this are not fitted from their own windows: they take the
# weight interpolated between the nearest fitted bins, or the nearest fitted bin's weight
# past either end. Increase this to require more evidence before trusting a bin.
# Default: 10
min_windows_per_bin: 10
```

**BAM filtering:** unmapped, secondary, and supplementary records are skipped automatically. Reads with mapping quality `<= min_mapq` are dropped before contributing to coverage. The default `min_mapq: 0` reproduces `samtools depth` defaults; set it to 20 if your downstream `gen-reads` runs assume mapq-filtered coverage.

**Large genomes:** The tool walks the BAM once and allocates a per-contig depth array (`u32` per reference base) only for contigs that receive at least one record. Peak memory for hg38 at 30× is on the order of all-contigs-summed depth arrays — see the HPC section below for sizing.

**Long reads:** Set `window_size` to approximately your typical read length (e.g. 5000–50000 for ONT or PacBio). Larger windows mean fewer windows per contig and faster model building. Non-overlapping windows (`window_stride: window_size`) are recommended so each observation is independent.

**If the model has no effect in gen-reads:** If every GC% bin in your reference has fewer observations than `min_windows_per_bin`, all weights will be neutral (1.0) and no bias will be applied. Soft-masked (lowercase) reference bases are counted as the bases they are. A BAM that covers only part of the reference, such as an exome or a panel, is refused when the median window has no coverage: set `bed_file` to the regions it covers. This is logged as a warning. To diagnose, lower `min_windows_per_bin` or increase the region covered (use a larger reference or remove the BED restriction).

## Using the model in gen-reads

```yaml
# in gen_reads_config.yml
gc_bias_model: /path/to/gc_bias_model.json.gz

# optional: inflate fragment count to compensate for low-weight regions
# true (default): total coverage stays close to the requested depth
# false: coverage in low-GC regions will be lower than requested
gc_bias_normalize_coverage: true
```
