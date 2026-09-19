# Building multiple BAM-derived models in one pass
If you need both a fragment length model and a GC bias model from the same BAM, `eidolon gen-bam-models` walks the BAM once and feeds both observers from the single pass — half the I/O of running the two per-tool commands back-to-back.

```bash
$ eidolon gen-bam-models -c gen_bam_models_config.yml
```

Copy `template_config/gen_bam_models.yml` and fill in the sections for whichever models you want:

```yaml
bam_file: /path/to/aligned.bam
min_mapq: 0   # walker-level filter; FragLengthObserver enforces MAPQ > 10 internally

frag_length:
  output_file: /path/to/frag_length_model.json.gz
  overwrite_output: false
  min_reads: 2

gc_bias:
  reference: /path/to/reference.fa
  output_file: /path/to/gc_bias_model.json.gz
  overwrite_output: false
  bed_file: .
  window_size: 100
  window_stride: 100
  min_windows_per_bin: 10
```

Both sections are optional but at least one must be present. The output models are interchangeable with what `gen-frag-length-model` and `gen-gc-bias-model` would produce on their own, so the resulting files plug into `gen-reads` exactly the same way.

**One filter for the shared walk.** The top-level `min_mapq` is the only MAPQ knob exposed to the unified walker, and it gates *both* observers. The standalone tools differ: `gen-frag-length-model` hard-codes its own MAPQ > 10 cutoff (via `BamWalkFilter::for_frag_length()`), while `gen-gc-bias-model` honors the `min_mapq` you pass it. If you want the unified runner's frag-length output to match the standalone command byte-for-byte, set `min_mapq: 10`; if you want the GC bias output to match a standalone run with `min_mapq: 20`, set `min_mapq: 20`. When the two policies need to differ, run `gen-bam-models` twice — one config per output — and the per-observer BAM iteration savings still beat running the standalone commands.
