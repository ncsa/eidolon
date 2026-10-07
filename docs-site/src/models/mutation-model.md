# Generating a Mutation Model
`eidolon` now includes the ability to read in real data and learn the parameters to reproduce it in a simulation run. So far, this is limited to the mutation model, using the `eidolon gen-mut-model` subcommand.

```bash
$ eidolon gen-mut-model -c gen_mut_model_config.yml
```

The inputs are a single-sample VCF and the reference FASTA the VCF was called against. `eidolon` computes statistics for indels and SNPs and builds the trinucleotide model of SNP generation, as in the original `eidolon`. The output is a gzipped JSON model file that can be passed directly to `gen-reads` via its `mutation_model` config key.

Copy `template_config/gen_mut_model_template.yml` to a directory of your choosing and fill in the fields:

```yaml
# required: path to the reference FASTA used to call the VCF
reference: /path/to/reference.fa

# required: path to a single-sample VCF file
# The VCF must include FORMAT and SAMPLE columns, and GT must be present in FORMAT.
# QUAL=. is accepted and treated as 0.
vcf_file: /path/to/variants.vcf

# required: output path for the generated model; should end in .json.gz
output_file: /path/to/output_model.json.gz

# optional: restrict model learning to regions in this BED file
# set to . (dot) to use the entire reference (default)
bed_file: .

# optional: set to true to overwrite an existing output file (default: false)
overwrite_output: false
```

## VCF requirements
- Single sample only — multi-sample VCFs are not yet supported (tracked in #412)
- Each variant record must have `GT` in the `FORMAT` column; `eidolon` hard-errors if GT is missing
- `QUAL=.` is accepted and treated as quality score 0
- Insertion and deletion lengths are each fitted from the VCF. If one type is absent (an insertions-only VCF, say), that type's length distribution is the default one and a warning says so; the type that was observed keeps its fitted lengths. With no indels at all, both use the defaults.
- Each SNP's and indel's `REF` must match the reference at its position. If any record disagrees, `gen-mut-model` errors out, giving the number of SNPs and indels affected and the first few positions. A mismatch usually means the VCF was called against a different build, or lifted over badly. If the VCF does belong to this reference, drop the mismatching records deliberately and fit on the result:
  ```
  bcftools norm -f reference.fa --check-ref x in.vcf.gz -Oz -o checked.vcf.gz
  ```
- A SNP at a contig edge has no trinucleotide context. It is left out of the model entirely (not counted toward the mutation rate or the SNP/indel split), with a warning giving the count. If **every** SNP is at an edge, `gen-mut-model` errors out rather than writing a model with no SNP context data.
- A `bed_file` must name at least one contig in the reference. If none match, the error lists both sets of names; a `chr1` versus `1` mismatch is the usual cause.

Caveats: Only one sample can be read at this point (#412). Currently, high-mutation regions and common variants features from Python NEAT are not yet implemented (#413).

Try it out and let us know if you run into any issues!

We will continue to improve this process in the future. If you have suggestions for parsing names to facilitate this, let us know in the Issues section!
