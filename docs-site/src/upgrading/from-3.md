# Upgrading from 3.x

v4.0.0 has three changes that can stop a run which completed under 3.x. Each fails before
any output is written and says what to change.

| What you see | Cause | What to do |
|---|---|---|
| A config error naming an unknown key | A key the subcommand does not read. 3.x ignored it and left that setting at its default. | Remove the key, or rename it to the suggested one. The error lists every accepted key, and `template_config/` has a template per subcommand. |
| `rneat: command not found` | The deprecated `rneat` alias is removed. | Run `eidolon` instead. Same subcommands, same flags. |
| `gen-mut-model`: "The VCF record at … has REF …, but the reference has … there" | The training VCF disagrees with the reference: a different build, a liftover problem, or a few bad records. | Check that the VCF belongs to this reference. If it does, drop the mismatching records and fit on the result (below). |

## Unknown config keys

Every top-level key must be one the subcommand reads; so must the keys inside `adapters:`
(`gen-reads`) and `frag_length:` / `gc_bias:` (`gen-bam-models`). A key whose value is `.` is
checked too.

A rejected key never did anything. In 3.x it was ignored, so the run used that setting's
default without warning: `tumor_mutation_model:`, meant as `tumor_model:`, simulated the tumor
with the germline model. **If a 3.x config now fails, the 3.x results it produced were made
without that setting.** Re-run anything that matters with the corrected key.

For a NEAT config, `gen-reads` names the eidolon equivalent where one exists (`include_vcf` →
`input_vcf`, `error_model` → `sequence_error_model`, …) and otherwise says the option is not
supported. Keys are not translated automatically.

## `gen-mut-model` and REF mismatches

Each SNP's and indel's REF must match the reference at its position. The fit stops at the
first record that does not, naming it and both bases. If the VCF does belong to this reference
and only a few records disagree, drop them explicitly:

```bash
bcftools norm -f reference.fa --check-ref x in.vcf.gz -Oz -o checked.vcf.gz
```

If many records disagree, the VCF was most likely called against a different build. Fit
against the matching reference instead.

Two related changes alter a model fitted from the same VCF:

- A SNP at a contig edge (no base on one side) is left out of the mutation rate. 3.x counted it.
- A `bed_file` that names no contig in the reference stops with both lists of contig names,
  usually a `chr1` versus `1` mismatch, instead of "Unknown error".

## Output that changes without any action

These are fixes, not compatibility breaks, but they change output for the inputs they affect.
Re-baseline anything measured on such inputs under 3.x:

- **Soft-masked (lowercase) references** now simulate exactly like uppercase ones.
- **`input_vcf` runs that combine a literal deletion of 50 bp or more with a `<DUP>`, `<DEL>` or
  `<CNV>`** no longer over-produce reads past the deletion.
- **`gen-cancer-reads` with the bundled COSMIC models** places somatic SNVs slightly differently;
  the models were rebuilt.
- **`gen-seq-error-model` with `max_reads`** samples the whole FASTQ, so a capped fit gives a
  different model and takes longer.
- **`filter-reads` output** is always gzip-compressed, so a plain `x.fastq` produces
  `x_filter.fastq.gz`.
- **Golden BAM records** gain an `NM` tag.

See `CHANGELOG.md` for the full list.
