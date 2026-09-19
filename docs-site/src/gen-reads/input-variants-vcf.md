# Input Variants VCF
You can supply a VCF of variants to force into the simulation:

```yaml
input_vcf: /path/to/variants.vcf.gz
```

`eidolon` will place every variant from the VCF into the corresponding position in the simulated reads and the output VCF. Random variants are still generated at `mutation_rate` for all positions not covered by the input VCF; set `mutation_rate: 0.0` to disable random variants entirely and output only the provided set of variants.

## Requirements

- The VCF must be single-sample (one sample column).
- Every record must include `GT` in the FORMAT field. `eidolon` uses the genotype to determine whether to apply the variant to all reads covering the position (homozygous, e.g. `1/1`) or only a probabilistic subset (heterozygous, e.g. `0/1`). Records without `GT` are rejected.
- Contig names must match the short names derived from the reference FASTA (text after `>` up to the first whitespace character). Variants on unrecognised contigs are skipped with a warning.
- Both `.vcf` and `.vcf.gz` files are accepted.

## Supported variant types

| Type | Condition | Handled |
|------|-----------|---------|
| SNP | REF and ALT both single base | Yes |
| Insertion | single-base REF, multi-base ALT | Yes |
| Deletion | multi-base REF, single-base ALT | Yes |
| Symbolic SV | ALT is `<DEL>` / `<DUP>` / `<CNV>` / `<INS>` / `<INV>` / breakend / other `<TAG>` | Yes — see "Symbolic / structural variants" below |
| Literal complex | multi-base REF **and** multi-base ALT (literal bases) | **No** — skipped with warning |

## Symbolic / structural variants

Symbolic ALTs (VCF 4.2 §1.4) are accepted and round-tripped to the output VCF verbatim, with `INFO/END`, `INFO/SVLEN`, and `INFO/CN` preserved. As of v1.10, `eidolon` can also generate symbolic SVs *de novo* from a learned model — opt in by setting `sv_rate_scale: 1.0` (or higher) in your gen-reads YAML; see "De novo SV generation" below.

| SV | Effect on read depth | Effect on read sequence |
|----|----------------------|-------------------------|
| `<DEL>` | hom → ×0, het → ×(ploidy−1)/ploidy; or ×CN/ploidy if `INFO/CN` is set | none — bases in the deleted span just stop producing reads |
| `<DUP>` | hom → ×2, het → ×(ploidy+1)/ploidy; or ×CN/ploidy if `INFO/CN` is set | none — extra reads come from the forward-strand reference |
| `<CNV>` | ×CN/ploidy when `INFO/CN` is set; otherwise warned and passed through with no depth change | none |
| `<INS>` | none (insertion at a single anchor base — no span to modulate) | none — the inserted sequence is not synthesized |
| `<INV>` | none in this release | **not modeled** — reads come from the forward-strand reference, not the inverted sequence |
| Breakends, unknown `<TAG>` | none | none — round-tripped only |

DEL anchor convention: for `<DEL>`, POS is the unaffected base immediately before the deletion (per VCF 4.2), so the modulated span is `[POS+1, END]` in 1-based coords. For `<DUP>` / `<CNV>` / `<INV>`, POS is the first base of the affected region, so the span is `[POS, END]`.

When an SV zeroes out coverage (hom `<DEL>` or `INFO/CN=0`), the mutation rate over the same span is also zeroed so de-novo SNPs don't pollute the output VCF with variants that never appear in reads.

## Current caveats

- *Multi-allelic records*: only the first ALT allele is used; additional alleles are silently ignored. Split multi-allelic records with `bcftools norm -m -` before passing to `eidolon`. (This is true for both literal and symbolic ALTs — `<DEL>,<DUP>` on the same line is treated as the first ALT only.)
- *REF allele verification*: `eidolon` does not check that the REF field matches the reference sequence at that position. Mismatches will produce biologically incorrect output without any warning.
- *Literal complex variants*: records whose REF and ALT are both multi-base strings (and the ALT is literal bases, not a `<TAG>`) are skipped with a logged warning and do not appear in the output.
- *`<INV>` interior is forward-strand*: inversion **junction** reads are emitted (the chimeric signal callers use to detect the inversion), but reads from the *interior* of an inverted span are still transcribed from the forward-strand reference. The `<INV>` record round-trips to the output VCF.
- *Breakends*: `<BND>` junction reads are generated (chimeric reads spanning the mate locus; validated against Manta — see the ACCESS report), and the record round-trips to the output VCF.

## De novo SV generation

`eidolon` can sample SVs directly from a learned `SvModel` rather than relying on user-supplied records: `<DEL>`, `<DUP>`, `<CNV>`, `<INV>` and literal insertions are drawn per contig, and `<BND>` translocations are drawn separately as mate pairs across two contigs. **`<BND>` therefore needs a reference with at least two contigs** — on a single-contig reference the correct yield is zero rather than a same-contig junction, so a chr22-only run reports `truth BND: 0` by construction. Generation is **off by default** (opt-in via `sv_rate_scale`) so v1.9 pipelines remain unchanged.

To enable, add a single line to your gen-reads YAML:

```yaml
sv_rate_scale: 1.0
```

- `0.0` (the default) disables de novo SV generation entirely — only `input_vcf` records flow through.
- `1.0` reproduces the rate from the trained model.
- Larger values scale the rate proportionally for stress testing.

When enabled, `eidolon` consults `MutationModel.sv_model`:

1. **If you trained your own model** with `eidolon gen-mut-model -c <yaml>` against an SV-rich VCF (e.g. a gnomAD-SV slice), the model file includes a fitted SV component covering type / length / copy-number / homozygous frequencies. Pass it via `mutation_model: /path/to/your.json.gz` in the gen-reads YAML.
2. **If you don't supply a model**, gen-reads loads the bundled default. The default carries **approximate** gnomAD-SV v2.1 parameters — useful for kicking the tires, but not a substitute for retraining on data that matches your downstream use case. See `eidolon-core/src/models/sv_model_defaults.rs` for the parameter sources.

Per-contig sampling: Poisson count from `per_base_rate × contig_len × sv_rate_scale` → weighted type pick → log-normal length (rejected outside `[50bp, contig_len / 4]`) → uniform anchor with overlap and N-gap rejection → `INFO/CN` draw for `<CNV>` → Bernoulli for genotype.

De novo records merge with `input_vcf` SVs and flow through the same depth-modulation path, so simulated coverage and the golden VCF round-trip behave identically for both sources. `compare-vcfs` skips both into the `skipped.symbolic` bucket — symbolic ALTs aren't byte-comparable.

Caveats:
- The bundled default is a literature-derived approximation, not a refit on the actual gnomAD VCFs. Don't use it for distribution-faithful benchmarking — retrain.
- `<INV>` and breakends are not generated (the read content for either isn't modeled yet).
- A trained model's `sv_model` field is `None` if the training VCF lacked sufficient SV observations (< 2 per type after filtering). Loading that model with `sv_rate_scale > 0` is harmless — generation just no-ops.
