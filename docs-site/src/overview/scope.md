# Scope: germline is general, somatic is human

These two halves of `eidolon` have deliberately different reach, and it is worth knowing which
one you are using.

**Germline simulation is species-agnostic.** `gen-reads` needs a reference and, optionally,
models trained from your own data — nothing in that path assumes a human genome. It has been
exercised on H1N1 (13.5 kb, eight contigs), *E. coli*, yeast, soy, and human references from a
single chromosome up to whole-genome GRCh38. If you are simulating germline variation for any
organism, that is a supported use.

**Somatic / cancer simulation is specialized to human.** The bundled mutational models are
derived from human cancer cohorts — COSMIC signatures and PCAWG structural-variant
distributions — and the per-tissue models (BRCA, skin, lung) are human tissue types. SV size
and type distributions, mutation rates, and signature weights all carry that provenance.

`gen-cancer-reads` will *run* against a non-human reference and produce output, but the model
driving it has no meaning there: a yeast genome given a human pan-cancer SV spectrum produces
numbers, not biology. **We do not claim cancer simulation for non-human references**, and the
validation campaigns behind those models are human-only. Use another organism's reference for
cancer work only if you have trained your own models against matched data, and treat the
bundled ones as human-specific.

This split is also why the two halves are validated differently: germline correctness is
checked across deliberately varied contig shapes and sizes, while somatic realism is measured
against real human sequencing data.
