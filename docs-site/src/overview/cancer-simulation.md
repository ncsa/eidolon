# Cancer simulation

`eidolon` simulates tumor / normal sequencing data end-to-end. The
`eidolon gen-cancer-reads -c <config.yml>` subcommand runs two `gen-reads` passes —
one normal-genotype, one tumor — over the same reference and merges them at a
configurable purity into a single "tumor biopsy" FASTQ that downstream somatic
callers (Mutect2, Strelka, Manta, …) consume directly, plus an origin-tagged truth
VCF (`INFO/EIDOLON_ORIGIN ∈ {germline, somatic, shared}`) for scoring. It also
generates the foundational structural-variant types cancer SVs depend on — `<BND>`
translocations, `<INV>` inversions, and de novo `<INS>` — and ships bundled
pan-cancer and per-tissue (BRCA / skin / lung) models plus Docker-based benchmark
pipelines.

**See [`docs/cancer_howto.md`](../guides/cancer_howto.md) for a copy-paste guide** with
worked examples, output reference, benchmarking, and model training. Design
rationale and calibration caveats live in
[`docs/cancer_simulator.md`](../guides/cancer_simulator.md).
