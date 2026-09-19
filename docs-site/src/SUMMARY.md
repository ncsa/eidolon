# Summary

[The eidolon project](index.md)

---

# Overview

- [Scope: germline is general, somatic is human](overview/scope.md)
- [Cancer simulation](overview/cancer-simulation.md)
- [How `eidolon` compares to NEAT](overview/comparison-with-neat.md)

# Getting started

- [Prerequisites](getting-started/installation.md)

# Generating reads

- [Fastq Output](gen-reads/fastq-output.md)
- [BAM Output](gen-reads/bam-output.md)
- [Targeted Generation with a BED File](gen-reads/targeted-generation.md)
- [Custom Mutation Rate Regions with a BED File](gen-reads/mutation-rate-regions.md)
- [Input Variants VCF](gen-reads/input-variants-vcf.md)
- [FASTQ Shuffling](gen-reads/fastq-shuffling.md)
- [3′ Adapter Readthrough](gen-reads/adapter-readthrough.md)
- [Parallel Processing](gen-reads/parallel-processing.md)
- [Reading Bed Data](gen-reads/reading-bed-data.md)

# Cancer simulation

- [Cancer simulation how-to](guides/cancer_howto.md)
- [Cancer simulator design](guides/cancer_simulator.md)

# Models

- [Generating a Mutation Model](models/mutation-model.md)
- [Generating a Sequencing Error Model](models/sequencing-error-model.md)
- [Generating a GC Bias Model](models/gc-bias-model.md)
- [Generating a Fragment Length Model](models/fragment-length-model.md)
- [Building multiple BAM-derived models in one pass](models/bam-models.md)
- [Model-builder resource baseline](guides/model_builder_baseline.md)
- [PCAWG SV measurements](guides/pcawg_sv_measurement.md)

# Other subcommands

- [Filtering Your Data](subcommands/filter-reads.md)
- [Comparing a caller's VCF against the golden VCF](subcommands/compare-vcfs.md)

# HPC and validation

- [Running on HPC](hpc/running-on-hpc.md)
- [Running eidolon on an HPC cluster](guides/hpc_guide.md)
- [Re-validation / regression protocol](guides/regression_protocol.md)

# Reference

- [Versioning and the public API](reference/versioning.md)
- [SV support matrix](guides/sv_support_matrix.md)
- [Current Benchmarks](reference/benchmarks.md)

# Upgrading

- [Upgrading from 2.0.0](upgrading/from-2-0-0.md)

---

[Project documents not in this site](appendix/other-documents.md)
