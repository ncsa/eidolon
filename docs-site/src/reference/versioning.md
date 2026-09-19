# Versioning and the public API

As of **v3.0.0**, eidolon follows [Semantic Versioning 2.0.0](https://semver.org): a
MAJOR bump means something you may depend on changed incompatibly, a MINOR bump adds
functionality compatibly, and a PATCH bump is fixes only. Releases before v3.0.0 were
versioned less strictly — notably v1.11.0 and v1.12.0 changed the FASTQ read-name format
in minor bumps.

SemVer requires a project to say what its public API is. For eidolon:

## Public API — a change here means a MAJOR bump
- Names and semantics of emitted VCF INFO tags, the VCF sample-column name, and the
  `FILTER` / `FORMAT` fields eidolon writes
- FASTQ/BAM **read-name (QNAME) format** — the `EIDOLON_generated_` /
  `EIDOLON_chimeric_` prefixes and the positional fields encoded after them
- CLI subcommand names, flags, and configuration-YAML keys
- Model-file compatibility, **in the backward direction only** — a model built by an
  older eidolon stays readable by a newer one. See *Model files* below for what that
  does and does not cover

## Not public API — may change in a MINOR or PATCH release
- The `eidolon-core` Rust library surface; it exists to serve the binary
- Log output, progress reporting, and human-readable messages
- The exact simulated *content* for a given seed — reads are a random draw, so sampler
  and model changes legitimately alter output while preserving the format

## Versioned separately from eidolon itself
- The `compare-vcfs` JSON report carries its own `schema_version` field (currently `1.3.0`),
  bumped only when a backward-incompatible field change lands. Read that field rather than
  inferring the report's shape from the eidolon version — the two move independently, so a
  minor eidolon release can carry a new schema, and a major one need not

If you parse eidolon's output, the first list is what you are relying on, and a major
version bump is your signal to check this file before upgrading.

## Model files

"Compatible" covers one direction, and neither direction is currently checked automatically.

- **Old model, newer eidolon — supported.** Fields added since the model was written
  deserialize to their shipped defaults, so the model keeps its previous behavior for
  anything it does not carry, rather than failing to load.
- **New model, older eidolon — NOT supported, and NOT detected.** Unrecognized fields are
  ignored silently. An older binary will read a newer model, discard what it does not
  understand, and simulate using its own defaults with no warning that it did so. Do not move
  model files backwards across releases.
- **Model files carry no version stamp.** A model cannot tell you which eidolon wrote it, so
  neither direction above can be verified — by you or by eidolon. Tracked in #708.

The practical rule: rebuild models with the eidolon you intend to run, and treat a model file
as belonging to the version that produced it.
