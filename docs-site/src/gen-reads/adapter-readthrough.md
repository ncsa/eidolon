# 3′ Adapter Readthrough
Real Illumina libraries whose insert is shorter than the read length "read through" the
fragment into the 3′ sequencing adapter — the read's tail is adapter sequence, not genome.
`eidolon` can simulate this so short-insert data looks realistic and exercises downstream
adapter-trimming QC. **Disabled by default** — output is byte-identical to prior versions
when the `adapters` block is omitted.

```yaml
adapters:
  enabled: true
  preset: truseq        # truseq | nextera | custom
  overrun: polyg        # polyg | random -- what follows the end of the adapter construct
  # for preset: custom, give explicit 5′→3′ uppercase-ACGT sequences:
  # r1: AGATCGGAAGAGC...
  # r2: AGATCGGAAGAGC...
```

When enabled, any fragment whose insert is shorter than `read_len` is kept (instead of
being resampled away) and the read is padded to `read_len` at its 3′ end with adapter
sequence — R1 gets `r1`, R2 gets `r2` (appended after R2 is reverse-complemented, matching
real read orientation). Adapter bases carry the same quality/error models as genomic bases,
are marked soft-clipped (`S`) in the golden BAM, and never introduce variants. **How much
readthrough you get is controlled by how many inserts fall below `read_len`** — i.e. by
`fragment_mean` / `fragment_st_dev` relative to `read_len`. With a long-insert setting
(`fragment_mean` ≫ `read_len`) essentially no reads reach the adapter and enabling this
changes nothing.

Each preset is the whole construct a read passes through after the insert: the adapter stem,
then a placeholder sample index, then the P7 (R1) or P5 (R2) flow-cell binding sequence. A read
longer still has run off the end of the molecule, and `overrun` sets what fills the rest:

- `polyg` (default) — G, which is how two-color instruments (NovaSeq, NextSeq, iSeq) call a
  cycle with no signal. This is the poly-G tail seen in real short-insert data.
- `random` — random bases, closer to four-color instruments (MiSeq, HiSeq).

Fill bases take quality scores and are soft-clipped like adapter bases, but take no
substitution errors, since there is no template base to misread. The adapter is never
repeated: no real molecule carries a second copy of it (#865).

Presets:
- `truseq` — Illumina TruSeq (`AGATCGGAAGAGC…`), the most common. Single-index adapter 1
  (`ATCACG`).
- `nextera` — Nextera / transposase Mosaic End (`CTGTCTCTTATACACATCT`), indexes N701/S501.
- `custom` — supply your own `r1` / `r2`. Give the full sequence each read runs into,
  including your own index and P5/P7 if you want them; `overrun` fill starts where it ends.

Not compatible with `long_reads: true` (adapter readthrough is a short-read phenomenon).

## When to enable it
- **Validating adapter-trimming pipelines** — generate reads with known adapter content to
  confirm fastp / cutadapt / Trim Galore settings actually detect and remove it.
- **Short-insert libraries** — small-RNA / cfDNA / degraded-DNA style data where inserts sit
  near or below the read length and readthrough is expected; the simulated reads then match
  what the sequencer really produces.
- **Aligner soft-clip benchmarking** — untrimmed adapter tails should be soft-clipped by the
  aligner; use this to check that behaviour end-to-end.
- **Realistic FASTQ→VCF tests** — exercise the whole trim → align → call path on data that
  carries adapters instead of idealized adapter-free reads.

To keep short inserts as **genomic** reads *without* modeling adapters (e.g. to study the
coverage behaviour of a short-insert library on its own), set `keep_short_fragments: true`
instead of enabling `adapters` — short fragments are then emitted as insert-length reads
with no adapter padding.

**Validated behaviour** (paired germline benchmark: chr22, 30×, TruSeq, 3 reps/arm,
BWA-MEM2 → GATK HaplotypeCaller vs the golden VCF): adapter readthrough is **detectable**
(fastp adapter content and aligner soft-clip fraction rise only when it is enabled) and
**fully callable** — SNP and indel recall/precision/F1 are unchanged within replication
noise whether the adapters are fastp-trimmed or left for BWA-MEM2 to soft-clip. An
adapter-free short-insert control showed the only fidelity difference versus a long-insert
baseline is the reduced effective coverage of short inserts (mate overlap), not the adapters.
That benchmark predates #865, when reads past the adapter stem repeated it instead of running
into the index, P5/P7 and poly-G; it has not been re-run on the current sequence.
