# SV polish roadmap — rebuilding confidence from the reads up

**Status: plan, not results.** Nothing here is verified until the gate it describes has been
run and its output recorded. Written 2026-08-14 after job 21072620.

## Why this exists

Two gaps, found together and easily confused:

1. **Scale.** Every input-VCF fidelity cell in `docs/sv_support_matrix.md` was measured on
   H1N1 — 8 contigs, 1–2.3 kb each, largest event ~1.2 kb. The whole-genome campaigns plant
   events up to **9.1 Mb**. Fidelity at real scale has never been checked; only *de novo*
   generation has whole-genome evidence, and only as scored by callers.
2. **Depth of evidence.** "A caller found it" is the only read-level evidence for most types.
   That is real evidence — a caller cannot detect what is not in the reads — but it is
   end-to-end and opaque. It cannot distinguish *correct* from *detectable*, and it says
   nothing about whether the artifacts a researcher actually inspects are interpretable.

The second is the one this roadmap addresses first, because a scale sweep over machinery we
have not verified in detail just produces more numbers of unknown meaning.

Related and tracked separately: [#537](https://github.com/ncsa/eidolon/issues/537) — subclonal
CCF never reaches SVs.

## Confidence ranking

Ordering is by **strength of evidence and simplicity of mechanism**, not by recall. A recall
figure is one observation; a simple mechanism with two independent confirmations is a stronger
position.

### Reading in (`input_vcf`)

| rank | type | why here |
|---|---|---|
| 1 | **literal DEL** | simplest possible mechanism — a reference slice minus bases. `D=49` vs 20 baseline |
| 2 | **`<DEL>` symbolic** | depth 0.54 het / 0.00 hom exact; whole-genome recall 0.83–0.86. Known ~8% dosage bias ([#499](https://github.com/ncsa/eidolon/issues/499)) |
| 3 | **`<DUP>`** | 1.61 / 2.17 measured; fragment-spanning separately verified. Same #499 bias; one residual (missing reads would be invisible) |
| 4 | **`<INV>`** | interior 1.00, junction dip explained and bounded; homozygous-junction test exists. Dip realism never compared to real data |
| 5 | **`<CNV>`** | semantics confirmed (CN total, `CN/ploidy`, GT ignored). Depth-only — no junction signature exists to cross-check against |
| 6 | **BND inter/intra** | 30 chimeric reads; whole-genome geometry clean (50 parsed, 0 unpaired, 0 mispaired). Ranked below CNV **only** for soak time — the geometry was wrong until v3.1.0 |
| 7 | **literal INS ≤150 bp** | works, but the ceiling is a read length |
| — | **BND + inserted seq** | ❌ [#498](https://github.com/ncsa/eidolon/issues/498) insert dropped from reads |
| — | **`<INS>` symbolic** | ❌ [#500](https://github.com/ncsa/eidolon/issues/500) silent no-op |
| — | **BND unpaired `A.`** | ❌ [#500](https://github.com/ncsa/eidolon/issues/500) destroys coverage |
| — | **literal INS ≳200 bp** | ❌ [#516](https://github.com/ncsa/eidolon/issues/516) head only |

### De novo

| rank | type | why here |
|---|---|---|
| 1 | **DEL** | recall 0.828 / 0.862 across two independent callers |
| 2 | **DUP** | 0.875 / 0.875 — identical across callers, which is worth a glance in its own right |
| 3 | **INV** | recall 1.000 / 1.000; precision artifact understood and attributable to our representation |
| 4 | **BND (cancer)** | 0.920 Manta, rates reproduced from PCAWG. Delly's 0.440 is a caller limitation per the pipeline's own guide |
| 5 | **CNV** | 4 of 6 by direction — a denominator of 6 is thin |
| 6 | **INS** | ~3 draws per 30× genome, ~1 surviving the #516 cap. Not validatable from single replicates |

## The three gates

Each SV type passes three gates before it is called done. A gate is a **falsifiable
assertion**, and every gate needs a case where it must *not* fire.

### Gate 1 — FASTQ: the reads carry the signal

The reads are the ground truth of a simulator. Everything else is derived.

| type | expected signal in the reads |
|---|---|
| DEL | reads from the affected haplotype spanning the junction carry the novel left→right adjacency; **no** read carries interior sequence from that haplotype |
| DUP | reads spanning the novel tandem junction carry end-of-copy → start-of-copy adjacency |
| INV | junction reads at **both** breakpoints carry the reverse-complemented adjacency |
| CNV | depth only — no junction signature exists, which is itself worth asserting |
| BND | chimeric reads carry piece A + piece B in the orientation the ALT bracket form declares |
| INS | reads carry novel bases from the **middle**, not only the head (the #516 lesson) |

Must-not-fire for every row: a no-variant control run must show none of it.

### Gate 2 — BAM: a researcher can interpret it

Two distinct artifacts, and both matter:

- **The golden BAM** is eidolon's own claim about where each read came from. Its CIGARs must
  reflect the variant (a `D` op spanning a deletion, an `I` op for an insertion), and read
  names must mark chimeras.
- **The realigned BAM** — the FASTQ put through a real aligner — is what a researcher actually
  works from. It must carry the signatures callers depend on:

| type | signature a caller needs |
|---|---|
| DEL | discordant pairs with insert size ≈ nominal + L; split reads clipped at both breakpoints; depth drop |
| DUP | outward-facing (everted) pairs; split reads at the junction; depth rise |
| INV | same-orientation (FF/RR) pairs; split reads at both junctions; depth flat over the interior |
| BND | split reads with `SA` tags pointing at the mate locus; cross-contig discordant pairs |
| INS | soft clips accumulating at the insertion point; assemblable novel sequence |

**This gate has never been run.** We have only ever asked whether a caller found something,
never whether the evidence it needs is present and well-formed. That is the same gap #516 hid
in for three campaigns.

### Gate 3 — VCF: the record is correct and spec-conformant

Per type: `SVTYPE`, `END`, `SVLEN` sign and magnitude, ALT form, anchor-base convention, `GT`,
and for BND both `MATEID` pairing and bracket orientation. Cross-check against what Manta and
Delly emit for the same event class — the truth should be representable in the same terms it
will be scored against.

## Order of work

One type at a time, in the ranking order above, **starting with the most confident**. Starting
with a type expected to pass validates the gate methodology itself; a surprise there is
information about the harness, not the SV. Starting with a known-broken type would leave a
failure ambiguous between "the SV is broken" and "the gate is wrong".

```
DEL → DUP → INV → BND → CNV → INS
```

Each type is done when all three gates pass **and** the gate tests are shown non-vacuous by
mutation. The broken cells (#498, #500, #516) are expected to fail their gates; that is the
point, and the gate is what will define "fixed" for them.

## Phase 1 — realistic sizes in realistic contexts

Only after a type has passed its gates on a fast fixture:

1. **Size sweep: 1 kb / 100 kb / 1 Mb.** The first question this answers is whether
   [#499](https://github.com/ncsa/eidolon/issues/499)'s ~8% over-delivery is **size-dependent**.
   If it is a boundary effect it should scale as roughly `fragment_length / event_length` and
   nearly vanish at 1 Mb; if it holds at 8% it is a multiplier bug. Either answer localises a
   mechanism the issue currently lacks, and the measurement is the same either way.
2. **Repeat context.** H1N1 has no segmental duplications, no centromere, no N-gaps. Real SVs
   are overwhelmingly repeat-mediated. Plant the same event inside a segdup, inside a simple
   repeat, and in unique sequence, and compare.
3. **N-gaps.** What happens when an event overlaps an assembly gap — is it planted, refused, or
   silently mangled?

## Phase 2 — scorers at scale

Only once Phase 1 is clean. Whole-genome campaigns with truvari + Manta + Delly, as today, but
now with the ability to attribute a low recall to a caller rather than guess. The existing
selftest/decoy calibration stays as the floor.

## Prerequisites and open decisions

- **No aligner is installed on the workstation** (samtools only). Gate 2 needs one. Installing
  `bwa-mem2` or `minimap2` locally keeps the fast loop; otherwise every Gate 2 run costs a
  Delta round trip. This is a development-environment change only — it does not touch the
  shipped artifact or the conda recipe.
- **A ≥2-contig real reference** is required for BND (chr22 alone cannot plant one — de novo
  BND is inter-chromosomal by design). chr20 + chr21 is the natural pair.
- **BED-restricted windows** should keep Phase 1 near the H1N1 loop rather than the 1174
  core-hour whole-genome loop.

## What this roadmap does not cover

- Germline SV realism. De novo BND for human germline stays ⛔ not recommended regardless.
- [#537](https://github.com/ncsa/eidolon/issues/537), subclonal CCF for SVs — a cancer-model
  gap on its own axis; the gates here are dosage-agnostic.
- Caller *tuning*. We measure what stock callers do; making them do better is not our problem.
