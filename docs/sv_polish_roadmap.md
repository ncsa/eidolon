# SV polish roadmap — rebuilding confidence from the reads up

**Status: plan, not results.** Nothing here is verified until the gate it describes has been
run and its output recorded. Written 2026-08-14 after job 21072620.

**Updated 2026-08-29 (v3.2.0).** Every cell this document listed as broken has since been
fixed and validated at genome scale — Delta jobs 21575385 and 21603825, chr20+21+22, 30x,
purity 0.6, both exit 0. The ranking below is kept because the *reasoning* about strength of
evidence still holds; the verdicts are updated. A stale roadmap is worse than no roadmap: for
several days this file told readers that deletions were only partially realized while the
gates were reporting 76 of 76.

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
| 2 | **`<DEL>` symbolic** | depth 0.54 het / 0.00 hom exact; whole-genome recall 0.83–0.86. The ~8% dosage bias once attributed to [#499](https://github.com/ncsa/eidolon/issues/499) does **not** reproduce at realistic event size — 0.98 +/- 1.9% on a 10 kb event over a 1 Mb reference; the original was a 1.2 kb event on a 13.5 kb fixture where depth noise alone is ~13% sigma. The residual mechanism is [#584](https://github.com/ncsa/eidolon/issues/584) |
| 3 | **`<DUP>`** | 1.61 / 2.17 measured on H1N1; fragment-spanning separately verified. Those figures carry [#584](https://github.com/ncsa/eidolon/issues/584) — junction reads are emitted **on top of** coverage-multiplied reads, an excess independent of event length and therefore worst on short events. Accurate at 10 kb, biased high at ~1 kb |
| 4 | **`<INV>`** | interior 1.00, junction dip explained and bounded; homozygous-junction test exists. Dip realism never compared to real data |
| 5 | **`<CNV>`** | semantics confirmed (CN total, `CN/ploidy`, GT ignored). Depth-only — no junction signature exists to cross-check against |
| 6 | **BND inter/intra** | 30 chimeric reads; whole-genome geometry clean (50 parsed, 0 unpaired, 0 mispaired). Ranked below CNV **only** for soak time — the geometry was wrong until v3.1.0 |
| 7 | **literal INS ≤150 bp** | works, but the ceiling is a read length |
| — | **BND + inserted seq** | ✅ fixed ([#498](https://github.com/ncsa/eidolon/issues/498)) — spliced between the reference pieces; insert comes out of the fragment budget, not added to it |
| — | **`<INS>` symbolic** | ✅ fixed ([#500](https://github.com/ncsa/eidolon/issues/500)) — realized with synthesised novel sequence from the same sampler the de novo path uses |
| — | **BND unpaired `A.`** | ⚠ **rejected**, not supported ([#500](https://github.com/ncsa/eidolon/issues/500)) — no longer destroys coverage, but a legitimate VCF 4.2 call eidolon refuses ([#623](https://github.com/ncsa/eidolon/issues/623)) |
| — | **literal INS ≳200 bp** | ✅ fixed ([#516](https://github.com/ncsa/eidolon/issues/516)) — 38/38 verified in reads, 62–2127 bp, Delta job 21382756 |
| — | **literal INS in the golden BAM** | ✅ fixed ([#589](https://github.com/ncsa/eidolon/issues/589)) — 11/11 present in the reads, 9 longer than a 151 bp read, largest 968 bp (job 21603825) |
| — | **literal DEL ≳ read length** | ✅ fixed ([#590](https://github.com/ncsa/eidolon/issues/590)) — deletions live on the haplotype; **76/76 junction sequence present, 76/76 coverage removed**, flank ratios 0.33–0.77 at purity 0.6 (job 21603825) |

### De novo

| rank | type | why here |
|---|---|---|
| 1 | **DEL** | recall 0.828 / 0.862 across two independent callers |
| 2 | **DUP** | 0.875 / 0.875 — identical across callers, which is worth a glance in its own right |
| 3 | **INV** | recall 1.000 / 1.000; precision artifact understood and attributable to our representation |
| 4 | **BND (cancer)** | 0.920 Manta, rates reproduced from PCAWG. Delly's 0.440 is a caller limitation per the pipeline's own guide |
| 5 | **CNV** | 4 of 6 by direction — a denominator of 6 is thin |
| 6 | **INS** | The #516 cap is gone; every draw is now realized. Still sparse per replicate — chr22 at `SV_RATE_SCALE=30` yields 4, chr20+21+22 yields 11, scale 100 yields 38. Manta recovers 2/11 all-calls, a caller ceiling with the reads verified to carry all 11 |

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

### Gate 2a — the realigned BAM: a caller can interpret it

The FASTQ put through a real aligner is what a researcher actually works from. It must carry the
signatures callers depend on:

| type | signature a caller needs |
|---|---|
| DEL | discordant pairs with insert size ≈ nominal + L; split reads clipped at both breakpoints; depth drop |
| DUP | outward-facing (everted) pairs; split reads at the junction; depth rise |
| INV | same-orientation (FF/RR) pairs; split reads at both junctions; depth flat over the interior |
| BND | split reads with `SA` tags pointing at the mate locus; cross-contig discordant pairs |
| INS | soft clips accumulating at the insertion point; assemblable novel sequence |

**Status: DONE for DEL, DUP, INV, BND** (`eidolon/tests/gate2_realigned_*.rs`, PRs #542–#545).
Each passes, each was mutation-verified, and each produced a finding beyond the pass — depth-only
checking is insufficient (DEL), a measured mechanism for #499 (DUP), the reads cleared of the
precision artifact (INV), and undocumented copy-number semantics (BND, #546). CNV and INS remain.

Before these, we had only ever asked whether a caller found something, never whether the evidence
it needs is present and well-formed — the same gap #516 hid in for three campaigns.

Requires `bwa-mem2`, so these tests are `#[ignore]`d and CI never runs them. Run deliberately:
`conda activate aln && cargo test --test gate2_realigned_del -- --ignored --nocapture`.

### Gate 2b — the golden BAM agrees with the FASTQ

**Its own phase, not a missing piece of 2a** ([#548](https://github.com/ncsa/eidolon/issues/548)).

The golden BAM is eidolon's claim about where each read came from. It and the FASTQ are produced
by the same generation but written by different code (`fastq_tools.rs`, `bam_writer.rs`), and
nothing asserts they agree — the shape that let `sv_model.rs` and `runner.rs` disagree about BND
geometry for eight releases.

Separated because **it needs no aligner.** Consistency between two eidolon outputs is checkable
from the outputs alone, so unlike 2a it can be an ordinary CI test that never gets skipped. The
strongest assertion needs no variant at all: for each QNAME the golden BAM's SEQ must equal the
FASTQ's (reverse-complemented for reverse-strand records). Then CIGARs reflecting the variant
(`D` for a deletion, `I` for an insertion), and chimeric read marking.

Deliberately ranked **above** the remaining 2a types: cheaper, no aligner, runs in CI, and covers
every SV type at once rather than one per gate.

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
mutation. The broken cells (#498, #500, and now #589/#590 — #516 is fixed) are expected to fail their gates; that is the
point, and the gate is what will define "fixed" for them.

## Phase 1 — realistic sizes in realistic contexts

Only after a type has passed its gates on a fast fixture:

1. ~~**Size sweep: 1 kb / 100 kb / 1 Mb.**~~ **RUN, 2026-08-22, during the
   release/fragment-placement investigation.** It is a boundary effect, and it does nearly
   vanish at realistic size: 1200bp event / ~400bp flanks (H1N1) needed a 1.13x correction;
   the SAME 1200bp event with megabase flanks (a real chr22 window) needed only 1.07x; at a
   realistic 100kb event size with megabase flanks, 1.003x -- see
   `eidolon::tests::sv_support_matrix::depth_modulation_is_accurate_at_realistic_scale`,
   which pins this at 1Mb scale as a permanent regression guard, and
   `docs/claude_engineering_audit.md` §5.6's 2026-08-22 addendum for the full measurement.
   Two mechanisms turned out to contribute, not one: #499's own (chimeric junction reads
   landing outside the coverage-multiplied budget) plus a second, independent one found by
   this investigation (fragment placement legitimately extending across a coverage-
   multiplier boundary, correct and necessary for removing an artificial dead zone there,
   but costing some of the segment's own declared depth to redistribution). Both are real,
   both vanish at scale, and #499 was closed on this basis. **Practical takeaway: H1N1
   cannot host an event "much larger than fragment length" with wide flanks at the same
   time -- it is 2280bp total -- so any coverage-multiplier depth check on it is a fast
   mechanism check, never a precision one. Prefer events large relative to fragment length,
   and genuinely large references, whenever depth-multiplier accuracy is what's being
   validated; small references (viral genomes, small scaffolds) are fine for mechanism
   testing but should not be used to validate depth precision.**
2. **Repeat context.** H1N1 has no segmental duplications, no centromere, no N-gaps. Real SVs
   are overwhelmingly repeat-mediated. Plant the same event inside a segdup, inside a simple
   repeat, and in unique sequence, and compare.
3. **N-gaps.** What happens when an event overlaps an assembly gap — is it planted, refused, or
   silently mangled? **Partially answered, 2026-08-23.** The *read-generation* half is now
   covered and guarded: `regions_of_interest` is built from `get_non_n_regions()`, so N bases
   are excluded from generation, and fragment placement must not extend a fragment's end
   across a gap even though it deliberately extends across chunk and coverage-multiplier
   boundaries. That distinction was gotten wrong once already (review of the
   fragment-placement branch: 103 of 6000 reads carried gap sequence, against 0 before that
   branch) and is now pinned by
   `eidolon::tests::sv_support_matrix::fragments_do_not_extend_across_an_assembly_gap`,
   shown non-vacuous by mutation. **Still open:** what happens to an *SV event* whose span
   overlaps a gap — whether it is planted, refused, or silently mangled — which is a
   different code path (`sv_modulation_range` / the SV samplers) and untested.

4. **Targeted (BED) coverage semantics.** Measured 2026-08-23 while validating the
   fragment-placement rewrite: on an exome-scale BED (436 targets, mean width 178 bp,
   real chr22 sequence) a requested `coverage: 60` delivers ~24.5x on target, because
   40.5% of read bases legitimately spill past targets narrower than the fragment
   length. Documented as a caveat in the README's `target_bed` section; making
   `coverage` mean on-target depth automatically is
   [#578](https://github.com/ncsa/eidolon/issues/578). Note the design trap recorded
   there: isolated exome targets lose ~59% of their spill while *adjacent* segments (an
   SV interval and its baseline neighbours) exchange spill and lose ~0-7%, so a uniform
   inflation factor would break SV depth modulation.

## Phase 2 — scorers at scale

Only once Phase 1 is clean. Whole-genome campaigns with truvari + Manta + Delly, as today, but
now with the ability to attribute a low recall to a caller rather than guess. The existing
selftest/decoy calibration stays as the floor.

## Prerequisites and open decisions

- **The local toolchain already covers Gates 1–3.** `bwa-mem2` lives in the `aln` conda
  environment, `samtools`/`bcftools`/`hap.py`/`som.py` in `hap_py_env`, and samtools 1.22 in
  `samtools122`. Check `conda env list` before concluding a tool is missing — several are
  installed only inside an environment. What is **not** local: truvari, Manta, Delly, so
  Phase 2 still needs Delta.
- **A ≥2-contig real reference** is required for BND (chr22 alone cannot plant one — de novo
  BND is inter-chromosomal by design). chr20 + chr21 is the natural pair.
- **BED-restricted windows** should keep Phase 1 near the H1N1 loop rather than the 1174
  core-hour whole-genome loop.

## What this roadmap does not cover

- Germline SV realism. De novo BND for human germline stays ⛔ not recommended regardless.
- [#537](https://github.com/ncsa/eidolon/issues/537), subclonal CCF for SVs — a cancer-model
  gap on its own axis; the gates here are dosage-agnostic.
- Caller *tuning*. We measure what stock callers do; making them do better is not our problem.
