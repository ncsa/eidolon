# SV support matrix — what eidolon can reproduce, and from where

**Measured 2026-08-02** on the H1N1 fixture (8 contigs), 60× diploid, `ploidy: 2`,
heterozygous unless noted. Every row is an observation, not a design intent. Cells
marked ⚠ or ✗ are known gaps, not omissions.

Two separate questions, and they are **not** required to agree:

- **Input** — you supply the record via `input_vcf`. The contract is fidelity: what you
  put in comes out, in the truth VCF *and* in the reads.
- **De novo** — eidolon samples the record from a model. The contract is realism: the
  distribution matches the source corpus.

Input fidelity is the higher priority. A researcher modelling a specific variant needs
to know it reached the output; a mismatched *distribution* is a modelling argument, but
a variant that silently fails to appear is a broken tool.

## 2026-08-16 characterization addendum

The focused Gate 2b checks and the first Delta performance/BAM measurements are now complete.
These are baselines for the current `develop` build, not pooled recall estimates.

### Post-#537 subclonal SV smoke (develop `4a3ea4a`)

Delta jobs `21281985` and `21282371` exercised the merged subclone-SV implementation on
multi-contig *S. cerevisiae*. The higher-rate run produced 12 SV records (7 at CCF 0.2,
1 at CCF 0.5, and 4 clonal at CCF 1.0); all 12 carried `EIDOLON_CCF` and
`EIDOLON_VAF`. Three reciprocal translocation pairs shared CCFs of 0.2, 1.0, and 0.2.
After bwa-mem2 alignment, 199 chimeric alignments were retained, 100% of 2,426,389 reads
mapped, 69 were supplementary, and 99.99% were properly paired. This validates CCF
stamping, reciprocal-BND pairing, and read-level junction survival; it is not a caller
recall estimate.

The full-GRCh38 post-merge stress run (job `21282504`, 30×, purity 0.6, nominal
`SV_RATE_SCALE=1.0`) produced 140 tagged SV records: 31 at CCF 0.2, 43 at CCF 0.5,
and 66 at CCF 1.0. Its merged BAM contained 588,995,182 reads, all mapped, with
2,837 chimeric alignments, 2,216 supplementary records, 99.99% proper pairing, and
71 GB compressed size. This confirms the subclonal SV path at primary-assembly scale;
it is still a simulation/read-level result rather than a pooled caller-recall estimate.

### Local performance and BAM baseline

The FASTQ-only benchmark (Delta job `21188169`, three replicates per row) remained consistent
with the earlier chr22 measurement after the BAM changes:

| genome | eidolon median | NEAT 4 median | eidolon RSS | NEAT 4 RSS |
|---|---:|---:|---:|---:|
| *E. coli* K-12 | 14.63 s | 94.36 s | 82.1 MB | 204.9 MB |
| *S. cerevisiae* S288C | 22.85 s | 230.11 s | 63.8 MB | 182.2 MB |
| GRCh38 chr22 | 62.01 s | 751.47 s | 251.3 MB | 1530.6 MB |

The yeast thread sweep was 22.85 s (1 thread), 11.15 s (4), 7.18 s (16), and 7.68 s (64).
This benchmark does not produce BAMs. A separate chr22 BAM timing smoke produced 7,780,196
records, including 60 chimeric records, in 3:05.39 at 230.5 MB peak RSS; `samtools quickcheck`
passed and the BAM was 843 MB. The result is a BAM-path baseline, not a before/after comparison
against an older binary.

A full-primary-assembly BAM-only run (Delta job `21262840`, GRCh38 at 30×, 128 requested
threads) completed successfully in 23:56.16. It wrote 583,088,206 records to a 62 GB BAM;
`samtools quickcheck` passed. Process-level `/usr/bin/time` reported 13.6 GiB peak RSS and
924% CPU utilization (about 9.2 cores), showing that this path does not scale to all 128
requested threads. Slurm reported a much larger 111 GiB `MaxRSS`; that accounting value is
retained as a cluster-side note, while the process-level measurement is the comparable RSS
figure. This run produced no input SVs, so it characterizes full-genome BAM throughput and
resource use rather than BND recovery.

### First multi-contig GRCh38 SV smoke

Delta job `21170801` (30×, purity 0.6, `SV_RATE_SCALE=1.0`, Manta only) produced 137 somatic
SV records: 29 DEL, 32 DUP, 15 INV, 2 INS, 52 BND, and 7 CNV. All 52 BND records parsed, and
the scorer controls passed (truth-vs-truth recall 1.000; shifted-truth decoy recall 0.000).
There were no intra-contig BND spans because this smoke's BND set was inter-contig, so BNDspan
scoring was correctly skipped.

| Manta scorer | TP | FN | FP | recall | precision |
|---|---:|---:|---:|---:|---:|
| overall scoreable | 67 | 11 | 21 | 0.859 | 0.761 |
| DEL | 27 | 2 | 3 | 0.931 | 0.900 |
| DUP | 26 | 6 | 6 | 0.813 | 0.813 |
| INV | 14 | 1 | 12 | 0.933 | 0.538 |
| BND | 40 | 12 | 52 | 0.769 | 0.435 |
| INS | 0 | 2 | 0 | 0.000 | — |

These are PASS-only caller measurements from one replicate, not simulator correctness verdicts.
Most importantly, the independent INS probe found **2/2 planted insertions in the reads** even
though Manta emitted no PASS INS calls. Five de novo insertions longer than the 151 bp read
length were dropped as expected under #516. The run retained 30 GB normal and 70 GB tumor BAMs
with `PRUNE_BAM=0` for follow-up inspection.

Follow-up characterization is complete. The initial smoke is retained above as a historical
single-run baseline; the follow-up observations below provide the caller comparison and the
nominal-rate Manta repeatability check. These are still separate from the larger pooled campaign
reported in the whole-genome tier below.

### Follow-up multi-contig GRCh38 observations

Delta job `21226562` reran the same 30×/0.6/`SV_RATE_SCALE=1.0` design with both Manta and Delly
on one simulated tumor/normal pair. The truth again contained 137 somatic records (including
52 BNDs), all 52 BNDs parsed with zero unpaired or mispaired junctions, and all scorer controls
passed. Manta recovered 42/52 BNDs (recall 0.808, precision 0.457); Delly recovered 20/52
(recall 0.385, precision 1.000). Overall scoreable recall/precision were 0.795/0.747 for Manta
and 0.821/0.660 for Delly. The independent insertion probe found both planted insertions in
the reads (2/2); the zero `manta_INS` score is therefore a PASS-filter limitation, not evidence
that the simulator omitted the insertions.

Delta job `21249416` supplied the third nominal-rate Manta observation. It reproduced the same
truth counts and clean BND geometry (52/52 parsed; 0 unpaired or mispaired), with scorer
self-tests and decoys passing. Manta scored 67 TP, 11 FN, and 20 FP overall (recall 0.859,
precision 0.770), and 42 TP, 10 FN, and 50 FP for BND (recall 0.808, precision 0.457). Across
the three Manta observations, BND recall was approximately 0.77–0.81 while aggregate recall
remained about 0.86. The consistent truth-side checks and the caller-to-caller difference in
BND recovery support a caller/representation limitation rather than a simulator defect.

---

## Input VCF → output

| variant | truth VCF | read-level evidence | status |
|---|---|---|---|
| `<DEL>` symbolic | preserved | het **0.54** vs 0.50 expected; hom **0.00** exact | ⚠ [#499](https://github.com/ncsa/eidolon/issues/499) ~8% high |
| `<DUP>` symbolic | preserved | het **1.61** vs 1.50; hom **2.17** vs 2.00 | ⚠ [#499](https://github.com/ncsa/eidolon/issues/499) ~8% high |
| `<CNV>` symbolic | preserved | CN 0/2 exact; CN 1,3,4,6 all **~8% high** | ⚠ [#499](https://github.com/ncsa/eidolon/issues/499) |
| `<INV>` symbolic | preserved | deep interior **1.00**; junction-proximal **0.74–0.82** | ✅ (see note) |
| `<INS>` symbolic | preserved | **identical to no-variant control** | ❌ [#500](https://github.com/ncsa/eidolon/issues/500) silent no-op |
| BND, inter-chromosomal | both records preserved | 30 chimeric reads | ✅ |
| BND, intra-chromosomal | both records preserved | 30 chimeric reads | ✅ |
| BND + inserted sequence | preserved **with the insert** | **0 of 30** chimeric reads carry it | ❌ [#498](https://github.com/ncsa/eidolon/issues/498) |
| BND, mate contig absent | — | — | ✅ hard error (correct) |
| BND, single/unpaired (`A.`) | preserved | 0 chimeric; **depth 42.4 vs 72.0** | ❌ [#500](https://github.com/ncsa/eidolon/issues/500) destroys coverage |
| Fragment spanning a whole DUP | — | 60–800 bp blocks all match the derived haplotype | ✅ measured; earlier ❌ was noise |
| Two junctions within one fragment | all 4 preserved | 60 chimeric; depth matches derived haplotype | ✅ |
| Literal INS, ≤ ~150 bp | preserved | inserted bases present in reads | ✅ |
| Literal INS, ≳ 200 bp | preserved **with full SVLEN** | interior and tail both present (600 bp: 25 head / 19 middle / 16 tail) | ✅ [#516](https://github.com/ncsa/eidolon/issues/516) fixed |
| Literal DEL (`TTTT…`→`A`) | preserved | **D=49** vs 20 baseline | ✅ |

### Untested cells: none remain

Every row above is measured. Read-level evidence is depth against a **separate no-variant
control run**, or BAM CIGAR `I`/`D` counts against that control's sequencing-error
background (15 `I`, 20 `D` in the probed window), or chimeric-read counts.

### The confirmed failures

**[#498] BND with inserted sequence.** VCF 4.2 allows a breakend ALT to carry novel bases
inserted at the junction — common in real rearrangements, since NHEJ frequently leaves
sequence at the breakpoint. `parse_bnd_alt` accepts it, but `get_bnd_pieces` rebuilds the
read from reference pieces only, so the insert never reaches the output. The truth VCF
keeps it. **A benchmark built from that data asserts an insertion the reads do not
contain.**

**[#500] Symbolic `<INS>` is a silent no-op.** `SVTYPE=INS;SVLEN=60` is preserved in the
truth VCF, and the reads are behaviourally identical to the no-variant control — same `I`
count, same `D` count, same depth. A 60 bp insertion was requested and nothing inserted.
The de novo path *does* synthesise novel sequence for INS, so the capability exists; it is
simply not reached from `input_vcf`. That last sentence was an assertion when first written
and is now **measured**: jobs 20885875/20892075 planted 4 insertions of 296–2500 bp whose
novel bases are present in the truth VCF, at a rate consistent with the model's
`Ins = 0.0279` (4 observed against ~5.2 expected-visible). It was nearly retracted as false
in 2026-08 on the strength of `truth INS: 0` in the SV harness — which turned out to be
measuring the harness's own selector, not the generator. See §"De novo generation" below and
`docs/claude_engineering_audit.md` §5.3.

**[#516] Large insertions were only PARTIALLY realized — FIXED.** Kept here because the
defect's shape is the most instructive one in this repo: the head was always present while the
interior collapsed to zero, so any probe near the insertion's start said it worked, which is
where a casual check looks. The truth VCF meanwhile declared the full `SVLEN`.

Measured on the H1N1 fixture at `read_len=100`, `fragment_mean=250`, `GT 0/1`, counting reads
containing a 30-mer from the head / middle / far end of the inserted sequence. Identical
fixture and config before and after, so the columns are directly comparable:

| insert | head (was → now) | middle (was → now) | tail (was → now) |
|---|---|---|---|
| 50 bp | 16 → 22 | 13 → 15 | 8 → 15 |
| 150 bp | 21 → 19 | 5 → **20** | **0 → 22** |
| 200 bp | 25 → 20 | **0 → 15** | **0 → 15** |
| 600 bp | 27 → 25 | **0 → 19** | **0 → 16** |

The fix samples fragments in *altered-haplotype* coordinates, where an insertion has real
width, so fragments begin inside it. `eidolon/tests/long_insertion_rework.rs` pins eight
criteria on eidolon's own FASTQ, and `eidolon/tests/gate2_realigned_ins.rs` pins the same
insertion through **bwa-mem2**: a homozygous 600 bp insertion at 60x yields 5/5 interior
probes against 0/5 in a same-seed no-variant control, 96 breakpoint clips against 0, and 51
reads unmapped-with-mapped-mate against 0.

At genome scale, de novo on a 5 Mb chr22 window (default germline model, `sv_rate_scale=8`):
**427 of 440** insertions ≥200 bp have their middle sequence in the reads, with zero
fallbacks. The remaining 13 are thinly covered rather than absent — mean interior depth 3.96
against 6.77 for the rest, still 68% of their length covered. Interior coverage of a
heterozygous insertion runs at ~half the flanking reference (7.0 vs 13.0 reads per 30-mer),
which is the correct ratio, not a shortfall.

**What is NOT yet established:** Delta job 21378484 (chr22, `SV_RATE_SCALE=30`, COSMIC
pancancer somatic) verified 4 of 4 planted insertions at 5/5 interior probes including one of
2127 bp — real evidence, but a denominator of four. `verify_planted_ins` now refuses to let a
run that small read as validation. A larger-denominator campaign is outstanding.

**The historical mechanism**, kept because it explains why the defect was invisible:

Fragments and read windows were chosen purely in *reference* offsets
(`cover_dataset`, `generate_fragments.rs:323`), and an insertion has zero reference width, so no
read window could ever *begin* inside one. Reads are assembled per-read by walking a reference
slice and expanding variants inline, bounded by `bases_written < read_length`
(`fastq_tools.rs:429`) with a `break 'outer` at `:555` that discards the rest of the insert
silently — no log, no counter. Hence:

> novel bases visible in one read = `min(L, read_length − anchor_offset − 1)`

**At most `read_length − 1` inserted bases are ever realized, at any declared SVLEN.** Probe
visibility follows exactly: head needs `anchor_offset ≤ 69` (no `L` term), middle needs
`anchor_offset ≤ 84 − L/2` (impossible for `L ≥ 170`), tail needs `≤ 99 − L` (impossible for
`L ≥ 100`) — which reproduces all fifteen measured cells.

The head count is **size-independent**: with the heterozygous coin removed (`GT 1/1`) it is
42/41/**50/50/50** across the five sizes. An earlier version of this note claimed it "rises with
size" and read meaning into 16 → 27; that was RNG drift, since a larger insert fires
`break 'outer` sooner and consumes fewer sequencing-error draws.

Consequence, measured: SV validation campaign 20925151 planted 22 de novo somatic insertions of
61–2155 bp and Manta detected exactly one — the 61 bp event, and only at `MinSomaticScore`.
That is not primarily a caller limitation. Manta documents a *fully assembled* ceiling of
"approximately twice the read-pair fragment size" and states it also reports very large
insertions from a breakend signature alone, as `IMPRECISE <INS>` with
`LEFT_SVINSSEQ`/`RIGHT_SVINSSEQ`. Most of the planted set was inside that range, and Manta
reported nothing within 2 kb of any of them — because the evidence for the declared length was
not in the reads. **INS recall figures from campaigns before the fix are invalid above ~a
read length and must not be re-quoted; DEL/DUP/INV/BND/CNV from those runs are unaffected.**
Job 21378484 is the first post-fix counter-example: Manta's single INS call was an exact
position-and-size match to a planted 93 bp insertion, dropped only by `--passonly` as
`MinSomaticScore`, while the three larger insertions it did not call were all confirmed
present in the reads at 5/5 interior probes. Evidence present, caller declined — which is a
different finding from evidence absent, and the read-level probe is what separates them.

**The interim `read_len - 1` cap on de novo insertions has been REMOVED.** While #516 was
open, `gen-reads` refused to plant a de novo insertion whose novel sequence exceeded what a
read could carry, rather than emit a truth record the reads could not support. That truncated
the model's own distribution — the bundled default puts `Ins` at log-normal (5.7, 1.0),
median ~299 bp, so at `read_len=150` roughly three quarters of INS draws were refused.

With insertions sampled in altered-haplotype coordinates the cap is unnecessary, and its
removal is pinned by `eidolon/tests/multi_sv_integration.rs`, which now asserts the opposite
of what it once did: a de novo insertion MAY declare more novel bases than `read_len - 1`.
Reinstating the cap, or unwiring the sampler, fails that assertion while every other
assertion in the file stays green — which is the specific regression it exists to catch.
De novo insertion sizes now follow the model rather than the read length; measured largest
realized on a 5 Mb chr22 window is 9,093 bp.

Pinned by `runner.rs`'s `unrealizable_insertions_are_dropped_at_exactly_the_boundary` and
`the_insertion_cap_does_not_touch_other_variant_types` (a must-not-fire: dropping large
DELs/DUPs/INVs would be far worse than the bug being worked around), plus an end-to-end
assertion in `multi_sv_integration.rs` that no de novo INS in the truth VCF exceeds
`read_len - 1`. Removing the call site makes that fail with eleven oversized records.

This is an interim measure. It removes the false-truth problem, not the capability gap: eidolon
still cannot simulate insertions longer than a read.

**[#500] A single breakend (`A.`) destroys coverage.** VCF 4.2 allows an unresolved
partner. `parse_bnd_alt` returns no mate, and the result is worse than being ignored:
depth falls **72.0 → 42.4** (−41%) with **zero** chimeric reads. The truth declares a
breakend, the reads contain no junction, and a depth caller sees a partial deletion no
record describes.

**Fragment spanning a whole DUP — retracted, it was never broken.** This cell previously read
❌, citing a three-piece stitch (`left + block + right`) that `get_dup_pieces` cannot express,
and "4 of 30 reads match no haplotype". All three parts are withdrawn:

- **Measured false.** `short_dup_spanned_by_a_fragment_still_matches_the_derived_haplotype`
  sweeps DUP blocks of 60/100/150/200/400/800 bp against a ~200 bp fragment — every one
  spannable by a fragment — and all match the derived haplotype at the ≥90% threshold used
  throughout that file.
- **The 13% was the error model.** `matches_with_one_small_indel` exists because a purely
  positional comparison frameshifts on sequencing-error indels; its own note says that without
  it "a correct implementation looks ~9% broken". 4 of 30 = 13% sits inside that band, so the
  original measurement almost certainly predates that tolerance.
- **The issue reference was wrong.** #474 is closed and concerns INV/DUP anchor-base
  conventions, not stitching — it never covered this.

Residual, and why this is ✅ rather than settled: the test asserts the reads that exist are
correct, not that none are missing. Fragments needing three pieces being silently *dropped*
would dip coverage across a short DUP and leave the test green. A depth comparison against a
no-variant control would close that.

### `<INV>` — investigated, and the *test* was wrong, not the code

The first measurement of this cell read **0.63**, which looked like a balanced inversion
losing a third of its coverage. It was flagged rather than filed, and that was the right
call: diagnosis showed the fixture was at fault.

Against a no-variant control run of the identical config:

```
300bp inversion, 250bp fragments (H1N1_NA:300-600)
  inside inversion      0.62
  flank left            0.81
  flank right           0.76
  far away              1.00     <- coverage is fine elsewhere, so the loss is real

1200bp inversion, 250bp fragments (H1N1_PB2:500-1700)
  deep interior         1.00     <- NO loss
  near left junction    0.82
  near right junction   0.74
  outside               0.99
```

**The inverted sequence is covered at full depth.** The depletion is junction-proximal,
extending roughly one fragment length around each breakpoint. A 300 bp inversion is
*smaller than a 250 bp fragment*, so nearly every fragment overlapping it touches a
junction and the entire event sits inside the dip zone — which is what produced 0.63.

The remaining junction dip (0.74–0.82 over ~180 bp) is consistent with the mechanism: a
fragment spanning a junction becomes a chimeric read placed at one locus, so its coverage
contribution to the other side is not counted. Real aligners soft-clip and place such
reads at one side too, so some dip is expected. **How much dip is realistic has not been
compared against real data** — that is unquantified, not verified.

Practical consequence for anyone writing SV tests: **an inversion shorter than the
fragment length is not a clean fixture.** Its coverage behaviour is dominated by junction
effects, not by the inversion.

### `<CNV>` — semantics confirmed, but a systematic ~8% bias found ([#499](https://github.com/ncsa/eidolon/issues/499))

The first reading of this cell (1.67 for `CN=4`) matched no sensible interpretation. It
was measured against an in-run control span on a 300 bp event — the same flawed method
that produced the bogus `<INV>` result. Re-measured properly (1200 bp event, deep
interior only, against a separate no-variant run):

```
CN   multiplier  expected  observed  obs/exp
 0      0.00       0.00      0.00     exact
 1      0.50       0.50      0.54     1.08
 2      1.00       1.00      1.00     exact
 3      1.50       1.50      1.61     1.07
 4      2.00       2.00      2.17     1.09
 6      3.00       3.00      3.28     1.09
```

**Semantics answered:** `CN` is the **total** copy number, multiplier `CN/ploidy`, and
genotype is ignored when `CN` is present. That is consistent and defensible.

**But the delivery is ~8% high** whenever the multiplier is neither 0 nor 1 — and it is
not CNV-specific. The same bias appears in DEL and DUP, both genotypes:

```
DEL het  0.54 vs 0.50   (1.079)      DUP het  1.61 vs 1.50   (1.073)
DEL hom  0.00 vs 0.00   exact        DUP hom  2.17 vs 2.00   (1.084)
```

`coverage_multiplier_for` returns the correct values; the gap is between the requested
multiplier and the depth actually delivered. Mechanism not diagnosed — see #499.

> **Update 2026-08-23 — mechanism diagnosed, and it is size-dependent.** #499 is closed.
> Two effects contribute, both boundary effects, both shrinking as the event grows
> relative to the fragment length: chimeric junction reads landing outside the
> coverage-multiplied budget (#499's own, excess ~`1/event_length`), and fragments
> legitimately extending across the event's coverage-multiplier boundary into its
> neighbours (introduced with the fragment-placement rewrite, and the thing that removes
> an artificial dead zone at every such boundary). See the size guidance below.

### How big does an SV need to be for its declared depth to be accurate?

**Rule of thumb: depth multipliers land within ~1% for events ≳10× the fragment length.
Below ~5×, expect several percent of under-delivery. On a small contig it is worse,
because the flanking baseline is narrow too.**

Measured, het `<DUP>` (declared 1.50×), `read_len: 100`, `fragment_mean: 250`,
3 seeds averaged except where noted:

| event | ÷ fragment_mean | flanks | delivered / declared |
|---|---|---|---|
| 1,200 bp (H1N1_PB2, 2,280 bp contig) | 4.8× | ~400 bp | **88%** |
| 1,200 bp (real chr22 window) | 4.8× | megabases | **93%** |
| 1,500–4,000 bp (synthetic 2 Mb) | 6–16× | megabases | **99–101%** |
| 100,000 bp (real chr22 window) | 400× | megabases | **99.7%** |

Read this conservatively: the 1,200 bp (93%) against 1,500 bp (~100%) difference is
within seed-to-seed noise at that scale, so the useful signal is the trend, not any
single row. What the table does establish is that both ends are real — a
fragment-length-scale event genuinely loses several percent, and a realistic event
genuinely does not.

**Practical guidance:**

- **Validating depth-multiplier accuracy** (DEL/DUP/CNV dosage, copy-ratio callers): use
  events ≳10× fragment length on a reference with wide flanking sequence. The
  `depth_modulation_is_accurate_at_realistic_scale` test does exactly this and is the
  reference measurement.
- **Small references** (viral genomes, small scaffolds, the H1N1 fixture) cannot host an
  event much larger than the fragment length *and* leave wide flanks — H1N1 is 2,280 bp
  total. They remain excellent for fast mechanism checks and are used that way
  throughout this file, but they should not be used to validate depth *precision*.
- **Small SVs on any reference** carry the same few-percent under-delivery. If you are
  simulating sub-kilobase events and depend on their exact dosage, measure it rather
  than assuming the declared multiplier.

---

## Whole-genome tier — de novo, scored by callers

Everything above is **H1N1**: 8 contigs of 1–2.3 kb, largest event ~1.2 kb. That fixture cannot
say whether an SV behaves at the sizes anyone simulates. This section is the other tier —
GRCh38 at 30×, purity 0.6, `SV_RATE_SCALE=1.0`, scored by Manta 1.6.0 and Delly against truvari.

**It answers a different question.** The tables above are *input fidelity* — what you supplied
came out. This is *de novo generation as an independent tool can recover it*. Neither
substitutes for the other, and the gap between them is the subject of
`docs/sv_polish_roadmap.md`.

> ⚠ **Four replicates, pooled — against a target of ~8.** `aggregate_sv_reps.sh` pools summed
> TP/FN/FP (not a mean of per-replicate ratios, which would weight a replicate with 3 events the
> same as one with 40). At `SV_RATE_SCALE=1.0` GRCh38 yields ~25 translocations per run, so ~8
> replicates is the design target. Per-replicate truth counts vary widely — BND 50–76, INV 9–20,
> DUP 18–34 — so single-replicate figures were never going to be stable.
>
> ⚠ **Every recall below is PASS-only.** truvari is run with `--passonly`, so a truth event whose
> only matching call was non-PASS is scored FN. Measured instance: `chr10:45850637` was called by
> Manta with *identical* POS/END/SVLEN and scored a DUP false negative because its FILTER was
> `MinSomaticScore`. These recalls are therefore **underestimates by an unmeasured amount** —
> see [#541](https://github.com/ncsa/eidolon/issues/541).

**Array 21072620** (2026-08-13, seeds 9–12, GRCh38 30x, purity 0.6, `SV_RATE_SCALE=1.0`),
pooled over 4 replicates:

| scorer | reps | TP | FN | FP | recall | precision |
|---|---|---|---|---|---|---|
| `manta_overall` | 4 | 255 | 30 | 84 | 0.895 | 0.752 |
| `manta_DEL` | 4 | 105 | 11 | 13 | 0.905 | 0.890 |
| `manta_DUP` | 4 | 91 | 14 | 15 | 0.867 | 0.858 |
| `manta_INV` | 4 | 59 | 3 | 56 | 0.952 | 0.513 |
| `manta_BND` | 4 | 214 | 28 | 230 | 0.884 | 0.482 |
| `manta_INS` | **2** | 0 | 2 | 0 | 0.000 | — |
| `delly_overall` | 4 | 253 | 32 | 114 | 0.888 | 0.689 |
| `delly_DEL` | 4 | 104 | 12 | 44 | 0.897 | 0.703 |
| `delly_DUP` | 4 | 91 | 14 | 17 | 0.867 | 0.843 |
| `delly_INV` | 4 | 58 | 4 | 52 | 0.935 | 0.527 |
| `delly_BND` | 4 | 108 | 134 | 0 | 0.446 | **1.000** |
| `delly_INS` | **1** | 0 | 1 | 1 | 0.000 | — |

The INS rows ran in 1 and 2 of 4 replicates respectively and are **not comparable** with the
rest; the aggregator says so itself rather than pooling them silently.

### DUP: both callers miss the identical set, and one of them is ours

`manta_DUP` and `delly_DUP` agree to the record — TP=91, FN=14 — differing only in FP. Two
independent callers missing the same events points at a property of the events, not the callers.
Examined for one replicate (FN = 4):

| missed DUP | size | N content | called by Manta? |
|---|---|---|---|
| chr14:66171649 | 16 kb | 0.0% | nothing |
| chr11:39805201 | 43 kb | 0.0% | nothing |
| chr10:45850637 | 122 kb | 0.0% | **yes, exactly — filtered by `--passonly`** |
| chr15:19940438 | 1.7 Mb | 5.2% | nothing |

So it is **not** a size floor (1.7 Mb is trivially detectable by depth) and **not** assembly
gaps (three are 0.0% N; SV placement already requires a ±200 bp N-free window at both
breakpoints, `sv_model.rs` `ALIGNABLE_WINDOW`, #224). One of the four is a filter artifact
(#541). The remaining three — including a megabase-scale duplication with nothing called
anywhere near it — are unexplained, and diagnosing them needs reads, so `PRUNE_BAM=0` on a
future run.

Calibration controls passed: truth-vs-truth recall 1.000 with **all** 71 scoreable records
scored, decoy (shifted truth) recall 0.000; same for BND across all 50. So the matching
configuration is not loose enough to accept anything.

### BND geometry is fixed, confirmed at genome scale

This is the first whole-genome confirmation of the v3.1.0 fix. `BND recall=0.000` was the
symptom of a truth VCF describing a rearrangement the reads did not carry; it now reads:

```
BND geometry: 50 record(s), 50 parsed into a form, 0 unparsable
  t[p[  direct/deletion-like    9      t]p]  head-to-head            18
  [p[t  tail-to-tail           14      ]p]t  direct/duplication-like  9
  reciprocity: 0 unpaired, 0 mispaired
```

All four bracket forms present, every record parsed, **0 unpaired and 0 mispaired** — the
#451-era failure (truth emitted unmatchable by construction) is gone, and Manta independently
recovers 46 of the 50 records from the reads alone.

### INV precision 0.500 is ours, not the callers' — and Gate 2 now shows the reads are clean

Both callers report TP=9, FP=9 — *identically*. The pipeline predicts this in its own pre-flight:
16 of the truth's junctions are inversion-oriented (`t]p]` or `[p[t`), and a caller's breakends
for those convert to `<INV>` via Manta's `convertInversion.py`, landing as INV false positives.
Two independent callers arriving at the same number is the evidence that it is a representation
artifact. **INV recall 1.000 is real; INV precision from this tier is not quotable.**

Gate 2 (`eidolon/tests/gate2_realigned_inv.rs`) closes the remaining inference. Realigning a
homozygous 1.2 kb inversion with bwa-mem2 gives **123 breakpoint-clipped reads against 0 in the
control**, **98 same-orientation (FF/RR) pairs against 0**, **0 everted (RF) pairs** — so nothing
suggests a duplication — and interior depth of 1.026x the flank, i.e. balanced. The reads carry
an unambiguous inversion signature, so the false positives are not attributable to them. That was
previously an argument from two callers agreeing; it is now a measurement.

### INS: one planted, one verified present, zero found

`truth INS: 1` (65 bp) was the #516 cap working as designed, and the drop log proved it rather
than implying it: 2 insertions dropped (chr2, chr19) + 1 kept = 3 draws, against ~3.5 expected
from `Ins = 0.0279` over ~130 SVs, of which ~0.9 were expected to survive the 150 bp cap.

**That cap no longer exists** (#516 fixed, 2026-08-24). Every draw is now realized, so a
run's INS count is the model's own rate rather than the rate times a survival fraction —
which is why the numbers in this section cannot be compared to a post-fix run.

`INS read support: 1 of 1` — 7 reads carry a 30-mer from the **middle** of the inserted
sequence. This is the first campaign to check that the reads contained what the truth declared
*before* the BAMs were pruned, which is exactly how #516 survived three campaigns.

> ⚠ **But the check itself is not trustworthy yet.** Job 21076622, same array, reports
> `INS read support: 0 of 1` — which turned out to measure nothing. Its probe file is well
> formed (`chr3 92124128 56bp` + a valid 30-mer), yet no per-insertion line was emitted on
> stdout or stderr: `count_probe_hits` produced no output, the reporting loop ran zero
> iterations, and the `0` is an untouched initialiser printed against an independently counted
> `n_probes`. Worse, `unsupported` stayed `0`, so the failure gate never fired and the
> replicate archived as a clean result ([#540](https://github.com/ncsa/eidolon/issues/540)).
>
> The `1 of 1` above came through the same code path. It did print its detail line, so it is
> probably a true pass, but **treated the #516 cap as unverified in production** until #540 is
> fixed and both replicates are re-checked. Moot for the cap itself, which is gone; the #540
> concern about a replicate archiving a clean result without printing its detail line stands.

Neither caller found it. Manta's only INS call was a false positive on another chromosome
(`chr11:99348945`, 58 bp, non-PASS); nothing was called within 2 kb of the planted event. Delly
emitted no INS records at all. So the miss is **verified rather than inferred** — the reads
demonstrably carried the insertion. With n=1 this does not constrain a recall rate (95% CI
≈ [0, 0.98]); why it was missed is unresolved, and separating "somatic-score threshold at ~7
supporting reads" from "the reads lack a signature the aligner can present" needs `PRUNE_BAM=0`
and a look at the CIGARs.

**A 30× GRCh38 replicate yields ~3 INS draws and ~1 after the cap.** INS cannot be validated at
caller level from single replicates at any runtime — it needs pooled replicates or an
INS-enriched model. Across the ~8 replicates the aggregator targets that is still only ~8
planted insertions, so an INS-enriched model is likely the only route to a usable denominator.

### What this tier does NOT establish

- **Input fidelity at scale.** Every cell above this section is still H1N1-only.
- **Dosage.** Detection is not dosage; the ~8% over-delivery (#499) was measured on a 1.2 kb
  event and has never been checked at these sizes.
- **Subclonal behaviour.** **Resolved for de-novo SVs in #537.** Supported DEL/DUP/INV/INS/CNV
  records now carry `EIDOLON_CCF`/`EIDOLON_VAF` and use CCF-adjusted depth and junction
  evidence; reciprocal translocation BND mates share one CCF. Delta yeast validation
  exercised CCF 0.2/0.5/1.0 and retained 199 chimeric alignments. This section's older
  nominal-rate campaign predates that fix; the full-GRCh38 post-fix run above is the
  current release-candidate stress test.
- **Small events.** The planted set skews very large (DUP 3.4 Mb, DEL 9.1 Mb, INV 2.0 Mb).
  Multi-megabase events are easy to detect by depth; these recall figures say little about the
  1–50 kb range where callers actually struggle.

---

## De novo generation

| type | source | status |
|---|---|---|
| DEL / DUP / INV / CNV | PCAWG (cancer) or gnomAD-SV (germline) | ✅ sampled, symbolic ALT |
| INS | same | ✅ sampled, **literal ALT** — see below |
| BND — **cancer** | PCAWG `TRA` per-donor mean, `Bnd = 0.2323`, reproduced exactly | ✅ inter-chromosomal only |
| BND — **human germline** | `Bnd = 0.1943`, a heuristic never derived from data | ⛔ **not recommended** |
| BND — intra-chromosomal | — | ⛔ not generated by design |

**BND de novo means inter-chromosomal translocation only.** PCAWG's `TRA` class is 100%
inter-chromosomal (68,547/68,547, `docs/pcawg_sv_measurement.md` M1), and a same-contig
junction is a DEL/DUP/INV by its orientation — already sampled at their own rates.

**De novo INS is a LITERAL record, and it must carry `SVTYPE=INS`.** Unlike every other
sampled type, an insertion is emitted as `REF=<anchor>`, `ALT=<anchor><novel bases>` rather
than as symbolic `<INS>` — that routes the novel sequence through gen-reads' literal-insertion
machinery, so reads spanning the locus actually carry the bases, and it matches what Manta
emits for a *resolved* insertion. `END` equals `POS` (an insertion consumes no reference
span) and `SVLEN` is positive and equal to `ALT.len() - REF.len()`.

The tag is not cosmetic. Without it the record is indistinguishable from a small-indel-model
insertion of the same size, so **nothing downstream can attribute it** — a length threshold
cannot separate them, because the `Ins` length fit puts ~17% of structural insertions below
100 bp while the indel tail reaches ~93 bp. Measured consequence: SV validation jobs 20885875
and 20892075 planted 4 insertions between them (296, 1057, 1492, 2500 bp) and scored none,
both reporting `truth INS: 0`, because `sv_pipeline.sbatch` selects SVs on
`ALT[0]~"<" || INFO/SVTYPE!="."`. Every SV validation run before this had the same hole.
Pinned by `sampled_ins_carry_svtype_and_svlen_so_they_are_separable_from_indels` and, across
the sampler/writer boundary, `sampled_ins_svtype_survives_the_vcf_writer`.

**De novo SVs are opt-in**: `sv_rate_scale` defaults to `0.0`. For human germline they are
not recommended at all — constitutional balanced translocations occur in roughly 1 in 500
live births, so a normal genome carries ~zero. The machinery stays for cancer (rates are
real and verified) and is expected to matter for plants.

---

## Where input and de novo deliberately differ

| | input | de novo |
|---|---|---|
| intra-chromosomal BND | accepted, reads generated | never produced |
| BND with inserted sequence | accepted but insert dropped (#498) | never produced |
| arbitrary contig pairs | whatever you supply | length-weighted ([#495](https://github.com/ncsa/eidolon/issues/495)) |

These are not inconsistencies to reconcile. Input is "reproduce what I gave you"; de novo
is "sample what the corpus says". A researcher wanting a specific rearrangement should use
`input_vcf` — that path preserves IDs, POS and ALT verbatim and generates junction reads.

---

## How this was measured

```bash
# per cell: one variant, input_vcf, sv_rate_scale 0.0, produce_bam
eidolon gen-reads -c cell.yml
samtools index o.bam
samtools depth -a -r CONTIG:START-END o.bam | awk '{s+=$3;n++} END{print s/n}'   # vs a control span
zcat o_r1.fastq.gz | paste - - - - | awk '$1 ~ /EIDOLON_chimeric/' | grep -c PROBE
```

**One methodological trap, hit while producing this table.** Do not grep a FASTQ for a
sequence probe without separating the sequence from the quality string: `paste - - - -`
puts them on one line, and Phred+33 quality characters overlap the DNA alphabet (`G` is
Q38, so a run of 12 `G`s is common in *quality*). A homopolymer probe reported 71 false
matches that way. Filter to the sequence field, and prefer a non-homopolymer probe. This
is the same class of error as the FASTA reader accepting a FASTQ because `>` is Q29.
