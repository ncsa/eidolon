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
| Literal INS, ≳ 200 bp | preserved **with full SVLEN** | **head only** — interior and far end in 0 reads | ❌ [#516](https://github.com/ncsa/eidolon/issues/516) partial |
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

**[#516] Large insertions are only PARTIALLY realized, which is worse than not at all.** An
insertion longer than roughly one read is spliced into the haplotype, but only its first
~100–150 bases are ever sequenced. Measured on the H1N1 fixture at `read_len=100`,
`fragment_mean=250`, counting reads containing a 30-mer from the head / middle / far end of the
inserted sequence:

| insert | head | middle | tail |
|---|---|---|---|
| 50 bp | 16 | 13 | 8 |
| 150 bp | 21 | 5 | 0 |
| 200 bp | 25 | **0** | **0** |
| 600 bp | 27 | **0** | **0** |

The head is always present while the interior collapses to zero — so any probe near the
insertion's start says it works, which is where a casual check looks. The truth VCF meanwhile
declares the full `SVLEN`.

**Located mechanism.** Fragments and read windows are chosen purely in *reference* offsets
(`cover_dataset`, `generate_fragments.rs:323`), and an insertion has zero reference width, so no
read window can ever *begin* inside one. Reads are assembled per-read by walking a reference
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
not in the reads. **INS recall figures are invalid above ~a read length; DEL/DUP/INV/BND/CNV
are unaffected.**

Suspected mechanism, same family as #498: fragments are placed against *reference*
coordinates and the insert is spliced in afterwards, so a fragment can reach into the start of
an insertion but none ever *begins* inside it — nothing samples the interior.

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
