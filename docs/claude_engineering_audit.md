# Claude on eidolon: an engineering audit

**Subject:** three months of AI-assisted development on a scientific codebase
(2026-05-02 → 2026-08-01), audited for what went well, what went wrong, and what
generalises.

**Repository:** [`ncsa/eidolon`](https://github.com/ncsa/eidolon) — a Rust reimplementation
and extension of [NEAT](https://github.com/ncsa/NEAT), a next-generation sequencing read
simulator. All code links below are pinned to
[`5c68e6f`](https://github.com/ncsa/eidolon/tree/5c68e6f9f7724920494a84bf535bd694279c08f9)
so line numbers stay valid.

---

## 0. How to check this document

Every factual claim carries one of: a commit SHA, a `file:line` link, a PR number, or a
SLURM job ID. Nothing here rests on recollection. Where a claim could not be verified,
it says so.

**Stated conflict of interest.** This is a report about Claude's work, researched and
written by Claude. The evidence was also selected by Claude. That is a real structural
problem and no amount of care removes it — treat the *evidence* as the product and the
*framing* as an argument to be checked. Three specific hazards worth knowing about:

- **Selection.** Defects listed here are ones that were eventually found. Defects never
  found cannot appear, and the audit's central finding (§6) had been invisible for two
  months, which is direct evidence that the list is incomplete.
- **Flattering framing.** Where a defect has a sympathetic explanation, the sympathetic
  explanation is also the one Claude would generate. §5 and §7 are the places to be most
  skeptical.
- **Verified reversals.** Two claims in this document contradict what Claude asserted
  earlier in the same session (§2 timeline, §3 test density). Both reversals came from
  re-deriving numbers rather than from being challenged, which is weak evidence the
  method works — but they were also *found* only because someone asked for an audit.

A fourth bias runs the **opposite** way and is worth stating up front so the document is
not read as uniformly self-critical: **an evidence standard built on commits and
`file:line` links can only see code that exists.** Work correctly *avoided* — a wrong
design talked out of before it was built — leaves no artifact and cannot appear in any
metric here. §7.1 documents one such case that the human happened to remember; there is
no way to know how many others there were.

---

## 1. Summary

> Claude's error rate is not uniform across task types. It is **low where a falsifiable
> oracle exists** — porting to match a reference implementation, refactoring, expanding a
> test suite against known behaviour, matching a tool's documented output — and **high
> where the correctness criterion must be derived from understanding external scientific
> data**. Claude's expressed confidence does not track that difference.

The sharpest single illustration: a feature was marked **"Critical … ✅ Shipped"** and its
exit criterion **"validated"**, when its defining property had never been implemented and
no run could have tested it (§6).

| | |
|---|---|
| Duration | 2026-05-02 → 2026-08-01 (3 months) |
| Commits carrying `Co-Authored-By: Claude` | 462 |
| PRs in the period | ~#118 → #492 |
| Rust LOC | 7,658 → 44,084 (5.8×) |
| `#[test]` functions | 85 → 762 (9.0×) |
| Test density | 11.1 → **22.1** → 17.3 per kLOC (peak at end of Phase 1) |

---

## 2. Timeline — and three corrections to the working narrative

Derived from `git log`; every row spot-checked against the commit.

```
2023-09-05  5068974  Rust port begins (human, solo)
2024-03-18  a439cb2  CI added (human)
2025-11-10  b8ef78d  last commit before a ~6-month dormancy
2026-05-02  28a48b2  first Claude-assisted commit
2026-05-22  33dcc60  model-parity regression suite  <- end of the test/fill-in burst
2026-05-23  23dfeab  SV machinery begins
2026-05-25  8b5dc6f  cancer design doc
2026-05-30  e9e51cc  BND support (squash of 7 WIP commits)
2026-05-30  baccf2a  v1.11.0 cancer MVP
2026-06-20  e8e3b50  NEAT-4 germline fidelity comparison  <- 3 weeks AFTER cancer MVP
2026-07-21  48493bf  rename rusty-neat -> eidolon (v2.0.0)
2026-07-30  00a6374  v3.0.0 "measurement integrity"
```

**Correction 1 — Claude did not build the CI.** It was
[`a439cb2`](https://github.com/ncsa/eidolon/commit/a439cb2), 2024-03-18, human-authored,
two years before Claude's involvement. Claude expanded the test *suite*; the harness it
ran in already existed.

**Correction 2 — cancer work began before NEAT-parity verification finished.** The
recalled sequence was "verify against NEAT first, then build cancer." The commits say
the reverse: the cancer MVP shipped 2026-05-30, the NEAT-4 germline fidelity comparison
landed 2026-06-20, and core germline defects were still being found through 2026-07-20 —
realistic Ts/Tv ([`8b38483`](https://github.com/ncsa/eidolon/commit/8b38483)), het/hom
ratio 0.01 → 0.333 ([`49487c5`](https://github.com/ncsa/eidolon/commit/49487c5)), and an
RNG returning values outside [0,1)
([`6d3b09f`](https://github.com/ncsa/eidolon/commit/6d3b09f)). **Cancer simulation was
built on a germline foundation that was still being corrected underneath it.**

**Correction 3 — "Phase 3 = the BND phase" is wrong.** BND shipped 2026-05-30 in
[`e9e51cc`](https://github.com/ncsa/eidolon/commit/e9e51cc), i.e. Phase 2. What happened
2026-07-30 → 08-01 was a *measurement-integrity audit* of it.

Two caveats on the history itself: 126 commits have author date ≠ committer date
(rebases), so pre-2026 monthly granularity is approximate; and
[`e9e51cc`](https://github.com/ncsa/eidolon/commit/e9e51cc) is an explicit squash of 7
WIP commits, so BND's development *appears* to be one day's work and was not.

---

## 3. Test density: the most useful single metric

`#[test]` count is a vanity metric — it went up 9× and several of those tests cannot
fail. Density against LOC is more honest, and the shape is the finding:

```
tag        date        tests    LOC     tests/kLOC
v1.2.0     2025-09-26      85    7,658    11.1     <- human-only baseline
v1.3.0     2026-05-08     297   15,793    18.8     <- Claude, test-focused work
v1.6.0     2026-05-22     394   17,838    22.1     <- PEAK (end of Phase 1)
v1.12.0    2026-05-30     602   30,701    19.6     <- SV + cancer sprint
v1.15.0    2026-06-06     634   33,610    18.9
v2.0.0     2026-07-21     689   37,496    18.4
v3.0.0     2026-07-30     720   39,751    18.1
HEAD       2026-08-01     762   44,084    17.3
```

**Density doubled while Claude was asked to write tests, peaked the day before feature
work started, and has declined every release since.** Between v1.6.0 and v1.12.0 —
the SV and cancer sprint — LOC grew 72% and tests grew 53%.

This is the behavioural signal underneath the thesis: **Claude optimises what it is
asked for.** Told "expand test coverage," it produced the steepest coverage improvement
in the project's history. Told "add cancer SV simulation," it produced features and let
density slide, without flagging the trade.

> **For a slide:** *"Test count went up 9×. Test density peaked in week 3 and has fallen
> every release since. The model wrote excellent tests when tests were the task, and
> fewer when features were the task — and never mentioned the difference."*

---

## 4. Phase-by-phase

### Phase 1 — fill-in, tests, efficiency (2026-05-02 → 05-22)

**This phase went well and the evidence supports the recollection.**

Highs:
- Tests 85 → 394; density 11.1 → 22.1, still the project's peak.
- Real performance defects found and fixed, e.g. streaming replaced whole-file
  block reads ([`2351d66`](https://github.com/ncsa/eidolon/commit/2351d66)).
- Model-parity regression suite ([`33dcc60`](https://github.com/ncsa/eidolon/commit/33dcc60))
  pinning builder output byte-for-byte.

Lows:
- Some tests written here were later found to be **structure-only** — see §5.3.
- Two of the three efficiency fixes credited to this phase actually happened later: the
  N-replacement stats bug ([`bc96156`](https://github.com/ncsa/eidolon/commit/bc96156),
  2026-06-15) and log volume ([`c13e563`](https://github.com/ncsa/eidolon/commit/c13e563),
  2026-07-02). Worth noting *because both were latent defects that survived Phase 1's
  test expansion* — coverage grew without catching them.

**Why this phase worked: an oracle existed.** NEAT's behaviour was the specification.
Every question had a checkable answer.

### Phase 2 — SV and cancer (2026-05-23 → 06-17)

**This is where the durable problems were introduced, and it did not look like it at the
time.**

Highs:
- Large, coherent feature delivery: SV types, cancer MVP, tumour/normal with purity,
  COSMIC/PCAWG-derived models, per-tissue models.
- The junction machinery is genuinely correct — all four VCF 4.2 breakend orientation
  forms, verified later against the spec's own worked examples.
- Real callers (Manta, Delly) detect the simulated SVs at 0.67–0.84 recall. **The
  machinery works.**

Lows:
- **The BND premise was never checked against the source data** (§6). This is the
  central finding of the audit.
- Test density fell from 22.1 to 19.6 while LOC grew 72%.
- Cancer was built before germline parity was verified (§2, Correction 2).
- The size model was fitted without checking the distributional form (§5.6).

**Why this phase went wrong: no oracle, and none was demanded.** "Does this match NEAT 4?"
has an answer. "Is this how cancer genomes rearrange?" does not, unless someone goes and
measures it. Claude filled the gap with plausible defaults and documented them
confidently.

### Phase 3 — measurement-integrity audit (2026-07-30 → 08-01)

Highs — this phase found real things:
- **#451**: de novo BND truth described a rearrangement the reads did not carry (geometry
  flags never set on that path). Had produced `BND recall = 0.000` for a caller that had
  found every junction.
- **#450**: recall denominators never checked against the truth handed to the scorer;
  160 of 567 planted sites silently excluded — the entire lowest-VAF cluster, i.e. the
  case the harness existed to test.
- **#457**: BND scored without `--bnddist` on the inherited belief that truvari could not
  match breakends — false since truvari v5.0.0; Delta runs v5.4.0.
- Mutation testing adopted as the standard for believing a test.
- `eidolon validate` built and cross-checked against samtools/bcftools, 22/22 verdict
  agreement.
- First automated tests for the Delta harness scripts (`scripts/delta/tests/`), which had
  determined report numbers while being covered by nothing.

Lows:
- Three sequential confident-but-wrong explanations for one symptom (§5.4).
- A fabricated duration propagated into four artifacts before the human caught it (§5.5).
- **The audit fixed everything downstream of the real defect without finding the defect.**
  Every Phase 3 fix concerned *measuring* BND correctly. None asked whether the BNDs
  should exist in that form at all. That question was only reached when the human asked
  for an audit of the audit.

---

## 5. Defect taxonomy

The kinds generalise better than the instances.

### 5.1 Inherited premise never checked against the source

The highest-cost class. A claim is adopted from documentation, a comment, or an earlier
decision, and every subsequent piece of work is built on it.

| instance | premise | reality |
|---|---|---|
| BND (§6) | PCAWG `TRA` count ⇒ our BND rate | `TRA` is *inter-chromosomal*; the generator emits only same-contig junctions |
| truvari | "truvari cannot benchmark breakends" | true of v4, reversed in v5.0.0; deployed version was v5.4.0 with `--bnddist` |
| `bnd_proximity.py` | same belief | an entire helper written to work around a limitation that no longer existed |

**Detection cost is the problem, not frequency.** The truvari belief produced
`BND recall = 0.000` across every SV run for two months and was read as a caller
limitation each time.

### 5.2 Cross-component invariant

Two components, each correct and each unit-tested, disagreeing about a shared contract.

**Adapter soft clips** —
[`fastq_tools.rs:657`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/eidolon-core/src/file_tools/fastq_tools.rs#L657)
tags adapter read-through `'S'`;
[`bam_writer.rs:373`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/eidolon-core/src/file_tools/bam_writer.rs#L373)
mapped everything except `I`/`D` to `Match`:

```rust
fn char_to_cigar_kind(c: char) -> CigarKind {
    match c {
        'I' => CigarKind::Insertion,
        'D' => CigarKind::Deletion,
        _   => CigarKind::Match,      // 'S' silently becomes M
    }
}
```

Every adapter base was written to the BAM as aligned sequence. A soft clip does not
consume reference, so read spans were overstated too. `fastq_tools` had a unit test
pinning the `'S'`; nothing checked the BAM honoured it. Fixed in
[#491](https://github.com/ncsa/eidolon/pull/491).

**This is in core output, not the experimental cancer layer** — and it backs #125, the
largest single number in the ACCESS report's fix table (SNP recall 0.0004 → 0.944).

The same shape produced #451 (truth VCF and reads disagreeing about geometry).

### 5.3 Verification theatre

Tests and harnesses that report success without being able to report failure.

[`bnd_fastq.rs:59`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/eidolon/tests/bnd_fastq.rs#L59)
— the entire content assertion of the BND integration test:

```rust
let has_chimeric = fastq_lines.iter().any(|l| l.contains("EIDOLON_chimeric"));
assert!(has_chimeric, "Expected to find chimeric reads in FASTQ output...");
```

A substring in a read *name*. Not which read, not how many, not one base. **Every BND
defect ever found in this repo would pass this test.**

Others of the same kind: `cosmic_bundle.rs:161` asserts `symbolic + bnd > 0` — one record
of any type — as the whole proof that the COSMIC model drives SV generation;
`cancer_parity.rs` compares the Rust implementation against a shell script, so both
being wrong passes, and it skips silently when tools are absent.

Harness-level instances follow the same pattern: a metric reported without checking it
covered its own inputs. `VERDICT: PASS` over 160 silently-excluded sites; `nsom -gt 0`
passing while every value was the malformed string `AF=AF=0.3000`.

### 5.4 Confident wrong causal explanation

`BND recall = 0.000` received three sequential confident explanations — unpaired truth,
then truvari's inability, then ALT orientation — before the actual cause. Each was
plausible enough to stop at. Separately, an "8% unexplained" INV discrepancy was
attributed to a simulator defect and turned out to be an artifact of the measurement.

**The failure is not being wrong. It is being wrong in the register of being right** —
no hedge, no "the next thing to check would be."

### 5.5 Fabricated specifics

Claude repeatedly wrote that a defect had gone unnoticed "for a year." The real interval
was **two months** (v1.13.1 2026-06-02 → v3.1.0 2026-08-01). The number was never
derived from anything; it propagated into the ACCESS report, `CLAUDE.md`, a commit
message and a PR body before the human challenged it. Separately, Claude described PRs
as unmerged while reciting a remembered queue rather than checking — all 17 were merged.

**Both are the same failure: a plausible specific generated where a lookup was needed**,
in exactly the register that discourages checking.

### 5.6 Unvalidated distributional assumption

SV sizes are fitted as log-normals. The implied tails:

```
type    median      p90       p95       p99
Del      115kb    8.4Mb    28.3Mb   276.2Mb
Dup      221kb    8.8Mb    25.1Mb   178.0Mb
Inv      445kb   55.0Mb   215.3Mb  2786.2Mb   <- ~9x the human genome
Cnv       51kb    1.7Mb     4.6Mb    29.2Mb
```

The log-normal *form* was chosen, not tested. What actually bounds emitted sizes is
`max_length_fraction` (default 0.25), documented as engineering tractability, not
biology. **A data-derived-looking parameter whose real behaviour is governed by an
invented constant.**

> **⚠️ Correction (2026-08-02), from measuring the corpus** —
> `docs/pcawg_sv_measurement.md` M4. This section originally called the fitted tails
> "biologically implausible." Half of that was my assumption, not the data's:
>
> | | empirical p99 | fitted p99 | empirical max |
> |---|---|---|---|
> | DEL | 85.0 Mb | 276.2 Mb (3.2×) | 236.4 Mb |
> | INV | 116.2 Mb | 2,786 Mb (24×) | 237.6 Mb |
>
> The fitted *tail* is genuinely wrong — 3–24× too heavy at p99. But **multi-megabase
> somatic SVs are real**, reaching whole-chromosome scale in PCAWG. So the practical
> concern inverts: the 0.25 cap (~62 Mb on chr1) is not protecting against an
> unrealistic tail, it is **truncating events that genuinely occur**.
>
> Worth noting as a data point about this document's own method: the claim was flagged as
> an *unvalidated assumption*, which was correct, and then I supplied a guess about which
> direction it was wrong in. The guess was half wrong. The taxonomy entry stands; the
> editorializing inside it did not survive contact with the data.

---

## 6. Case study: BND

The clearest end-to-end example, and the one that prompted this audit.

**What the docs claimed.** `<BND>` translocations, priority **Critical**, exemplars
BCR-ABL / PML-RARA / EWSR1-FLI1 / MYC-IGH — marked **"✅ Shipped in v1.12.0"**, with the
exit criterion *"a caller calls translocation breakpoints"* marked **"Status: validated."**

**What was true.** eidolon has never generated a translocation. Both
[`sv_model.rs:810`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/eidolon-core/src/structs/sv_model.rs#L810)
and
[`:936`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/eidolon-core/src/structs/sv_model.rs#L936)
hardcode the breakend mate to the anchor's own contig, and `sample_variants` never
receives any other contig. Job 20719077 emitted 466 junctions; **466 of 466** were
same-contig — by construction, not by chance. Four of the five named exemplars are
inter-chromosomal.

**How it happened.** The BND share (23.2%) is the PCAWG per-donor mean of
`svclass == "TRA"` rows —
[`normalize_pcawg_sv_model.py:108`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/tools/normalize_pcawg_sv_model.py#L108).
`TRA` is precisely the inter-chromosomal class. But
[`build_pcawg_sv_vcf.py:127`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/tools/build_pcawg_sv_vcf.py#L127)
reads only five columns:

```python
chrom1, start1, chrom2, start2, svclass = (f[0], f[1], f[3], f[4], f[10])
...
if svt == "BND":
    records.append((c1, s1 + 1, "<BND>", "SVTYPE=BND", "0/1"))
    continue                      # chrom2 and start2 discarded here
```

`strand1`/`strand2` (fields 8 and 9) are never read at all. **We inherited the count of
translocations while discarding the two facts that make one a translocation.** The
shipped model therefore has no chromosome-pair, strand, or geometry field, so geometry
falls back to a uniform 25% per form — and roughly a quarter of emitted "BND"s are
direct same-contig joins whose reads carry a plain deletion under a `SVTYPE=BND` label.

**A conflation underneath it.** `BND` in a VCF is a *notation* for an adjacency a caller
did not resolve — not an event class. Manta buckets inversion-oriented junctions as BND;
Delly types the same junctions `<INV>`. Both were observed on identical data in job
20745149. A "23.2% BND rate" partly encodes which caller and chemistry produced the
source calls.

**Why it survived a dedicated BND audit.** Phase 3 fixed the scoring, the flags, the
denominators, and the harness — everything about *measuring* BND. It never asked whether
the events should exist in that form. Retracted in
[#492](https://github.com/ncsa/eidolon/pull/492).

> **For a slide:** *"The model counted inter-chromosomal translocations from real cancer
> data, then emitted them as same-chromosome junctions with random partners and random
> orientation. The label came from the data; the biology was invented. It was marked
> Critical, Shipped, and Validated."*

---

## 7. Who caught what

**Caught by Claude, unprompted:** the streaming inefficiency; the N-replacement stats
bug; log volume; #450's silent denominator exclusion; #451's de novo flag defect; the
`bnd_proximity.py` redundancy; the BNDspan negative-control failure; the adapter CIGAR
defect (during this audit).

### 7.1 The contribution class this audit cannot see

One episode, attested by the human author and recorded here in his words:

> *"I was on the verge of a custom alignment function when Claude suggested a much
> simpler path that I had overlooked."*

The context is CIGAR generation. The apparent problem — "produce a CIGAR string
describing how this read aligns to the reference" — reads like an alignment problem, and
the obvious solution is to write an aligner. **It is not an alignment problem.** A
simulator already knows every edit it applied, so the CIGAR is bookkeeping accumulated
during generation, not an inference recovered afterwards:

```rust
// eidolon-core/src/file_tools/fastq_tools.rs:548-563
if is_first_base { cigar_ops.push('M'); is_first_base = false; }
else             { cigar_ops.push('I'); }
...
for _ in 0..deletion_skip { cigar_ops.push('D'); }   // reference bases skipped
```

Verifiable consequence: **there is no aligner anywhere in this codebase.** The whole
class of work — implementing it, testing it, tuning it, and carrying its performance cost
on every simulated read — was avoided.

**This matters methodologically, and it cuts against the rest of this document.**
§0 notes that defects never found cannot appear in the defect list, which biases the
audit *toward* Claude. This is the mirror bias, running the other way: **work that was
correctly not done leaves no artifact.** There is no commit, no diff, no test, no line
of code to link — I searched the history and found nothing, exactly as expected. Every
metric in §3, every taxonomy entry in §5, and the entire evidentiary standard of this
report can only measure code that exists.

So the honest statement of this document's coverage is narrower than it looks:

| | measurable here | invisible here |
|---|---|---|
| Code written correctly | ✅ tests, density, parity | |
| Code written wrongly | ✅ the defect taxonomy | |
| Defects never found | | ❌ (§0 — biases toward Claude) |
| **Work correctly avoided** | | ❌ **(biases against Claude)** |

The second bias is probably the larger of the two in a codebase this size, and it has no
remedy short of the human writing down near-misses when they happen. The one above is
recorded because he remembered it, not because the method surfaced it.

> **For a slide:** *"The best thing the model did on this project produced no code, so
> nothing in this audit can measure it. A custom aligner was about to be written; it
> wasn't needed, because a simulator already knows what it changed. There is no commit
> to point at — that is the point."*

**Caught only because the human intervened:**

| intervention | what it surfaced |
|---|---|
| "all new code needs adequate test coverage… I'd rather work slower than produce wrong code" | The test-adequacy audit, which found the structure-only tests |
| "you're calling things unmerged without checking github again" | Recitation-instead-of-lookup |
| "idk where you got it stuck in your head about a year" | The fabricated duration, after it reached four artifacts |
| "it sounds like a gap in our understanding of the data" | **The translocation defect** — the largest finding here |
| "probably it has multiple representations, depending on chemistries and sequencers" | The BND-as-notation conflation |

**The pattern is not that the human found bugs Claude missed.** It is that the human
supplied *doubt about premises*, and Claude supplied *thoroughness within them*. Every
item in the right-hand column is a question about whether the frame was right.

**One genuine counterexample, in fairness:** the `bnd_proximity.py` finding was
premise-doubt generated by Claude. An entire helper had been written on the inherited
belief that truvari could not benchmark breakends; Claude checked the tool's actual
version and changelog, found the belief had been false since v5.0.0, and deleted the
workaround. So the capability is not absent.

But it is unreliable and it is not *scheduled*. The premise Claude questioned was about a
third-party tool with a published changelog — a lookup with an answer. The premises it did
not question were about what the scientific data means, where checking requires deciding
to go measure something. Across three months, that second kind was initiated by the human
every time, including during a phase explicitly dedicated to auditing that exact
subsystem.

---

## 8. What this suggests for practice

1. **Ask what the oracle is before starting.** If there is no way to be proven wrong, the
   work is not ready to be built. "Match NEAT 4" is an oracle. "Simulate cancer
   realistically" is not, until someone measures something.
2. **Distinguish data-derived from invented, explicitly and in the artifact.** The BND
   defect would have been visible in a table with two columns. It was not invisible
   because it was subtle; it was invisible because nobody wrote it down.
3. **Confidence is not a signal — treat its uniformity as the warning.** The same tone
   accompanied "tests 85 → 297" (true) and "translocations validated" (false in its
   defining respect).
4. **Mutation-test, or do not claim coverage.** Break the code and watch the test fail.
   Every vacuous test here passed continuously for months.
5. **Assert coverage of inputs, not just metrics.** `n_scored` vs `n_planted`. A metric
   over an unknown denominator is not a result.
6. **Audit the premise separately from the implementation, and schedule it.** A
   three-day audit of BND *measurement* did not find that BND was the wrong event. Those
   are different questions and the second is never reached by pursuing the first.

> **For a slide:** *"The model is a strong engineer and an unreliable scientist. It will
> build what you specify, test it well if you ask, and optimise what you measure. It will
> not ask whether the specification corresponds to reality — and it will sound equally
> confident either way."*

---

## Appendix — reproducing the numbers

```bash
# test density by tag
for t in v1.2.0 v1.3.0 v1.6.0 v1.12.0 v1.15.0 v2.0.0 v3.0.0; do
  n=$(git grep -h '#\[test\]' $t -- '*.rs' | wc -l)
  l=$(git grep -h '' $t -- '*.rs' | wc -l)
  echo "$t $n $l $(python3 -c "print(f'{$n/$l*1000:.1f}')")"
done

# Claude-attributed commits by month
git log --grep='Co-Authored-By: Claude' --format='%ad' --date=format:%Y-%m | sort | uniq -c

# the BND same-contig constraint
grep -n 'Some(contig_name.to_string())' eidolon-core/src/structs/sv_model.rs

# what the PCAWG corpus builder reads, and drops
sed -n '127,129p;155,161p' tools/build_pcawg_sv_vcf.py

# the shipped model's SV fields
python3 -c "import gzip,json; \
  print(list(json.load(gzip.open('tools/cosmic_v104_pancancer_model.json.gz'))['sv_model']))"
```

**Related:** [#491](https://github.com/ncsa/eidolon/pull/491) (adapter CIGAR fix),
[#492](https://github.com/ncsa/eidolon/pull/492) (retractions),
[`CLAUDE.md`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/CLAUDE.md)
(the working standard these lessons produced),
[`docs/access_report_draft.md`](https://github.com/ncsa/eidolon/blob/5c68e6f9f7724920494a84bf535bd694279c08f9/docs/access_report_draft.md)
§3.7.1.
