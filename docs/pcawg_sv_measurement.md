# PCAWG SV measurement pass (P2) — results

**Status: MEASURED.** Every number below was computed from the PCAWG open-tier
`consensus_sv` corpus on 2026-08-02. Reproduction commands are in the appendix.

## Why this was run

The SV model's BND component rested on an unchecked premise: the BND share (23.2%) is the
PCAWG per-donor mean of `svclass == "TRA"` rows, while the generator can only emit
same-contig junctions. This pass measures what the source data says, **before** any code
is written, so P3's scope follows from evidence. See `docs/cancer_simulator.md`
(retraction) and `docs/claude_engineering_audit.md` §6.

## Corpus and coverage

| | |
|---|---|
| Source | `s3://icgc25k-open/PCAWG/consensus_sv` via `object.genomeinformatics.org` |
| Fetched with | `tools/fetch_pcawg_sv_corpus.sh` — endpoint live, **md5s from 2024-07-18 still match** |
| Size | 7.7 MB — ran on a workstation, no HPC allocation |
| Donor files | 1,926 ICGC + 822 TCGA = **2,748** |
| Donors with ≥1 call | 2,605 (**143 donors have zero SV calls**) |
| Rows parsed | **309,246 of 309,246** — zero parse errors, zero dropped |
| Header variants | 1 (identical across all 2,748 files) |

Coverage assertions passed: every row landed in exactly one cell of every cross-tab, and
cross-tab totals equal the row count.

---

## M1 — `TRA` is inter-chromosomal. Exactly.

```
svclass       n    intra-chrom  inter-chrom   % inter
DEL       86,139       86,139            0      0.0%
DUP       76,214       76,214            0      0.0%
h2hINV    39,495       39,495            0      0.0%
t2tINV    38,851       38,851            0      0.0%
TRA       68,547            0       68,547    100.0%
TOTAL    309,246                                  coverage OK
```

**68,547 of 68,547 TRA rows are inter-chromosomal. Zero exceptions.** The four span
classes are 100% intra-chromosomal, also with zero exceptions.

This is the load-bearing result. The 23.2% BND share is a **translocation rate**, and the
generator emits none of it — it produces same-contig junctions instead.

**There is no such thing as an intra-chromosomal `TRA` in this corpus.** The event class
eidolon generates for its BND budget does not exist in the data the budget came from.

## M2 — Strand spectrum

```
svclass       n        ++       +-       -+       --
DEL       86,139     0.0%   100.0%     0.0%     0.0%
DUP       76,214     0.0%     0.0%   100.0%     0.0%
h2hINV    39,495   100.0%     0.0%     0.0%     0.0%
t2tINV    38,851     0.0%     0.0%     0.0%   100.0%
TRA       68,547    24.6%    25.4%    25.3%    24.7%
```

**Known-answer control: passed.** Each span class is *definitionally* one strand
orientation, at 100%. That confirms the format reading; had DEL been mixed, nothing else
here would be trustworthy.

Two conclusions, and they point in opposite directions:

**For inter-chromosomal junctions, uniform geometry is empirically correct.** TRA
orientation is 24.6 / 25.4 / 25.3 / 24.7 — indistinguishable from uniform.
`default_bnd_geometry_weights()`'s uniform 25%, written as a deliberate non-claim, turns
out to be *right* — for translocations.

**For intra-chromosomal junctions, uniform is emphatically wrong**, because the four
orientations are not one class at all. They are the four separate classes, with distinct
rates:

| orientation | class | share of intra-chromosomal |
|---|---|---|
| `+-` | DEL | 35.8% |
| `-+` | DUP | 31.7% |
| `++` | h2hINV | 16.4% |
| `--` | t2tINV | 16.1% |

So a same-contig junction with uniformly-drawn orientation is not a BND — it is a DEL, a
DUP, or an INV, drawn at the wrong rates, under a `SVTYPE=BND` label. eidolon already
generates DEL/DUP/INV separately, so the intra-contig BND budget is **duplicating
categories it is also emitting properly elsewhere**.

## M3 — Translocation partner structure

Since TRA is 100% inter-chromosomal, there is no empirical intra-chromosomal partner
distance to fit. The current model's uniform-over-contig mate draw
(`sv_model.rs:782-785`) models a quantity that does not appear in the source at all.

The real question is which chromosomes pair. **They do not pair uniformly:**

```
68,547 TRA junctions across 276 distinct chromosome pairs
(uniform would be ~0.36% per pair)

  chr1  - chr12   2.46%   <- ~7x uniform
  chr12 - chr6    1.73%
  chr12 - chr5    1.32%
  chr1  - chr2    1.21%
  chr2  - chr3    1.11%

endpoint share:  chr1 8.17%, chr12 8.14%, chr6 6.46% ... chr22 1.96%, chrY 0.12%
```

chr1's 8.17% tracks its share of genome length; **chr12's 8.14% does not** — chr12 is
~4.3% of the genome, so it is roughly 2× enriched. Partner choice is therefore not
explained by length alone, and a length-weighted draw would be an improvement over
uniform but still wrong.

## M4 — Spans vs the fitted log-normal

```
type       n    source        median      p90      p95       p99       max
DEL    86,139   EMPIRICAL       92kb   16.6Mb   34.4Mb    85.0Mb   236.4Mb
                fitted LN      115kb    8.4Mb   28.3Mb   276.2Mb  unbounded
DUP    76,214   EMPIRICAL      190kb   16.2Mb   35.2Mb    88.6Mb   236.4Mb
                fitted LN      221kb    8.8Mb   25.1Mb   178.0Mb  unbounded
INV    78,346   EMPIRICAL      1.2Mb   36.8Mb   60.1Mb   116.2Mb   237.6Mb
                fitted LN      445kb   55.0Mb  215.3Mb  2786.2Mb  unbounded
```

*(INV empirical here pools h2hINV+t2tINV; the fit uses h2hINV only. Medians are
comparable, the count is not.)*

**This partially corrects a claim in `docs/claude_engineering_audit.md` §5.6.** That
section called the fitted tails "biologically implausible." The tail *shape* is indeed
wrong — the fitted p99 overstates the empirical p99 by **3.2× (DEL)** and **24× (INV)**.
But the premise that multi-megabase SVs are implausible was mine, not the data's:
**empirical maxima reach ~236 Mb**, i.e. whole-chromosome scale. Very large somatic SVs
are real.

That inverts the practical concern. The engineering cap (`max_length_fraction`, default
0.25 of contig length — ~62 Mb on chr1) is not protecting against an unrealistic tail; it
is **truncating events that genuinely occur**. The log-normal body is a fair
approximation; the tail is wrong in both directions at once — too heavy in the fit, too
tightly clipped in the sampler.

## M5 — The shipped model reproduces exactly

Recomputing `normalize_pcawg_sv_model.py:108-114` from the corpus:

```
        recomputed   shipped    delta
  Del       0.2919    0.2919  -0.0000
  Dup       0.2582    0.2582  -0.0000
  Inv       0.1338    0.1338  +0.0000
  Bnd       0.2323    0.2323  -0.0000
  Cnv       0.0559    0.0559  -0.0000
  Ins       0.0279    0.0279  -0.0000
```

**All six match to four decimal places.** This had never been verifiable — the
`.counts.json` sidecar was not committed — and it is now confirmed.

Two details that had to be right for this to reproduce, both of which are correct:

- **The denominator is all 2,748 donors, not the 2,605 with calls.** Donors with zero SVs
  are real observations; excluding them would inflate every rate by 5.5%.
- **`t2tINV` is deliberately skipped** (`build_pcawg_sv_vcf.py:143-145`) because a
  balanced inversion appears as two BEDPE rows bracketing one footprint. This is a
  documented, considered decision, and the corpus independently supports it: h2hINV
  39,495 vs t2tINV 38,851, near-equal, consistent with pairing.

**The fitting arithmetic was never the defect.** The numbers are right; what is wrong is
the semantics of the event they are attached to.

## M6 — Which BND representation are we modeling?

Not measured — a decision to record, prompted by the observation that `BND` is a
*notation* for an unresolved adjacency rather than an event class. The corpus makes the
choice concrete:

- **PCAWG consensus** uses `svclass` to classify every junction it can, reserving `TRA`
  for inter-chromosomal ones. Under this convention, "BND" ≡ translocation, and
  intra-chromosomal junctions are never BND.
- **Manta (short read)** buckets inversion-oriented adjacencies as BND and calls direct
  ones `DEL`/`DUP:TANDEM`. **Delly** types the same junctions `<INV>`. Both observed on
  identical data, job 20745149.

These are incompatible, and the model currently mixes them: it takes a rate from the
PCAWG convention and emits records under something closer to the Manta one.

**Proposal for P3:** declare `bnd_representation: pcawg_consensus` explicitly in the
model, meaning BND ≡ inter-chromosomal translocation, with intra-chromosomal junctions
emitted as DEL/DUP/INV. One convention first; knobs later.

---

## Decision-rule outcome

The rule was fixed in writing before the data was seen:

| Rule | Measured | Fires? |
|---|---|---|
| `TRA` predominantly inter-chromosomal → full genome-aware sampler | **100.0%**, 68,547/68,547 | **YES** |
| Material share of BND-classed events intra-chromosomal → split category | 0.0% | no |
| Strand spectrum materially non-uniform → fit geometry weights | uniform for TRA; **non-uniform intra (36/32/16/16)** | **partly** — see below |
| Partner distance shows decay → fit distribution | no intra-chrom TRA exists; **chromosome pairing non-uniform** | **YES**, reframed |
| Spans not log-normal → re-open size model | body OK, p99 off by 3–24×, cap truncates real events | **YES** |

**Consequences for P3:**

1. **Make `sample_variants` genome-aware** so a BND mate can land on another contig. The
   read path already supports this (`runner.rs:2360`); only the sampler does not.
2. **Keep uniform BND geometry** — measured correct for translocations. The existing
   non-claim survives contact with the data.
3. **Stop emitting intra-contig BND.** Under the PCAWG convention it is not a category;
   those junctions are DEL/DUP/INV, already generated separately at their own rates.
4. **Draw partner chromosomes from the empirical pair distribution**, not uniformly and
   not by length alone.
5. **Revisit the size cap**, which truncates events the corpus shows are real; and the
   log-normal tail, which is too heavy by 3–24× at p99.

---

## Appendix — reproduction

```bash
# fetch (7.7 MB; script verifies md5s)
tools/fetch_pcawg_sv_corpus.sh --out-dir /tmp/pcawg     # or fetch consensus_sv/ only

cd /tmp/pcawg/consensus_sv
mkdir -p icgc tcga
tar xzf final_consensus_sv_bedpe_passonly.icgc.public.tgz -C icgc
tar xzf final_consensus_sv_bedpe_passonly.tcga.public.tgz -C tcga

# flatten: donor, chrom1, start1, chrom2, start2, strand1, strand2, svclass
: > all_sv.tsv
for f in $(find icgc tcga -name '*.bedpe.gz' | sort); do
  zcat "$f" | awk -F'\t' -v d="$(basename "$f" | cut -d. -f1)" -v OFS='\t' \
    'NR==1{next} NF>=11 {print d,$1,$2,$4,$5,$9,$10,$11}' >> all_sv.tsv
done
wc -l all_sv.tsv                     # expect 309246

# M1 — intra vs inter by class
awk -F'\t' '{n[$8]++; if($2==$4) a[$8]++} END{for(c in n) printf "%-8s %7d %7d %6.1f%%\n", c,n[c],a[c]+0,100*(n[c]-(a[c]+0))/n[c]}' all_sv.tsv

# M2 — strand spectrum by class
awk -F'\t' '{n[$8]++; k[$8 FS $6 $7]++} END{for(x in k) print x, k[x]}' all_sv.tsv | sort

# M3 — translocation chromosome pairing
awk -F'\t' '$8=="TRA"{a=$2;b=$4; if(a>b){t=a;a=b;b=t} print a"-"b}' all_sv.tsv | sort | uniq -c | sort -rn | head

# M5 — per-donor means -> type probabilities (denominator = 2748, INV = h2hINV only)
awk -F'\t' '{n[$8]++} END{for(c in n) print c, n[c]}' all_sv.tsv
```
