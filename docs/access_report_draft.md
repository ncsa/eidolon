# Benchmarking eidolon: Enhancing Performance of a Genome Sequencing Simulator for Cancer Genomics with Rust

**ACCESS allocation report — DRAFT**
Allocation: `bhrd-delta-cpu` (NCSA Delta). Core-hours used to date: ~2,500 of ~397,500.

---

## 1. Project summary

`eidolon` (formerly `rusty-neat` / `rneat`; renamed at v2.0.0) is a Rust port and
extension of NEAT, a next-generation-sequencing read simulator. It generates FASTQ, a golden BAM, and a truth VCF whose statistical
properties match real data, and it adds a native tumor/normal cancer workflow
(`gen-cancer-reads`) with structural variants and an origin-tagged somatic truth
set. This project used Delta to (1) verify eidolon produces high-quality data on a
small test set, (2) compare eidolon to its predecessor NEAT 4 at scale on both
performance and fidelity, and (3) validate the cancer workflow end-to-end through
a real somatic-variant-calling pipeline.

A central outcome is methodological: **benchmarking eidolon against NEAT 4 and
against GATK/Mutect2 surfaced — and we fixed — six real simulator defects** that
unit tests had not caught. The verification process is itself a result — most
strikingly for the sixth defect (§3.11), a malformed-FASTQ bug that was invisible
to local read-count checks and only surfaced when a real aligner ran the data at
scale on Delta.

---

## 2. Methods

All runs on Delta CPU nodes (dual-socket AMD EPYC, 128 cores). Pipeline tools:
BWA-MEM2 (alignment), GATK HaplotypeCaller (germline calling), GATK Mutect2 +
FilterMutectCalls (somatic calling), hap.py / som.py (scoring against truth),
NEAT 4.6.1 (reference simulator). Every run archives its result artifacts and
SLURM resource usage (core-hours) to durable project storage; `collect_report.sh`
aggregates them.

Reference genomes: chr22 (GRCh38) for the verification tier; E. coli, yeast for
performance scaling. Default 30× coverage, 151 bp paired-end, ploidy 2 unless
noted.

---

## 3. Results (Phase 1 — complete)

### 3.1 Simulator data quality (chr22, 30×)

| Property | Result |
|---|---|
| Determinism (FASTQ + VCF, 1 vs 8 threads) | **byte-identical (PASS)** |
| Mean coverage (covered bases) | 29.96× (requested 30×) |
| Breadth | 0.77 (chr22 N-gaps correctly left uncovered) |
| Heterozygous allele fraction | 0.486 (clean diploid binomial ~0.5) |
| Het/hom ratio | **~2.0** (human-realistic; a default-model bug that inflated this to ~95 was found and fixed during this validation — see note) |
| Ts/Tv | **2.35** (human-realistic) |

Determinism holding byte-for-byte across thread counts is a guarantee NEAT 4 does
not make.

**Het/hom ratio (bug found and fixed).** The germline runs initially showed a
heterozygous/homozygous variant-site ratio of ~95 (truth and recovered tracking
together — e.g. 94.7 / 95.5 on the input-variety set), far above a real
resequenced human's ~1.5–2.0. Investigation traced this to the bundled default
mutation model: its `homozygous_frequency` was **0.01**, so each variant site is
made homozygous-alt with probability 0.01, forcing a het/hom ratio of ~99. That
value was a faithful port of the NEAT lineage's own default (NEAT2 used 0.010,
NEAT4 0.001) — a long-standing upstream default that is unrealistic for any
diploid, where a resequenced sample is ~1/3 homozygous-alt. The generative
sampling itself is correct: `homozygous_frequency` is fit as (homozygous
sites)/(total sites) and consumed directly as P(homozygous), with no inversion —
a properly *fit* model already yields a realistic ratio. The defect was purely
the unrealistic default parameter. It has been corrected to **0.333** (het/hom
≈ 2.0), confirmed by re-running the germline pipeline on chr22 with the
corrected default (42,997 variant sites → het/hom 2.00, P(hom) 0.333) and locked
by a regression test. Because the caller recovers whatever zygosity is
simulated (truth and query ratios matched throughout), recall, precision and
Ts/Tv are unaffected by this correction — the fix improves realism of the
absolute ratio, not the validation conclusions. A population / common-variant
(LD) structured zygosity model (#316) remains a separate, optional realism
enhancement; it is *not* required for a realistic het/hom ratio.

### 3.2 Performance vs NEAT 4 (single thread)

| Genome | eidolon wall | NEAT wall | **Speedup** | eidolon RSS | NEAT RSS | **Memory** |
|---|---|---|---|---|---|---|
| E. coli (4.5 MB) | 9.27 s | 96.85 s | **10.4×** | 40 MB | 281 MB | **7.0×** |
| yeast (11 MB) | 20.68 s | 231.1 s | **11.2×** | 53 MB | 168 MB | **3.2×** |
| chr22 (49 MB) | 60.94 s | 810.3 s | **13.3×** | 227 MB | 1525 MB | **6.7×** |

eidolon is **~10–13× faster** than NEAT 4 and uses 3–7× less, flatter memory (both
tools single-threaded, 10× coverage, n=3 median on one exclusive node).

**Revision (v1.19.1).** An earlier draft reported ~2.7–2.9×. Those eidolon timings
were taken while eidolon defaulted to `--log-level trace`, whose per-base debug I/O
dominated its runtime (issue #340); NEAT's figures were unaffected. With the
default corrected to `info`, eidolon's true single-thread advantage is **~10–13×**.
Memory was never affected (logging is CPU/I/O-bound), so the 3–7× memory advantage
is unchanged. Both tools ran at minimal logging — eidolon `warn`, NEAT its default
`info`, which has no per-read logging and is timing-equivalent to `warn` (verified
in NEAT's source; #346 makes NEAT `--log-level WARNING` explicit going forward).

**Multicore behavior (v1.19.1).** The same fix also reverses eidolon's apparent
thread *regression*. On the multi-contig yeast genome, in-process threading now
**scales** — 1 thread 20.7 s → 4 threads 9.5 s → 16 threads 5.3 s (**~3.9×**) —
whereas the pre-fix build got *slower* with more threads (77 → 116 s). The
regression was per-thread trace-log I/O contention (every worker writing the
multi-GB log), not the compute engine. This recovery is characterized fully in
**§3.6.2** — threading scales up to the genome's contig count (~6–9×) — but the
validated whole-genome path remains multi-process region-sharding (§3.6.1/§3.6.2),
which parallelizes sub-contig and packs near-linearly across a node, well beyond
what threading a single process can reach.

### 3.3 Germline fidelity vs NEAT 4 (chr22, 30×, identical pipeline)

Caller recovery of each simulator's truth set through the same GATK pipeline:

| | eidolon SNP | NEAT SNP | eidolon INDEL | NEAT INDEL |
|---|---|---|---|---|
| Recall | 0.989 | 0.989 | 0.982 | 0.969 |
| Precision | 0.9996 | 0.991 | 0.989 | 0.972 |

**Statistically equivalent**, with eidolon marginally cleaner on precision. eidolon
faithfully reproduces NEAT 4's germline behavior.

### 3.4 Cancer end-to-end (chr22, 30×, purity 0.6)

Mutect2 tumor/normal recovery of eidolon's somatic truth (after fixes + truth/query
normalization):

| | Recall | Precision |
|---|---|---|
| Somatic SNV | 0.925 | 0.872 |
| Somatic indel | 0.900 | 0.844 |

eidolon's simulated tumor/normal data is recovered well by a standard somatic
pipeline. **Scope:** SNV/indel only (Mutect2 does not call SVs), single
purity/coverage point — both expanded in Phase 2.

### 3.5 Defects found and fixed through verification

Benchmarking surfaced six real bugs invisible to unit tests; each was diagnosed,
fixed, and re-verified end-to-end. A second audit cycle in August 2026 surfaced eleven
more, of a different kind — defects in the *measuring* layer rather than the simulator —
and those are tabulated after the first six.

| # | Defect | Effect when fixed |
|---|---|---|
| #280 | Default mutation model near-random (uniform substitution matrix) | Ts/Tv 0.72 → 2.35 (NEAT parity) |
| #285 | Alt alleles written only to forward-strand reads | somatic SNV recall 0.33 → 0.94 |
| #287 | Indels garbled on reverse-strand reads | indel precision 0.07 → 0.49 (read-level) |
| #289 | Coverage regression introduced by #287 (R2 deletion truncation) | restored −14.5% lost read pairs |
| #290 | som.py indel normalization off by default (scoring) | cancer indel recall 0.60 → 0.90 |
| #125 | Adapter path wrote a malformed FASTQ record (quality one char longer than the sequence) for zero-length inserts; bwa-mem2 halts at the first such record | short-insert adapter alignment 2,986 → 3.89M reads (silent truncation eliminated); SNP recall 0.0004 → 0.944 |

**August 2026 — a second cycle, and a different class of defect.** The six above are
simulator bugs found by benchmarking. The eleven below were found by auditing the
*measuring* layer, and every one of them had been reporting a plausible number rather
than failing:

| # | Defect | What it was reporting instead |
|---|---|---|
| #451 | de novo BND truth described a rearrangement the reads did not carry (flags never set for the de novo path) | `BND recall = 0.000` for a caller that had found every junction to within 1–2 bp |
| #450 | subclonal-VAF harness excluded the sites it exists to test (`bcftools call` discards a hom-ref alt and its read count) | `VERDICT: PASS` with clean bias/MAE over 407 of 567 planted sites |
| #457 | per-type SV scoring filtered the truth but not the query | `manta_BND FP=43` against `manta_overall FP=24` — a subset with more false positives than the whole |
| #457 | BND scored without `--bnddist`, on the inherited belief truvari could not match breakends | `0.000`, across every SV run to date |
| #444 | `compare-vcfs` wrote VCFs whose own header did not declare their INFO tags | opaque `[E::bcf_hdr_read]` failures far downstream |
| — | `build_bnd_spans` was **called 250 lines above its definition** — `command not found` (127) swallowed by `\|\|` | "no intra-contig junctions?" — impossible on a single-contig reference |
| — | recall denominators never checked against the truth handed to the scorer | 4 of 152 per-type checks silently scored a subset (a DEL reported 0.800 where the full denominator gives 0.727) |
| — | a stale binary ran with no version check and no log line | a full set of plausible numbers testing none of the fix it was submitted to verify |
| — | POS meant "last unaffected base" for `<DEL>` and "first affected base" for `<DUP>`/`<CNV>`/`<INV>` | `END−POS == \|SVLEN\|` for DEL and `\|SVLEN\|−1` for the rest, in the same file |
| — | `<INV>` and `<DUP>` junction geometry had **no tests at all** | four mutations, including inverting the junction dispatch, passed the entire suite |
| #458 | only one of four breakend geometries was ever emitted | BND was not a separable category from INV — our junctions converted to `<INV>` and scored as false positives |

The pattern is uniform and worth stating plainly: **each reported a metric without
asserting it had measured what it was given.** A record count passed while every value
was the malformed string `AF=AF=0.3000`; a recall was computed over a denominator nobody
had checked; a scorer was blamed for a defect in its input. Several were caught only
because a test was written to fail — reverting the fix and watching it break — rather
than to pass.

All six of the original defects affect every paired-end run, not just cancer (the sixth
affects any short-insert adapter run). The somatic-SNV journey (0.33 → 0.94) and cancer-indel journey
(0.20 → 0.90) trace the cumulative impact. Defect #125 is detailed in §3.11: it is
the clearest case of verification-as-result — local `zcat`/`wc` read counts were
correct (line counts are unaffected by a too-long quality string), so the bug was
only exposed by pushing the data through a real aligner at scale.

### 3.6 Multicore scaling and the whole-genome strategy

> **Update (v1.19.1):** the wall-clock figures in this section were measured on
> builds that defaulted to `--log-level trace`, whose per-base debug I/O we later
> found dominated the runtime (issue #340). They are therefore **~3–4× too slow**
> and are kept here for the investigation narrative. See **§3.6.1** for the
> corrected simulation-only timing on the fixed build.

Benchmarking exposed a scaling characteristic that reshaped the whole-genome
plan. On Delta's dual-socket EPYC nodes, eidolon's in-process (rayon) threading
initially *regressed* — yeast 30× took 218 s on 1 thread but 278 s at 2 and
337 s at 16, monotonically slower with more threads. A systematic isolation
sweep (thread count × {glibc, mimalloc, jemalloc} × {default, numactl
--interleave}, 2 reps) showed *swapping* the allocator or NUMA policy made no
difference — but that pointed at allocation *volume* rather than allocator
*choice*. Eliminating the read loop's per-base heap allocation removed the
at-scale degradation: 16-thread dropped from 337 s to **249 s (−26 %)**, now
within ~14 % of single-thread, and the curve recovers with thread count instead
of diverging.

A residual 1→2-thread step (+28 %, ~218→280 s) remains; it is allocator- and
NUMA-independent (identical under mimalloc/jemalloc), i.e. a memory-bandwidth /
cross-socket floor, not something further code removes. Net: in-process
threading is now **break-even — no longer harmful, but no net speedup** over
single-thread, which remains the fastest per-worker mode.

The decisive experiment for the whole-genome strategy was a **procs-per-node
sweep**: K independent single-threaded eidolon *processes* concurrently on one
node. Processes **scale** — aggregate throughput rises ~3.5× from 1 to 64
(0.0118 → 0.0409 jobs/s) before the bandwidth ceiling plateaus it:

| processes/node | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|---|---|
| throughput (jobs/s) | 0.012 | 0.014 | 0.018 | 0.028 | 0.037 | 0.041 | **0.041** | 0.039 |
| per-process efficiency | 100% | 60% | 39% | 29% | 19% | 11% | 5% | 3% |

Separate processes share no thread pool, sequence-block, or counter state and no
allocator, yet still plateau at 32–64 — so the ceiling is genuine per-node
**memory bandwidth**, library-independent. But because each single-threaded
process runs at full speed, packing many per node and spreading across nodes
extracts that bandwidth far better than threading a single process does.

**Consequence — eidolon's HPC parallelism model is multi-process region-sharding,
not threading.** A genome is partitioned into anchor-windows (eidolon's
generation-time BED filter assigns each fragment to exactly one window in
absolute coordinates, so shards reconstruct a whole-genome run with no
double-counting or gaps), simulated as one single-threaded process per window
via a SLURM array spread across nodes, then concatenated.

**First full GRCh38 run — and the contention discovery.** The initial
whole-genome run (306 × 25 Mb shards, one single-threaded eidolon per shard,
submitted as a SLURM array on Delta's *shared* `cpu` partition with a 16-core
slice per shard) completed in **3 h 41 m** — far above the ~20–26 min the
procs-per-node sweep projected. Per-shard timing was sharply **bimodal**: median
2.5 min, but p90 131 min and a worst shard of **170 min** for an ordinary 25 Mb
window that runs in ~8 min on an unloaded node (confirmed by re-running that exact
window on an isolated node: 7 m 53 s). Grouping shard time by compute node exposed
the cause: the slow shards clustered on nodes holding only 1–2 of our shards,
while nodes that ran 12–17 of our shards each were uniformly fast — the *inverse*
of what our own packing would produce. The slowdown therefore came from **other
tenants' jobs** co-scheduled on those nodes: because eidolon is
memory-bandwidth-bound (above), a bandwidth-hungry neighbour starves it,
inflating an 8-min window 10–20×. On a shared partition this is uncontrollable
and makes the whole-genome wall-clock irreproducible.

**Fix — exclusive nodes, packed by us.** Requesting whole nodes (`--exclusive`)
and packing K shards per node ourselves removes every stranger: the only
remaining contention is our own K processes competing for the node's bandwidth.
This worked — a single isolated 25 Mb window dropped back to **7 m 53 s** — but
the path to a clean whole-genome number then exposed two further trade-offs that,
together with the shared-partition result, form the complete HPC story.

**Trade-off 1 — straggler ceiling vs. completeness.** The first exclusive sweep
paired packing with a *tight* per-task walltime (`10 + 6K` min) as
straggler-proofing. But packing slows each window (Trade-off 2), so the ceiling
silently *killed* the heavy windows — every K's slowest shard landed exactly on
its walltime (K=4: 34 min, K=8: 58 min, K=16: 103 min) and runs finished
incomplete (197–268 of 306 shards). Lesson: once `--exclusive` removes the
cross-tenant straggler, the tight ceiling is not only unnecessary, it is
*harmful*; the walltime must comfortably exceed the packed-window time, and
`--requeue` alone covers genuine node failure. (Fixed to `30 + 20K` min.)

**Trade-off 2 — packing density vs. per-shard speed.** Even with strangers gone,
eidolon is hard memory-bandwidth-bound, so our *own* co-packed processes slow each
other **superlinearly**. From the clean, complete (uncapped) runs, the heaviest
25 Mb window took ~8 min alone (K=1), **29 min at K=2, and 101 min at K=4** — the
node's bandwidth saturates almost immediately (by K=2), so packing more shards per
node buys no extra throughput, only steeply rising per-shard latency. (The
human-genome penalty is harsher than the earlier yeast procs-per-node sweep
implied: a 25 Mb window's working set drives far more memory traffic than yeast.)
Wall-clock compute is therefore minimized by *low* K (each shard near full speed)
spread across *many* nodes, not by dense packing.

**Trade-off 3 — compute speed vs. schedulability.** But low-K-many-nodes is
exactly the hardest thing to schedule: `--exclusive` claims a full 128-core node
even for a 1-core (K=1) task, so minimizing compute wall-clock means requesting
dozens of whole nodes, which on a busy shared cluster queues behind everyone
else (observed: 77- and 153-task low-K arrays stuck `PD (Resources/Priority)`
with no estimable start). The bandwidth-optimal configuration is the
scheduling-pessimal one.

These three trade-offs resolve into one rule and three regimes:

| Approach | Compute speed | Schedulability | Reproducible? | Best when |
|---|---|---|---|---|
| Shared partition, any packing | unpredictable (3 h 41 m; 8-min windows → 170 min under co-tenants) | instant | no | never the right choice |
| Exclusive, **low** K, many nodes | fastest (~8 min/window) | poor on a busy cluster | yes | cluster is idle / reserved nodes |
| Exclusive, **moderate** K, fewer nodes | good | schedules quickly | yes | busy shared cluster (pragmatic) |

**Recommendation for HPC users.** Always run on **exclusive nodes** — a shared
per-core slice exposes a bandwidth-bound simulator to unbounded cross-tenant
contention and can inflate wall-clock 10–20× with no warning. Then **match
packing density to node availability**: pack *light* (1–2 shards/node) across as
many nodes as you can get when the cluster is idle or you hold a reservation
(minimizes compute wall-clock); pack *heavier* (8–16/node) when nodes are scarce
(fewer whole-node grabs → schedules sooner, at the cost of slower per-shard time).
Use a generous per-task walltime plus `--requeue`; never a tight ceiling. Tooling
(`make_shard_beds.sh`, `genome_array.sbatch` [exclusive packing],
`run_genome_reps.sh` [K-sweep × reps], `aggregate_reps.sh`, `merge_shards.sh`,
`tune_procs_per_node.sbatch`) is implemented and validated.

**Measured results (each run complete, 306/306 shards).** Two
completeness-verified exclusive runs against the shared baseline:

| Run | Whole-genome wall | Heaviest window | Median window | Effective nodes |
|---|---|---|---|---|
| Shared partition (16-core slice) | 3 h 41 m | 170 min | 2.4 min | ~varies (contended) |
| Exclusive, **K=2** | **2 h 30 m** | **29 min** | 0.8 min | ~30 (in waves) |
| Exclusive, **K=4** | 12 h 32 m | 101 min | 0.4 min | ~10 (starved) |

Two results stand out. First, **exclusive K=2 beat the shared run on every axis** —
faster wall-clock (2 h 30 m vs 3 h 41 m) and, more durably, a **6× tighter
per-window ceiling (29 min vs 170 min)** with no contention outliers; that ceiling
improvement is reproducible and scheduling-independent, where the wall-clock is
not. Second, **K=4's 12.5 h wall is almost entirely scheduling latency, not
compute**: its heaviest window was 101 min, but its 77 whole-node tasks trickled
through only ~10 concurrently-available exclusive nodes (vs ~30 for K=2, which was
submitted first and got nodes first). The K2/K4 wall gap is therefore
node-availability luck, not packing — Trade-off 3 made literal. With full node
availability, K=2's compute floor is a single ~29-min wave (153 nodes), i.e. a
whole human genome at 30× in **~30 min** when the nodes are actually there; the
2 h 30 m reflects getting only ~30 of them, in waves, on a busy day.

The honest framing for the report: eidolon's advantage is **single-thread
efficiency + low, flat memory**, and HPC throughput is reclaimed by **process-
level sharding** — a defensible, mechanistically-grounded story rather than a
(false) claim of linear thread scaling. Since each shard runs single-threaded on
a (relatively slow) Delta core, ongoing single-thread optimization (reducing
per-read allocations) compounds directly across the whole genome.

### 3.6.1 Post-fix simulation-only timing (v1.19.1)

The scaling investigation above ran on builds whose default log level was
`trace`; eidolon's per-base debug events fire on the order of `coverage × read_len
× reference_bp`, so a default run wrote a multi-gigabyte `.neat.log` and spent
most of its wall-clock on log I/O. Fixing the default to `info` and removing the
per-base debug calls (issue #340, v1.19.1) is **output-preserving** — the Tier-1
regression suite re-confirmed every fidelity and determinism metric unchanged —
but it is a large *performance* change. On the seeded regression benchmark it cut
wall-clock **3.4× on E. coli (33.0 → 9.79 s) and 4.2× on chr22 (255.0 → 61.37 s)**
at flat peak RSS, and the perf baseline was re-frozen on the fixed build. The
figures in §3.6 are therefore pre-fix and conservative, and the whole-genome
strategy was re-measured.

**Simulation-only timing (fixed build).** We re-ran the region-sharded
whole-genome pipeline in generation-only mode (FASTQ + truth VCF, no downstream
alignment or calling) on two whole genomes, capturing per-window wall-clock from
each shard's `/usr/bin/time`:

| Genome | Size | Windows | Coverage | Packing | per-window median / max | Total simulation compute |
|---|---|---|---|---|---|---|
| *Glycine max* (Wm82.a4) | ~1.01 Gb | 369 × 10 Mb | 30× | K=4 | 3.8 s / 49 s | **1.50 core-hours** |
| GRCh38 (primary) | ~3.1 Gb | 306 × 25 Mb | 30× | K=16 | 18.3 s / 128 s | **4.87 core-hours** |

Both runs completed every window (369/369, 306/306). "Total simulation compute"
is the single-core sum of all window wall-clocks — the reproducible,
node-availability-independent figure users can plan around: **a whole human
genome at 30× is ≈4.9 core-hours of simulation; a whole soybean genome ≈1.5.**

**How much the fix changed the picture.** The GRCh38 run is directly comparable
to the pre-fix numbers above — same reference, same 25 Mb windows. Pre-fix, an
isolated 25 Mb window took ~8 min (K=1) and packed heavy windows blew out to
29 min (K=2) / 101 min (K=4), with the whole genome at 2 h 30 m. On the fixed
build the **heaviest 25 Mb window is 128 s even at K=16 packing** — i.e. faster,
under 16-way contention, than the old single-process window was in isolation.
This points to the old superlinear packing penalty being substantially **trace-log
I/O contention** (each packed process writing a multi-GB log to Lustre) rather
than pure memory bandwidth: with logging gone, dense packing is cheap.

**Revised HPC guidance.** This simplifies the three-regime recommendation above.
Because per-window cost is now seconds-to-~2-minutes even at high packing, the
"low-K-for-speed vs high-K-for-scheduling" tension largely dissolves: **dense
packing (K≈16) is the sensible default** — fewest whole nodes (306 windows →
~20 node-tasks, a single scheduling wave) at negligible per-window latency cost.
Both runs sat far inside a generous walltime (human heaviest window 128 s against
a 110-min ceiling). The residual memory-bandwidth effect is real (per-window time
still rises with packing) but now second-order; a fixed-build K-sweep would
quantify the now-shallow curve. Exclusive nodes remain mandatory — the
shared-partition contention finding (§3.6) is independent of logging.

*Framing:* the core-hour totals are hard numbers (summed per-window wall-clock);
end-to-end wall-clock is still dominated by exclusive-node scheduling latency, so
we report the reproducible compute figure rather than a node-availability-dependent
wall time. All numbers are simulation-only — the eidolon contribution — deliberately
excluding downstream alignment/calling.

### 3.6.2 Threading and packing re-characterized (fixed v1.19.1 build)

Both scaling curves in §3.6 were measured under the old `trace` default. Re-running
them on the fixed build changes the *mechanism* — though not the whole-genome
*strategy*, which it puts on firmer ground.

**In-process threading scales — bounded by contig count.** The apparent thread
*regression* of §3.6 was per-thread trace-log I/O. On the fixed build a single
process scales with threads, up to the genome's contig count (eidolon's unit of
parallelism is the contig):

| threads | yeast 30× (16 contigs) | soy 10× (100s of contigs) |
|---|---|---|
| 1 | 44.9 s (1.0×) | 1237 s (1.0×) |
| 2 | 26.2 s (1.7×) | 1228 s (1.0×) |
| 8 | 9.8 s (4.6×) | 414 s (3.0×) |
| 16 | 7.2 s (**6.2×**) | 231 s (5.4×) |
| 64 | 7.5 s (6.0×) | 139 s (**8.9×**) |

yeast (16 contigs) peaks at ~6× at 16 threads and flattens beyond — no contigs left
to split. soy (hundreds of contigs) keeps climbing to ~9× at 64. Two real,
non-logging limits remain: scaling is **sub-linear**, and on a large working set the
**1→2-thread step is a wash** (soy 1.01×, the cross-socket bandwidth floor §3.6
noted), while a cache-resident genome takes it cleanly (yeast 1.7×).

**Multi-process packing is near-linear to 64/node.** §3.6's procs-per-node
efficiency collapse (100% → 5% at 64, "bandwidth saturates by K≈2") was per-process
log I/O. On the fixed build, K independent single-thread processes pack
near-linearly, and a follow-up on a **real 25 Mb shard-window** confirms it holds at
production window size — in fact *better* than the small ecoli workload:

| procs/node | 1 | 16 | 64 | 128 |
|---|---|---|---|---|
| ecoli (4.5 Mb) efficiency | 100% | 90% | 90% | 82% |
| **25 Mb window efficiency** | 100% | 100% | **99%** | 75% |

The 25 Mb window's per-process time is flat at ~93 s from 1 to 64 processes, rising
to 110 s only at 128 — a heavier per-process compute load amortizes fixed per-node
overhead better, so it packs *better* than ecoli. **64/node (one per physical core)
is the efficient sweet spot**, dropping only when hyperthreads are oversubscribed at
128.

**Strategy — region-sharding wins, on corrected and stronger grounds.** The
whole-genome recommendation (multi-process region-sharding, §3.6.1) stands, but the
reasoning is now: (1) region-sharding parallelizes at *sub-contig* window
granularity — arbitrary parallelism, and the *only* option for few-/single-contig
inputs (chr22) where threading does nothing; (2) independent shards pack
near-linearly to ~64/node (~99 % on a real 25 Mb window), whereas threading one
process is capped at contig count and reaches only ~6–9×. So many thin shards ≫ few
fat-threaded processes, and packing is far cheaper than §3.6's pre-fix "K≈2" claim —
64 shards/node (one per physical core) is efficient, and a whole human genome (306
windows) is then just ~5 node-tasks. In-process threading remains a useful
*secondary* lever for running a multi-contig genome as a single job (~6–9×), but
sharding is the throughput path.

### 3.7 Structural-variant validation (chr22, Manta + truvari)

The cancer workflow's structural variants (DEL / DUP / INV / BND / CNV) had only
ever been validated through Mutect2, which scores SNVs and indels — never through
an SV-aware caller. We closed that gap: simulate an SV-rich tumor/normal pair
(60×, purity 0.8), call somatic SVs with **Manta** (tumor/normal mode), and score
against eidolon's SV truth with **truvari**.

On the classes Manta and truvari can score (DEL/DUP/INV), eidolon's variants —
**including large events up to 1 Mb** — are recovered at **0.88–1.00 recall
(32/35 = 91%)** in this single run.

**Replicated (n=18 = 3 tissue models × 3 seeds × 2 callers).** Because that single
run rests on few events per type, we replicated SV recovery at higher SV count.
The replicated means are **DEL 0.84 ± 0.11, DUP 0.85 ± 0.11, INV 0.92 ± 0.09**
(Manta ≈ Delly throughout). Replication **corrected the single run's optimism** —
in particular `INV 1.000` was a small-n artifact; the replicated mean is 0.92 with
a 0.75–1.00 range. The qualitative conclusion is unchanged and now carries honest
confidence intervals: DEL/DUP/INV recover at ~0.84–0.92, tool-independently. The
strongest evidence, though, looks *below* the caller, at the
reads themselves: eidolon does not edit the reference; it emits breakpoint signal
through a dedicated chimeric-read pass that stitches flanks across each junction,
and that signal is present and correctly placed for **every** class — including
where the caller fails to recover it.

- **CNV — two ways.** Manta cannot call copy number. (i) *By depth:* a `CN=4`
  region showed 110.5× vs a 60.4× flank = **1.83×**, matching the expected
  `(0.8·4 + 0.2·2)/2 = 1.8×` at purity 0.8. (ii) *By a copy-number caller:* GATK
  `DenoiseReadCounts`→`ModelSegments`→`CallCopyRatioSegments` recovered **4/7 truth
  CNVs with the correct gain/loss direction** (the misses are the small/low-amplitude
  het CNVs — kb-scale events near the read-depth segmentation limit). Both confirm
  eidolon's CNVs carry correct, caller-detectable copy-ratio signal.
- **Somatic specificity (no leak).** A homozygous somatic deletion was depleted
  5× in the tumor (60×→12×) while the normal stayed at full depth — the SV is real
  in the reads and correctly tumor-restricted.
- **BND signal is in the reads even when uncalled.** Every truth-BND locus
  inspected carried **8–73 discordant/split reads** — the same range as the
  *detected* BNDs. Manta calls only ~half of them, re-representing them as
  `DUP:TANDEM`/`DEL` at the correct coordinates.

  > **⚠️ Retraction (2026-08-01).** This bullet previously concluded: *"the misses
  > are a *caller* limitation — translocations are the hardest short-read class —
  > not absent signal."* The first and last clauses stand: the signal is present,
  > and the reads are correct. **The attribution does not.** eidolon has never
  > generated a translocation — the de novo sampler hardcodes the breakend mate to
  > the same contig (`sv_model.rs:810`, `:936`), so all 466 junctions in job 20719077
  > were intra-chromosomal by construction. The explanation invoked a class of event
  > that was not in the data.
  >
  > The measured cause is our own representation. We type every junction `SVTYPE=BND`
  > regardless of geometry, but a direct same-contig adjacency *is* a deletion or a
  > tandem duplication, and Manta calls it as one. Job 20745149 confirmed this
  > exactly — of 932 truth BND records, **all 452 direct-form records were false
  > negatives and none appeared in `tp-base`**:
  >
  > ```
  > tp-base: t[p[ 0    t]p] 188  [p[t 166  ]p]t 0     -> direct 0,   inv-oriented 354
  > fn:      t[p[ 226  t]p] 68   [p[t 58   ]p]t 226   -> direct 452, inv-oriented 126
  > ```
  >
  > Manta's recall over junctions it can represent as breakends is **354/480 = 0.738**
  > (head-to-head 0.734, tail-to-tail 0.741 — no bias within the class). The
  > "~half" figure was a ceiling set by us, not a caller limitation. See
  > `docs/cancer_simulator.md` for the full retraction and its cause.

Per-type recovery (Manta + truvari, 60×/0.8):

| SV type | Truth | Recovery | Bounded by |
|---|---|---|---|
| DEL | 17 | **0.882** (15/17) | — (incl. 376 / 645 kb) |
| DUP | 14 | **0.929** (13/14) | — (incl. 1 Mb) |
| INV | 4 | **1.000** (4/4) | — |
| BND | 22 | **superseded — see §3.7.1** | *the claim in this row was wrong* |
| CNV | 7 | depth 1.83× + GATK caller 4/7 by direction | Manta does not call CNV |

Overall figures across all 64 truth SVs (recall 0.500, precision 0.561) are bounded
by BND and CNV. The CNV bound is a genuine tool constraint (Manta does not call copy
number). **The BND bound was not**: it was a simulator defect compounded by two
scoring errors, and the paragraphs below that attribute it to truvari are retracted.
§3.7.1 gives the corrected account.

**Confirmed by an independent second caller.** Delly (read-pair/split/depth-based,
no assembly) recovers the same SVs at the *same* rates — **DEL 0.882, DUP 0.929,
INV 1.000, identical to Manta** — with slightly higher precision (0.667 vs 0.561),
and likewise scores 0/22 BND and 0/7 CNV. Two independent callers agreeing this
closely is strong evidence eidolon's SV data is not tuned to one tool.

The inference we drew from the shared BND zero was wrong, and instructively so. Two
callers agreeing told us the cause was common to both — we concluded that meant the
*scorer*, when the other thing common to both is the *input*. It was the input. See
§3.7.1. The CNV half of the claim stands: neither caller emits truvari-matchable
copy-number calls.

**Per-tissue models shift the spectrum correctly.** Swapping the pan-cancer model
for tissue-specific COSMIC SV models reproduces each tissue's dominant SV type in
the truth set, tracking the model's type probabilities, while per-type *recovery*
stays tissue-independent:

| Tissue | Dominant truth type (✓ model) | DEL / DUP / INV recall | BND frac | Overall |
|---|---|---|---|---|
| BRCA | **DUP** (37, > DEL 27) — model Dup 0.32 | 0.89 / 0.94 / 1.00 (0.94) | 32% | 0.585 |
| skin | **BND**-enriched (45) — model Bnd 0.31 | 0.87 / 0.94 / 0.80 (0.89) | 47% | 0.432 |
| lung | **DEL** (35) — model Del 0.35 | 0.86 / 0.93 / 1.00 (0.90) | 27% | 0.614 |

Each tissue's largest bucket matches its model's largest type probability.
Scoreable-class recall holds at **~0.89–0.94 regardless of tissue**; a tissue's
*overall* recall simply tracks how much of its spectrum is BND (skin 47% → 0.43,
lung 27% → 0.61) — a scorer limitation, not detection quality.

**Recall vs. tumor purity — robust, with a coverage-model caveat at the extreme.**
Holding model and total coverage fixed (60×) and sweeping purity, DEL+DUP+INV recall
is **0.91 / 0.94 / 0.94** at purity 0.3 / 0.5 / 0.7, then drops to **0.17 at 0.9**.
The collapse is *not* a detection failure: `cancer_simulate` splits a fixed coverage
budget by purity (`normal = (1−purity)·total`), so at purity 0.9 the matched normal
falls to **6×**, and Manta — unable to score somatic confidence against so thin a
normal — filters 37 of 68 calls as `MinSomaticScore` (which `--passonly` then drops).
Recall is stable as long as the matched normal stays ≳18× (purity ≤0.7 here); the
extreme-purity drop is a simulation coverage-model artifact (#315: decouple
matched-normal depth from purity), not SV generation. Real tumor/normal pairs
sequence the normal at an independent depth, avoiding this.

**Confirmed by a controlled re-run (§4).** To check that the collapse is this
normal-starvation artifact and not a purity-handling defect, we re-ran the sweep
on the SNV/Mutect2 path with the matched normal **pinned at 30×** (`total =
30/(1−purity)`, so tumor depth varies while the normal stays healthy). Somatic
SNV recall then **holds across purity — 0.87 / 0.97 / 0.97 / 0.97** at 0.3 / 0.5 /
0.7 / 0.9 — with no collapse at the extreme; it in fact *rises* with purity and
plateaus at the detection ceiling for purity ≥ 0.7, dipping at 0.3 only because
the pinned-normal design leaves less tumor depth (13×) there. This confirms the
extreme-purity drop is matched-normal starvation common to any tumor/normal
caller (Manta or Mutect2 losing the germline reference it subtracts) — the #315
coverage-model artifact — not eidolon's variant generation. The same-caller
fixed-budget SNV arm (the direct "before" to this "after") makes it airtight:
splitting a fixed 60× budget, SNV recall is **0.88 / 0.97 / 0.96** at purity
0.3 / 0.5 / 0.7 and then **collapses to 0.003 at 0.9** (normal 6×) — a near-total
loss, *even more severe than the SV path's 0.17*, because Mutect2 is more
normal-dependent than Manta. Same caller, same metric, only the budget policy
differs: the collapse is entirely matched-normal starvation.

This is the first end-to-end validation of eidolon's structural-variant output, and
the read-level signal is correct across all classes and sizes. (Detection needs
adequate depth and SV count: a 30×/0.6 chr22 run yields too few somatic SVs to
measure; 60×/0.8 gives a scoreable set.)

**Per-tissue SV spectra reproduce at scale (§4).** Exercising the three bundled
per-tissue SV models (BRCA / skin / lung) on chr22 at the default SV rate yields
only ~3 somatic SVs per run — far too few to express a type distribution. Forcing
the SV rate up to ~200–280 events per tissue makes the signatures clear: **skin is
the most BND-enriched (59 %), lung the most DEL-leaning (40 %), and BRCA the most
DUP-enriched (24 %, ~2× the others)** — the last exactly the expected BRCA1/2
tandem-duplicator phenotype, and all three consistent with the PCAWG-derived models
(#202 / #237). Two caveats: differentiation requires enough events (a scale
requirement, not a defect — a small target won't show it), and BND is elevated
across all tissues, so the tissue signal lives in the *relative* enrichments rather
than the absolute mode.

> **⚠️ Correction (2026-08-01).** This paragraph previously called the BND elevation
> "a generation tendency, and unscored downstream since truvari does not benchmark
> breakends." Both halves were wrong. Truvari has matched breakends since v5.0.0 via
> `-B/--bnddist` and Delta runs v5.4.0 — see §3.7.1, where the same inherited belief is
> retracted. And the elevation is not a tendency but a direct inheritance: the BND share
> is the PCAWG per-donor mean of `svclass == "TRA"` rows
> (`normalize_pcawg_sv_model.py:108-114`), i.e. the **inter-chromosomal translocation
> rate** — while the generator emits only same-contig junctions. The per-tissue BND
> column therefore reports each tissue's translocation burden attached to events that are
> not translocations. The *relative* enrichments across tissues may still be meaningful;
> the absolute BND fraction is not. See the retraction in `docs/cancer_simulator.md`.

**Honest caveats — tracked as realism enhancements (epic #311).** The simulated SVs
are *idealized*: clean cuts in unique sequence, lacking the microhomology, imprecise
breakpoints, and repeat / segmental-duplication context that make real somatic SVs
hard to call. The 0.88–1.00 DEL/DUP/INV recall is therefore an **upper bound**
relative to real, repeat-embedded variants (#312). Two further items: eidolon encodes
tandem-dup junctions as generic BND pairs rather than the canonical `DUP` that
callers and truth sets expect (#313), and SV-scale insertions did not appear in any
run (`INS: 0`) despite a non-zero model probability (#314). None is a defect in the
validated read generation — all are realism / interoperability improvements.

### 3.7.1 BND: measuring our own defect (v1.13.1 → v3.1.0)

`BND recall = 0.000` stood in this report from the first SV validation until
2026-08-01, attributed across successive drafts to truvari's inability to match
breakends. The defect dates to v1.13.1 (2026-06-02). **All of those
attributions were wrong**, and the sequence of wrong answers is itself the finding.

**What was actually true.** eidolon's read generator dispatches on two flags that were
only ever assigned by the input-VCF parser, so every *de novo* breakend left them at
their defaults — VCF 4.2 case 4, a DIRECT join — while the ALT written to the truth VCF
said `t]p]`, case 2, reverse-complemented. **The truth VCF described a rearrangement the
reads did not carry**, from v1.13.1 (2026-06-02) until v3.1.0 (2026-08-01). And a direct intra-contig join is not a breakend at
all: joining A-left to B-right *is* a deletion. eidolon had been planting DELs and DUPs
and labelling them BND.

Manta had been right the whole time. On Delta job 20675480 it placed DEL/DUP candidates
within **1–2 bp of all seven** planted junctions and emitted no breakend at any of them —
exactly what direct joins should produce. The reads were never at fault either: 7–15
discordant pairs and 3–15 split reads at every one of the 14 breakends.

**Three claims retracted:**

1. *"truvari cannot benchmark breakends."* False since truvari v5.0.0, which added
   `-B/--bnddist`; Delta runs v5.4.0. The harness simply never passed the flag. This was
   inherited from an early draft and repeated without being checked against the tool's
   changelog.
2. *"Two callers agreeing on zero proves the scorer is at fault."* Two callers agreeing
   localises the cause to something common to both — which is as easily the **input** as
   the scorer. It was the input.
3. *"BND signal is present at read level, so the simulator is fine."* The signal was
   present, but it was the signal for a different event than the truth declared.

**Corrected measurements (August 2026, jobs 20684366 / 20699004 / 20716772):**

| metric | before | after | replication |
|---|---|---|---|
| `manta_BND` recall | 0.000 | **1.000** (14/14) | 2 runs |
| `manta_INV` recall | 0.000 | **0.833** (5/6) | 2 runs |
| `manta_overall` recall | 0.625 | **0.833** | 2 runs |

`manta_INV` was itself a representation artifact: Manta stopped emitting `<INV>` around
v1.5 and reports inversions as breakend pairs, shipping `convertInversion.py` to convert
them back. The harness now runs it, and scores BND against the **raw** calls and INV
against the **converted** ones — one file per question, because the converter consumes
the breakends it converts.

**A residual that is ours, quantified.** `manta_INV precision = 0.312` (FP = 11), of
which **7 are our own junctions**: eidolon emitted only the head-to-head geometry, which
*is* the inversion signature, so `convertInversion.py` converts our junctions into `<INV>`
records that then score as false positives against a truth that never contained them.
Confirmed by matching each FP's `END` against the truth's junction mates — 7 land on
junctions, 4 on real inversions (Manta redundancy). Excluding our own artifact, precision
is **5/9 ≈ 0.556**.

**Status of this section.** BND geometry sampling was implemented (#458 item 2) and the
post-geometry campaign has now run — jobs 20719077 and 20745149, soy whole assembly, 466
junctions against the 7 the numbers above rest on. **The numbers above are superseded.**

> **⚠ Read these as MECHANISM verification, not DATA verification.** Both jobs ran at
> `SV_RATE_SCALE=40` — forty times the model's nominal rate — chosen to get enough
> junctions for a four-way geometry distribution to be readable. `sv_pipeline.sbatch`
> defaults to `1.0`, so this was a deliberate override and not the pipeline's normal
> setting.
>
> What that does and does not affect:
> - **Unaffected:** the geometry mix (226/256/224/226) and reciprocity. Geometry is drawn
>   per junction, independently of how many junctions there are.
> - **Affected, and unquantified:** every recall/precision figure, including
>   `BNDinv = 0.738`. A genome carrying 40x the realistic SV density is not the genome a
>   caller is tuned for — junctions sit closer together, overlap rejection saturates
>   (visible as `requested N ... only placed M` warnings), and the realized SV spectrum
>   can drift from the model's. Whether caller recall moves with SV density has not been
>   measured here.
>
> A data-verification run means `SV_RATE_SCALE=1.0` on a real genome, at whatever scale
> is needed to accumulate junctions — not extrapolation from a rate-inflated one. That
> run has not been done, and until it is, §3.7's recall numbers characterise the
> mechanism working, not the caller's performance on realistic input.

Measured on the real truth, all four VCF 4.2 forms present and indistinguishable from the
model's declared uniform weights (χ² ≈ 1.5, df = 2, p ≈ 0.47), with 0 unpaired and 0
mispaired records across 932:

```
t[p[ 226 (24.2%)   t]p] 256 (27.5%)   [p[t 224 (24.0%)   ]p]t 226 (24.2%)
```

The headline correction: **`manta_BND recall = 0.380` is not a caller result.** Callers
bucket only inversion-oriented adjacencies as BND; a direct same-contig adjacency is a
deletion or tandem duplication and is called as one. All 452 direct-form truth records
were false negatives and none appeared in `tp-base`. Manta's recall over the junctions it
can represent as breakends is **354/480 = 0.738** (head-to-head 0.734, tail-to-tail 0.741).
The harness asserts truth reciprocity before scoring. It also used to report the
inversion-oriented subset separately as `BNDinv`; **that split was retired in #507** and the
paragraph above should be read as history. It was a consequence of the de novo sampler
hardcoding the breakend mate to the anchor's own contig — within one chromosome a direct
adjacency genuinely *is* a deletion or tandem duplication. Since P3a, BND means an
inter-chromosomal translocation, where no deletion or duplication interpretation exists and a
caller emits BND for all four geometries. The two figures then collapse (job 20824321: BND
0.935 vs BNDinv 0.929; pooled over campaign 20925151, manta 0.868 vs 0.876 and delly 0.417 vs
0.412), so `manta_BND` is simply the BND number. `BNDinv` precision was additionally an
artifact — the truth was restricted while the query was not, so every unmatched call became a
false positive (#512).

**§3.7 is still not submittable, for a deeper reason than the one previously given.**
eidolon had never generated a translocation at the time of writing — the de novo sampler
hardcoded the breakend mate to the same contig. **Superseded: P3a made de novo BND
inter-chromosomal**, so the paragraph below describes the state that motivated §3.7's caveat
rather than the current one. The caveat itself stands until a translocation-era run is
scored. "BND recall" here means **same-contig junction recall**, and the
BND category itself is under audit: its 23.2% share is inherited from PCAWG's
inter-chromosomal `TRA` count, and `BND` in a VCF is a representation convention rather
than an event class (Manta and Delly typed identical junctions differently in job
20745149). See the retraction in `docs/cancer_simulator.md`. A PCAWG measurement pass
precedes any further BND claims here.

### 3.8 At-scale validation (chr1–3 subset, ~690 Mb)

To confirm the chr22 metrics are not overfit to one chromosome, both pipelines were
re-run on a **chr1–chr3 subset of GRCh38** — ~690 Mb, ~13× chr22, entirely different
sequence. Every metric held or slightly improved:

| Pipeline | Metric | chr1–3 | chr22 |
|---|---|---|---|
| Germline (eidolon) | SNP recall / precision | **0.990 / 0.9997** | 0.989 / 0.9996 |
| Germline (eidolon) | INDEL recall / precision | **0.988 / 0.990** | 0.982 / 0.989 |
| Cancer (Mutect2 T/N) | somatic SNV recall / precision | **0.944 / 0.910** | 0.925 / 0.872 |
| Cancer (Mutect2 T/N) | somatic INDEL recall / precision | **0.908 / 0.885** | 0.900 / 0.844 |

Germline Ts/Tv held at 2.34 (truth = query). eidolon's fidelity is therefore not an
artifact of chr22 — it holds at 13× the scale on different sequence. Whole-genome
runs were not pursued (no divergence to investigate); the Mutect2 tumor/normal step
is the practical scale ceiling (~3.2 h on chr1–3, and ≫48 h projected genome-wide
without interval scatter).

*Tooling note:* hap.py/som.py hardcodes UCSC chr-prefixing of input VCFs, so an
Ensembl-named reference (GRCh38: `1/2/3`) needs a chr-prefixed copy for scoring —
now handled automatically by the cancer pipeline.

### 3.9 Cross-caller coverage and mutational-signature fidelity

To avoid judging eidolon's cancer features by a single caller, second callers were
added across classes on the existing Delta harnesses — **Delly** (somatic SV, vs
Manta), **Strelka2** (somatic SNV/indel, vs Mutect2), and **GATK somatic-CNV** (vs
the depth check of §3.7) — plus a **mutational-signature** check
(SigProfilerAssignment). Cross-caller *agreement* is the strongest anti-overfit
signal; the signature check probes something variant-recall cannot: does eidolon
reproduce the COSMIC mutational *signature* its tumor model encodes?

**We found that at the signature level, it did not** Fitting
COSMIC signatures to 6,225 simulated somatic SNVs (chr1–3) assigned them entirely
to flat, clock-like signatures (SBS5 45 %, SBS3 31 %, SBS41 17 %, SBS12 7 %), with
**no SBS1** — the near-universal C>T-at-CpG signature present in essentially every
real tumor. The 96-context spectrum confirms it: the four CpG `[C>T]G` contexts are
the *lowest* C>T contexts (16–23 counts vs ~65 average).

**Root cause (verified both ends).** The bundled COSMIC model was faithful — its
per-context substitution model encodes strong CpG C>T enrichment (conditional
0.78–0.88 at CpG vs 0.39–0.62 elsewhere). Digging into the code, we found a limitation
that we had introduced early in eidolon's lifecycle, but had not noticed on normal, small
data runs. The code places mutations at a **context-independent rate** and conditions
only the *alt allele* on trinucleotide context. The realized spectrum is therefore
*(genome trinucleotide frequency) × (conditional alt)*, not the COSMIC signature;
because CpG dinucleotides are ~10× genome-depleted, the CpG C>T peak never forms. 
COSMIC signatures are **rate patterns** (mutations-per-context), which a
uniform-placement model cannot represent.

**Context** 
When we dug into the legacy code, likely introduced early in eidolon's lifecycle as a 
stopgap until the more complex version could be ported. That time turned out to be now,
as it was limiting out ability to model cancer genetics effectively.

The code could eidolon faithfully reproduce the substitution-*type*
distribution and transition bias (Ts/Tv 2.34, §3.1/3.3) and recovers variants well
across independent callers and at scale (§3.3/3.4/3.7/3.8) — but simulated somatic
SNVs did not reconstruct the input COSMIC *signature* under signature extraction.
NEAT 2.x and 4.x both weight SNP *placement* by trinucleotide context (`sample_trinucs` /
`init_trinucBias`), but the Rust port applied context only to the *alt allele* and
placed positions uniformly (the germline default model shared the gap). The fix was to
complete the context port, which had been neglected and missed in previous tests. Since this was
exactly the class of gap we added the signature check to find, this bug also shows why running
at scale and vetting the statistics was an essential development step.

> **Update (v1.20.0, #372).** Context-weighted placement has since been restored: SNP
> positions are now drawn with probability proportional to each site's fitted
> trinucleotide mutability (matching NEAT 2/4), on by default; context-flat models keep
> uniform placement byte-identically. Re-validated on SEQC2 HCC1395, the SBS-96 cosine of
> real-vs-simulated somatic SNVs rose from **0.72 to 0.99**. Signature-*extraction* was then
> re-run on the current build (chr1–3, 6,263 somatic SNVs, matching this section's ~6,225):
> SigProfilerAssignment now recovers **SBS1 (391) and APOBEC SBS2/13 (709)** — both absent in
> the run above — in a sensible pan-cancer blend. **#320 is resolved/closed.** The chr22 spot-check
> (SBS1 6/366) was small-N noise; at the proper SNV count the recovery is unambiguous.

**Cross-caller results (#317).** The broadened coverage reinforces the anti-overfit
picture across independent callers. **SVs:** Delly reproduces Manta's per-type recall
*exactly* (DEL 0.882 / DUP 0.929 / INV 1.000; §3.7). **SNV/indel:** VarScan2 — a
different statistical model from Mutect2 — calls eidolon's somatic variants at high
**precision** (0.93 SNV / 0.94 indel, vs Mutect2's 0.87 / 0.84): the calls it makes
*match the truth*, so the data is not a Mutect2-specific artifact. Its lower recall
(0.63 / 0.53 vs Mutect2's 0.93 / 0.90) reflects VarScan2's more conservative model
and the high-confidence (`Somatic.hc`) filter, not a simulator property — Mutect2
(sensitive) and VarScan2 (precise) bracket the truth as expected for two callers.
(Strelka2 was dropped — its 2018 build SIGSEGVs on Delta's stack. GATK somatic-CNV
recovered **4/7 CNVs by direction** — no-PoN, since its panel step is a Spark tool
incompatible with the env's Java 25 — corroborating the depth check, §3.7.)

### 3.10 Reproducibility, determinism, and sharding correctness

**Reproducibility** is used here in the standard computational sense — the same inputs
(reference file, config, seed) run through the same code produce the same result
([ACM](https://arxiv.org/pdf/2402.07530); [The Turing Way](https://book.the-turing-way.org/reproducible-research/overview/overview-definitions/)).
A dedicated harness (`run_order_independence.sbatch`) changes exactly one variable at a
time at a fixed seed and compares the header-less, sorted truth-VCF body by md5 (the
contig-name check compares the `(POS,REF,ALT)` multiset, since renaming legitimately
relabels `CHROM`). Run on the soybean-cyst-nematode assembly (`GCA_040805935.1`, ~145 Mb
and highly fragmented — 1,178 shard windows, 157,945 variants), the invariances that
*must* hold, do:

| Check | Result | Evidence |
|---|---|---|
| Determinism (rerun, 1 thread) | **PASS** | identical (`5bcbaf5…`) |
| Thread-invariance (1 thread vs 8, set + count) | **PASS** | identical (`5bcbaf5…`) — parallelism never perturbs results |
| Contig-name invariance (rename every contig, order fixed) | **PASS** | identical `(POS,REF,ALT)` (`431120a…`) + count (157,945); only `CHROM` relabels |
| Shard-order independence (fwd vs reversed merge) | **PASS** | identical (`44caf2b…`) |
| Shard disjointness (1,178 windows) | **PASS** | 0 duplicate `CHROM:POS:REF:ALT` keys |

**Reproducibility horizon.** eidolon's variant realization is a function of the seed,
config, each contig's sequence, and contig *order* — not of contig *name*, thread count,
or shard-merge order. The RNG derives one child stream per contig from its *index*
(`derive_child(contig_idx)`), so renaming a contig (order preserved) yields the identical
variant set, merely relabeled, and thread scheduling never perturbs it (a dedicated
contig-name-invariance check confirms this). Re-ordering contigs in the reference is a
*different input file* — a different checksum, hence a different artifact under any of the
definitions above — and produces a different realization, which is the expected behavior
of a deterministic program on changed input, not a reproducibility gap. Making output
invariant to contig re-ordering would be an optional *canonicalization* (key each stream
by name rather than index); it is out of scope and accepted-and-documented, since no
workflow requires reproducing a run from a genome with different contig ordering.

**Sharding.** The region-sharded whole-genome path draws each window from an independent
per-`(contig, chunk)` seed, so shard-order independence and disjointness both hold and the
merge never depends on reproducing a monolithic RNG stream. A sharded run is therefore
intentionally distinct from a single monolithic run — here 157,945 variants monolithic vs
157,962 across the 1,178 shards; the ~17-variant difference is the independent per-window
seeds, and the disjointness check (0 duplicate keys) confirms it is not double-counting.
That independence is exactly what makes the parallel path correct. Net: eidolon is
reproducible and parallelism-invariant where it matters, and its HPC sharding is verifiably
correct.

### 3.11 Adapter readthrough validation (chr22, 30×, TruSeq) — and a bug only Delta caught

`eidolon` added optional Illumina 3′ adapter readthrough (#125): when a fragment's
insert is shorter than the read length, the read is padded at its 3′ end with adapter
sequence, exactly as a real sequencer produces. We used Delta both to **confirm the new
feature is callable** and — as it turned out — to **catch a correctness bug that no local
check had surfaced**.

**Design.** A four-arm matrix on the germline variant-calling pipeline
(eidolon → BWA-MEM2 → GATK HaplotypeCaller → hap.py), chr22 at 30×, TruSeq preset, three
replicates per arm, paired by seed so the truth set is identical across arms. The arms are
built to *isolate the adapter effect from the insert-size effect* — necessary because the
no-adapter baseline rejects short fragments, so a naïve off-vs-on comparison confounds the
two:

| arm | inserts | adapter |
|---|---|---|
| `off` | long (short fragments rejected) | none — baseline |
| `short_ctrl` | short (kept) | none — genomic reads |
| `on_raw` | short | readthrough, aligned raw |
| `on_trim` | short | readthrough, fastp-trimmed |

**The bug Delta caught.** The first at-scale run collapsed: the adapter arms recovered
almost no variants (SNP recall ≈ 0.0004) and their BAMs held ~2,986 reads instead of ~3.9M.
This was *not* reproducible from read counts alone — `zcat | wc` reported the full read
count locally, and eidolon's own logs reported success. The cause was a malformed FASTQ record
(quality string one character longer than the sequence) emitted for degenerate zero-length
inserts; a too-long quality line leaves the 4-lines-per-record count intact but makes
**bwa-mem2's parser stop at the first offending record**, silently truncating alignment.
The defect was therefore invisible to every local sanity check and only manifested when the
data passed through a real aligner at production scale — precisely the kind of failure the
Delta verification pipeline exists to surface. Once diagnosed (bisected to zero-insert
fragments, 5,513 malformed records in a 3.9M-read run), the fix was a one-line guard plus a
regression test asserting `seq.len() == qual.len()` for every emitted record. It is
catalogued as defect #125 in §3.5.

**Result after the fix.** Adapter readthrough is realistic and fully callable. The realism
signals move only when the feature is enabled, and the fidelity contrasts show the adapters
themselves cost nothing:

| realism signal | off | on_raw | on_trim |
|---|---|---|---|
| fastp adapter reads | 5,913 | 2,425,876 | 2,425,876 |
| soft-clip fraction | 0.0004 | 0.0804 | 0.0006 (recovered by trim) |
| mean insert size | 214.6 | 178.5 | 178.5 |
| mean read length | 151 | 151 | 139 (adapter trimmed) |

| fidelity contrast (Δ = A − B) | isolates | snp_recall | indel_recall |
|---|---|---|---|
| `on_trim` − `short_ctrl` | pure adapter effect (same inserts) | +0.0003 | +0.0018 |
| `on_raw` − `on_trim` | adapter handling (soft-clip vs trim) | −0.0002 | −0.0009 |
| `short_ctrl` − `off` | short-insert coverage effect (no adapters) | −0.0263 | −0.0334 |

Adding adapter readthrough to short-insert reads costs no fidelity (within replication
noise), whether the adapters are fastp-trimmed or left for BWA-MEM2 to soft-clip. The only
fidelity difference versus the long-insert baseline is the reduced effective coverage of
short-insert libraries (mate overlap) — present in `short_ctrl` with adapters entirely
absent, so it is not attributable to readthrough.

**Infrastructure note.** These runs also hardened the harness against Delta's shared
filesystem: after a Lustre OST dropped mid-run and wedged jobs in uninterruptible I/O wait,
the germline pipeline was moved to node-local NVMe staging (heavy FASTQ/BAM I/O on the
compute node's local disk, only small artifacts copied back to Lustre), making it immune to
transient OST outages; a fastp thread-count livelock at high core counts was also capped.

### 3.12 Subclonal heterogeneity — intra-tumor VAF spectrum (#405)

NEAT and pre-#405 eidolon model a tumor with a single global `purity` scalar, so every
somatic variant collapses to essentially one effective VAF — no intra-tumor
heterogeneity. eidolon now models a tumor as a **mixture of subclones at distinct
cancer-cell fractions (CCF)**: each somatic variant's observed alt fraction composes as
`purity × dosage × CCF` (three orthogonal axes — normal contamination, per-copy
multiplicity, subclonal prevalence). The architecture can be authored inline, **fitted
from real data** (PyClone-VI / PCAWG cluster tables), or **replayed** from a real somatic
VCF at its observed VAF. Every somatic truth record carries the intended CCF (`INFO/EIDOLON_CCF`)
and the intended observed VAF (`INFO/EIDOLON_VAF`) as ground truth.

> **Correction (2026-07-31).** The figures in this section were, from #405 landing
> (2026-07-24) until #450 was fixed a week later, computed over **407 of 567 planted
> sites** — an entire CCF cluster was silently
> excluded, because the harness measured observed VAF through `bcftools call`, which
> discards a hom-ref call together with its read count (#450). The exclusion was
> VAF-dependent, so it flattered the result in the direction that mattered. Re-measured
> on Delta job 20675479 with the genotype-calling step removed: **543 of 567 sites**,
> three CCF clusters present where two were, `bias = −0.0072`, `MAE = 0.0291`. The
> clusters land where `purity × dosage × CCF` predicts (0.40 / 0.16 / 0.04 at purity
> 0.8), which the code cannot fake. A 4.2 % residual remains unscored, distributed
> 7/9/8 across the three clusters — reported per stratum by the harness rather than
> requiring two tables to be differenced by hand.

**Real-data validation (Delta).** A round-trip the unit tests cannot reach: simulate,
align the merged reads with bwa-mem2, and measure the observed VAF at each somatic site
against the intended `EIDOLON_VAF`. Two complementary runs — a *controlled synthetic*
architecture and a *real tumor's empirical* spectrum:

| Metric | Generative — soy scaffold | Reproductive — HCC1395 chr1 |
|---|---|---|
| Input VAF spectrum | 3 synthetic CCF clusters (≈0.04 / 0.16 / 0.40) | real tumor's empirical VAF, all deciles 0.1–1.0 |
| Coverage / purity | 200× / 0.8 | 200× / 0.9 |
| Somatic sites (compared) | 567 planted (387 scored) | 2,948 (2,698 scored) |
| Bias, mean(observed − intended) | −0.003 | −0.004 |
| MAE (vs ~200× noise floor) | 0.026 | 0.024 |
| Pearson r | 0.95 | **0.99** |

The **generative** run (SEQC2-agnostic, single soybean scaffold, three subclones)
confirms the `purity × dosage × CCF` composition is unbiased and accurate to the sampling
limit; per-site Pearson is a spread-dependent summary for tight discrete clusters, so
fidelity is gated on bias + MAE-vs-noise-floor rather than r.

The **reproductive** run replays SEQC2 **HCC1395** (a triple-negative breast cancer cell
line) — its high-confidence somatic SNVs honored at their observed VAF over GRCh38 chr1.
This is a real tumor's full *continuous* VAF distribution (2,948 sites spanning every
decile), so Pearson is fully meaningful: the aligned reads reproduce it at **r = 0.99**,
unbiased (−0.004), with per-decile MAE uniform at 0.009–0.027 across the whole 0.1–1.0
range. Together these show eidolon reproduces both a *designed* subclonal architecture and
a *real tumor's* empirical VAF spectrum — a complexity axis (intra-tumor heterogeneity) the
predecessor could not represent at all. Harnesses: `scripts/delta/run_subclonal_vaf_validation.sh`
(generative) and `scripts/delta/run_hcc1395_reproductive.sh` (reproductive).

---

## 4. Phase 2 — robustness at scale (planned, key deliverables)

Phase 1 established correctness and performance on a focused test set (chr22 and
small genomes). Phase 2 closes the four gaps that Phase 1 deliberately did not
cover, with an explicit goal of **avoiding over-tuning to chr22**: exercise eidolon
on substantial and *varied* inputs to confirm robustness and surface any
remaining defects of the kind already found.

1. **Whole-genome scale — done at chr1–3 (§3.8).** Germline and cancer pipelines
   re-run on a ~690 Mb chr1–3 subset of GRCh38 (~13× chr22): all metrics held or
   slightly improved (germline SNP 0.990 / indel 0.988; cancer somatic SNV 0.944 /
   indel 0.908), confirming nothing was overfit to chr22. Full whole-genome not
   pursued — no divergence to investigate, and Mutect2 T/N is the scale ceiling
   (≫48 h genome-wide without interval scatter).

2. **Multicore scaling / HPC tuning — COMPLETE (§3.6).** The thread regression
   was characterized (memory-bandwidth-bound, allocator- and NUMA-independent) and
   the strategy resolved and *measured*: multi-process region-sharding on exclusive
   nodes. The full GRCh38 sweep produced the shared-vs-exclusive contrast (3 h 41 m
   contended → 2 h 30 m clean at K=2, with a 6× tighter per-window ceiling), the
   superlinear packing curve (8→29→101 min/window at K=1/2/4), and the three-regime
   HPC-usage recommendation. The only follow-on is opportunistic: a low-K run on a
   *reservation* to demonstrate the ~30 min compute floor when whole nodes are
   actually available.

3. **Exercise the new cancer features deeply — SV validation DONE (§3.7).**
   Structural-variant realism was validated downstream with Manta + truvari: eidolon's
   DEL/DUP/INV (to 1 Mb) recover at 0.88–1.00 recall, CNV confirmed by depth,
   somatic specificity confirmed — the first end-to-end check of the SV machinery.
   The BND half of this claim, and the "no simulator defect found" that accompanied
   it, are **retracted**: a de novo breakend defect had been misread as a scorer
   limitation (§3.7.1). The post-geometry campaign has since run (jobs 20719077 /
   20745149): the interim `manta_BND recall = 1.000` on 7 junctions is superseded by
   **0.738 over 480 inversion-oriented records** across 466 junctions, and the
   aggregate 0.380 is floored by our own SVTYPE choice rather than by the caller.
   **A second, deeper retraction applies to this item:** the BND work cannot be
   called DONE at all, because eidolon has never generated a translocation — the de
   novo sampler places every breakend mate on the same contig. See the retraction in
   `docs/cancer_simulator.md`; a PCAWG measurement pass precedes any further BND claim.
   Remaining sub-items: exercise the **per-tissue cancer models** (BRCA / skin /
   lung) to confirm tissue-specific SV spectra downstream, and sweep the
   **purity-mixing** machinery across purity/coverage.

4. **Input variety and parameter sweeps.** Run across a range of genome sizes and
   compositions (bacterial, fungal, invertebrate, mammalian) and across cancer
   parameters (purity, coverage, tumor mutation rate) to ensure eidolon handles
   diverse inputs rather than a single tuned case.

5. **Simulation-only timing for the community.** Benchmark *just the eidolon
   simulation stage* (FASTQ + truth VCF generation, no downstream alignment or
   variant calling) across genome sizes, coverages, and the sharded whole-genome
   configuration, to publish concrete real-world wall-clock / throughput numbers
   users can plan around (e.g. "whole human genome at 30× in N minutes on K
   nodes"). This isolates the simulator's own performance — eidolon's core
   contribution — from pipeline overhead, and is the headline deliverable for
   communicating eidolon's HPC performance to the community. (The sharded GRCh38
   sweep in item 2 is the first such measurement at whole-genome scale: a complete
   human genome at 30× in 2 h 30 m on ~30 exclusive nodes, ~30 min compute-bound.)

---

## 5. Phase 3 — stretch goals (time permitting)

- **Long-read simulation (ONT / PacBio-style)** (tracked: #319). Validate eidolon against the
  long-read paradigm: simulate single-molecule reads (kb-scale lengths and the
  higher, indel-dominated, homopolymer-dependent error profile characteristic of
  nanopore / HiFi data, unpaired), align with **minimap2**, and score recall with
  long-read callers — **Clair3 / DeepVariant** for SNV/indel and **Sniffles2 /
  cuteSV** for SVs (truvari for SV scoring, as in §3.7). Because a single long read
  spans an SV breakpoint, this is the natural complement to the short-read SV
  validation: it stresses the SV machinery from the other side and would exercise
  the breakpoint-realism work in epic #311 directly. *Prerequisite:* a long-read
  mode in eidolon (read-length distribution + long-read error model) — today eidolon
  targets short paired-end Illumina data, so this is feature work gated ahead of
  the validation, mirroring the short-read harnesses already built
  (`germline_e2e` / `cancer_pipeline` / `sv_pipeline`).
- **Plant genomes and polyploidy.** Exercise large, repetitive plant genomes and
  scope tuning for polyploid simulation (eidolon's `ploidy` parameter and the
  allele-dosage work tracked separately), where structural variation and high
  copy number stress the simulator differently than human data.

---

## 6. Resource usage

| | Core-hours |
|---|---|
| Documented production runs (9, archived) | **491** |
| Total charged incl. exploratory / failed / timed-out runs | ~2,500 |
| Remaining allocation | ~397,500 |

The 9 archived production runs (3 baseline, 1 benchmark, 4 cancer, 1 germline)
total 491 core-hours. The gap to the ~2,500 charged is dominated by a single
timed-out benchmark that ran `--exclusive` (charging all 128 cores for 8 h ≈
1,000 core-hours before being redesigned) plus failed/iterated runs during
bug-fixing. Phase 1 was deliberately economical (chr22-scale); Phase 2's
whole-genome runs and tuning sweeps are the primary consumers of the remaining
allocation, where the binding constraint is per-run wall-clock, not core-hours.

**Note (v1.19.1):** these charged core-hours were real, but they *overstate* the
compute eidolon actually needs — most runs predate the logging fix (#340) and ran
3–4× slower than necessary under the old `trace` default. On the fixed build the
figures are strikingly small: a whole soybean genome at 30× is **~1.5 core-hours**
of simulation and a whole human genome **~4.9** (§3.6.1), and the tuning sweeps
that motivated much of the "expensive" framing are no longer needed. Future
Phase-2 work costs a fraction of the Phase-1 charge rate; the allocation is not a
constraint.

---

_Run-level metrics and provenance for every job (version, git commit, reference,
core-hours) are archived under the project results directory and assembled by
`scripts/delta/collect_report.sh`. Note the archive spans **two roots** because the
project was renamed `rusty-neat`/`rneat` → `eidolon` partway through (v2.0.0): Phase 1
runs — everything backing §3.1–§3.11 — are under
`/projects/bhrd/jallen17/rneat-access-results/`, while runs from v2.0.0 onward are under
`/projects/bhrd/jallen17/eidolon-access-results/`. Each job directory records its own
version and git commit, so provenance is per-run regardless of which root holds it._
