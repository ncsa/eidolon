# Running eidolon on an HPC cluster

This guide is for users who want to run `eidolon` at scale — from a bacterial genome
on one node up to a whole human (or larger) genome across many nodes — on a SLURM
cluster. It distills what we learned validating eidolon on NCSA Delta (see
`docs/access_report_draft.md` §3.6 for the underlying measurements) into
cluster-agnostic guidance.

The SLURM scripts referenced here live in `scripts/delta/`. They carry Delta
specifics (module names, an `--account`, conda env names); adapt those three
things to your site and the rest is portable.

---

## 1. The one-paragraph mental model

eidolon is **single-thread efficient**, with **low, flat memory**. In-process (rayon)
threading does scale — but only up to the genome's *contig count* (eidolon's unit of
parallelism), so it caps at ~6–9× and does nothing for a single-contig input. The
way to use a big machine is therefore **many single-threaded processes**, each
simulating a different region of the genome, packed onto nodes: this parallelizes
at *sub-contig* granularity and, on current builds, packs **near-linearly** — many
more independent shards than threading one process can ever reach. This is *region
sharding*, and it is the HPC story.

Three practical consequences follow, and the rest of this guide is just their
detail:

1. **Run one thread per process** (`num_threads: 1`) and pack many per node.
   (Threading a single process is a fine *secondary* lever for a multi-contig
   genome run as one job — ~6–9× — but sharding gives far more.)
2. **Run on exclusive nodes.** On a shared node a co-tenant's bandwidth-hungry job
   can inflate per-window time many-fold with no warning; a whole node to yourself
   is the only way to get reproducible wall-clock (see §5).
3. **Pack densely.** Independent single-thread shards pack at ~99% efficiency to
   ~64 per node (one per physical core) and ~75% at 128 — so match K to core count,
   not the old "keep it low" rule (see §5).

---

## 2. Before anything at scale: the log-level footgun

**Set `--log-level info` (the default as of v1.19.1) and never run a large job at
`debug` or `trace`.**

eidolon's per-base `debug!`/`trace!` events fire on the order of
`coverage × read_length × reference_bp` times. At `trace`, a single chromosome-22
run wrote a multi-gigabyte `.neat.log` and spent *most of its wall-clock on log
I/O*. Before v1.19.1 the default was `trace`; if you are on an older build, pass
`--log-level info` explicitly on every run. On any version, only turn up verbosity
for a *small* debugging run, never a production one.

```bash
eidolon --log-level info gen-reads -c sim.yml    # good for production
# eidolon --log-level debug ...                  # ONLY on a tiny debugging input
```

---

## 3. Quickstart: a small genome on one node

For anything that fits comfortably in one node's memory and finishes in a
reasonable single-process wall-clock (bacteria, fungi, small invertebrates), you
do not need sharding at all. One single-threaded process:

```yaml
# sim.yml
reference: /path/to/genome.fa
read_len: 250
coverage: 30
ploidy: 2
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
produce_fastq: true
produce_vcf: true
produce_bam: false
overwrite_output: true
output_dir: /scratch/$USER/run1
output_filename: sample
rng_seed: my run 1
num_threads: 1
```

```bash
#SBATCH --partition=<your-cpu-partition>
#SBATCH --nodes=1
#SBATCH --exclusive          # or a generous --cpus-per-task + --mem if you can't get exclusive
#SBATCH --time=01:00:00
eidolon --log-level info gen-reads -c sim.yml
```

If you want more throughput on that node, do not raise `num_threads` — instead run
several of these processes concurrently on different regions (see the sharding
approach below), or several independent samples at once.

**Memory:** eidolon's resident memory is dominated by the reference contigs it loads,
not by coverage. Reference loading is `target_bed`-aware and loads only the contigs
a run touches, so peak RSS per process ≈ *the size of the largest contig that
process simulates* plus a small overhead. Budget accordingly: a 25 Mb window on
human chr1 still loads all of chr1 (~250 MB); a whole small genome loads all of it.

---

## 4. Whole-genome scale: region-sharded array jobs

The genome is partitioned into contiguous **anchor-windows**. eidolon's `target_bed`
is a *generation-time* filter in absolute reference coordinates that anchors each
fragment to exactly one window (reads may extend past a window edge; ownership is
by anchor point), so the union of all shards reconstructs a whole-genome run **with
no double-counting and no gaps**, and the per-shard outputs **merge by simple
concatenation**. Each window is one single-threaded eidolon process; a SLURM array
spreads them across nodes.

The tooling (in `scripts/delta/`):

| Script | Role |
|---|---|
| `make_shard_beds.sh` | partition `REF.fa.fai` into `shard_NNNN.bed` windows; prints the shard count |
| `genome_array.sbatch` | SLURM array; each task owns one **exclusive node** and runs `SHARDS_PER_NODE` windows on it concurrently. Idempotent + `--requeue`. |
| `run_genome_reps.sh` | driver: sweep packing density `K` and/or run replicates with distinct seeds |
| `aggregate_reps.sh` | fold rep manifests into mean ± sd wall-clock tables |
| `merge_shards.sh` | concatenate per-shard FASTQ + `bcftools concat|sort` the truth VCF into one dataset |

**Do a fail-fast smoke test first** (see §7) so you never commit a big node
request to a config that dies in the first shard.

Then the full run:

```bash
# 1. index + partition (25 Mb windows → ~124 shards for GRCh38 primary assembly)
samtools faidx REF.fa
N=$(bash scripts/delta/make_shard_beds.sh REF.fa.fai 25000000 $SCRATCH/shards)

# 2. choose packing density K and derive node-tasks (K=64 packs ~99%; see §5)
K=64; T=$(( (N + K - 1) / K ))

# 3. submit the exclusive-node array (%20 caps concurrent nodes — tune to availability)
SHARD_DIR=$SCRATCH/shards OUTROOT=$SCRATCH/wg_run REFERENCE=REF.fa SHARDS_PER_NODE=$K \
  sbatch --array=0-$((T-1))%20 scripts/delta/genome_array.sbatch

# 4. after the array finishes, merge (on a compute node — it's I/O-heavy)
SHARD_OUTROOT=$SCRATCH/wg_run sbatch scripts/delta/merge_shards.sh
```

For a **simulation-only timing** run (no downstream alignment/calling — the number
users actually want to plan around) this is the whole job: steps 1–3 produce the
FASTQ + truth VCF and the array's `time.txt` files hold the per-window wall-clock.

---

## 5. Choosing packing density (K) and whether to request a reservation

Good news: on current builds this is barely a decision. Independent single-thread
shards pack **near-linearly** — measured on a real 25 Mb human window (30×), one
exclusive node runs:

| shards/node (K) | per-window time | per-process efficiency |
|---|---|---|
| 1 | ~94 s | 100% |
| 16 | ~93 s | 100% |
| 64 | ~94 s | **99%** |
| 128 | ~110 s | 75% |

Per-window time is essentially flat to **64 shards/node** (one per physical core),
dropping only when you oversubscribe hyperthreads at 128. **Default to K=64.** A
whole human genome (306 × 25 Mb windows) is then just **~5 node-tasks**, each
finishing in ~2 min of compute.

Two knobs:
- **K=128** packs the node fullest → fewest node-tasks (~3 for a human genome) at
  ~75% efficiency. Use it when whole nodes are scarce.
- **Lower K** only if a genome's largest contig is big and RAM is tight — peak RSS
  per node ≈ K × (largest contig loaded), since each shard loads its window's contig.

**Do you need a reservation? Almost never.** A whole genome needs only ~3–5
exclusive nodes at K=64–128, which schedules easily even on a busy cluster — you're
not trying to grab 150 nodes for a low-K run. Submit with a `%`-throttle matched to
node availability and let it drain; with no hard deadline, a reservation buys
nothing.

**Still run exclusive.** Dense packing only works because you hold the whole node.
A shared per-core slice re-exposes eidolon to cross-tenant memory-bandwidth theft,
which is unbounded and makes wall-clock unreproducible — that finding is real and
independent of packing density.

**Cost:** a whole-genome 30× simulation is only ~1.5–5 core-hours of actual compute
(soybean ~1.5, human ~4.9). Trivial against any allocation — the binding constraint
is wall-clock / node availability, not credits.

---

## 6. Reproducibility scope

eidolon is deterministic and parallelism-invariant *where it matters*, with one
documented boundary:

- **Determinism:** same seed + same reference + same params → byte-identical output.
- **Thread-invariant:** 1 thread vs N threads → byte-identical. Parallelism never
  perturbs results.
- **Shard-order independent + disjoint:** the sharded whole-genome path gives the
  same merged result regardless of shard completion order, with zero duplicate
  variants across windows. This is why sharding is safe.
- **Contig-order sensitive (by design, monolithic runs only):** a single
  *monolithic* (non-sharded) run consumes one global RNG stream in contig order, so
  reversing contig order in the FASTA changes the realization (even the variant
  count). **Consequence:** to reproduce a monolithic run exactly, use the
  *byte-identical* reference (same contig order), not merely the same genome. The
  sharded path is immune to this — each shard seeds independently from an
  absolute-coordinate window. (Tracked as #322: optional per-contig seeding to make
  even monolithic output contig-order-independent.)

Practical guidance: **for reproducible large runs, prefer the sharded path and
record your seed(s) + the exact reference file** (or its checksum).

---

## 7. Fail-fast before you commit big node requests

A large array that dies on a bad config or a malformed-output bug wastes both your
time and your queue priority. Two cheap guards:

1. **Smoke-test one window first.** `scripts/delta/wg_smoke.sh` runs a single small
   window as a normal (non-exclusive, few-minute) job and then *validates the output
   the way a downstream aligner would* — it checks that every FASTQ record has
   `len(seq) == len(qual)` and that reads and truth variants were actually produced.
   This is the exact class of defect (a malformed FASTQ that `wc -l` can't see but
   bwa-mem2 silently truncates on) that cost us a full at-scale run during adapter
   validation. Run it whenever you change the build, the config, or the reference:

   ```bash
   REFERENCE=$SCRATCH/neat_data/genome.fa bash scripts/delta/wg_smoke.sh
   ```

2. **Throttle the first real array** (`--array=...%N` with a small `N`) so a
   systematic problem surfaces on the first few nodes, not all of them at once.

`genome_array.sbatch` is idempotent (a completed shard is skipped on retry) and
uses `--requeue`, so a task killed by a sick node reruns elsewhere and only the
unfinished windows in its batch redo — you never lose completed work.

---

## 8. Shared-filesystem hygiene (Lustre / GPFS)

Large simulation runs move a lot of FASTQ/BAM. On a shared parallel filesystem:

- **Stage heavy I/O on node-local disk** (`$SLURM_TMPDIR` / local NVMe) when your
  site provides it, and copy only the small final artifacts back. A dropped Lustre
  OST otherwise wedges jobs in uninterruptible I/O wait. (`germline_e2e.sbatch`
  demonstrates this pattern.)
- **Keep manifests and bookkeeping on durable storage, never scratch** — scratch can
  throw transient EIO or hit quota mid-run, and you don't want a driver half-submit
  because a tiny write failed. (`run_genome_reps.sh` writes its manifest to `$HOME`.)
- **Prune intermediates.** Whole-genome 30× FASTQ is ~35–40 GB per run; delete BAMs
  and FASTQ after scoring if you only need the metrics (`PRUNE=1` in the harnesses).
- **Watch your scratch quota** — on Delta it is a *project* quota, not per-user, so a
  teammate's data counts against yours.

---

## 9. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Run is far slower than expected, huge `.neat.log` | `--log-level trace/debug` | use `--log-level info` (default ≥ v1.19.1) |
| Whole-genome wall-clock wildly variable, some windows 10–20× slow | shared-partition co-tenant contention | request `--exclusive` nodes |
| Array tasks stuck `PD (Resources/Priority)` forever | asking for more whole nodes than are free | raise K (toward 64–128 → fewer node-tasks), or lower the `%`-throttle |
| Packed windows all hit the walltime and finish incomplete | walltime ceiling below packed-window time | raise `--time` (never a tight ceiling once exclusive) |
| Downstream aligner recovers almost no reads despite full read count | malformed FASTQ record | run `wg_smoke.sh`'s integrity check; upgrade to ≥ v1.19.0 |
| Can't reproduce a previous run's variants | different reference contig order | use the byte-identical reference; or use the sharded path |
| Peak RSS higher than expected | a shard's largest contig is large | lower K so fewer big contigs load per node simultaneously |

---

*Measurements in this guide come from the eidolon validation campaign on NCSA Delta;
see `docs/access_report_draft.md` (§3.6 for scaling, §3.10 for reproducibility,
§3.11 for the malformed-FASTQ story) and `docs/regression_protocol.md` for how the
performance and fidelity baselines are frozen and gated.*
