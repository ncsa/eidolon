# Parallel Processing
`eidolon gen-reads` processes contigs in parallel by default using rayon's work-stealing thread pool. Each contig is an independent unit of work — variant generation, fragment sampling, and FASTQ/BAM writing all happen concurrently across contigs — so references with many contigs scale well across cores.

Output is **byte-identical regardless of `num_threads`** (the same seed always produces the same reads in the same order), so you can change the thread count freely without affecting results.

A note on scaling: read generation is largely **memory-bandwidth bound**, so on references dominated by one or a few large chromosomes the wall-time gain from extra cores is limited — the bottleneck is memory throughput, not CPU. An experimental sub-contig **chunk size** knob (below) can split a single chromosome across cores, but it is **disabled by default** because in our benchmarks it did not improve wall time and occasionally regressed it.

## Thread count

By default `eidolon` uses all available logical cores. You can cap the thread count with the `num_threads` config key:

```yaml
# use 4 threads instead of all available cores
num_threads: 4
```
or disable parallelism entirely (useful for debugging or reproducibility testing)
```yaml
# disable parallelism entirely (useful for debugging or reproducibility testing)
num_threads: 1
```

Omit `num_threads` (or set it to `.`) to restore the default all-cores behaviour.

## Thread count and hardware — fewer threads can be faster

Read generation moves a lot of data relative to the arithmetic it does, so it is
largely **memory-bandwidth bound**. On a typical desktop or laptop (which has
only a couple of memory channels), a single `eidolon` thread can already saturate
much of the available memory bandwidth — so adding cores yields little speedup
and, past a point, can even run *slower* as threads contend for the memory bus.
In our desktop benchmarks, large references ran about as fast (sometimes faster)
on **1–4 threads** as on all 8.

Practical guidance:

- **Desktop / laptop:** try `num_threads: 1` (or a small number) and compare — it
  is often as fast or faster than all-cores for a single large run, and leaves
  cores free for other work. If you are simulating **many samples**, running
  several single-threaded `eidolon` jobs in parallel typically beats one
  many-threaded job.
- **HPC nodes** with many memory channels (and multiple sockets) have far more
  aggregate bandwidth, so higher `num_threads` scales better there.
- Output is byte-identical regardless of `num_threads`, so it is safe to tune the
  thread count purely for speed on your hardware.

## Chunk size (experimental, opt-in)

By default each contig is one unit of work. Setting `chunk_size` splits contigs
into sub-contig chunks so a single large chromosome can be worked by several
cores at once. It is **off by default**: read generation is memory-bandwidth
bound, so chunking did not improve wall time in our benchmarks (and added a
little overhead). It is kept for CPU-bound or very-many-core scenarios where it
may help. The size is in base pairs and independent of the thread count, so
output stays byte-identical regardless of `num_threads`.

```yaml
# default — disabled (one chunk per whole contig); omit or set to . or 0
chunk_size: 0

# opt in: fixed chunk size in base pairs
chunk_size: 5000000
```

*When might it help?* Only when read generation is **CPU-bound** rather than
memory-bandwidth bound — the opposite of what we measured on a typical desktop.
Consider trying it when **both** are true:

1. Your reference is dominated by **one or a few large contigs**, so the default
   per-contig parallelism leaves most cores idle (e.g. a single-chromosome
   assembly, or a genome where one chromosome dwarfs the rest).
2. You have **memory bandwidth to spare relative to cores** — a multi-socket or
   many-memory-channel HPC node rather than a commodity desktop — and/or you run
   compute-heavier settings (GC-bias-weighted coverage, long-read mode, very high
   coverage) where per-read CPU work dominates memory traffic.

It is **not** worth enabling on a typical workstation, or on references with many
contigs (those already parallelize per contig). Always benchmark on your own
hardware: start with `chunk_size ≈ longest_contig_bp / (4 × cores)` (a few chunks
per core), compare wall time against the default, and keep it only if it is
actually faster.

## BAM output is fully parallel

`eidolon` uses a per-contig temp-file strategy: each contig worker writes its alignment records to a private temporary BAM body file, then a single concatenation pass assembles them in reference order into the final coordinate-sorted BAM.

## Reproducibility

Each contig's random number generator is derived deterministically from the parent seed and the contig's position in the reference, so output is identical across runs with the same seed even when the number of threads changes.
