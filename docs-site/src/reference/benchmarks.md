# Current Benchmarks
We ran some benchmarks on `eidolon` on my home desktop, a very basic Linux desktop, using data pulled from public resources on the internet (mostly ncbi). The "yeast" is brewer's yeast. These are the results:

```text
[16:26:39] Build complete.
eidolon Benchmark Report
Run:        20260510_162639
Machine:    pop-os  |  8 logical CPUs
Coverage:   10x    |  Read length: 151 bp
Binary:     /home/joshfactorial/code/eidolon/target/release/eidolon
──────────────────────────────────────────────────────────────────────
SECTION 1: Single-ended read generation
Genome          Size(MB)   Wall time  Peak RSS(MB)      CPU%
──────────────────────────────────────────────────────────────────────
[16:26:40] Single-ended: ecoli  (4.6 MB)
ecoli                4.6     0:14.83         416.9       99%
[16:26:54]   Done: wall=0:14.83  peak_rss=416.9 MB  cpu=99%
[16:26:54] Single-ended: pneumonia  (2.1 MB)
pneumonia            2.1     0:05.39         197.0       99%
[16:27:00]   Done: wall=0:05.39  peak_rss=197.0 MB  cpu=99%
[16:27:00] Single-ended: yeast  (12 MB)
yeast                 12     0:10.94         981.7      352%
[16:27:11]   Done: wall=0:10.94  peak_rss=981.7 MB  cpu=352%
[16:27:11] Single-ended: c_elegans  (97 MB)
c_elegans             97     7:35.23        8591.2      387%
[16:34:46]   Done: wall=7:35.23  peak_rss=8591.2 MB  cpu=387%
──────────────────────────────────────────────────────────────────────
Section 1 — Passed: 4 / 4  |  Failed: 0
──────────────────────────────────────────────────────────────────────
SECTION 2: Paired vs single comparison  (genome: yeast, 10x)
  Paired-end uses fragment_mean=400, fragment_st_dev=50
  Same reference, same coverage, same seed — only the read-generation mode differs.
Mode               Wall time  Peak RSS(MB)      CPU%
──────────────────────────────────────────────────────────────────────
[16:34:46] Paired comparison — single-ended pass
single-ended         0:10.48         970.4      365%
[16:34:57]   Single done: wall=0:10.48  peak_rss=970.4 MB  cpu=365%
[16:34:57] Paired comparison — paired-ended pass
paired-ended         0:13.30        1126.8      354%
[16:35:10]   Paired done: wall=0:13.30  peak_rss=1126.8 MB  cpu=354%
──────────────────────────────────────────────────────────────────────
  Paired/single wall-time ratio: 1.27x  (10.48s → 13.30s)
  Paired/single peak-RSS ratio:  1.16x  (970.4 MB → 1126.8 MB)
──────────────────────────────────────────────────────────────────────
```
