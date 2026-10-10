# Running on HPC

`eidolon` runs as a single process and fits naturally onto a single HPC compute node. No MPI, distributed computing, or special environment setup is required. The notes below cover resource budgets for whole-genome human-scale runs.

### gen-reads

`gen-reads` is already multi-threaded via rayon (see [Parallel Processing](../gen-reads/parallel-processing.md) above). Set `num_threads` in your config to match your CPU allocation:

```yaml
num_threads: 16
```

Memory scales with `num_threads × largest_contig_size`, not total genome size — each thread processes one contig at a time. For hg38, chr1 is ~249 Mbp, so 16 simultaneous threads need roughly 16 GB for sequence data plus overhead. Budget additional scratch space for per-contig temp files before the final assembly: roughly 3× the expected output size.

| Genome | Threads | Recommended RAM | Scratch space | Typical wall time |
|--------|---------|-----------------|---------------|-------------------|
| Bacterial (~5 Mbp) | 4 | 4 GB | 5 GB | < 1 min |
| Human hg38, 10× SE | 16 | 32 GB | 100 GB | ~20 min |
| Human hg38, 30× PE | 16 | 48 GB | 300 GB | ~60 min |

Example SLURM header:

```bash
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --tmp=300G
```

### gen-gc-bias-model

Single-threaded; walks the BAM once and accumulates per-base reference coverage as a `Vec<u32>` per observed contig. Peak memory is roughly `4 bytes × total reference bases that received at least one record`. For a full human genome (~3.2 Gbp) this is ~13 GB. For a single chromosome or a targeted/exome BAM it is much smaller.

If memory is tight on whole-genome runs, split the BAM by chromosome and run `gen-gc-bias-model` once per chromosome on a representative subset that spans your GC range of interest, then average or combine the resulting models.

Recommended SLURM header (whole-genome BAM):

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
```

### gen-mut-model

Single-threaded; processes the VCF one chromosome at a time. A full human genome WGS VCF (~10 GB) typically completes in 30–60 minutes with peak RSS under 4 GB.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:30:00
```

### gen-seq-error-model

Single-threaded; streams through the FASTQ. For a full genome FASTQ (~600 M reads), expect 30–60 minutes. `max_reads` caps how many records the fit uses, drawn uniformly from the whole file rather than taken from its start. A FASTQ is written in flowcell order, so its first records all come from one tile of one lane. A capped run still reads every record, once to count them and once to sample, so it shortens the fit but not the pass over the file:

```yaml
max_reads: 5000000   # a uniform sample of ~5 M records; 0 uses every record
```

Measured on Delta with eidolon v4.0.0, one CPU, on a 29 GB gzipped HG002 FASTQ: a 5 M-record cap took 11 minutes and 53 MB; fitting every record took 14 minutes. Most of the time is the pass over the file, which is why the cap saves little.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:30:00   # capped or not, every record is read
```

### gen-frag-length-model

Single-threaded; streams TLEN fields from the BAM. A full human genome BAM at 30× typically completes in under 15 minutes with peak RSS under 2 GB.

```bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
```

### Environment notes

- The `eidolon` binary has no runtime dependencies beyond a standard C library (glibc), which is present on all Linux HPC systems.
- No module loads or conda environments are required — copy the release binary to your scratch or project directory and run it directly.
- If you compile from source on the cluster, ensure `cmake` is available (`module load cmake` on most systems) before running `cargo build --release`.
- **You do not need a Rust module from your administrators.** If the site has no `rust`
  module, `rustup` installs a toolchain into `$HOME/.cargo` as an ordinary user, and it
  persists across jobs:
  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
  source "$HOME/.cargo/env"
  ```
- If you do build on the cluster, send the artifacts to scratch:
  `export CARGO_TARGET_DIR="$SCRATCH/cargo-target/eidolon"`. The toolchain itself is
  modest, but a release build's `target/` directory runs to several GB of crate artifacts
  against a ~5 MB binary, and an HPC home quota is usually small enough that this is what
  fills it.
