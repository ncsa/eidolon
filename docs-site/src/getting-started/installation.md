# Installing eidolon

Three routes. The first two need no Rust toolchain at all.

## Bioconda

```
conda install -c bioconda eidolon
```

A prebuilt binary with dependencies handled. This is the easiest route on a workstation
and works on a cluster wherever you already have conda.

## A release binary

Each release publishes binaries on the
[releases page](https://github.com/ncsa/eidolon/releases). Pick the one matching your
system:

| asset | for |
|---|---|
| `eidolon-x86_64-unknown-linux-gnu` | most Linux distributions |
| `eidolon-x86_64-unknown-linux-gnu-rhel8` | RHEL 8 and derivatives, where the general Linux build's glibc is too new |
| `eidolon-aarch64-unknown-linux-gnu-rhel8` | 64-bit ARM on RHEL 8 |
| `eidolon-x86_64-apple-darwin` | macOS on Intel |
| `eidolon-x86_64-pc-windows-msvc.exe` | Windows |

The binary is self-contained — it needs nothing beyond a standard C library, so there is
no environment to set up. Download it, make it executable, and run it. Let us know if one
of these does not work on your system.

## Building from source

You will need the Rust toolchain including `cargo`; see the
[cargo installation guide](https://doc.rust-lang.org/cargo/getting-started/installation.html).
You will also need `git`.

```bash
git clone git@github.com:ncsa/eidolon.git
cd eidolon
cargo build --release
```

Compilation may stop on a missing system dependency — `cmake` is the usual one, since some
of the crates `eidolon` depends on need it. On Debian/Ubuntu that is
`sudo apt install cmake`; on RHEL/Rocky, `sudo dnf install cmake`. There may be others
depending on your system. Drop a comment on an issue if you need specific help.

Build from source when you are testing a change: a released binary is, by construction,
not the code you are working on. See [Running on HPC](../hpc/running-on-hpc.md) for
building on a cluster.

```bash
$ eidolon --help
Usage: eidolon [OPTIONS] [SUB-COMMAND]

SUB-COMMANDS:
  gen-reads              Generates reads for an input dataset
  filter-reads           Filters the output of gen-reads
  gen-mut-model          Generates a mutation model from real VCF data
  gen-seq-error-model    Generates a sequencing error model from real FASTQ data
  gen-frag-length-model  Generates a fragment length model from a BAM or SAM file
  gen-gc-bias-model      Generates a GC bias model from a reference FASTA and aligned BAM
  gen-bam-models         Builds multiple models (frag-length, GC bias) from one BAM in a single pass
  compare-vcfs           Compares a NEAT-simulated golden VCF against a downstream variant-caller VCF
  gen-cancer-reads       Simulates a tumor/normal mixture (two gen-reads passes + merge)
  validate               Checks an emitted FASTQ/VCF/BAM against what its consumers accept
  compare-af             Per-allele AF correlation between a truth VCF and a simulated one
  help                   Print this message or the help of the given subcommand(s)

Options:
      --log-level [<log_level>]  Verbosity of the written .neat.log. The on-screen log is always
                                 info. [default: info] [possible values: trace, debug, info,
                                 warn, error, off]
      --log-dest <log_dest>      Sets the log destination (full path with full filename) for the written log
  -h, --help                     Print help
```

`validate` and `compare-af` take file paths and flags directly rather than a
configuration YAML — both are diagnostic/measurement runs invoked ad hoc and from
harness scripts, where writing a YAML per check is friction.

To check options for a subcommand:

```bash
$ eidolon gen-reads --help
Generates reads for an input dataset

Usage: eidolon gen-reads [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

To run filter reads, check the help menu.

```bash
$ eidolon filter-reads --help
Filters the output of gen-reads

Usage: eidolon filter-reads --configuration-yaml <configuration_yaml>

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

To run gen-mut-model:

```bash
$ eidolon gen-mut-model --help
Generates a mutation model from input VCF data

Usage: eidolon gen-mut-model [OPTIONS]

Options:
  -c, --configuration-yaml <configuration_yaml>  Path to configuration file.
  -h, --help                                     Print help
```

Use the help menu to see the available options and leave an issue if you find something bad happening.

To compile and run `eidolon` yourself, besides the rust toolchain, you will need `git` installed for your operating system. You will then need to git clone and cd into the repo directory. From your home directory in Linux the process might look something like:

```bash
~/$ git clone git@github.com:ncsa/eidolon.git
~/$ cd eidolon
~/eidolon/$
```

For Windows and Mac users, try the binary packaged with the latest version of `eidolon`.

Once in the repo, you can build the program either in debug (default) or release mode. The main difference is how much info it gives you if there is an error. Release mode also has some optimizations to run it faster.

```bash
~/eidolon/$ cargo build --release
```

If you prefer to run the package directly without using the binary, you can also use
```bash
~/eidolon/$ cargo run -- gen-reads -c my_config.yml
```

Rust will download any required packages. Compiling Rust code is the slowest part of the process. The final binary will be built and the program run immediately after in the second case. To run the program manually, from the repo main dir, run

```bash
~/eidolon/$ ./target/release/eidolon -h
```

`eidolon` uses a configuration file to read values it needs for the run. A command line execution might look like this:

```bash
~/eidolon/$ ./target/release/eidolon -c /path/to/filled/in/config.yml
```
If you record the output in the logs of Seed string to regenerate these exact results: XXXXXXX, you should be able to use that string as input with rng_seed and reproduce your results.
