# Prerequisites

The easiest way to install `eidolon` is via [Bioconda](https://bioconda.github.io/):

```
conda install -c bioconda eidolon
```

This pulls a prebuilt binary with all dependencies handled — no Rust toolchain
required. If you prefer to build from source or grab a release binary, read on.

You will need to install the rust toolchain to compile `eidolon`, including `cargo`. Check the cargo documentation for instructions (https://doc.rust-lang.org/cargo/getting-started/installation.html). Alternatively, you can try one of the binaries on the release page. Select the one that matches your system and let us know if you run into errors. During compilation, you may run into errors, such as cmake not found. Some of the packages `eidolon` uses have these dependencies. For Debian/Ubuntu this should be a simple `sudo apt install cmake` and for RHEL/Rocky type distros this should be `sudo dnf install cmake`. There may be some other requirements. Drop a comment if you need specific help.

Download the executable in the release (current version 3.2.0).

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
