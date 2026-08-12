#!/usr/bin/env bash
# One-time environment setup for eidolon on Delta (NCSA ACCESS).
# Run this interactively before submitting any SLURM jobs.
#
# What this does:
#   1. Installs rustup + a stable toolchain into $HOME/.cargo (persists across jobs)
#   2. Creates a conda env with NEAT 4 (the Python predecessor) for benchmarking
#   3. Creates a conda env (bioinf) with tools not available as Delta modules:
#        bcftools, bwa-mem2, gatk4
#      Note: samtools and htslib (tabix/bgzip) ARE available as Delta modules
#        (samtools/1.22-cce19.0.0 and htslib/1.22-gcc13.3.1) and are loaded
#        by the SLURM jobs directly. bwa-mem2 and gatk are not available.
#   4. Pulls the hap.py Apptainer image used by cancer_pipeline.sbatch
#
# Usage (from the repo root):
#   bash scripts/delta/setup.sh
#
# Re-running is safe — each step skips if already done.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Resolve Delta paths ($SCRATCH isn't exported on Delta). setup.sh runs
# interactively, so $0 is reliable here (unlike the spooled sbatch jobs).
source "$(dirname "$0")/lib_report.sh"
SIF_DIR="${SIF_DIR:-$SCRATCH/sif}"          # where to cache Apptainer .sif files
CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$SCRATCH/cargo-target/eidolon}"
CONDA_ENV_NAME="neat4"
BIOINF_ENV_NAME="bioinf"

echo "=== eidolon Delta setup ==="
echo "Repo:       $REPO_ROOT"
echo "SIF dir:    $SIF_DIR"
echo "Cargo target: $CARGO_TARGET_DIR"

# ── 1. Rust toolchain + eidolon binary ────────────────────────────────────
# Check for a module first; fall back to rustup if absent.
if module load rust 2>/dev/null; then
    echo "[1/4] Rust module loaded: $(rustc --version)"
elif command -v cargo >/dev/null 2>&1; then
    echo "[1/4] Rust already installed via rustup: $(rustc --version)"
else
    echo "[1/4] Installing rustup (no rust module found)..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
    source "$HOME/.cargo/env"
    echo "      $(rustc --version)"
fi
[[ -f "$HOME/.cargo/env" ]] && source "$HOME/.cargo/env"

# Redirect build artifacts to $SCRATCH — a full release build generates
# 10-15 GB of crate artifacts that would fill a typical HPC home quota.
export CARGO_TARGET_DIR
mkdir -p "$CARGO_TARGET_DIR"
echo "[1/4] Building eidolon release binary (artifacts → \$SCRATCH)..."
cd "$REPO_ROOT"
# Say what is about to be built, and refuse to build a checkout that is behind its
# remote. A partial or interrupted `git pull` leaves the tree on older code and setup
# happily builds it: job 20682989 ran a 3.0.1 binary against a 3.1.0 checkout, spent
# 25 core-hours, and produced a complete set of plausible numbers that tested none of
# the fix it was submitted to verify. Nothing in the old output could have revealed
# that, because setup never reported what it built.
echo "      checkout: $(git describe --tags --always --dirty 2>/dev/null || echo unknown)"
if git rev-parse --abbrev-ref --symbolic-full-name '@{u}' >/dev/null 2>&1; then
    upstream="$(git rev-parse --abbrev-ref --symbolic-full-name '@{u}')"
    if git fetch -q origin 2>/dev/null; then
        behind="$(git rev-list --count "HEAD..$upstream" 2>/dev/null || echo 0)"
        if [[ "${behind:-0}" -gt 0 ]]; then
            echo "      ERROR: checkout is ${behind} commit(s) behind ${upstream}." >&2
            echo "      Pull before building, or the binary will not contain the code you" >&2
            echo "      think it does. Set ALLOW_BEHIND=1 to build anyway." >&2
            [[ "${ALLOW_BEHIND:-0}" == "1" ]] || exit 1
        fi
    else
        echo "      WARNING: could not fetch ${upstream} — cannot confirm this checkout" >&2
        echo "      is current. If a pull failed earlier, you may be building stale code." >&2
    fi
fi
# LINK WITH gcc, NOT the Cray `cc` wrapper.
#
# Delta's site-wide `default` module set loads PrgEnv-gnu together with
# craype-accel-nvidia80, which tells the Cray compiler driver to target NVIDIA A100s. `cc`
# then composes its link line from pkg-config's `virtual:world`, which pulls in a CUDA
# runtime plus cray-mpich, cray-libsci, cray-dsmml and libfabric. eidolon uses NONE of
# those — it is pure Rust plus zlib-ng/libdeflate via cmake.
#
# After the RHEL/PE upgrade that world referenced cray-sdk-cudatoolkit-25.3_11.8, whose .pc
# file is absent, so `cc` refused to link ANYTHING — including a hello-world, and including
# every dependency build script. Every build failed with "Package 'cray-sdk-cudatoolkit...',
# required by 'virtual:world', not found" (2026-08-12).
#
# Confirmed working 2026-08-12. Unloading cudatoolkit ALONE was measured and still fails — it
# removes the package while leaving craype-accel-nvidia80, which is what demands it. Unloading
# craype-accel-nvidia80 instead is untested, and would not persist regardless: `default`
# reloads at every login. This export is immune to that. rustc already does the real linking
# with its own bundled lld (-fuse-ld=lld); `cc` was only ever supplying a driver front-end
# that insisted on resolving a GPU stack we do not use. gcc-native/13.2 is in the default
# set, so `gcc` is present.
#
# Override CARGO_LINKER=cc to go back to the wrapper if a future PE actually needs it.
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER="${CARGO_LINKER:-gcc}"
echo "      Linker: $CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER (see the comment above; Cray cc pulls in CUDA)"
cargo build --release 2>&1 | tail -5
echo "      Binary: $CARGO_TARGET_DIR/release/eidolon"
# Report the version actually produced. This is the line whose absence cost a run.
built_ver="$("$CARGO_TARGET_DIR/release/eidolon" --version 2>/dev/null || echo unknown)"
repo_ver="$(awk -F\' '/^version = /{print $2; exit}' Cargo.toml 2>/dev/null || echo unknown)"
echo "      Built:  ${built_ver}   (Cargo.toml says ${repo_ver})"
if [[ "$built_ver" != *"$repo_ver"* ]]; then
    echo "      ERROR: built binary reports '${built_ver}' but Cargo.toml says ${repo_ver}." >&2
    exit 1
fi

# ── 2. NEAT 4 conda env ─────────────────────────────────────────────────
setup_conda   # load the conda module + bootstrap `conda` (Delta: miniforge3-python)

# NEAT 4 + bcftools. NEAT shells out to `bcftools sort` at runtime (its
# environment.yml pins bcftools==1.22*), but the bioconda `neat` package does
# NOT pull it — without bcftools NEAT dies with FileNotFoundError: 'bcftools'.
# Used by benchmark.sbatch (throughput) and germline_e2e.sbatch (NEAT fidelity
# arm). If the env already exists, RETROFIT bcftools (a pre-existing neat4 from
# before this was added wouldn't otherwise get it).
if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "[2/4] Conda env '$CONDA_ENV_NAME' exists — ensuring bcftools is present..."
    conda run -n "$CONDA_ENV_NAME" bcftools --version >/dev/null 2>&1 \
        || conda install -y -n "$CONDA_ENV_NAME" -c conda-forge -c bioconda bcftools
else
    echo "[2/4] Creating conda env for NEAT 4 (bioconda)..."
    conda create -y -n "$CONDA_ENV_NAME" -c conda-forge -c bioconda neat bcftools
fi
# `neat` has no --version flag; verify the tools resolve instead.
echo "      neat4: $(conda run -n "$CONDA_ENV_NAME" which neat 2>/dev/null || echo 'neat MISSING')"
echo "             $(conda run -n "$CONDA_ENV_NAME" bcftools --version 2>/dev/null | head -1 || echo 'bcftools MISSING')"

# ── 3. bioinf conda env (bwa-mem2, gatk4, bcftools) ────────────────────
# samtools and htslib (tabix/bgzip) are available as Delta modules.
# bcftools, bwa-mem2, and gatk4 are not — install via bioconda.
# mimalloc/jemalloc are for the HPC tuning sweep (tune_sweep.sbatch tests them
# via LD_PRELOAD to isolate allocator lock contention on Delta's EPYC); numactl
# (a system tool on Delta) covers the NUMA lever. Harmless extras in this env.
if conda env list | grep -q "^${BIOINF_ENV_NAME} "; then
    echo "[3/4] Conda env '$BIOINF_ENV_NAME' exists — ensuring tuning allocators + seqkit present..."
    conda run -n "$BIOINF_ENV_NAME" sh -c 'ls "$CONDA_PREFIX"/lib/libmimalloc.so* >/dev/null 2>&1' \
        || conda install -y -n "$BIOINF_ENV_NAME" -c conda-forge mimalloc jemalloc
    conda run -n "$BIOINF_ENV_NAME" which seqkit >/dev/null 2>&1 \
        || conda install -y -n "$BIOINF_ENV_NAME" -c bioconda seqkit
    # fastp drives the adapter-readthrough QC/trim arm (germline_e2e MEASURE_REALISM/TRIM, #125).
    conda run -n "$BIOINF_ENV_NAME" which fastp >/dev/null 2>&1 \
        || conda install -y -n "$BIOINF_ENV_NAME" -c bioconda -c conda-forge fastp
else
    echo "[3/4] Creating bioinf conda env (bcftools, bwa-mem2, gatk4, seqkit, fastp, mimalloc, jemalloc)..."
    conda create -y -n "$BIOINF_ENV_NAME" \
        -c bioconda -c conda-forge \
        bcftools bwa-mem2 gatk4 seqkit fastp mimalloc jemalloc
    echo "      Tools installed:"
    # `|| true`: head -1 closes the pipe after one line, so the tool gets
    # SIGPIPE (exit 141) and `set -o pipefail` would otherwise abort setup.
    conda run -n "$BIOINF_ENV_NAME" bcftools  --version 2>&1 | head -1 || true
    conda run -n "$BIOINF_ENV_NAME" bwa-mem2  version   2>&1 | head -1 || true
    conda run -n "$BIOINF_ENV_NAME" gatk      --version 2>&1 | head -1 || true
fi

# ── 3b. Somatic-caller coverage envs (cross-validation) ──────────────────
# Independent second callers to confirm eidolon's cancer features aren't tuned to a
# single tool, plus a mutational-signature fidelity check. Each gets its own env
# (strelka2 is python2, sigprofiler python3, delly a static binary), so they don't
# clash with bioinf or each other. GATK somatic CNV needs NO new env (gatk4 is in
# bioinf); Manta + truvari already have their own envs (see sv_pipeline.sbatch).
#   delly       — somatic SV, 2nd caller vs Manta (esp. the BND recovery question)
#   strelka     — somatic SNV/indel, 2nd caller vs Mutect2
#   sigprofiler — SigProfilerAssignment, COSMIC-signature fidelity of the model
# NOTE: sigprofiler downloads a reference genome on first use
# (SigProfilerMatrixGenerator); the signature_check helper installs GRCh38 once.
#   varscan     — somatic SNV/indel, robust 2nd caller vs Mutect2 (Strelka2 2.9.10
#                 SIGSEGVs on Delta's stack, so it's the working cross-check)
for env_spec in \
    "delly:-c bioconda -c conda-forge delly" \
    "strelka:-c bioconda -c conda-forge strelka" \
    "varscan:-c bioconda -c conda-forge varscan" \
    "sigprofiler:-c bioconda -c conda-forge sigprofilerassignment"; do
    env_name="${env_spec%%:*}"; env_pkgs="${env_spec#*:}"
    if conda env list | grep -q "^${env_name} "; then
        echo "[3b/4] Conda env '$env_name' already present — skipping."
    else
        echo "[3b/4] Creating conda env '$env_name' ($env_pkgs)..."
        conda create -y -n "$env_name" $env_pkgs
    fi
done

# ── 4. hap.py Apptainer image ───────────────────────────────────────────
mkdir -p "$SIF_DIR"
HAPPY_SIF="$SIF_DIR/happy_v0.3.12.sif"
if [[ -f "$HAPPY_SIF" ]]; then
    echo "[4/4] hap.py SIF already present: $HAPPY_SIF"
else
    echo "[4/4] Pulling hap.py image (this takes a few minutes)..."
    apptainer pull "$HAPPY_SIF" docker://jmcdani20/hap.py:v0.3.12
    echo "      Saved to: $HAPPY_SIF"
fi

cat <<EOF

Setup complete.

Next steps:
  1. Place reference genome(s) and input data in \$SCRATCH
  2. ACCESS account is preset to bhrd-delta-cpu in the *.sbatch files.
     Override per-submit if needed:  sbatch --account=<other> scripts/delta/<job>.sbatch
  3. (optional) Point results at durable storage for the ACCESS final report —
     defaults to \${WORK:-\$HOME}/eidolon-access-results (NOT scratch, which is purged):
       export RESULTS_DIR=/projects/<account>/eidolon-access-results
  4. Submit jobs:
       sbatch scripts/delta/baseline_capture.sbatch   # simulator baseline
       sbatch scripts/delta/germline_e2e.sbatch       # SNP/indel recall
       sbatch scripts/delta/benchmark.sbatch          # NEAT-vs-eidolon throughput
       sbatch scripts/delta/cancer_pipeline.sbatch    # somatic SNV/indel recall
  5. Build the ACCESS final-report document once runs finish:
       bash scripts/delta/collect_report.sh           # -> \$RESULTS_DIR/REPORT.md

Each job persists its result artifacts + a provenance/resource manifest under
\$RESULTS_DIR (off scratch); collect_report.sh aggregates them (incl. total
core-hours ≈ SUs) into REPORT.md.
EOF