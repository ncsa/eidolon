#!/usr/bin/env bash
# Shared reporting/persistence helpers for the Delta SLURM jobs.
#
# Source this from each *.sbatch right after REPO_ROOT is set:
#     source "$REPO_ROOT/scripts/delta/lib_report.sh"
# (NOT via $(dirname "$0") — under sbatch $0 is a spooled copy). Sourcing also
# resolves the Delta filesystem paths below, so do it before any $SCRATCH use.
# Call archive_run() once near the end of the job.
#
# WHY THIS EXISTS: the jobs do their heavy work in $SCRATCH, which Delta PURGES
# (files untouched ~30 days are deleted). For the ACCESS final report we need
# outcomes preserved on DURABLE storage, consolidated, with per-run provenance
# and resource (core-hour / SU) accounting. archive_run() copies the small
# report artifacts off scratch and writes a run_manifest.tsv; collect_report.sh
# later aggregates all runs into a single REPORT.md.

# ── Delta filesystem paths ───────────────────────────────────────────────
# Delta does NOT export $SCRATCH/$WORK to jobs, and these scripts run under
# `set -u`, so resolve them here (this file is sourced before any $SCRATCH use).
# Layout (confirmed): scratch = /scratch/<project>/$USER (purged ~30 days);
# durable project space = /projects/<project>/$USER.
# Override by exporting SCRATCH / RESULTS_DIR, or by setting ACCESS_PROJECT.
: "${USER:=$(id -un)}"
: "${ACCESS_PROJECT:=bhrd}"
: "${SCRATCH:=/scratch/$ACCESS_PROJECT/$USER}"
if [[ -n "${SLURM_JOB_ID:-}" && ! -d "$SCRATCH" ]]; then
    echo "ERROR: scratch dir '$SCRATCH' not found — export SCRATCH=... or set ACCESS_PROJECT=..." >&2
    exit 1
fi

# Durable results root for the ACCESS final report — the project filesystem
# (persists), NEVER scratch. Override with RESULTS_DIR=...
RESULTS_DIR="${RESULTS_DIR:-/projects/$ACCESS_PROJECT/$USER/eidolon-access-results}"

# ── conda activation (Delta) ─────────────────────────────────────────────
# On Delta the miniforge module puts conda's legacy `activate` on PATH but NOT
# `conda` itself, so `source activate` bootstraps the base env + shell
# functions. Override the module with CONDA_MODULE=...; no-op if conda is
# already available. The `set +u` is required: conda's activate scripts read
# unbound vars (e.g. $PS1) and would abort under our `set -u`.
setup_conda() {
    # 1. Make the `conda` command available (module + legacy `source activate`).
    if ! command -v conda >/dev/null 2>&1; then
        module load "${CONDA_MODULE:-miniforge3-python}" 2>/dev/null || true
        set +u
        command -v conda >/dev/null 2>&1 || source activate 2>/dev/null || true
        set -u
    fi
    command -v conda >/dev/null 2>&1 || {
        echo "ERROR: conda unavailable after 'module load ${CONDA_MODULE:-miniforge3-python}' + 'source activate'." >&2
        echo "       Set CONDA_MODULE=<your module> or initialize conda before running." >&2
        return 1
    }
    # 2. Install the activate hook. `module load`/`source activate` expose the
    #    `conda` command but NOT the shell hook, so `conda activate <env>` errors
    #    "Run 'conda init' before 'conda activate'" (seen on Delta — the
    #    `|| source activate` fallback in conda_activate masked it). Sourcing
    #    conda.sh is the scriptable equivalent of `conda init` and makes
    #    `conda activate` work cleanly.
    local base; base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "$base" && -f "$base/etc/profile.d/conda.sh" ]]; then
        set +u; source "$base/etc/profile.d/conda.sh"; set -u
    fi
}

# set-u-safe environment activation (same unbound-var reason as above).
conda_activate() {
    set +u
    conda activate "$1" || source activate "$1"
    local rc=$?
    set -u
    return "$rc"
}

# Emit five space-separated resource values for the current SLURM job:
#   elapsed_s alloc_cpus alloc_nodes maxrss_kb core_hours
# Best-effort: prints zeros when not under SLURM or sacct is unavailable, so the
# caller can always `read` exactly five fields.
_resource_values() {
    local jid="${SLURM_JOB_ID:-}"
    if [[ -z "$jid" ]] || ! command -v sacct >/dev/null 2>&1; then
        echo "0 0 0 0 0.0"; return
    fi
    # Allocation line carries ElapsedRaw/AllocCPUS/AllocNodes.
    local elapsed cpus nodes
    read -r elapsed cpus nodes < <(
        sacct -j "$jid" --noheader --parsable2 \
            --format=ElapsedRaw,AllocCPUS,AllocNodes 2>/dev/null \
        | head -1 | awk -F'|' '{print $1+0, $2+0, $3+0}')
    : "${elapsed:=0}"; : "${cpus:=0}"; : "${nodes:=0}"
    # MaxRSS (KB). archive_run is called at job END but while the job is still
    # RUNNING, so sacct hasn't finalized the .batch step's MaxRSS yet -> it comes
    # back empty (every Delta run reported 0). sstat reports a *running* step's
    # live RSS, so try it first (.batch then .0), and fall back to sacct. awk
    # normalizes the K/M/G/T suffix to KB and takes the max across rows.
    local max_rss_awk='
        { v=$1; if (v=="" || v=="-") next
          u=substr(v,length(v),1); n=v
          if      (u=="K") n=substr(v,1,length(v)-1)
          else if (u=="M") n=substr(v,1,length(v)-1)*1024
          else if (u=="G") n=substr(v,1,length(v)-1)*1048576
          else if (u=="T") n=substr(v,1,length(v)-1)*1073741824
          if (n+0>max) max=n+0 }
        END { printf "%.0f", max+0 }'
    local rss_kb=0
    if command -v sstat >/dev/null 2>&1; then
        rss_kb=$(sstat -a -j "$jid" --noheader --parsable2 --format=MaxRSS 2>/dev/null | awk "$max_rss_awk")
    fi
    if [[ -z "$rss_kb" || "$rss_kb" == 0 ]]; then
        rss_kb=$(sacct -j "$jid" --noheader --parsable2 --format=MaxRSS 2>/dev/null | awk "$max_rss_awk")
    fi
    # Core-hours = elapsed_hours * allocated CPUs (≈ Delta CPU SUs).
    local core_hours
    core_hours=$(awk -v e="$elapsed" -v c="$cpus" 'BEGIN{printf "%.3f", (e/3600.0)*c}')
    echo "${elapsed:-0} ${cpus:-0} ${nodes:-0} ${rss_kb:-0} ${core_hours:-0.0}"
}

# archive_run <kind> <source_outdir> [artifact-file ...]
# Copies the named (small) report artifacts off scratch into the durable
# results dir and writes run_manifest.tsv with provenance + resource usage.
# Large data (FASTQ/BAM) is intentionally left in scratch — only report-
# relevant files (csv/tsv/json/logs) are persisted.
archive_run() {
    local kind="$1" outdir="$2"; shift 2
    local jid="${SLURM_JOB_ID:-local}"
    local dest="$RESULTS_DIR/$kind/job_${jid}"
    mkdir -p "$dest"

    local f
    for f in "$@"; do
        [[ -e "$outdir/$f" ]] && cp -f "$outdir/$f" "$dest/" 2>/dev/null || true
    done
    # SLURM stdout/err (written to the submit dir as <jobname>_<jobid>.out/.err).
    local base="${SLURM_JOB_NAME:-job}_${jid}"
    [[ -f "${base}.out" ]] && cp -f "${base}.out" "$dest/" 2>/dev/null || true
    [[ -f "${base}.err" ]] && cp -f "${base}.err" "$dest/" 2>/dev/null || true

    local ver git_desc
    ver="$("${EIDOLON_BIN:-eidolon}" --version 2>/dev/null || echo unknown)"
    git_desc="$(git -C "${REPO_ROOT:-.}" describe --tags --always --dirty 2>/dev/null || echo unknown)"
    local elapsed_s alloc_cpus alloc_nodes maxrss_kb core_hours
    read -r elapsed_s alloc_cpus alloc_nodes maxrss_kb core_hours < <(_resource_values)

    # Tab-separated key/value manifest — robust to parse with awk, no jq needed.
    printf '%s\t%s\n' \
        kind          "$kind" \
        date_utc      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        slurm_job_id  "$jid" \
        eidolon_version "$ver" \
        git           "$git_desc" \
        reference     "$(basename "${REFERENCE:-NA}")" \
        elapsed_s     "${elapsed_s:-0}" \
        alloc_cpus    "${alloc_cpus:-0}" \
        alloc_nodes   "${alloc_nodes:-0}" \
        maxrss_kb     "${maxrss_kb:-0}" \
        core_hours    "${core_hours:-0.0}" \
        artifacts     "$*" \
        > "$dest/run_manifest.tsv"

    echo "[archive] $kind -> $dest  (core_hours=${core_hours:-0.0}, files: $*)"
}

# ── OUTDIR provenance guard ─────────────────────────────────────────────────
# Reusing an OUTDIR across a rebuild silently mixes artifacts from two code
# versions: fresh caller VCFs scored against a truth VCF written by the OLD
# binary. The per-script "simulation regenerated — invalidating stale downstream
# outputs" logic only fires when the simulation RE-RUNS, so nothing catches the
# inverse (new binary, skipped simulation). That is a wrong answer with no
# warning, which is the failure mode this repo keeps hitting.
#
# check_outdir_version sets EIDOLON_VERSION and FORCE_RESIM (0/1).
# Call it before the simulation stage; call stamp_outdir_version after a
# successful simulation.
check_outdir_version() {
    local bin="$1" outdir="$2"
    EIDOLON_VERSION="$("$bin" --version 2>/dev/null | tr -d '\r\n')"
    if [[ -z "$EIDOLON_VERSION" ]]; then
        echo "ERROR: could not read '$bin --version' — is it built?" >&2
        return 1
    fi
    FORCE_RESIM=0
    local stamp="$outdir/.eidolon_version"
    if [[ ! -f "$stamp" ]]; then
        # No stamp: either a fresh OUTDIR, or one predating this guard. If it
        # already holds simulation outputs we cannot know what produced them —
        # say so rather than assume, but don't force an expensive redo of a
        # legacy directory the user may know is fine.
        if compgen -G "$outdir"/*_merged_truth.vcf.gz >/dev/null 2>&1; then
            echo "WARNING: $outdir has simulation outputs but no .eidolon_version stamp," >&2
            echo "         so their provenance is unknown. If they predate your last" >&2
            echo "         rebuild, the truth VCF and the caller inputs come from" >&2
            echo "         different code. Use a fresh OUTDIR if in doubt." >&2
        fi
        return 0
    fi
    local prev
    prev="$(tr -d '\r\n' < "$stamp")"
    [[ "$prev" == "$EIDOLON_VERSION" ]] && return 0
    if [[ "${ALLOW_VERSION_MISMATCH:-0}" == "1" ]]; then
        echo "WARNING: $outdir was built by '$prev' but this binary is" >&2
        echo "         '$EIDOLON_VERSION'; ALLOW_VERSION_MISMATCH=1 — reusing anyway." >&2
        return 0
    fi
    echo "  OUTDIR was built by '$prev' but this binary is '$EIDOLON_VERSION'."
    echo "  Regenerating the simulation so the truth VCF and the caller inputs come"
    echo "  from the same code. Set ALLOW_VERSION_MISMATCH=1 to reuse as-is."
    FORCE_RESIM=1
    return 0
}

# Gate emitted artifacts through `eidolon validate` before anything consumes them.
#
# This replaces accreting more ad-hoc greps. The reason it exists: the malformed
# `AF=AF=0.3000` INFO value sailed past an `nsom > 0` guard because that guard counted
# RECORDS, not content — and bcftools would never have complained either, since it
# silently converts a type-mismatched value to `.`. A record count cannot see that; a
# validator that knows the declared type can.
#
# ERROR findings are fatal: producing numbers from an artifact a downstream tool will
# reject, or has silently emptied, is precisely the failure this whole harness keeps
# hitting. WARNINGs are printed and do not stop the run, because by construction nothing
# downstream rejects them.
#
# Skips (loudly) if the binary predates the subcommand, so an older checkout still runs.
validate_artifacts() {  # <eidolon-bin> <label> <file>...
    local bin="$1" label="$2"; shift 2
    [[ $# -gt 0 ]] || return 0
    if [[ ! -x "$bin" ]]; then
        echo "  WARNING: $bin is not executable — skipping $label validation." >&2
        return 0
    fi
    if ! "$bin" validate --help >/dev/null 2>&1; then
        echo "  WARNING: this eidolon build has no \`validate\` subcommand — skipping" >&2
        echo "    $label validation. Rebuild to enable the artifact gate." >&2
        return 0
    fi
    local present=()
    local f
    for f in "$@"; do
        [[ -f "$f" ]] && present+=("$f")
    done
    if [[ ${#present[@]} -eq 0 ]]; then
        echo "  WARNING: none of the $label artifacts exist yet — nothing validated." >&2
        return 0
    fi
    echo "  Validating $label (${#present[@]} artifact(s))..."
    local rc=0
    "$bin" validate "${present[@]}" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        echo "ERROR: $label failed validation (see the findings above). Each one names the" >&2
        echo "  tool and operation that will reject it. Refusing to compute results from an" >&2
        echo "  artifact a consumer will not accept." >&2
        return 1
    fi
    return 0
}

stamp_outdir_version() {
    printf '%s\n' "$EIDOLON_VERSION" > "$1/.eidolon_version"
}

# ── OUTDIR concurrency lock ─────────────────────────────────────────────────
# Two jobs pointed at one OUTDIR silently interleave their writes — the same
# FASTQ, BAM, truth-VCF and caller-output paths — and still produce numbers,
# derived from corrupted intermediates. Observed with SLURM jobs 20635663 and
# 20636020, which both simulated into $SCRATCH/t6_sv_v3 concurrently; neither
# reported anything wrong.
#
# The lock is held for the life of the job because fd 9 stays open. flock is
# released by the kernel when the holder exits, so a killed job leaves NO stale
# lock to clean up — the lockfile persisting on disk is harmless.
lock_outdir() {
    local outdir="$1" lockfile="$1/.lock"
    mkdir -p "$outdir"
    # <> not > : opening for write would TRUNCATE the holder's identity before we
    # even try to acquire, so a rejected job couldn't report who has it.
    if ! exec 9<>"$lockfile"; then
        echo "ERROR: cannot open lock file $lockfile" >&2
        return 1
    fi
    if ! flock -n 9; then
        local holder
        holder="$(tr -d '\r' < "$lockfile" 2>/dev/null | head -1)"
        echo "ERROR: another job is already running in $outdir" >&2
        echo "       held by: ${holder:-<unknown>}" >&2
        echo "  Concurrent runs sharing an OUTDIR interleave their writes and yield" >&2
        echo "  results computed from corrupted intermediates — with no warning." >&2
        echo "  Use a different OUTDIR, or wait for that job to finish." >&2
        return 1
    fi
    printf 'SLURM job %s on %s since %s\n' \
        "${SLURM_JOB_ID:-manual}" "$(hostname)" "$(date -Is 2>/dev/null || date)" \
        > "$lockfile"
    return 0
}

# Serialize the shared bwa-mem2 index build. Two jobs that both see no index will
# both run `bwa-mem2 index` against the SAME path and interleave writes, poisoning
# it for both — the 0-byte .bwt.2bit.64 state the call sites already guard against
# after the fact. Blocking lock: the second job waits, then finds the index present
# and skips. Degrades to the unlocked check if the reference dir isn't writable.
index_reference_locked() {
    local ref="$1"
    # BRACES ARE LOAD-BEARING. `exec` with redirections and no command applies them to the
    # CURRENT SHELL, permanently — so the bare `exec 8<>lock 2>/dev/null` this replaces sent
    # every subsequent >&2 in the job to /dev/null. Job 21532276 lost all three read-evidence
    # gates' diagnostics and its whole failure gate that way: .err stopped at 127 bytes the
    # moment this function ran (immediately before bwa), while the job continued for another
    # ten minutes and exited 1. The three failing loci it named existed nowhere on disk.
    # Grouping scopes the 2>/dev/null to the group while fd 8 still persists after it.
    if { exec 8<>"${ref}.index.lock"; } 2>/dev/null; then
        flock 8
        if [[ ! -s "${ref}.bwt.2bit.64" ]] && [[ ! -s "${ref}.bwt" ]]; then
            echo "  Indexing reference for BWA-MEM2 (one-time, lock held)..."
            bwa-mem2 index "$ref" 2>&1 | tail -3
        fi
        flock -u 8
        exec 8>&-
    else
        echo "  NOTE: cannot create ${ref}.index.lock (read-only dir?) — indexing" >&2
        echo "        unserialized; do not run concurrent jobs against a fresh reference." >&2
        if [[ ! -s "${ref}.bwt.2bit.64" ]] && [[ ! -s "${ref}.bwt" ]]; then
            echo "  Indexing reference for BWA-MEM2 (one-time)..."
            bwa-mem2 index "$ref" 2>&1 | tail -3
        fi
    fi
    [[ -f "${ref}.fai" ]] || samtools faidx "$ref"
}
