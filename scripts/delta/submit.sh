#!/usr/bin/env bash
# Wrap a one-liner in sbatch with this project's defaults.
#
# WHY THIS EXISTS: the submit rule in CLAUDE.md says anything that might take more than five
# minutes gets submitted rather than run on a login node, which caps at 30 minutes. A rule is
# only followed if it is cheap, and wrapping by hand means recalling eight flags -- account,
# partition, nodes, tasks, memory, time, and a log path -- to run one `samtools` command. That
# friction is what gets things run on the head node instead.
#
# USAGE
#   scripts/delta/submit.sh 'samtools fastq -F 0x900 in.bam | gzip -c > out.fastq.gz'
#   NAME=bam2fq TIME=08:00:00 MEM=16G scripts/delta/submit.sh '...'
#   CPUS=8 scripts/delta/submit.sh 'bwa-mem2 mem -t 8 ...'
#
# The command is passed as ONE argument, quoted. It runs under bash in the submitting
# directory, so relative paths behave as they do interactively.
#
# For a SCRIPT, do not use this -- give the script its own `#SBATCH` directives and `sbatch`
# it directly. This is for one-liners only.

set -euo pipefail

[[ $# -ge 1 && -n "${1:-}" ]] || {
    sed -n '2,20p' "$0" | sed 's/^# \?//' >&2
    exit 2
}

CMD="$*"
NAME="${NAME:-eidolon-oneliner}"
ACCOUNT="${ACCOUNT:-bhrd-delta-cpu}"
PARTITION="${PARTITION:-cpu}"
CPUS="${CPUS:-1}"
MEM="${MEM:-4G}"
TIME="${TIME:-04:00:00}"
# Deterministic and findable. %j is the job id, so concurrent submissions do not collide.
LOGDIR="${LOGDIR:-${SCRATCH:-$PWD}}"
LOG="$LOGDIR/${NAME}-%j.log"

echo "submitting: $CMD"
echo "  account=$ACCOUNT partition=$PARTITION cpus=$CPUS mem=$MEM time=$TIME"
echo "  log: $LOG"

# -o alone captures BOTH streams: SLURM merges stderr into --output unless --error is given
# separately. No redirect inside the --wrap.
out="$(sbatch --parsable \
    -A "$ACCOUNT" -p "$PARTITION" \
    -N 1 -n 1 -c "$CPUS" --mem="$MEM" -t "$TIME" \
    -J "$NAME" -o "$LOG" \
    --wrap="$CMD")"
jobid="${out%%;*}"

echo "submitted job $jobid"
echo
# Capturing output is not enough. A killed job leaves partial output that reads like a short
# successful run -- which is how a 30-minute kill went unnoticed until the file turned out to
# be empty. The exit state is what says otherwise.
echo "when it finishes, CHECK THE EXIT STATE, not just the log:"
echo "  sacct -j $jobid --format=JobID,State,ExitCode,Elapsed,MaxRSS"
echo "  cat ${LOG/\%j/$jobid}"
