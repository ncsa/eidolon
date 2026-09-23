#!/usr/bin/env bash
# Tests for submit.sh. No SLURM here, so `sbatch` is stubbed and the assertions are on the
# command line it is handed -- which is the whole content of this script.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT:-$HERE/../submit.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$WORK/bin"
cat > "$WORK/bin/sbatch" <<'STUB'
#!/usr/bin/env bash
# Record the argv we were given, then answer like `sbatch --parsable`.
printf '%s\n' "$@" > "$SBATCH_ARGV"
echo "4242;delta"
STUB
chmod +x "$WORK/bin/sbatch"
export PATH="$WORK/bin:$PATH"
export SBATCH_ARGV="$WORK/argv.txt"
export SCRATCH="$WORK"

if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label from to; do
        [[ -n "$label" ]] || continue
        cp "$SCRIPT" "$WORK/mutant.sh"
        FROM="$from" TO="$to" perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant.sh"
        if cmp -s "$SCRIPT" "$WORK/mutant.sh"; then
            printf '  ERROR    %-48s mutation did not apply\n' "$label"; survived=$((survived+1)); continue
        fi
        if SCRIPT="$WORK/mutant.sh" bash "$0" >/dev/null 2>&1; then
            printf '  SURVIVED %-48s <- nothing caught this\n' "$label"; survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTS'
account default changed@ACCOUNT="${ACCOUNT:-bhrd-delta-cpu}"@ACCOUNT="${ACCOUNT:-wrong-account}"
env override ignored for time@TIME="${TIME:-04:00:00}"@TIME="04:00:00"
--parsable dropped, job id unparseable@out="$(sbatch --parsable \@out="$(sbatch \
stderr split into a separate file@    -J "$NAME" -o "$LOG" \@    -J "$NAME" -o "$LOG" -e "$LOG.err" \
command not passed through@--wrap="$CMD")"@--wrap="true")"
no-argument usage accepted@[[ $# -ge 1 && -n "${1:-}" ]] || {@[[ 1 ]] || {
MUTS
    echo
    [[ "$survived" -eq 0 ]] && { echo "all mutations caught"; exit 0; } || { echo "$survived survived"; exit 1; }
fi

PASS=0; FAIL=0
MIN_ASSERTIONS=15
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
# argv is one token per line. An exact line match suits standalone tokens; a flag/value PAIR
# has to be matched against the joined line, since matching the value alone is ambiguous --
# "1" appears as the value of -N, -n AND -c, so asserting it proves nothing about which.
argv_has()  { grep -qxF -- "$3" "$SBATCH_ARGV" && ok "$1" || bad "$1" "argv line: $3" "$(tr '\n' ' ' < "$SBATCH_ARGV")"; }
argv_pair() { local j; j="$(tr '\n' ' ' < "$SBATCH_ARGV")"
              case "$j" in *"$3"*) ok "$1";; *) bad "$1" "argv contains: $3" "$j";; esac; }

echo "=== defaults ==="
OUT="$(bash "$SCRIPT" 'echo hello' 2>&1)"
argv_has "account defaults to bhrd-delta-cpu" "" "bhrd-delta-cpu"
argv_has "partition defaults to cpu"          "" "cpu"
argv_pair "one cpu by default"                "" "-c 1 "
argv_has "memory defaults to 4G"              "" "--mem=4G"
argv_has "walltime defaults to 4h"            "" "04:00:00"
argv_has "the command is passed verbatim"     "" "--wrap=echo hello"
argv_has "--parsable so the job id can be read" "" "--parsable"
# -o alone must carry both streams; a separate -e would split stderr away from it.
if grep -qxF -- "-e" "$SBATCH_ARGV"; then
    bad "stderr is NOT split into a separate file" "no -e flag" "$(tr '\n' ' ' < "$SBATCH_ARGV")"
else
    ok "stderr is NOT split into a separate file"
fi
has "reports the parsed job id" "$OUT" "submitted job 4242"
# The point of the reminder: a killed job's log looks like a short successful one.
has "tells you to check the exit state" "$OUT" "sacct -j 4242"
has "and resolves %j in the log path it prints" "$OUT" "eidolon-oneliner-4242.log"

echo "=== overrides ==="
OUT2="$(NAME=bam2fq TIME=08:00:00 MEM=16G CPUS=8 bash "$SCRIPT" 'samtools fastq in.bam' 2>&1)"
argv_has "TIME override honoured" "" "08:00:00"
argv_has "MEM override honoured"  "" "--mem=16G"
argv_pair "CPUS override honoured" "" "-c 8 "
argv_has "NAME override reaches the job name" "" "bam2fq"

echo "=== refuses to submit nothing ==="
if bash "$SCRIPT" >/dev/null 2>&1; then
    bad "no argument is an error" "non-zero exit" "exit 0"
else
    ok "no argument is an error"
fi

echo
printf 'assertions: %d passed, %d failed\n' "$PASS" "$FAIL"
if [[ $((PASS + FAIL)) -lt $MIN_ASSERTIONS ]]; then
    echo "FAIL: only $((PASS + FAIL)) assertions ran, expected at least $MIN_ASSERTIONS" >&2
    exit 1
fi
[[ "$FAIL" -eq 0 ]] || exit 1
