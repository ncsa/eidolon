#!/usr/bin/env bash
# Put a pair of FASTQs into random order, with bounded memory, keeping mates in register (#870).
#
# WHY: eidolon writes reads in genome order, and aligning the genome-ordered CHM13 simulation
# needed 166 GB where the real library, with more pairs, fit in 128 GB (#861, #870). Aligning
# the same reads shuffled says whether order is the cause. seqkit shuffle cannot do this at
# genome scale: it holds every record in memory, and its low-memory two-pass mode is FASTA only.
#
# HOW: each pair becomes one tab-separated line (8 FASTQ lines), prefixed with a seeded random
# key; `sort` orders on the key using disk, not memory, beyond SORT_MEM; the lines are split back
# into two gzipped FASTQs. Deterministic for a given SEED and input.
#
# CHECKS, all hard failures:
#   - every input record has 4 lines on each mate (a pair that does not split into 8 fields);
#   - R1 and R2 names agree for every pair after the shuffle (same name, /1 against /2);
#   - pairs out == pairs in. The input passes through process substitutions, whose failure
#     set -e cannot see, so this count is what establishes that nothing was dropped.
#
# USAGE
#   IN_DIR=$SCRATCH/raw_pairs/<tag> OUT_DIR=$SCRATCH/raw_pairs/<tag>_shuf \
#     sbatch scripts/delta/shuffle_pairs.sbatch
# Reads IN_DIR/r1.fastq.gz and r2.fastq.gz and writes the same names into OUT_DIR, which is
# where stage_raw_pairs.sbatch looks for local reads (TAG=<tag>_shuf).

#SBATCH --job-name=eidolon-shuffle
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=shuffle-%j.log

set -euo pipefail

IN_DIR="${IN_DIR:?set IN_DIR to a directory holding r1.fastq.gz and r2.fastq.gz}"
OUT_DIR="${OUT_DIR:?set OUT_DIR}"
SEED="${SEED:-11}"
SORT_MEM="${SORT_MEM:-40G}"
# sort's temp files hold the whole library uncompressed (~170 GB for 256.6M pairs of 2x125), so
# they are compressed by default. SORT_COMPRESS= (empty) turns that off when scratch has room.
SORT_COMPRESS="${SORT_COMPRESS-gzip}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

in1="$IN_DIR/r1.fastq.gz"; in2="$IN_DIR/r2.fastq.gz"
for f in "$in1" "$in2"; do [[ -s "$f" ]] || { echo "FATAL: missing $f" >&2; exit 1; }; done
[[ "$(realpath -m "$IN_DIR")" != "$(realpath -m "$OUT_DIR")" ]] || { echo "FATAL: OUT_DIR must differ from IN_DIR" >&2; exit 1; }
mkdir -p "$OUT_DIR/sort_tmp"
out1="$OUT_DIR/r1.fastq.gz"; out2="$OUT_DIR/r2.fastq.gz"
if [[ -s "$out1" || -s "$out2" ]]; then
    echo "FATAL: $OUT_DIR already holds reads; delete them to reshuffle" >&2; exit 1
fi
echo "=== shuffle_pairs (#870): $IN_DIR -> $OUT_DIR, seed $SEED ==="

paste <(zcat "$in1" | paste - - - -) <(zcat "$in2" | paste - - - -) \
    | awk -F'\t' -v seed="$SEED" -v cnt="$OUT_DIR/pairs_in.txt" '
        BEGIN { srand(seed) }
        NF != 8 { printf "FATAL: input pair %d splits into %d fields, not 8\n", NR, NF > "/dev/stderr"; bad = 1; exit 1 }
        { printf "%.15f\t%s\n", rand(), $0 }
        END { if (!bad) print NR > cnt }' \
    | LC_ALL=C sort -t $'\t' -k1,1 -S "$SORT_MEM" -T "$OUT_DIR/sort_tmp" --parallel="$THREADS" \
        ${SORT_COMPRESS:+--compress-program="$SORT_COMPRESS"} \
    | awk -F'\t' -v o1="$out1.part" -v o2="$out2.part" -v cnt="$OUT_DIR/pairs_out.txt" '
        {
            n1 = $2; n2 = $6; sub(/\/1$/, "", n1); sub(/\/2$/, "", n2); sub(/ .*/, "", n1); sub(/ .*/, "", n2)
            if (n1 != n2) { printf "FATAL: mates out of register at output pair %d: %s vs %s\n", NR, $2, $6 > "/dev/stderr"; bad = 1; exit 1 }
            print $2 "\n" $3 "\n" $4 "\n" $5 | ("gzip -1 > " o1)
            print $6 "\n" $7 "\n" $8 "\n" $9 | ("gzip -1 > " o2)
        }
        END { if (!bad) { close("gzip -1 > " o1); close("gzip -1 > " o2); print NR > cnt } }'

n_in=$(cat "$OUT_DIR/pairs_in.txt" 2>/dev/null || echo missing)
n_out=$(cat "$OUT_DIR/pairs_out.txt" 2>/dev/null || echo missing)
echo "pairs in: $n_in, pairs out: $n_out"
[[ "$n_in" =~ ^[0-9]+$ && "$n_in" -gt 0 && "$n_in" == "$n_out" ]] \
    || { echo "FATAL: pair counts disagree or are missing; the shuffle did not complete" >&2; exit 1; }
mv "$out1.part" "$out1"
mv "$out2.part" "$out2"
rm -rf "$OUT_DIR/sort_tmp"
echo "DONE. Align with: TAG=$(basename "$OUT_DIR") sbatch scripts/delta/stage_raw_pairs.sbatch"
