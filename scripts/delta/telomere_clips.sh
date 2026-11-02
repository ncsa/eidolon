#!/usr/bin/env bash
# What are the clipped bases on reads in the telomeric arrays? (#861)
#
# telomere_ends.sbatch found real array reads clipped at 0.40 (median over 46 ends) against
# 0.02 simulated, with control windows agreeing at 1-2%. Candidate causes leave different marks:
#   - reads running off the assembly's end: the clip is on the side toward the chromosome tip,
#     and the alignment reaches the tip;
#   - variant repeats the assembly lacks: the clipped bases are telomere-like but not canonical
#     (TCAGGG, TGAGGG, TTGGGG, TTCGGG, CTAGGG and reverse complements);
#   - sequencing error: clipped bases at low quality, with no particular sequence;
#   - adapter: the TruSeq stem in the clip.
# For each library this reports, over primary non-duplicate reads overlapping the arrays with
# a soft clip of at least MIN_CLIP bp: which side the clip is on, whether the alignment touches
# the tip, what the clipped bases are made of, and their quality against the same reads'
# aligned bases. It prints a few clipped sequences to look at.
#
# Reads arrays.tsv and regions.tsv from a telomere_ends.sbatch run (ENDS_DIR). Small regions,
# so it runs in seconds; the #SBATCH lines are for consistency with the rest of scripts/delta.
#
# USAGE
#   ENDS_DIR=$SCRATCH/telomere_861_pcrfree REAL_BAM=... SIM_BAM=... \
#     bash scripts/delta/telomere_clips.sh

#SBATCH --job-name=eidolon-teloclips
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=teloclips-%j.log

set -euo pipefail

module load samtools/1.22-cce19.0.0 2>/dev/null || true
command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH" >&2; exit 1; }
SCRATCH="${SCRATCH:-$PWD}"
ENDS_DIR="${ENDS_DIR:-$SCRATCH/telomere_861_pcrfree}"
REAL_BAM="${REAL_BAM:-$SCRATCH/raw_pairs/chm13_pcrfree/chm13_pcrfree.bam}"
SIM_BAM="${SIM_BAM:-$SCRATCH/raw_pairs/chm13_pcrfree_sim/chm13_pcrfree_sim.bam}"
MIN_CLIP="${MIN_CLIP:-10}"
TIP_BP="${TIP_BP:-10}"       # an alignment within this many bp of the end "touches the tip"
N_EXAMPLES="${N_EXAMPLES:-8}"

for f in "$ENDS_DIR/arrays.tsv" "$ENDS_DIR/regions.tsv" "$REAL_BAM" "$SIM_BAM"; do
    [[ -s "$f" ]] || { echo "FATAL: missing $f" >&2; exit 1; }
done

# Array regions with their end and the contig length, from telomere_ends' own outputs.
mapfile -t regions < <(awk 'NR == FNR { if ($2 != "ctrl") len[$1] = $5; next }
    $3 > 0 { if ($2 == "p") print $1 ":1-" $3; else print $1 ":" len[$1] - $3 + 1 "-" len[$1] }' \
    "$ENDS_DIR/regions.tsv" "$ENDS_DIR/arrays.tsv")
[[ "${#regions[@]}" -gt 0 ]] || { echo "FATAL: no arrays in $ENDS_DIR/arrays.tsv" >&2; exit 1; }
lens=$(awk '$2 != "ctrl" { print $1, $5 }' "$ENDS_DIR/regions.tsv" | sort -u)
echo "=== telomere_clips (#861): ${#regions[@]} arrays, clips >= $MIN_CLIP bp ==="

for lib in real sim; do
    bam_var="${lib^^}_BAM"
    # Collect first: a samtools failure must stop the script, not read as "no clipped reads".
    reads=$(samtools view -M -F 0xD04 "${!bam_var}" "${regions[@]}")
    awk -v lib="$lib" -v min="$MIN_CLIP" -v tip="$TIP_BP" -v nex="$N_EXAMPLES" '
        FNR == NR { len[$1] = $2; next }
        function qmean(q,   i, s) { s = 0; for (i = 1; i <= length(q); i++) s += ord[substr(q, i, 1)] - 33; return length(q) ? s / length(q) : 0 }
        function rotations(h, set,   i) { for (i = 0; i < 6; i++) set[substr(h, i + 1) substr(h, 1, i)] = 1 }
        BEGIN {
            for (i = 33; i < 127; i++) ord[sprintf("%c", i)] = i
            # Every 6-mer window of a clean repeat is a rotation of its hexamer, whatever the
            # phase, so windows are classified by rotation sets rather than counting whole
            # hexamers (which a clip cut mid-repeat would undercount).
            rotations("TTAGGG", CAN); rotations("CCCTAA", CAN)
            nv = split("TCAGGG TGAGGG TTGGGG TTCGGG CTAGGG CCCTGA CCCTCA CCCCAA CCCGAA CCCTAG", v, " ")
            for (i = 1; i <= nv; i++) rotations(v[i], VAR)
        }
        NF {
            n++
            c = $6; lead = (match(c, /^[0-9]+S/) ? substr(c, 1, RLENGTH - 1) + 0 : 0)
            trail = (match(c, /[0-9]+S$/) ? substr(c, RSTART, RLENGTH - 1) + 0 : 0)
            # Reference span from M/D/N/=/X ops.
            span = 0; x = c
            while (match(x, /[0-9]+[MDN=X]/)) { span += substr(x, RSTART, RLENGTH - 1); x = substr(x, RSTART + RLENGTH) }
            aln_end = $4 + span - 1
            pend = ($4 < len[$3] / 2)
            # Toward the tip: a left clip at a p end, a right clip at a q end.
            out = pend ? lead : trail; inn = pend ? trail : lead
            if (out < min && inn < min) next
            nc++
            touches = pend ? ($4 <= 1 + tip) : (aln_end >= len[$3] - tip)
            if (out >= min) { n_out++; if (touches) n_out_tip++ }
            if (inn >= min) n_in++
            seq = $10; q = $11; L = length(seq)
            clip = ""; cq = ""
            if (lead >= min) { clip = clip substr(seq, 1, lead); cq = cq substr(q, 1, lead) }
            if (trail >= min) { clip = clip substr(seq, L - trail + 1); cq = cq substr(q, L - trail + 1) }
            alq = qmean(substr(q, lead + 1, L - lead - trail))
            cb += length(clip)
            # Windows within each clip separately, so none straddles a left and a right clip.
            nseg = split((lead >= min ? substr(seq, 1, lead) : "") " " (trail >= min ? substr(seq, L - trail + 1) : ""), seg, " ")
            for (k = 1; k <= nseg; k++)
                for (i = 1; i + 5 <= length(seg[k]); i++) {
                    w = substr(seg[k], i, 6); nw++
                    if (w in CAN) can++; else if (w in VAR) var++
                }
            if (index(clip, "AGATCGGAAGAGC") || index(clip, "GCTCTTCCGATCT")) adapt++
            sq += qmean(cq) * length(cq); sa += alq
            if (ne < nex) ex[++ne] = sprintf("  %s:%d %s %s clip=%s", $3, $4, (pend ? "p" : "q"), $6, clip)
        }
        END {
            printf "\n%s: %d reads in the arrays; %d (%.3f) with a clip >= %d bp\n", lib, n, nc, n ? nc / n : 0, min
            if (!nc) exit
            printf "  clip toward the tip: %d (%.3f of clipped), of which the alignment reaches the tip: %d (%.3f)\n", n_out, n_out / nc, n_out_tip, n_out ? n_out_tip / n_out : 0
            printf "  clip away from the tip: %d (%.3f of clipped)\n", n_in, n_in / nc
            printf "  clipped bases: %d in %d 6-mer windows: canonical repeat %.3f, variant repeats %.3f, other %.3f\n", cb, nw, can / nw, var / nw, 1 - (can + var) / nw
            printf "  reads with the TruSeq stem in a clip: %d\n", adapt + 0
            printf "  mean base quality: clipped %.1f, aligned part of the same reads %.1f\n", sq / cb, sa / nc
            print "  examples:"; for (i = 1; i <= ne; i++) print ex[i]
        }' <(echo "$lens") - <<< "$reads"
done
