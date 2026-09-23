#!/usr/bin/env bash
# Measure homopolymer slippage from the READ BASES, not from the CIGAR.
#
# WHY THIS EXISTS (#692). On matched, untrimmed data, 41% of real reads overlapping a >= 21 bp
# homopolymer carry a soft clip of >= 20 bp; simulated reads manage 2.16%. That 19x gap
# accounts for half of `cand_per_mb`. The obvious way to ask what eidolon should generate
# there — read the `I`/`D` operations at those loci — cannot answer it: those are the events
# bwa-mem2 chose to GAP, and the ones that matter are the 41% it clipped instead. A
# measurement whose denominator excludes the phenomenon under study is not a measurement.
#
# This reads the run length out of each read's own sequence, anchored on flanking reference,
# and recovers the implied indel size behind each clip by finding where the clipped bases
# actually belong. `SEQ` in a BAM contains soft-clipped bases and is stored forward-strand, so
# the alignment is used only to FIND reads — never to interpret them.
#
# Usage:
#   scripts/delta/homopolymer_run_lengths.sh --bam <file> --reference <fa> \
#       --regions <bed> [--outdir <dir>]
#
# Env / flags:
#   PAD       400  reference fetched either side of an interval (flanks + clip search window)
#   FLANK     12   anchor length; ~4% of reads lose their anchor to sequencing error, and
#                  that loss is INDEPENDENT of whether the read clipped, which is the point
#   MIN_RUN   12   shortest homopolymer to measure
#   TAIL      15   bases of a clip's far end used to locate it (4^15 makes a chance hit ~0)
#   MIN_CLIP  20   clip size that counts, matching the realism panel's candidate threshold
#
# Compare two BAMs by running it twice against the same --regions. The real/simulated
# difference is the number; neither arm alone is one.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCI_AWK="${LOCI_AWK:-$HERE/homopolymer_loci.awk}"
MEASURE_AWK="${MEASURE_AWK:-$HERE/homopolymer_measure.awk}"
SUMMARISE_AWK="${SUMMARISE_AWK:-$HERE/homopolymer_summarise.awk}"

BAM=""; REFERENCE=""; REGIONS=""; OUTDIR=""
PAD="${PAD:-400}"; FLANK="${FLANK:-12}"; MIN_RUN="${MIN_RUN:-12}"
TAIL="${TAIL:-15}"; MIN_CLIP="${MIN_CLIP:-20}"

usage() {
    echo "usage: homopolymer_run_lengths.sh --bam <file> --reference <fa> --regions <bed> [--outdir <dir>]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)       BAM="${2:?--bam needs a value}"; shift 2;;
        --reference) REFERENCE="${2:?--reference needs a value}"; shift 2;;
        --regions)   REGIONS="${2:?--regions needs a value}"; shift 2;;
        --outdir)    OUTDIR="${2:?--outdir needs a value}"; shift 2;;
        -h|--help)   usage; exit 0;;
        *)           echo "unknown argument $1" >&2; usage; exit 1;;
    esac
done

[[ -n "$BAM" && -n "$REFERENCE" && -n "$REGIONS" ]] || { usage; exit 1; }
[[ -f "$BAM" ]]       || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: reference not found: $REFERENCE" >&2; exit 1; }
[[ -f "$REGIONS" ]]   || { echo "ERROR: regions BED not found: $REGIONS" >&2; exit 1; }
[[ -f "${REFERENCE}.fai" ]] || { echo "ERROR: reference is not indexed: ${REFERENCE}.fai (samtools faidx)" >&2; exit 1; }
# Assert the programs exist before anything reads their output. A missing awk file makes every
# downstream comparison compare empty against empty, which passes.
for f in "$LOCI_AWK" "$MEASURE_AWK" "$SUMMARISE_AWK"; do
    [[ -f "$f" ]] || { echo "ERROR: awk program not found: $f" >&2; exit 1; }
done
command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found. On Delta: module load samtools, or conda activate bioinf" >&2
    exit 1
}

OUTDIR="${OUTDIR:-$(mktemp -d)}"
mkdir -p "$OUTDIR"

echo "BAM:        $BAM"
echo "Reference:  $REFERENCE"
echo "Regions:    $REGIONS"
echo "Outdir:     $OUTDIR"
echo "Settings:   pad=$PAD flank=$FLANK min_run=$MIN_RUN tail=$TAIL min_clip=$MIN_CLIP"
echo ""

# ── 1. padded regions, clamped to the contigs ────────────────────────────────
#
# A region running off the end of a contig makes `samtools faidx` return a SHORTER sequence
# than asked for, which would silently shift every coordinate derived from it.
echo "[1/4] building padded windows ..."
awk -v pad="$PAD" -v OFS="" '
    NR == FNR { len[$1] = $2; next }
    /^#/ || NF < 3 { next }
    {
        if (!($1 in len)) { bad[$1]++; next }
        s = $2 + 1 - pad; if (s < 1) s = 1
        e = $3 + pad;     if (e > len[$1]) e = len[$1]
        if (e > s) print $1, ":", s, "-", e
    }
    END { for (c in bad) printf("  WARNING: %d interval(s) name contig %s, absent from the reference index\n", bad[c], c) > "/dev/stderr" }
' "${REFERENCE}.fai" "$REGIONS" > "$OUTDIR/windows.txt"
NWIN=$(wc -l < "$OUTDIR/windows.txt" | tr -d ' ')
echo "      $NWIN window(s)"
[[ "$NWIN" -gt 0 ]] || { echo "ERROR: no usable windows — do the BED's contig names match the reference?" >&2; exit 1; }

# ── 2. locus table ───────────────────────────────────────────────────────────
echo "[2/4] locating homopolymer runs ..."
samtools faidx "$REFERENCE" -r "$OUTDIR/windows.txt" 2>"$OUTDIR/faidx.err" \
  | awk -v pad="$PAD" -v flank="$FLANK" -v min_run="$MIN_RUN" -f "$LOCI_AWK" \
        > "$OUTDIR/loci.tsv" 2>"$OUTDIR/loci.log"
sed 's/^/      /' "$OUTDIR/loci.log"
NLOCI=$(wc -l < "$OUTDIR/loci.tsv" | tr -d ' ')
[[ "$NLOCI" -gt 0 ]] || {
    echo "ERROR: no loci survived. Zero loci and zero slippage look identical downstream, so" >&2
    echo "       this is a failure rather than a result. See the drop reasons above." >&2
    exit 1
}

# ── 3. per-locus measurement ─────────────────────────────────────────────────
echo "[3/4] measuring $NLOCI loci ..."
: > "$OUTDIR/measurements.tsv"
n=0
while IFS=$'\t' read -r contig rstart rlen rbase lflank rflank winstart winseq; do
    [[ -n "$contig" ]] || continue
    n=$((n + 1))
    rend=$((rstart + rlen - 1))
    # -F 0x904 matches what the realism panel counts: no unmapped, secondary or supplementary.
    # A supplementary record is hard-clipped, so its SEQ is not the read.
    samtools view -F 0x904 "$BAM" "$contig:$rstart-$rend" 2>/dev/null \
      | awk -v locus="$contig:$rstart" -v base="$rbase" -v lref="$rlen" \
            -v lflank="$lflank" -v rflank="$rflank" -v win="$winseq" \
            -v winstart="$winstart" -v tail="$TAIL" -v minclip="$MIN_CLIP" \
            -f "$MEASURE_AWK" >> "$OUTDIR/measurements.tsv"
done < "$OUTDIR/loci.tsv"
NMEAS=$(wc -l < "$OUTDIR/measurements.tsv" | tr -d ' ')
echo "      $NMEAS record(s) from $n locus queries"
[[ "$NMEAS" -gt 0 ]] || {
    echo "ERROR: no reads produced a record. Does the BAM cover these loci?" >&2
    exit 1
}

# ── 4. report ────────────────────────────────────────────────────────────────
echo "[4/4] summarising ..."
echo ""
awk -f "$SUMMARISE_AWK" "$OUTDIR/measurements.tsv" | tee "$OUTDIR/report.txt"
echo ""
echo "Raw records: $OUTDIR/measurements.tsv"
echo "Locus table: $OUTDIR/loci.tsv"
