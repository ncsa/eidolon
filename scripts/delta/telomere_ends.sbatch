#!/usr/bin/env bash
# Compare real and simulated CHM13 libraries at every chromosome end (#861).
#
# Both BAMs come from stage_raw_pairs.sbatch, aligned to the same reference by the same
# commands; see telomere_sim.sbatch for the chain. For each end of each contig this measures,
# per library:
#   - the telomeric array's length in the REFERENCE: counting 100 bp windows inward from the
#     end, the last window holding >= TEL_MIN copies of TTAGGG or CCCTAA before three windows
#     in a row that do not. Arrays carry degenerate stretches (CHM13 chr21 q: 8 copies at
#     3.7-3.8 kb, then full repeat again to 4.3 kb), so stopping at the first weak window
#     undercounts;
#   - depth in BIN-bp bins over the terminal END_BP, at all MAPQ and at MAPQ >= 20
#     (samtools depth, which skips unmapped, secondary, QC-fail and duplicate reads);
#   - over the array and over the rest of the terminal END_BP: reads, the MAPQ-0 fraction and
#     the clipped fraction (primary, non-duplicate records only);
#   - and, from the FASTQs with no mapping involved, reads carrying 4 consecutive TTAGGG or
#     CCCTAA per million reads.
# Depth is reported relative to each library's own median depth over interior control
# windows (CTRL_FRAC of the way along each contig), so the comparison does not depend on the
# two libraries having exactly the same coverage.
#
# MUST NOT FIRE. The control windows are measured the same way, and their MAPQ-0 and clip
# fractions should agree between libraries about as well as the realism panel says the genome
# does. A difference at the ends means something only next to that baseline.
#
# USAGE (after stage_raw_pairs has aligned both libraries):
#   sbatch scripts/delta/telomere_ends.sbatch
# Inputs and knobs are env vars, defaults below. Locally, the same file runs inline.

#SBATCH --job-name=eidolon-teloends
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=teloends-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
if [[ "${NO_ARCHIVE:-0}" != 1 ]]; then
    source "$REPO/scripts/delta/lib_report.sh"
    module load samtools/1.22-cce19.0.0 2>/dev/null || true
fi
SCRATCH="${SCRATCH:-$PWD}"
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/chm13/hs1.fa}"
REAL_BAM="${REAL_BAM:-$SCRATCH/raw_pairs/chm13_novaseq/chm13_novaseq.bam}"
SIM_BAM="${SIM_BAM:-$SCRATCH/raw_pairs/chm13_sim/chm13_sim.bam}"
OUT="${OUT:-$SCRATCH/telomere_861}"
END_BP="${END_BP:-20000}"
BIN="${BIN:-100}"
TEL_MIN="${TEL_MIN:-10}"          # of a possible 16 hexamers in a 100 bp window
MIN_CONTIG="${MIN_CONTIG:-1000000}"
EXCLUDE="${EXCLUDE:-chrY,chrM}"   # CHM13 is 46,XX; chrM is not a chromosome end
CTRL_FRAC="${CTRL_FRAC:-0.75}"

command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH" >&2; exit 1; }
for f in "$REFERENCE.fai" "$REAL_BAM" "$SIM_BAM"; do
    [[ -s "$f" ]] || { echo "FATAL: missing input: $f" >&2; exit 1; }
done
mkdir -p "$OUT"
echo "=== telomere_ends (#861) ==="
echo "real: $REAL_BAM"
echo "sim:  $SIM_BAM"
echo "ends: terminal $END_BP bp in $BIN bp bins; contigs >= $MIN_CONTIG bp, excluding $EXCLUDE"

# ── regions: contig, end (p|q|ctrl), 1-based start and stop, contig length ──────
awk -v min="$MIN_CONTIG" -v ex=",$EXCLUDE," -v w="$END_BP" -v f="$CTRL_FRAC" '
    $2 >= min && index(ex, "," $1 ",") == 0 {
        c = int($2 * f)
        print $1, "p", 1, w, $2
        print $1, "q", $2 - w + 1, $2, $2
        print $1, "ctrl", c, c + w - 1, $2
    }' OFS='\t' "$REFERENCE.fai" > "$OUT/regions.tsv"
n_regions=$(wc -l < "$OUT/regions.tsv")
[[ "$n_regions" -gt 0 ]] || { echo "FATAL: no contigs >= $MIN_CONTIG bp in $REFERENCE.fai" >&2; exit 1; }
echo "regions: $n_regions ($((n_regions / 3)) contigs)"

# Distance from the chromosome end (p, q) or from the window start (ctrl), 0-based.
DIST='function dist(pos) { return end == "p" ? pos - 1 : (end == "q" ? len - pos : pos - start) }'

# ── telomeric array length per end, from the reference ──────────────────────────
: > "$OUT/arrays.tsv"
while IFS=$'\t' read -r contig end start stop len; do
    [[ "$end" == ctrl ]] && continue
    samtools faidx "$REFERENCE" "$contig:$start-$stop" | awk -v end="$end" -v w="$BIN" -v min="$TEL_MIN" \
        -v contig="$contig" '
        NR > 1 { s = s toupper($0) }
        END {
            n = int(length(s) / w)
            arr = 0; miss = 0
            for (k = 0; k < n && miss < 3; k++) {
                # Window k counted from the chromosome end.
                off = (end == "p") ? k * w : length(s) - (k + 1) * w
                win = substr(s, off + 1, w)
                c = gsub(/TTAGGG/, "", win); win2 = substr(s, off + 1, w); c2 = gsub(/CCCTAA/, "", win2)
                if ((c > c2 ? c : c2) >= min) { arr = (k + 1) * w; miss = 0 } else miss++
            }
            print contig, end, arr
        }' OFS='\t' >> "$OUT/arrays.tsv"
done < "$OUT/regions.tsv"
echo "telomeric arrays in the reference (bp): $(awk '{ printf "%s%s:%s=%s", (NR > 1 ? " " : ""), $1, $2, $3 }' "$OUT/arrays.tsv")"

# ── per library: depth bins and read stats ──────────────────────────────────────
: > "$OUT/depth_bins.tsv"
: > "$OUT/read_stats.tsv"
for lib in real sim; do
    bam_var="${lib^^}_BAM"; bam="${!bam_var}"
    while IFS=$'\t' read -r contig end start stop len; do
        for q in 0 20; do
            # Collect first, then bin: a failure in samtools must stop the job, not yield no bins.
            depth=$(samtools depth -a -Q "$q" -r "$contig:$start-$stop" "$bam")
            awk -v lib="$lib" -v contig="$contig" -v end="$end" -v start="$start" -v len="$len" \
                -v q="$q" -v w="$BIN" -v want="$END_BP" "$DIST"'
                { b = int(dist($2) / w) * w; s[b] += $3; n[b]++; tot++ }
                END {
                    if (tot != want) { printf "FATAL: %s %s:%s depth has %d positions, expected %d\n", lib, contig, end, tot, want > "/dev/stderr"; exit 1 }
                    for (b in s) print lib, contig, end, q, b, s[b] / n[b]
                }' OFS='\t' <<< "$depth" >> "$OUT/depth_bins.tsv"
        done
    done < "$OUT/regions.tsv"

    # Reads over the array and over the rest of the terminal window (ctrl: the whole window).
    while IFS=$'\t' read -r contig end start stop len; do
        arr=0
        [[ "$end" != ctrl ]] && arr=$(awk -v c="$contig" -v e="$end" '$1 == c && $2 == e { print $3 }' "$OUT/arrays.tsv")
        for part in array flank; do
            if [[ "$end" == ctrl ]]; then
                [[ "$part" == array ]] && continue
                a=$start; b=$stop
            elif [[ "$part" == array ]]; then
                [[ "$arr" -gt 0 ]] || continue
                if [[ "$end" == p ]]; then a=1; b=$arr; else a=$((len - arr + 1)); b=$len; fi
            else
                if [[ "$end" == p ]]; then a=$((arr + 1)); b=$stop; else a=$start; b=$((len - arr)); fi
            fi
            reads=$(samtools view -F 0xD04 "$bam" "$contig:$a-$b")
            awk -v lib="$lib" -v contig="$contig" -v end="$end" -v part="$part" -v bp="$((b - a + 1))" '
                NF { n++; if ($5 == 0) z++; if ($6 ~ /[SH]/) c++ }
                END { print lib, contig, end, part, bp, n + 0, z + 0, c + 0 }' OFS='\t' <<< "$reads" >> "$OUT/read_stats.tsv"
        done
    done < "$OUT/regions.tsv"
done

# ── telomeric reads in the FASTQs, no mapping involved ──────────────────────────
: > "$OUT/telomeric_reads.tsv"
for lib in real sim; do
    bam_var="${lib^^}_BAM"; dir=$(dirname "${!bam_var}")
    fqs=("$dir/r1.fastq.gz" "$dir/r2.fastq.gz")
    if [[ ! -s "${fqs[0]}" || ! -s "${fqs[1]}" ]]; then
        echo "FATAL: no r1/r2.fastq.gz beside ${!bam_var}" >&2; exit 1
    fi
    counts=$(zcat "${fqs[@]}" | awk 'NR % 4 == 2 { n++; if (index($0, "TTAGGGTTAGGGTTAGGGTTAGGG") || index($0, "CCCTAACCCTAACCCTAACCCTAA")) t++ }
        END { print n + 0, t + 0 }')
    read -r n t <<< "$counts"
    # The pass must have read every read: the BAM's primary count is the same library.
    primary=$(samtools view -c -F 0x900 "${!bam_var}")
    [[ "$n" -eq "$primary" ]] || { echo "FATAL: $lib FASTQs hold $n reads, BAM has $primary primary" >&2; exit 1; }
    echo -e "$lib\t$n\t$t" >> "$OUT/telomeric_reads.tsv"
done

# ── summary ─────────────────────────────────────────────────────────────────────
awk -v w="$BIN" '
    FILENAME ~ /arrays/ { arr[$1 SUBSEP $2] = $3; next }
    FILENAME ~ /depth_bins/ {
        lib = $1; key = $2 SUBSEP $3; q = $4; b = $5; d = $6
        if ($3 == "ctrl") { if (q == 0) { cs[lib, $2] += d; cn[lib, $2]++ } ; next }
        a = arr[key]
        zone = (b < a) ? "array" : "flank"
        ds[lib, key, q, zone] += d; dn[lib, key, q, zone]++
        if (b < 500) { es[lib, key, q] += d; en[lib, key, q]++ }
        ends[key] = 1; libs[lib] = 1
        next
    }
    FILENAME ~ /read_stats/ {
        rs_n[$1, $2 SUBSEP $3, $4] = $6; rs_z[$1, $2 SUBSEP $3, $4] = $7; rs_c[$1, $2 SUBSEP $3, $4] = $8
        if ($3 == "ctrl") { cr_n[$1] += $6; cr_z[$1] += $7; cr_c[$1] += $8 }
        next
    }
    function median(v, m,   i, j, t) {
        for (i = 2; i <= m; i++) for (j = i; j > 1 && v[j - 1] > v[j]; j--) { t = v[j]; v[j] = v[j - 1]; v[j - 1] = t }
        return m % 2 ? v[(m + 1) / 2] : (v[m / 2] + v[m / 2 + 1]) / 2
    }
    function frac(a, b) { return b > 0 ? sprintf("%.3f", a / b) : "NA" }
    function rel(s, n, lib) { return n > 0 ? sprintf("%.2f", s / n / ctrl[lib]) : "NA" }
    END {
        for (lib in libs) {
            m = 0; delete v
            for (k in cn) { split(k, p, SUBSEP); if (p[1] == lib) v[++m] = cs[k] / cn[k] }
            ctrl[lib] = median(v, m)
            printf "%s: control depth median %.2f over %d windows; control MAPQ0 %s, clipped %s\n", lib, ctrl[lib], m, frac(cr_z[lib], cr_n[lib]), frac(cr_c[lib], cr_n[lib])
        }
        print ""
        print "Depth is relative to the library'"'"'s control median. edge = terminal 500 bp; arr = the array; q20 = MAPQ >= 20."
        printf "%-8s %-3s %6s | %-35s | %-35s\n", "contig", "end", "arr_bp", "real: edge arr arrq20 flank mq0 clip", "sim: edge arr arrq20 flank mq0 clip"
        n = 0
        for (key in ends) order[++n] = key
        for (i = 2; i <= n; i++) for (j = i; j > 1 && order[j - 1] > order[j]; j--) { t = order[j]; order[j] = order[j - 1]; order[j - 1] = t }
        for (i = 1; i <= n; i++) {
            key = order[i]; split(key, p, SUBSEP)
            line = sprintf("%-8s %-3s %6d", p[1], p[2], arr[key])
            for (li = 1; li <= 2; li++) {
                lib = li == 1 ? "real" : "sim"
                line = line sprintf(" | %5s %5s %5s %5s %5s %5s", rel(es[lib, key, 0], en[lib, key, 0], lib),
                    rel(ds[lib, key, 0, "array"], dn[lib, key, 0, "array"], lib),
                    rel(ds[lib, key, 20, "array"], dn[lib, key, 20, "array"], lib),
                    rel(ds[lib, key, 0, "flank"], dn[lib, key, 0, "flank"], lib),
                    frac(rs_z[lib, key, "array"], rs_n[lib, key, "array"]),
                    frac(rs_c[lib, key, "array"], rs_n[lib, key, "array"]))
            }
            print line
        }
    }' "$OUT/arrays.tsv" "$OUT/depth_bins.tsv" "$OUT/read_stats.tsv" | tee "$OUT/summary.txt"

echo
awk '{ printf "%s: %d reads, %d with 4x TTAGGG/CCCTAA = %.1f per million\n", $1, $2, $3, $2 ? $3 / $2 * 1e6 : 0 }' \
    "$OUT/telomeric_reads.tsv" | tee -a "$OUT/summary.txt"

if [[ "${NO_ARCHIVE:-0}" != 1 ]]; then
    archive_run "telomere_861/ends" "$OUT" summary.txt arrays.tsv depth_bins.tsv read_stats.tsv telomeric_reads.tsv regions.tsv
fi
echo "DONE."
