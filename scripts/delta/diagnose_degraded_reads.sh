#!/usr/bin/env bash
# Ask WHY a read ends up in the degraded population, using only the FASTQ.
#
# WHAT THIS IS FOR. #694's model splits reads into two populations and reproduces the library,
# but "9.80% of reads are degraded" is a statement about where we put the cut, not about the
# sequencer. Three mechanisms would each produce that shape, and they are distinguishable from
# the FASTQ alone because each predicts something different:
#
#   cluster-level optics/fluidics   -> degraded reads CLUSTER SPATIALLY (tiles, tile edges)
#   fragment length                 -> degraded R2 does NOT predict degraded R1 for the same
#                                      cluster. Fragments over ~500 nt yield low-quality R2
#                                      while R1 stays high quality (PMID 30814542), so mate
#                                      concordance stays near the base rate.
#   template sequence               -> degraded reads differ in BASE COMPOSITION (GC, runs)
#
# None of the three is assumed. Each section reports its denominator so a near-zero stratum is
# visible rather than silently dropped.
#
# WHAT IT CANNOT ANSWER. Fragment length itself needs alignment, so this measures the mate
# concordance the fragment-length hypothesis predicts, not the fragment lengths. A weak
# concordance is consistent with that hypothesis and does not establish it.
#
# USAGE
#   sbatch scripts/delta/diagnose_degraded_reads.sh
#   R1=/path/a_R1.fastq.gz R2=/path/a_R2.fastq.gz NAME=lib2 sbatch scripts/delta/diagnose_degraded_reads.sh

#SBATCH --job-name=eidolon-degdiag
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=degdiag-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"   # resolves $SCRATCH; Delta does not export it

DATA_DIR="${DATA_DIR:-/work/nvme/bhrd/jallen17/hg002}"
R1="${R1:-$DATA_DIR/hg002_R1.fastq.gz}"
R2="${R2:-$DATA_DIR/hg002_R2.fastq.gz}"
NAME="${NAME:-hg002}"
STRIDE="${STRIDE:-65}"
WINDOW="${WINDOW:-50}"          # tail window, matching the fitter
CUT="${CUT:-25}"                # tail cut, matching the fitter
OFFSET="${OFFSET:-33}"
OUT="${OUT:-degdiag_$NAME.txt}"

[[ -f "$R1" ]] || { echo "FATAL: no R1 at $R1" >&2; exit 1; }
[[ -f "$R2" ]] || { echo "FATAL: no R2 at $R2" >&2; exit 1; }

echo "=== diagnose_degraded_reads ($NAME) ==="
echo "repo:   $REPO ($(cd "$REPO" && git rev-parse --short HEAD))"
echo "R1:     $R1"
echo "R2:     $R2"
echo "stride: $STRIDE   window: last $WINDOW bases   cut: Q$CUT"
echo

# R2 is streamed through a coprocess-free `getline` so both mates are classified together.
# Both files are strided identically, which relies on Illumina's record order being paired;
# the header check below turns that assumption into an assertion rather than a hope.
awk_status=0
zcat -f -- "$R1" | awk -v r2cmd="zcat -f -- '$R2'" -v s="$STRIDE" -v w="$WINDOW" \
    -v cut="$CUT" -v off="$OFFSET" '
function tailmean(q,   i, n, sum) {
    n = length(q); if (n < w) return -1          # too short to classify
    sum = 0
    for (i = n - w + 1; i <= n; i++) sum += index(QCHARS, substr(q, i, 1)) - 1
    return sum / w
}
function gcfrac(seq,   t, n, g) {
    n = length(seq); if (n == 0) return -1
    t = seq; g = gsub(/[GCgc]/, "", t)
    return g / n
}
BEGIN {
    # index() is 1-based, so position 0 of QCHARS is Phred 0 after the -1 above.
    for (i = off; i < off + 64; i++) QCHARS = QCHARS sprintf("%c", i)
    HOMO = "AAAAAAAAAA|CCCCCCCCCC|GGGGGGGGGG|TTTTTTTTTT"   # 10-mers
}
{ line1 = $0
  getline line2; getline line3; getline line4
  rec++
  if (((rec - 1) % s) != 0) next
  n_pairs++

  # pull the matching R2 record
  if ((r2cmd | getline h2) <= 0) { print "FATAL: R2 ran out at pair " rec > "/dev/stderr"; exit 3 }
  r2rec++
  while (((r2rec - 1) % s) != 0) {
      for (k = 0; k < 3; k++) if ((r2cmd | getline junk) <= 0) { print "FATAL: R2 truncated" > "/dev/stderr"; exit 3 }
      if ((r2cmd | getline h2) <= 0) { print "FATAL: R2 ran out at pair " rec > "/dev/stderr"; exit 3 }
      r2rec++
  }
  if ((r2cmd | getline s2) <= 0 || (r2cmd | getline p2) <= 0 || (r2cmd | getline q2) <= 0) {
      print "FATAL: R2 record truncated at pair " rec > "/dev/stderr"; exit 3 }

  # Assertion, not assumption: the two files must be in the same order, or every mate-pair
  # number below is meaningless. Compare the coordinate part, before any /1 /2 or " 1:N:0:".
  split(line1, a1, /[ \/]/); split(h2, a2, /[ \/]/)
  if (a1[1] != a2[1]) {
      printf "FATAL: mate headers diverge at pair %d: %s vs %s\n", rec, a1[1], a2[1] > "/dev/stderr"
      exit 3
  }

  # header: @instrument:run:flowcell:lane:tile:x:y
  nf = split(a1[1], H, ":")
  if (nf < 7) { n_unparsed++; next }
  lane = H[4]; tile = H[5]; xc = H[6] + 0; yc = H[7] + 0

  t1 = tailmean(line4); t2 = tailmean(q2)
  if (t1 < 0 || t2 < 0) { n_unclassifiable++; next }
  n_class++
  d1 = (t1 < cut); d2 = (t2 < cut)
  if (d1) n_d1++
  if (d2) n_d2++
  if (d1 && d2) n_both++

  # 2. spatial
  tile_n[lane ":" tile]++; if (d1) tile_d[lane ":" tile]++
  yb = int(yc / 2000); yband_n[yb]++; if (d1) yband_d[yb]++
  xb = int(xc / 2000); xband_n[xb]++; if (d1) xband_d[xb]++

  # 4. composition, R1 only (R2 composition is the reverse strand of a different region)
  g = gcfrac(line2)
  if (g >= 0) {
      if (d1) { gc_d += g; gc_d_n++ } else { gc_h += g; gc_h_n++ }
      if (line2 ~ HOMO) { if (d1) homo_d++; else homo_h++ }
  }
  nn = line2; n_amb = gsub(/N/, "", nn)
  if (d1) { amb_d += n_amb } else { amb_h += n_amb }
}
END {
    if (n_class == 0) { print "FATAL: nothing was classified" > "/dev/stderr"; exit 4 }
    printf "== 1. headline rates ==\n"
    printf "  pairs sampled          %d  (of %d records, stride %d)\n", n_pairs, rec, s
    printf "  classifiable pairs     %d\n", n_class
    printf "  unclassifiable (short) %d\n", n_unclassifiable + 0
    printf "  unparsed headers       %d\n", n_unparsed + 0
    printf "  R1 degraded            %d (%.2f%%)\n", n_d1, 100.0 * n_d1 / n_class
    printf "  R2 degraded            %d (%.2f%%)\n", n_d2, 100.0 * n_d2 / n_class
    printf "\n"

    printf "== 2. mate concordance: does a bad R2 predict a bad R1? ==\n"
    printf "  the cluster-level hypothesis predicts a large lift; the fragment-length one\n"
    printf "  predicts almost none, since R1 stays high quality on long fragments\n"
    p1 = n_d1 / n_class; p2 = n_d2 / n_class
    exp_both = p1 * p2 * n_class
    printf "  P(R1 deg)              %.4f\n", p1
    printf "  P(R2 deg)              %.4f\n", p2
    printf "  both degraded          %d observed, %.1f expected if independent\n", n_both + 0, exp_both
    if (n_d2 > 0) printf "  P(R1 deg | R2 deg)     %.4f   lift %.2fx\n", n_both / n_d2, (n_both / n_d2) / p1
    printf "\n"

    printf "== 3. spatial: tiles ==\n"
    nt = 0; worst = 0; best = 1; sum_r = 0
    for (t in tile_n) {
        if (tile_n[t] < 200) { skipped++; continue }     # denominator floor, counted not hidden
        r = (tile_d[t] + 0) / tile_n[t]; nt++; sum_r += r
        if (r > worst) { worst = r; worst_t = t }
        if (r < best)  { best = r;  best_t = t }
        rr[nt] = r
    }
    if (nt == 0) { printf "  no tile had 200+ reads; %d tiles skipped\n", skipped + 0 }
    else {
        mean_r = sum_r / nt; ss = 0
        for (i = 1; i <= nt; i++) ss += (rr[i] - mean_r) ^ 2
        sd = (nt > 1) ? sqrt(ss / (nt - 1)) : 0
        printf "  tiles with 200+ reads  %d  (%d skipped below the floor)\n", nt, skipped + 0
        printf "  mean rate              %.4f   sd %.4f   cv %.3f\n", mean_r, sd, (mean_r > 0 ? sd / mean_r : 0)
        printf "  worst tile             %s at %.4f\n", worst_t, worst
        printf "  best  tile             %s at %.4f\n", best_t, best
        printf "  worst/best             %.2fx\n", (best > 0 ? worst / best : 0)
    }
    printf "\n"

    printf "== 4. spatial: position within tile (2000-unit bands) ==\n"
    printf "  y-band    n         rate\n"
    for (b = 0; b <= 40; b++) if (b in yband_n && yband_n[b] >= 200)
        printf "  %-8d  %-8d  %.4f\n", b * 2000, yband_n[b], (yband_d[b] + 0) / yband_n[b]
    printf "  x-band    n         rate\n"
    for (b = 0; b <= 40; b++) if (b in xband_n && xband_n[b] >= 200)
        printf "  %-8d  %-8d  %.4f\n", b * 2000, xband_n[b], (xband_d[b] + 0) / xband_n[b]
    printf "\n"

    printf "== 5. composition of degraded vs healthy R1 ==\n"
    if (gc_d_n > 0 && gc_h_n > 0) {
        printf "  mean GC   degraded %.4f (n %d)   healthy %.4f (n %d)   ratio %.3f\n",
               gc_d / gc_d_n, gc_d_n, gc_h / gc_h_n, gc_h_n, (gc_h / gc_h_n > 0 ? (gc_d / gc_d_n) / (gc_h / gc_h_n) : 0)
        hd = (homo_d + 0) / gc_d_n; hh = (homo_h + 0) / gc_h_n
        printf "  10-mer homopolymer  degraded %.4f (%d)   healthy %.4f (%d)   ratio %s\n",
               hd, homo_d + 0, hh, homo_h + 0, (hh > 0 ? sprintf("%.3f", hd / hh) : (hd > 0 ? "inf" : "n/a"))
        printf "  N bases per read    degraded %.4f   healthy %.4f\n",
               (amb_d + 0) / gc_d_n, (amb_h + 0) / gc_h_n
    } else printf "  one composition arm was empty: degraded n %d, healthy n %d\n", gc_d_n + 0, gc_h_n + 0
    printf "\n== done ==\n"
}' > "$OUT" || awk_status=$?

# The pass must be shown to have finished, outside anything the pass itself prints: a killed
# awk never reaches END, and a partial file reads like a short successful run.
if [[ "$awk_status" -ne 0 ]] || ! grep -q "^== done ==" "$OUT" 2>/dev/null; then
    echo "FATAL: the measurement pass did not complete (awk status $awk_status)." >&2
    echo "  Partial output is in $OUT; do not read numbers off it." >&2
    exit 1
fi

cat "$OUT"
echo
echo "Written to $OUT"
