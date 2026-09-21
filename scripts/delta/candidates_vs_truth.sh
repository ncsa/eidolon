#!/usr/bin/env bash
# Are the realism panel's real-side candidate breakpoints the DONOR'S OWN variants?
#
# WHY. The panel simulates with no variants on purpose (realism_panel.sbatch:528), so every
# artifact on the simulated side is one the simulator produced by itself. The real side has no
# such restriction: HG002 differs from GRCh38 at millions of small indels, and every one of
# them is a place where real reads disagree with the reference and simulated reads cannot.
#
# Job 22243396 classified its 550 real candidates as 3.8% repeat-like, 26.7% junction-like and
# 69.5% "minority-support, mappable, normal MAPQ" -- which is what a heterozygous small indel
# looks like. If most of those sit on known HG002 variants, then a large part of the headline
# cand_per_mb gap is a property of the COMPARISON rather than of the simulator, and the gap
# needs restating before it is quoted again.
#
# Read-only, seconds, no simulation. It reads the candidate dump the panel already wrote.
#
# COORDINATES. `real_candidates.tsv` column 2 is 0-BASED (reader.rs:149 converts noodles'
# 1-based alignment_start with `.get() - 1`); VCF POS is 1-based. The +1 below is that, and
# getting it wrong would shift every comparison by one base.
#
# WHY A WINDOW. A clip boundary sits near an indel, not on its POS -- aligners place the
# breakpoint a few bases either side depending on local sequence. Exact matching would
# undercount badly. W is the tolerance.
#
# THE CONTROL IS THE POINT. "70% are near a variant" means nothing without knowing what
# fraction of ARBITRARY positions in the same loci are near one. Variants are dense enough
# that a loose window hits often by chance. The shuffled row is that denominator.
#
# USAGE
#   bash scripts/delta/candidates_vs_truth.sh
#   OUTDIR=... HG002_VCF=... W=20 bash scripts/delta/candidates_vs_truth.sh

#SBATCH --job-name=eidolon-candtruth
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=candtruth-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"

OUTDIR="${OUTDIR:-$SCRATCH/realism_hg002_post734_b}"
CAND="${CAND:-$OUTDIR/real_candidates.tsv}"
REGIONS="${REGIONS:-$OUTDIR/regions.bed}"
VCF="${HG002_VCF:-$SCRATCH/neat_data/hg002/HG002_GRCh38_v4.2.1_benchmark.vcf.gz}"
W="${W:-10}"
SEED="${SEED:-22243396}"

for f in "$CAND" "$REGIONS" "$VCF"; do
    [[ -f "$f" ]] || { echo "FATAL: not a file: $f" >&2; exit 1; }
done
command -v bcftools >/dev/null || { echo "FATAL: bcftools not on PATH (conda activate bioinf)" >&2; exit 1; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

echo "=== candidates vs truth ==="
echo "candidates: $CAND"
echo "truth:      $VCF"
echo "regions:    $REGIONS"
echo "window:     +/-${W} bp"
echo

# Truth variants inside the measured loci only. Restricting to the regions keeps the shuffled
# control comparable -- both are drawn from the same span.
bcftools query -R "$REGIONS" -f '%CHROM\t%POS\n' "$VCF" > "$TMP/truth.tsv"
n_truth=$(wc -l < "$TMP/truth.tsv")
span=$(awk '{s += $3 - $2} END {print s+0}' "$REGIONS")
[[ "$n_truth" -gt 0 ]] || { echo "FATAL: no truth variants in the measured regions -- wrong VCF, or a contig-naming mismatch (chr1 vs 1)" >&2; exit 1; }
echo "truth variants in the measured $span bp: $n_truth  ($(awk -v n="$n_truth" -v s="$span" 'BEGIN{printf "1 per %.0f bp", s/n}'))"
echo

awk -F'\t' -v W="$W" -v seed="$SEED" -v span="$span" '
    # Truth positions, 1-based, hashed by contig:pos.
    FNR == NR { t[$1 ":" $2] = 1; next }
    FNR == 1  { next }                                  # candidate header
    {
        # 0-based candidate -> 1-based, then the window.
        p = $2 + 1
        cls = ($8 >= 0.5) ? "repeat-like" : (($7 >= 0.30) ? "junction-like" : "minority-support")
        n[cls]++; n["ALL"]++
        hit = 0
        for (d = -W; d <= W; d++) if (($1 ":" (p + d)) in t) { hit = 1; break }
        if (hit) { h[cls]++; h["ALL"]++ }
        # Remember the contig so the control can be drawn from the same loci.
        ctg[$1] = 1
    }
    END {
        printf "%-20s %8s %8s %8s\n", "class", "sites", "near", "pct"
        for (k in n) if (k != "ALL")
            printf "%-20s %8d %8d %7.1f%%\n", k, n[k], h[k]+0, 100*(h[k]+0)/n[k]
        printf "%-20s %8d %8d %7.1f%%\n", "ALL", n["ALL"], h["ALL"]+0, 100*(h["ALL"]+0)/n["ALL"]
        printf "\n(denominators are the site counts; a percentage over an unknown count is not a result)\n"
    }
' "$TMP/truth.tsv" "$CAND"

# ── the control ──────────────────────────────────────────────────────────────
# Same number of positions, drawn uniformly from the same loci. Without this the number
# above has nothing to sit against: in a region dense with variants, a +/-10 bp window hits
# often by chance.
n_cand=$(( $(wc -l < "$CAND") - 1 ))
awk -v n="$n_cand" -v seed="$SEED" 'BEGIN { srand(seed) }
    { ctg[NR] = $1; lo[NR] = $2; hi[NR] = $3; tot += $3 - $2; cum[NR] = tot; m = NR }
    END {
        for (i = 0; i < n; i++) {
            x = rand() * tot
            for (j = 1; j <= m; j++) if (x < cum[j]) break
            printf "%s\t%d\n", ctg[j], lo[j] + int(rand() * (hi[j] - lo[j])) + 1
        }
    }' "$REGIONS" > "$TMP/shuffled.tsv"

awk -F'\t' -v W="$W" '
    FNR == NR { t[$1 ":" $2] = 1; next }
    { n++; for (d = -W; d <= W; d++) if (($1 ":" ($2 + d)) in t) { h++; break } }
    END { printf "%-20s %8d %8d %7.1f%%   <- chance level\n", "shuffled control", n, h+0, 100*(h+0)/n }
' "$TMP/truth.tsv" "$TMP/shuffled.tsv"

echo
echo "READ IT THIS WAY. The minority-support row is the one that matters -- it is 69.5% of the"
echo "candidates and the panel's own note says it looks like a small indel. If it sits far above"
echo "the shuffled control, those candidates are the donor's variation, which the simulated arm"
echo "is constructed not to have, and the cand_per_mb gap is partly a comparison artifact."
