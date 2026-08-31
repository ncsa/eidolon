# Join candidate clip boundaries to the indels near them.
#
# Settles a question the two measurements cannot answer separately: is there a VARIANT at
# the loci where reads agree on a clip boundary, or none at all?
#
#   (a) an indel sits in the homopolymer and reads carrying it are clipped rather than
#       gapped  -> the fix is variant placement, issue #378
#   (b) no indel anywhere near  -> the locus is simply low-complexity and the aligner clips
#       reads whose ends land in the repeat with too few anchoring bases. Then neither #378
#       nor the error model is the fix: eidolon's reads are too close to the reference to be
#       clipped ANYWHERE, and that is the thing to change.
#
# Both produce clustered clips in homopolymers, so the homopolymer enrichment alone does not
# separate them. The presence and support of a nearby indel does.
#
# A SEPARATE FILE so the test runs the same program the job runs.
#
# Pass -v win=<bp> -v hf=<high frac> -v lf=<low frac>. Reads: indels, depth, candidates.
FILENAME ~ /indels/ { isup[$1 SUBSEP $2] = $3; ipos[++ni] = $1 SUBSEP $2; next }
FILENAME ~ /depth/  { dep[$1 SUBSEP $2] = $3; next }
# candidates: the realism panel's --dump-candidates TSV. Header line skipped by the
# non-numeric position, not by NR -- this is the third file, so NR is already large.
FILENAME ~ /candidates/ {
    if ($2 !~ /^[0-9]+$/) next
    cn++
    best = -1; bestf = -1
    for (d = 0; d <= win; d++) {
        for (s = -1; s <= 1; s += 2) {
            p = $2 + d * s
            k = $1 SUBSEP p
            if (k in isup) {
                dd = dep[k] + 0
                f = (dd > 0) ? isup[k] / dd : -1
                # Nearest wins; among equals the better-supported one.
                if (best < 0 || f > bestf) { best = d; bestf = f }
            }
        }
        if (best >= 0) break     # nearest first
    }
    if (best < 0)          { none++ }
    else if (bestf < 0)    { nodep++ }
    else if (bestf >= hf)  { high++; dsum += best }
    else if (bestf < lf)   { low++;  dsum += best }
    else                   { mid++;  dsum += best }
    next
}
END {
    print ""
    print "════════════════════════════════════════════════════════════════"
    print "CANDIDATE CLIP BOUNDARIES vs NEARBY INDELS"
    printf "  window: +/- %d bp\n", win
    print ""
    withi = high + mid + low + nodep
    printf "  %6d candidate sites\n", cn
    printf "  %6d (%5.1f%%) have an indel within the window\n", withi, (cn ? withi*100/cn : 0)
    printf "  %6d (%5.1f%%) have NONE\n", none, (cn ? none*100/cn : 0)
    print ""
    printf "  of those with an indel:  %d high-support, %d mid, %d low, %d without depth\n",
           high+0, mid+0, low+0, nodep+0
    if (withi > 0) printf "  mean distance to the nearest indel: %.1f bp\n", dsum / (high+mid+low)
    print ""
    print "HOW TO READ IT."
    print "  mostly HIGH  -> variants in homopolymers; clipped instead of gapped -> #378"
    print "  mostly LOW   -> PCR/sequencer slippage                -> the error model"
    print "  mostly NONE  -> no variant at all: the aligner clips in low-complexity sequence"
    print "                  on its own, and the gap is that simulated reads are too close"
    print "                  to the reference to be clipped anywhere. Neither #378 nor the"
    print "                  error model would close it."
    print "════════════════════════════════════════════════════════════════"
}
