# Summarise indel positions by reference homopolymer run length and support class.
#
# A SEPARATE FILE so the test runs the same program the job does. Reads four files,
# identified by EXACT PATH via -v, not by matching a substring of FILENAME.
#
# WHY EXACT. The first version matched `FILENAME ~ /ctx/`. Job 21671697 ran with
# OUTDIR=/scratch/.../indelctx_21671697, so bg.tsv's full path contained "ctx" and was read
# as context data -- the background never accumulated, bgtot stayed 0, and every enrichment
# printed 0.00x. The observed histogram was fine, so the table looked complete and the one
# column the job exists for was silently empty. The local test passed because `mktemp -d`
# produced a path without "ctx" in it: the fixture did not reproduce the production path.
#
# Pass -v f_ind= -v f_dep= -v f_ctx= -v f_bg= with the four paths, plus
# -v mx=<max run bin> -v hf=<high support fraction> -v lf=<low support fraction>.
BEGIN {
    if (f_ind == "" || f_dep == "" || f_ctx == "" || f_bg == "") {
        print "indel_context_summarise.awk: need -v f_ind= -v f_dep= -v f_ctx= -v f_bg=" > "/dev/stderr"
        exit 2
    }
}
FILENAME == f_ind { sup[$1 SUBSEP $2] = $3; next }
FILENAME == f_dep { dep[$1 SUBSEP $2] = $3; next }
FILENAME == f_ctx { hp[$1 SUBSEP $2]  = $3; next }
# Binned the SAME way as the observed side. Capping one and not the other makes the top row
# divide by an empty bin and report 0.00x for the longest runs -- exactly the rows this
# measurement exists to read.
FILENAME == f_bg  { b = $1 + 0; if (b > mx) b = mx; bgn[b] += $2; bgtot += $2; next }
END {
    for (k in sup) {
        # A position with no context row is not in a run we measured; treat it as run 1
        # rather than dropping it, so the denominator stays the full set of indels.
        h = (k in hp) ? hp[k] + 0 : 1
        if (h > mx) h = mx
        n[h]++; tot++
        d = dep[k] + 0
        if (d <= 0) { nodep++; continue }
        f = sup[k] / d
        if (f >= hf)     { hi[h]++; thi++ }
        else if (f < lf) { lo[h]++; tlo++ }
        else             { mid[h]++; tmid++ }
    }
    print ""
    print "════════════════════════════════════════════════════════════════"
    print "INDEL POSITIONS BY REFERENCE HOMOPOLYMER RUN LENGTH"
    print ""
    printf "%-7s %9s %9s %9s %9s %11s %12s\n", "run", "n", "high", "low", "mid", "bg_share", "enrichment"
    for (h = 1; h <= mx; h++) {
        if (!(h in n) && !(h in bgn)) continue
        bs = bgtot ? bgn[h] / bgtot : 0
        os = tot   ? n[h]   / tot   : 0
        printf "%-7s %9d %9d %9d %9d %11.6f %11.2fx\n", (h == mx ? ">=" mx : h ""), \
               n[h]+0, hi[h]+0, lo[h]+0, mid[h]+0, bs, (bs > 0 ? os / bs : 0)
    }
    print ""
    printf "totals: %d positions -- %d high, %d mid, %d low, %d without depth\n",
           tot, thi+0, tmid+0, tlo+0, nodep+0
    printf "background: %d reference bases (N runs excluded)\n", bgtot
    # Rule 4 applied to the denominator of the headline column. An empty background makes
    # every enrichment read 0.00x, which looks like a measurement and is the absence of one.
    if (bgtot == 0) {
        print ""
        print "FATAL: the background is empty, so every enrichment above is 0.00x by"
        print "       construction and none of them means anything. The run measured where"
        print "       indels are and not whether that is more often than chance."
        exit 1
    }
    print ""
    print "HOW TO READ IT. enrichment > 1 means indels favour that run length more than"
    print "the reference offers it. WHICH COLUMN carries the enrichment is the answer:"
    print "  high -> germline/somatic variants in homopolymers   -> issue #378, placement"
    print "  low  -> PCR/sequencer slippage                      -> the error model"
    print "  both -> both mechanisms; #378 first, and split the acceptance criterion"
    print ""
    print "This job MEASURES. It does not say the simulator is wrong, only where real"
    print "indels are and what kind they are."
    print "════════════════════════════════════════════════════════════════"
}
