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
FILENAME == f_ind { sup[$1 SUBSEP $2] = $3; ilen[$1 SUBSEP $2] = $4; next }
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
        # Length is banked per support class, because the two classes feed different
        # models: `low` is the sequencing-error length distribution
        # (SequencingErrorModel's ins/del_length_distribution), `high` is variant indel
        # size. Pooling them would hand the error model variant-sized events.
        L = ilen[k] + 0
        if (f >= hf)     { hi[h]++; thi++; hlen[L]++; }
        else if (f < lf) { lo[h]++; tlo++; llen[L]++; lxr[h SUBSEP (L<0?-L:L)]++; }
        else             { mid[h]++; tmid++ }
    }
    print ""
    print "════════════════════════════════════════════════════════════════"
    print "INDEL POSITIONS BY REFERENCE HOMOPOLYMER RUN LENGTH"
    print ""
    # THREE enrichment columns, because the pooled one feeds no model.
    #
    # Each is a share-over-share against the SAME background, differing only in which
    # subset of indels forms the numerator and its denominator:
    #
    #   enr_all  -- every indel position. Mixes slippage with variants, so it is a
    #               description of the data and NOT the input to anything.
    #   enr_low  -- slippage only. This is the sequencing-error curve (#661, fitted by #662).
    #   enr_high -- variants only. This is the placement curve (#378).
    #
    # Printing only the pooled column is what made this file misleading: on job 21674484 it
    # read 52.06x at runs >=10 where the error curve the model actually ships is 39.20x,
    # because 270 of those 777 positions were variants. Anyone checking a fitted curve
    # against that column would see a 33% shortfall and "correct" the fit toward a
    # variant-contaminated number. Each model's own column is now printed beside it.
    printf "%-7s %8s %8s %8s %8s %11s %10s %10s %10s\n", \
           "run", "n", "high", "low", "mid", "bg_share", "enr_all", "enr_low", "enr_high"
    for (h = 1; h <= mx; h++) {
        if (!(h in n) && !(h in bgn)) continue
        bs = bgtot ? bgn[h] / bgtot : 0
        # Each subset is normalized by ITS OWN total, not by `tot`. Dividing the low count
        # by the pooled total would report a share of the wrong population.
        e_all  = (bs > 0 && tot > 0) ? (n[h]  + 0) / tot / bs : 0
        e_low  = (bs > 0 && tlo > 0) ? (lo[h] + 0) / tlo / bs : 0
        e_high = (bs > 0 && thi > 0) ? (hi[h] + 0) / thi / bs : 0
        printf "%-7s %8d %8d %8d %8d %11.6f %9.2fx %9.2fx %9.2fx\n", (h == mx ? ">=" mx : h ""), \
               n[h]+0, hi[h]+0, lo[h]+0, mid[h]+0, bs, e_all, e_low, e_high
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
    print "WHICH ENRICHMENT COLUMN TO USE:"
    print "  enr_low  IS the sequencing-error context curve. Compare a gen-bam-models fit"
    print "           (#662) against THIS column, and against the shipped default in"
    print "           eidolon-core/src/models/sequencing_error_model.rs."
    print "  enr_high IS the variant placement curve for #378."
    print "  enr_all  is neither. It pools both mechanisms and is here for description"
    print "           only -- a fit that reproduces it has absorbed variants into the"
    print "           error model, which is the exact conflation #661 forbids."
    # ── indel LENGTH by support class ────────────────────────────────────────
    # The input to SequencingErrorModel's length distribution. Its shipped default is
    # [0.999, 0.001] over lengths [1, 2], carried over from NEAT2 and never measured.
    print ""
    print "════════════════════════════════════════════════════════════════"
    print "INDEL LENGTH BY SUPPORT CLASS  (+ = insertion, - = deletion)"
    print ""
    printf "%-8s %10s %9s %10s %9s\n", "length", "low(n)", "low(frac)", "high(n)", "high(frac)"
    # Union of observed lengths, printed in signed order so deletions and insertions
    # read as two arms around zero rather than as one pooled magnitude.
    for (L in llen) seen[L] = 1
    for (L in hlen) seen[L] = 1
    nL = 0
    for (L in seen) { nL++; ord[nL] = L + 0 }
    for (a = 2; a <= nL; a++) { v = ord[a]; b = a - 1
        while (b >= 1 && ord[b] > v) { ord[b+1] = ord[b]; b-- }
        ord[b+1] = v }
    for (a = 1; a <= nL; a++) {
        L = ord[a]
        printf "%-8s %10d %9.4f %10d %9.4f\n", (L > 0 ? "+" L : L ""), \
               llen[L]+0, (tlo ? llen[L]/tlo : 0), hlen[L]+0, (thi ? hlen[L]/thi : 0)
    }
    printf "\ntotals: %d low (slippage), %d high (variant)\n", tlo+0, thi+0
    # Rule 4 on this table's own denominator: a length column over zero events is not a
    # distribution, and printing 0.0000 for every row would look like one.
    if (tlo == 0) {
        print ""
        print "FATAL: no low-support indels, so the slippage length distribution above is"
        print "       empty by construction and cannot be fitted from this run."
        exit 1
    }
    print ""
    # ── slippage length CONDITIONED on homopolymer run length ────────────────
    # Marginals cannot answer whether long indel errors are slippage. #661 already makes
    # indel errors ~39x more likely inside a 10+-mer while their LENGTH is drawn
    # independently, so the simulator places long deletions inside homopolymers without
    # any evidence that real ones sit there. This is the table that decides it.
    print ""
    print "════════════════════════════════════════════════════════════════"
    print "SLIPPAGE INDEL LENGTH BY HOMOPOLYMER RUN LENGTH  (low support only)"
    print ""
    printf "%-6s", "run"
    for (b = 1; b <= 6; b++) printf "%8s", b
    printf "%8s%8s%8s%10s\n", "7-9", "10-19", ">=20", "n"
    for (h = 1; h <= mx; h++) {
        rn = 0
        for (b = 1; b <= 60; b++) rn += lxr[h SUBSEP b] + 0
        if (rn == 0) continue
        printf "%-6s", (h == mx ? ">=" mx : h "")
        for (b = 1; b <= 6; b++) printf "%8d", lxr[h SUBSEP b] + 0
        s79 = 0;  for (b = 7;  b <= 9;  b++) s79 += lxr[h SUBSEP b] + 0
        s19 = 0;  for (b = 10; b <= 19; b++) s19 += lxr[h SUBSEP b] + 0
        s20 = 0;  for (b = 20; b <= 60; b++) s20 += lxr[h SUBSEP b] + 0
        printf "%8d%8d%8d%10d\n", s79, s19, s20, rn
    }
    print ""
    print "HOW TO READ IT. If slippage length is a function of run length, the mass moves"
    print "right as run length grows and >=20 bp events sit only in long runs. If the >=20"
    print "column is populated at run 1-2, those events are not slippage -- they are"
    print "mapping artifacts or rare somatic events landing in the low-support class, and"
    print "indel-error length should NOT be conditioned on run length."
    print ""
    print "A deletion longer than its run cannot be slippage by construction: a 3-mer"
    print "cannot lose 10 bases of register. Cells right of the diagonal are that case."
    print ""
    print "WHAT THE low COLUMN IS FOR. It is the measured length distribution of"
    print "sequencing-error indels -- the input to SequencingErrorModel's"
    print "ins_length_distribution / del_length_distribution, whose shipped default is"
    print "[0.999, 0.001] over lengths [1, 2] inherited from NEAT2. The high column is"
    print "variant indel size and belongs to placement, not to the error model."
    print ""
    print "This job MEASURES. It does not say the simulator is wrong, only where real"
    print "indels are and what kind they are."
    print "════════════════════════════════════════════════════════════════"
}
