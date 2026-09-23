# Summarise the per-read records from homopolymer_measure.awk.
#
# Input:  the measurements TSV — `<locus>\tRUN\t<delta>`, `<locus>\tOFF\t<side>\t<offset>`,
#         `<locus>\tSKIP\t<reason>`.
# Output: a report on stdout.
#
# A SEPARATE FILE so the job and its test run the same program (see realism_summarise.awk for
# what happens otherwise).
#
# EVERY DENOMINATOR IS STATED. A slippage distribution over an unknown number of reads is not
# a result, and the skip reasons are the difference between "reads do not vary here" and "the
# anchor did not find them" (rule 4). They are printed even when zero.
#
# POSIX AWK ONLY.

BEGIN {
    FS = "\t"
    n_run = 0
    n_off = 0
    n_skip = 0
}

$2 == "RUN" {
    d = $3 + 0
    run[d]++
    n_run++
    loci_run[$1] = 1
    if (d != 0) n_run_nonzero++
    if (d < run_min || n_run == 1) run_min = d
    if (d > run_max || n_run == 1) run_max = d
    next
}

$2 == "OFF" {
    o = $4 + 0
    off[o]++
    n_off++
    side[$3]++
    a = (o < 0) ? -o : o
    if (a == 0)       off_band["0 (no indel implied)"]++
    else if (a <= 3)  off_band["1-3 (slippage)"]++
    else if (a <= 19) off_band["4-19"]++
    else              off_band["20+ (large indel)"]++
    next
}

$2 == "SKIP" {
    skip[$3]++
    n_skip++
    next
}

function pct(a, b) { return (b > 0) ? a * 100.0 / b : 0 }

END {
    nloci = 0
    for (k in loci_run) nloci++

    printf "════ homopolymer run lengths, read directly from the reads ════\n\n"
    printf "loci with at least one measured read : %d\n", nloci
    printf "run-length measurements              : %d\n", n_run
    printf "clip-offset measurements             : %d\n", n_off
    printf "skipped read/measurement pairs       : %d\n\n", n_skip

    if (n_run > 0) {
        printf "READ RUN LENGTH minus REFERENCE RUN LENGTH\n"
        printf "  %d of %d reads (%.1f%%) differ from the reference\n",
            n_run_nonzero, n_run, pct(n_run_nonzero, n_run)
        printf "  range %+d to %+d\n", run_min, run_max
        # Ascending, without gawk's asort: walk the observed range.
        for (d = run_min; d <= run_max; d++) {
            if (!(d in run)) continue
            printf "  %+4d  %7d  %5.1f%%  ", d, run[d], pct(run[d], n_run)
            bar = int(pct(run[d], n_run) / 2)
            for (i = 0; i < bar; i++) printf "#"
            printf "\n"
        }
        printf "\n"
    } else {
        printf "NO RUN-LENGTH MEASUREMENTS. Check the skip reasons below before reading\n"
        printf "anything else: an empty distribution and an anchor that never matched look\n"
        printf "identical in a histogram.\n\n"
    }

    if (n_off > 0) {
        printf "IMPLIED INDEL AT SOFT-CLIPPED READS (+ deletion, - insertion)\n"
        printf "  measured from %d clips", n_off
        for (s in side) printf "  [%s %d]", s, side[s]
        printf "\n"
        # Fixed band order, so an absent band prints as zero rather than vanishing.
        split("0 (no indel implied)|1-3 (slippage)|4-19|20+ (large indel)", bands, "|")
        for (i = 1; i <= 4; i++) {
            b = bands[i]
            printf "  %-22s %7d  %5.1f%%\n", b, off_band[b] + 0, pct(off_band[b] + 0, n_off)
        }
        printf "\n  WHAT THE BANDS MEAN\n"
        printf "    0     the clip would have mapped with no indel -- the aligner clipped\n"
        printf "          over ambiguity, not over an event\n"
        printf "    1-3   slippage the aligner declined to gap\n"
        printf "    20+   a genuine large indel; that is variant placement, not error\n\n"
    } else {
        printf "NO CLIP-OFFSET MEASUREMENTS.\n\n"
    }

    printf "SKIPS (why a read yielded no measurement)\n"
    if (n_skip == 0) printf "  none\n"
    for (k in skip) printf "  %-52s %7d  %5.1f%%\n", k, skip[k], pct(skip[k], n_skip + n_run + n_off)
    printf "\nThis job MEASURES. Nothing above is a pass or a fail.\n"
}
