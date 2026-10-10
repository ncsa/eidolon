# indel_probability from one BAM pass (#746).
#
# WHAT IS FITTED. For each base that takes a sequencing error, SequencingErrorModel makes it
# an indel with probability indel_probability x curve(run), and a substitution otherwise
# (sequencing_error_model.rs, generate_sequencing_error). curve is the homopolymer context
# curve (#661), normalized so its mean over the reference background is 1 on the sample it
# was fitted from. If errors fall on bases independently of run length, the share of error
# EVENTS that are indels is
#     E_indel / (E_indel + E_sub) = indel_probability x mean_bg(curve)
# so indel_probability = share / mean_bg(curve). The mean is computed here over THIS run's
# background with the shipped curve, and reported, rather than assumed to be 1.
#
# ERRORS ARE THE LOW-SUPPORT CLASS. Each indel or mismatch position is classed by its
# support over local depth with indel_context_summarise.awk's rule: >= hf is a variant,
# < lf is an error, anything between is ambiguous and counted in neither. Both kinds use
# the same rule and the same unit, reads carrying the alteration, so the ratio compares like
# with like and its denominator cancels. The per-base rates below are the sanity check on
# that, and they are the only place the denominator matters.
#
# in (vars): hf, lf, mx (run-length cap), curve (comma-separated, run 1..), f_ind, f_sub,
#            f_dep, f_bg, f_bases (aligned_bases), f_out (key/value results)
FILENAME == f_ind   { sup_i[$1 SUBSEP $2] = $3; next }
FILENAME == f_sub   { sup_s[$1 SUBSEP $2] = $3; next }
FILENAME == f_dep   { dep[$1 SUBSEP $2] = $3; next }
FILENAME == f_bg    { b = $1 + 0; if (b > mx) b = mx; bgn[b] += $2; bgtot += $2; next }
FILENAME == f_bases { if ($1 == "aligned_bases") bases = $2 + 0; next }
END {
    if (f_ind == "" || f_sub == "" || f_dep == "" || f_bg == "" || f_bases == "" || f_out == "") {
        print "FATAL: every input file must be named (f_ind f_sub f_dep f_bg f_bases f_out)" > "/dev/stderr"
        exit 2
    }
    nc = split(curve, cv, ",")
    for (h = 1; h <= mx; h++) cw[h] = cv[h <= nc ? h : nc] + 0

    for (k in sup_i) {
        d = dep[k] + 0; if (d <= 0) { inodep++; continue }
        f = sup_i[k] / d
        if (f >= hf) { ihi++; ivar += sup_i[k] }
        else if (f < lf) { ilo++; iev += sup_i[k] }
        else { imid++ }
    }
    for (k in sup_s) {
        d = dep[k] + 0; if (d <= 0) { snodep++; continue }
        f = sup_s[k] / d
        spool += sup_s[k]
        if (f >= hf) { shi++; svar += sup_s[k] }
        else if (f < lf) { slo++; sev += sup_s[k] }
        else { smid++ }
    }

    cm = 0
    for (h = 1; h <= mx; h++) cm += bgn[h] * cw[h]
    cmean = bgtot > 0 ? cm / bgtot : 0

    print ""
    print "════════════════════════════════════════════════════════════════"
    print "INDEL_PROBABILITY (#746): share of sequencing-error events that are indels"
    print ""
    printf "  %-34s %10s %10s %10s %10s\n", "", "error", "ambiguous", "variant", "no depth"
    printf "  %-34s %10d %10d %10d %10d\n", "indel positions", ilo, imid, ihi, inodep
    printf "  %-34s %10d %10d %10d %10d\n", "substitution positions", slo, smid, shi, snodep
    printf "  %-34s %10d %10s %10d\n", "indel events (reads)", iev, "", ivar
    printf "  %-34s %10d %10s %10d\n", "substitution events (reads)", sev, "", svar
    printf "  substitution events in every class: %d (error %d is the part that is not variants)\n", spool, sev
    print ""
    printf "  denominator: %d aligned bases\n", bases

    if (sev == 0 || iev == 0 || bases == 0) {
        print ""
        print "FATAL: nothing to take a ratio of -- error indel events " iev ", error substitution"
        print "       events " sev ", aligned bases " bases ". A ratio over an empty class is not a"
        print "       measurement. Check the regions, the read filter and the support thresholds."
        exit 1
    }
    share = iev / (iev + sev)
    fit = share / cmean
    printf "  indel error events per aligned base:        %.3g\n", iev / bases
    printf "  substitution error events per aligned base: %.3g\n", sev / bases
    printf "  observed indel share of error events:       %.6f\n", share
    printf "  background mean of the shipped curve:       %.6f  (1 on the sample it was fitted from)\n", cmean
    printf "  indel_probability = share / mean:           %.6f\n", fit
    print ""
    print "SANITY. #746 records an EXPECTATION, not a measurement: indel errors around 1e-6 to"
    print "1e-5 per base (the v3.4.0 CHANGELOG's ~1e-5) and substitutions around 1e-3. A rate far"
    print "outside it points first at the denominator or the support thresholds; at high depth the"
    print "error class admits more reads per position, so check the depth before the simulator."

    printf "key\tvalue\n" > f_out
    printf "indel_error_events\t%d\n", iev > f_out
    printf "sub_error_events\t%d\n", sev > f_out
    printf "indel_variant_events\t%d\n", ivar > f_out
    printf "sub_variant_events\t%d\n", svar > f_out
    printf "sub_events_pooled\t%d\n", spool > f_out
    printf "aligned_bases\t%d\n", bases > f_out
    printf "indel_error_rate\t%.8f\n", iev / bases > f_out
    printf "sub_error_rate\t%.8f\n", sev / bases > f_out
    printf "indel_share\t%.6f\n", share > f_out
    printf "curve_bg_mean\t%.6f\n", cmean > f_out
    printf "indel_probability\t%.6f\n", fit > f_out
}
