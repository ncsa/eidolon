# Summarise a realism panel TSV: per metric, the median across loci and the min-max range,
# for each of REAL and SIMULATED, plus the ratio between the medians.
#
# A SEPARATE FILE, not inlined in realism_panel.sbatch, so the job and its test run the exact
# same program. An earlier version lived inside the sbatch and the test extracted it with sed;
# the extraction silently produced nothing and the assertions passed against empty output.
# Same helper both sides, per CLAUDE.md.
#
# POSIX AWK ONLY. The first version used `vals[key][i]` — true multidimensional arrays, which
# are a gawk extension. Under mawk or BWK awk that is a syntax error, so the program emitted
# nothing at all and the job would have printed an empty table on Delta while exiting 0.
# Keys are SUBSEP-joined instead, which every awk supports.
#
# MEDIAN, NOT MEAN. Real dispersion varies several-fold between loci on one chromosome:
# measured on chr22, VMR reads 5.51, 6.87, 7.85, 8.88 at four loci and 36.10 at a fifth. The
# mean of those is 13.04 — a value no locus has, and a threshold derived from it would be
# ~66% high. The median lands in the cluster; the range beside it says how much of a gap
# belongs to the data and how much to whichever locus was picked.
#
# Usage:  awk -F'\t' -f realism_summarise.awk panel.tsv

NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    next
}

{
    lab = $1
    for (m in col) {
        if (m == "label" || m == "contig" || m == "start" || m == "end") continue
        k = lab SUBSEP m
        cnt[k]++
        vals[k, cnt[k]] = $(col[m]) + 0
    }
}

# Median of the cnt[k] values stored under key k. Sorts a local copy: an insertion sort is
# ample for the handful of loci a panel run measures, and avoids depending on gawk's asort.
function median(k, n,    i, j, t, a, mid) {
    for (i = 1; i <= n; i++) a[i] = vals[k, i]
    for (i = 1; i < n; i++)
        for (j = i + 1; j <= n; j++)
            if (a[j] < a[i]) { t = a[i]; a[i] = a[j]; a[j] = t }
    mid = a[int((n + 1) / 2)]
    lo = a[1]
    hi = a[n]
    return mid
}

END {
    split("cand_per_mb improper_pct clip_pct mapq0_pct depth_mean depth_vmr depth_excess depth_acf ins_sd ins_skew ins_p99", order, " ")
    printf "%-14s %26s %26s %10s\n", "metric", "REAL (median [min-max])", "SIM (median [min-max])", "gap"
    for (oi = 1; oi in order; oi++) {
        m = order[oi]
        printf "%-14s", m
        have_real = 0
        have_sim = 0
        for (li = 1; li <= 2; li++) {
            lab = (li == 1) ? "REAL" : "SIMULATED"
            k = lab SUBSEP m
            n = cnt[k]
            if (!n) { printf " %26s", "-"; continue }
            mid = median(k, n)
            if (li == 1) { real_med = mid; have_real = 1 } else { sim_med = mid; have_sim = 1 }
            printf " %26s", sprintf("%.4g [%.4g-%.4g]", mid, lo, hi)
        }
        # A ratio against zero is not a large number, it is a different statement: the
        # simulator produced NONE of this artifact. Printing a finite value there is how
        # "no background at all" gets read as "a small gap".
        # A ratio is only meaningful when both medians share a sign. Autocorrelation can go
        # negative — eidolon's reads -0.002 against real data's +0.8 — and 0.8 / -0.002
        # prints "-400x", which describes nothing. Those are reported as a difference
        # instead, which is what the two numbers actually differ by.
        if (have_real && have_sim && sim_med == 0 && real_med != 0) printf " %10s", "inf"
        else if (have_real && have_sim && (real_med < 0 || sim_med < 0))
            printf " %9.3f%s", real_med - sim_med, "d"
        else if (have_real && have_sim && sim_med != 0) printf " %9.1fx", real_med / sim_med
        else printf " %10s", "-"
        print ""
    }
}
