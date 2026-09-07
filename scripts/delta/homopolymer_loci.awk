# Turn padded reference windows into a homopolymer locus table.
#
# Input:  multi-FASTA from `samtools faidx -r <regions>`, headers `>contig:start-end`
#         (1-based inclusive, as samtools writes them).
# Output: one tab-separated locus per line —
#           contig  run_start  run_len  base  left_flank  right_flank  win_start  win_seq
#         where run_start and win_start are 1-based reference positions.
#
# A SEPARATE FILE, not inlined in the driver, so the job and its test run the same program.
# `realism_summarise.awk` records what happens otherwise: an extracted copy silently produced
# nothing and the assertions passed against empty output.
#
# POSIX AWK ONLY — no gawk extensions (no `asort`, no true multidimensional arrays). The
# harness runs under whatever awk Delta's environment supplies.
#
# Env (via -v):
#   pad      bases of padding on each side of the input interval (default 400)
#   flank    anchor length (default 12)
#   min_run  shortest homopolymer to report (default 12)
#
# WHY THE FLANK MUST BE UNIQUE IN THE WINDOW. The measurement locates the run in a read by
# finding the flank as a substring. If the flank occurs more than once in the surrounding
# reference, a read could anchor on the wrong copy and report a run length from somewhere
# else — silently, and biased in an unknown direction. Such loci are dropped and counted.

function count_occ(hay, needle,   n, p, s) {
    n = 0
    s = hay
    while ((p = index(s, needle)) > 0) {
        n++
        s = substr(s, p + 1)
    }
    return n
}

# Longest single-base run overlapping [lo, hi] (1-based, inclusive) of `s`.
# Results land in BEST_LEN / BEST_START / BEST_BASE because awk cannot return a tuple.
# N runs are skipped: a run of N is an assembly gap, not a homopolymer, and reads there
# carry no usable sequence.
function longest_run(s, lo, hi,   i, c, cur, curb, n, st) {
    BEST_LEN = 0
    BEST_START = 0
    BEST_BASE = ""
    cur = 0
    curb = ""
    n = length(s)
    for (i = 1; i <= n; i++) {
        c = substr(s, i, 1)
        if (c == curb) cur++
        else { curb = c; cur = 1 }
        if (curb == "N") continue
        st = i - cur + 1
        # Overlaps the core interval, and is the longest such run seen.
        if (cur > BEST_LEN && st <= hi && i >= lo) {
            BEST_LEN = cur
            BEST_START = st
            BEST_BASE = curb
        }
    }
}

function emit(   core_lo, core_hi, rs, re, lf, rf, run_abs) {
    if (name == "" || seq == "") return
    # Header `contig:start-end`; split on the LAST colon so a contig name containing one
    # (they exist in draft assemblies) does not silently shift every coordinate.
    p = length(name)
    while (p > 0 && substr(name, p, 1) != ":") p--
    if (p == 0) { skipped["header without coordinates"]++; return }
    contig = substr(name, 1, p - 1)
    coords = substr(name, p + 1)
    d = index(coords, "-")
    if (d == 0) { skipped["header without a range"]++; return }
    win_start = substr(coords, 1, d - 1) + 0
    win_end = substr(coords, d + 1) + 0

    # The core is the caller's original interval; the pad exists so flanks and the clip
    # search window have reference to work with.
    core_lo = pad + 1
    core_hi = length(seq) - pad
    if (core_hi < core_lo) { skipped["window smaller than its padding"]++; return }

    longest_run(seq, core_lo, core_hi)
    if (BEST_LEN < min_run) { skipped["no run >= min_run"]++; return }

    rs = BEST_START
    re = BEST_START + BEST_LEN - 1
    if (rs - flank < 1 || re + flank > length(seq)) { skipped["run too close to the window edge"]++; return }

    lf = substr(seq, rs - flank, flank)
    rf = substr(seq, re + 1, flank)
    # Anchors must be unambiguous within the window, or a read can latch onto the wrong copy.
    if (count_occ(seq, lf) != 1) { skipped["left flank not unique in the window"]++; return }
    if (count_occ(seq, rf) != 1) { skipped["right flank not unique in the window"]++; return }

    run_abs = win_start + rs - 1
    printf "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\n",
        contig, run_abs, BEST_LEN, BEST_BASE, lf, rf, win_start, seq
    kept++
}

BEGIN {
    if (pad == "") pad = 400
    if (flank == "") flank = 12
    if (min_run == "") min_run = 12
    name = ""
    seq = ""
    kept = 0
}

/^>/ {
    emit()
    name = substr($0, 2)
    seq = ""
    next
}

{ seq = seq toupper($0) }

END {
    emit()
    # Rule 4: a locus table is a denominator, and a denominator with silent losses is not
    # one. Every rejection reason is reported, on stderr so it cannot pollute the table.
    printf("loci kept: %d\n", kept) > "/dev/stderr"
    for (k in skipped) printf("  dropped (%s): %d\n", k, skipped[k]) > "/dev/stderr"
}
