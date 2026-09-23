# Measure homopolymer run length and implied indel size FROM THE READ BASES, ignoring the
# aligner's interpretation of them.
#
# Input:  SAM records (no header) for one locus, e.g. `samtools view -F 0x904 bam region`.
# Output: one tab-separated record per measurement —
#           <locus>  RUN     <read_run_len - reference_run_len>
#           <locus>  OFF     <3prime|5prime>  <implied indel: + deletion, - insertion>
#           <locus>  SKIP    <reason>
#
# WHY NOT THE CIGAR. The obvious way to measure slippage at a homopolymer is to read the
# `I`/`D` operations. That is biased by construction: those are the events bwa-mem2 chose to
# GAP, and at long homopolymers 41% of reads are soft-clipped instead (#692), so the events
# that matter are exactly the ones a CIGAR-based measurement cannot see. The indel-error
# length distribution now shipped in `sequencing_error_model.rs` was fitted that way and
# inherits the same blind spot.
#
# TWO FACTS MAKE THIS WORK. `SEQ` in a BAM contains soft-clipped bases, and it is stored
# forward-strand relative to the reference. So a clipped read still carries its full sequence,
# correctly oriented, and the alignment is used only to FIND the read — never to interpret it.
#
# MEASUREMENT 1 — run length. Locate the reference flank as an exact substring of the read,
# then count the homopolymer base from where it ends. The run must terminate inside the read,
# or the length is a lower bound rather than a measurement.
#
# MEASUREMENT 2 — implied indel. For a read soft-clipped at a long homopolymer, take the tail
# of the clipped sequence (past the run, so it is unique), find where it actually sits in the
# reference, and compare against where it would sit with no indel. That offset IS the indel
# size, recovered without the aligner's cooperation:
#
#   0        the clip would have mapped with no indel -- the aligner clipped over ambiguity
#   +/- 1-3  slippage the aligner declined to gap
#   +/- 20+  a genuine large indel (variant placement, #378)
#   absent   the clipped sequence is not downstream reference at all
#
# POSIX AWK ONLY. No gawk extensions.
#
# Env (via -v): locus base lref lflank rflank win winstart tail minclip

function count_occ(hay, needle,   n, p, s) {
    n = 0
    s = hay
    while ((p = index(s, needle)) > 0) {
        n++
        s = substr(s, p + 1)
    }
    return n
}

# Reference bases consumed: M, D, N, =, X. Insertions and clips consume query only.
function ref_span(cig,   i, c, num, span) {
    span = 0
    num = ""
    for (i = 1; i <= length(cig); i++) {
        c = substr(cig, i, 1)
        if (c >= "0" && c <= "9") num = num c
        else {
            if (c == "M" || c == "D" || c == "N" || c == "=" || c == "X") span += num + 0
            num = ""
        }
    }
    return span
}

function lead_clip(cig,   i, c, num) {
    num = ""
    for (i = 1; i <= length(cig); i++) {
        c = substr(cig, i, 1)
        if (c >= "0" && c <= "9") num = num c
        else return (c == "S") ? num + 0 : 0
    }
    return 0
}

function trail_clip(cig,   j, c, num) {
    if (substr(cig, length(cig), 1) != "S") return 0
    j = length(cig) - 1
    num = ""
    while (j >= 1) {
        c = substr(cig, j, 1)
        if (c < "0" || c > "9") break
        num = c num
        j--
    }
    return num + 0
}

function skip(reason) {
    printf "%s\tSKIP\t%s\n", locus, reason
}

BEGIN {
    if (tail == "") tail = 15
    if (minclip == "") minclip = 20
    FS = "\t"
}

{
    pos = $4 + 0
    cig = $6
    seq = toupper($10)
    if (cig == "*" || seq == "*" || seq == "") { skip("no sequence"); next }

    # ── Measurement 1: run length from the read's own bases ──────────────────
    #
    # A flank appearing twice in one read means the anchor is ambiguous for THIS read even
    # though it is unique in the reference window; take the measurement from the other side
    # rather than guess.
    got = 0
    if (count_occ(seq, lflank) == 1) {
        p = index(seq, lflank) + length(lflank)
        n = 0
        while (p + n <= length(seq) && substr(seq, p + n, 1) == base) n++
        if (p + n <= length(seq)) {
            printf "%s\tRUN\t%d\n", locus, n - lref
            got = 1
        } else {
            skip("run reaches the read end (left anchor)")
            got = 1
        }
    }
    if (!got && count_occ(seq, rflank) == 1) {
        p = index(seq, rflank) - 1
        n = 0
        while (p - n >= 1 && substr(seq, p - n, 1) == base) n++
        if (p - n >= 1) {
            printf "%s\tRUN\t%d\n", locus, n - lref
            got = 1
        } else {
            skip("run reaches the read start (right anchor)")
            got = 1
        }
    }
    if (!got) skip("neither flank anchors in this read")

    # ── Measurement 2: implied indel from a soft clip ────────────────────────
    span = ref_span(cig)
    endpos = pos + span - 1

    tc = trail_clip(cig)
    if (tc >= minclip && tc >= tail) {
        clip = substr(seq, length(seq) - tc + 1, tc)
        t = substr(clip, tc - tail + 1, tail)
        occ = count_occ(win, t)
        if (occ == 1) {
            # With no indel the clip would continue straight on from the alignment, so its
            # last `tail` bases would sit at endpos + tc - tail + 1.
            printf "%s\tOFF\t3prime\t%d\n", locus,
                (winstart + index(win, t) - 1) - (endpos + tc - tail + 1)
        } else if (occ == 0) skip("3' clip tail is not in the reference window")
        else skip("3' clip tail is ambiguous in the reference window")
    }

    lc = lead_clip(cig)
    if (lc >= minclip && lc >= tail) {
        t = substr(seq, 1, tail)
        occ = count_occ(win, t)
        if (occ == 1) {
            # Mirrored: a 5' clip would begin at pos - lc. Reference skipped between the clip
            # and the alignment pushes its true position LEFT, so the sign is inverted to keep
            # "+ means deletion" consistent with the 3' case.
            printf "%s\tOFF\t5prime\t%d\n", locus,
                (pos - lc) - (winstart + index(win, t) - 1)
        } else if (occ == 0) skip("5' clip tail is not in the reference window")
        else skip("5' clip tail is ambiguous in the reference window")
    }
}
