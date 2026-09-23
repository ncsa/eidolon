# CIGAR -> indel reference positions, with per-position support counts.
#
# A SEPARATE FILE, not inline in indel_context.sbatch, so the test exercises the same
# program the job runs. The reference-cursor semantics are the part most likely to be wrong
# and least likely to look wrong: getting them subtly off shifts every indel position by a
# few bases, which still produces a plausible table.
#
# The rule: M, D, N, = and X consume reference and advance the cursor. I, S and H do not.
# A deletion is recorded at the cursor BEFORE advancing past it, because that is where the
# junction is.
#
# LENGTH is emitted as a fourth column because the size of a slippage error is what the
# simulator has to reproduce, and it was being discarded here. `SequencingErrorModel`
# currently draws from [0.999, 0.001] over lengths [1, 2] -- NEAT2 placeholders that nobody
# has measured against real data. This column is the measurement.
#
# It is SIGNED: + for an insertion, - for a deletion. An inserted base and a deleted base
# at the same junction are different events, and an unsigned length would pool them.
#
# Where reads disagree at one junction, the MODAL signed length is reported. Tracked
# incrementally rather than by scanning the length table per position, so the cost stays
# linear in reads rather than positions x lengths.
#
# in:  <contig> <pos> <cigar>   (1-based POS, as SAM)
# out: <contig> <pos> <support> <modal signed length>
{
    pos = $2; cig = $3; n = ""
    for (i = 1; i <= length(cig); i++) {
        ch = substr(cig, i, 1)
        if (ch ~ /[0-9]/) { n = n ch; continue }
        len = n + 0; n = ""
        if (ch == "I" || ch == "D") {
            k = $1 SUBSEP pos
            c[k]++
            sl = (ch == "I") ? len : -len
            lc[k SUBSEP sl]++
            if (lc[k SUBSEP sl] > bestn[k]) { bestn[k] = lc[k SUBSEP sl]; bestl[k] = sl }
            if (ch == "D") pos += len
        } else if (ch ~ /[MN=X]/) {
            pos += len
        }
        # S and H are soft/hard clip: query-only, cursor unchanged.
    }
}
END { for (k in c) { split(k, a, SUBSEP); print a[1] "\t" a[2] "\t" c[k] "\t" bestl[k] } }
