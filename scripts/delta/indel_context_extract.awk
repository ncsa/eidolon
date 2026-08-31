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
# in:  <contig> <pos> <cigar>   (1-based POS, as SAM)
# out: <contig> <pos> <support>
{
    pos = $2; cig = $3; n = ""
    for (i = 1; i <= length(cig); i++) {
        ch = substr(cig, i, 1)
        if (ch ~ /[0-9]/) { n = n ch; continue }
        len = n + 0; n = ""
        if (ch == "I" || ch == "D") {
            c[$1 SUBSEP pos]++
            if (ch == "D") pos += len
        } else if (ch ~ /[MN=X]/) {
            pos += len
        }
        # S and H are soft/hard clip: query-only, cursor unchanged.
    }
}
END { for (k in c) { split(k, a, SUBSEP); print a[1] "\t" a[2] "\t" c[k] } }
