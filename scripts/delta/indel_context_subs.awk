# CIGAR + SEQ -> reference positions where reads disagree with the reference, with
# per-position support counts. The substitution half of the indel_probability measurement
# (#746); indel_context_extract.awk is the indel half, and both feed the same support rule.
#
# A SEPARATE FILE, not inline in indel_context.sbatch, so the test exercises the same
# program the job runs. The two cursors are the part most likely to be wrong and least
# likely to look wrong: off by a base, every mismatch lands on a neighbour, and the counts
# still look plausible.
#
# The rules, matching indel_context_extract.awk:
#   M, =, X  consume query AND reference: each aligned base is compared.
#   I, S     consume query only.
#   D, N     consume reference only.
#   H, P     consume neither.
# An N on either side is a no-call, not an error, and is skipped. A base outside every
# loaded reference segment is skipped too: there is nothing to compare it against.
#
# in:  first file  <contig> <0-based start> <sequence>   (one reference segment per line)
#      then        <contig> <pos> <cigar> <seq>           (1-based POS, as SAM)
# out: <contig> <pos> <support>   on stdout, and `aligned_bases <n>` to f_bases.
FILENAME == f_ref {
    nseg[$1]++
    sstart[$1, nseg[$1]] = $2 + 1                 # 1-based first base
    sseq[$1, nseg[$1]]   = toupper($3)
    next
}
function refbase(c, p,    i, off) {
    for (i = 1; i <= nseg[c]; i++) {
        off = p - sstart[c, i] + 1
        if (off >= 1 && off <= length(sseq[c, i])) return substr(sseq[c, i], off, 1)
    }
    return ""
}
{
    c = $1; rpos = $2; cig = $3; seq = toupper($4); q = 1; n = ""
    for (i = 1; i <= length(cig); i++) {
        ch = substr(cig, i, 1)
        if (ch ~ /[0-9]/) { n = n ch; continue }
        len = n + 0; n = ""
        if (ch == "M" || ch == "=" || ch == "X") {
            for (j = 0; j < len; j++) {
                rb = substr(seq, q + j, 1)
                refb = refbase(c, rpos + j)
                compare(c, rpos + j, rb, refb)
            }
            q += len; rpos += len
        } else if (ch == "I") {
            q += len
        } else if (ch == "S") {
            q += len
        } else if (ch == "D" || ch == "N") {
            rpos += len
        }
        # H and P: neither cursor moves.
    }
}
function compare(c, p, rb, refb,    k) {
    if (rb == "N" || refb == "N" || refb == "") return
    bases++
    if (rb != refb) { k = c SUBSEP p; mm[k]++ }
}
END {
    for (k in mm) { split(k, a, SUBSEP); print a[1] "\t" a[2] "\t" mm[k] }
    if (f_bases != "") print "aligned_bases\t" (bases + 0) > f_bases
}
