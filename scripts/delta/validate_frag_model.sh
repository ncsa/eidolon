#!/usr/bin/env bash
# Does a fragment model built from a real BAM reproduce that BAM's fragment distribution?
#
# Ground truth is measured with samtools + awk, NOT with eidolon, so this is a known-answer
# check rather than the code grading its own homework.
#
# COMPARE LIKE WITH LIKE. The builder trims outliers on purpose: a real BAM carries read
# pairs whose mates sit megabases apart on the same chromosome -- discordant and chimeric
# pairs, which are SV signal, not the library's fragment distribution. Measured on HCC1395
# normal chr20/21/22, leaving them in gives the "truth" a standard deviation of 147,396 and
# a skew of 228.5 over a support 61 million integers wide. Comparing a trimmed model
# against an untrimmed measurement measures the outliers, not the fit. So the truth is
# reported BOTH ways: raw, and restricted to the support the model actually covers, with
# the excluded mass named. The assertions use the restricted one.
#
# THE FILTER MUST MATCH THE BUILDER'S EXACTLY. `BamWalkFilter::for_frag_length()` takes
# paired, first-in-pair, mate mapped to the SAME reference, not secondary/supplementary,
# and MAPQ > 10 (the code is `if mq <= 10 { skip }`, so samtools needs -q 11, not -q 10).
# A difference in FILTERING would otherwise read as a difference in FITTING, which is the
# way this comparison is easiest to get quietly wrong.
#
# usage: validate_frag_model.sh <bam> <model.json.gz> [more_models.json.gz ...]
set -uo pipefail

BAM="${1:?usage: validate_frag_model.sh <bam> <model.json.gz> [...]}"; shift
[[ $# -ge 1 ]] || { echo "need at least one model file" >&2; exit 2; }
for t in samtools awk jq zcat; do
    command -v "$t" >/dev/null || { echo "FATAL: $t not found" >&2; exit 2; }
done

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

# ── moments + quantiles from a "<value> <weight>" table ─────────────────────
moments() {  # reads value/weight pairs on stdin, prints: mean sd skew p05 p50 p95 p99 bins gaps
    awk '
        { v[NR]=$1+0; w[NR]=$2+0; tot+=$2+0; if(NR==1||$1+0<lo) lo=$1+0; if(NR==1||$1+0>hi) hi=$1+0 }
        END {
            if (tot <= 0) { print "0 0 0 0 0 0 0 0 0"; exit }
            for (i=1;i<=NR;i++) mean += v[i]*w[i]; mean/=tot
            for (i=1;i<=NR;i++) var += (v[i]-mean)^2*w[i]; var/=tot
            sd = sqrt(var)
            if (sd > 0) { for (i=1;i<=NR;i++) sk += ((v[i]-mean)/sd)^3*w[i]; sk/=tot }
            # quantiles: walk in ascending value order (input is sorted by the caller)
            split("0.05 0.5 0.95 0.99", qs, " ")
            qi=1; acc=0
            for (i=1;i<=NR;i++) {
                acc += w[i]
                while (qi<=4 && acc >= tot*qs[qi]) { q[qi]=v[i]; qi++ }
            }
            for (j=qi;j<=4;j++) q[j]=v[NR]
            gaps = (hi-lo+1) - NR
            printf "%.4f %.4f %.4f %d %d %d %d %d %d\n", mean, sd, sk, q[1], q[2], q[3], q[4], NR, gaps
        }'
}

echo "Measuring ground truth from $BAM"
echo "  filter: -f 0x41 -F 0x90C -q 11 and RNEXT '=' (mirrors for_frag_length)"
samtools view -f 0x41 -F 0x90C -q 11 "$BAM" \
  | awk '$7=="=" && $9>0 { c[$9]++ } END { for (l in c) print l, c[l] }' \
  | sort -n > "$WORK/truth.tsv"
NPAIRS=$(awk '{s+=$2} END{print s+0}' "$WORK/truth.tsv")
# Rule 4: a metric over an unknown denominator is not a result.
[[ "${NPAIRS:-0}" -gt 0 ]] || { echo "FATAL: zero pairs passed the filter -- nothing to compare against" >&2; exit 1; }
read -r TMEAN TSD TSKEW TP05 TP50 TP95 TP99 TBINS TGAPS < <(moments < "$WORK/truth.tsv")
printf "  %'d pairs, %d distinct lengths, %d integer gaps\n" "$NPAIRS" "$TBINS" "$TGAPS"
# Not decoration: a wildly inflated sd here is the tell that the raw row is dominated by
# discordant pairs, and that only the restricted rows below are a comparison.
awk -v sd="$TSD" -v m="$TMEAN" 'BEGIN{ if (sd > 3*m)
    printf "  NOTE: raw sd (%.0f) far exceeds the mean (%.0f) -- discordant pairs dominate\n         the raw moments. Compare the restricted rows.\n", sd, m }'
echo

HDR=$(printf "%-24s %10s %9s %9s %6s %6s %6s %6s %7s %6s" \
      "" mean sd skew p05 p50 p95 p99 bins gaps)
echo "$HDR"; printf '%*s\n' "${#HDR}" '' | tr ' ' -
printf "%-24s %10.2f %9.2f %9.3f %6d %6d %6d %6d %7d %6d\n" \
       "REAL BAM (raw)" "$TMEAN" "$TSD" "$TSKEW" "$TP05" "$TP50" "$TP95" "$TP99" "$TBINS" "$TGAPS"

FAIL=0
declare -a REPORT
for MODEL in "$@"; do
    KIND=$(zcat "$MODEL" | jq -r 'if has("Normal") then "Normal" else "Discrete" end')
    if [[ "$KIND" == "Normal" ]]; then
        MEAN=$(zcat "$MODEL" | jq -r '.Normal.mean')
        SD=$(zcat "$MODEL" | jq -r '.Normal.st_dev')
        # A Normal has skew 0 and closed-form quantiles; state them rather than sampling.
        read -r MMEAN MSD MSKEW MP05 MP50 MP95 MP99 MBINS MGAPS < <(awk -v m="$MEAN" -v s="$SD" '
            BEGIN { printf "%.4f %.4f %.4f %d %d %d %d %d %d\n",
                    m, s, 0, m-1.6449*s, m, m+1.6449*s, m+2.3263*s, 0, 0 }')
    else
        # Stored weights are CUMULATIVE; recover per-bin mass before taking moments.
        zcat "$MODEL" | jq -r '.Discrete.distribution
              | [.values, .weights] | transpose | .[] | "\(.[0]) \(.[1])"' \
          | awk '{ v=$1; c=$2; print v, c-prev; prev=c }' > "$WORK/model.tsv"
        read -r MMEAN MSD MSKEW MP05 MP50 MP95 MP99 MBINS MGAPS < <(moments < "$WORK/model.tsv")
    fi
    # The truth, restricted to the support this model covers. Same data, same filter --
    # only the outliers the builder deliberately dropped are excluded, and how much mass
    # that was gets printed rather than assumed (rule 4).
    if [[ "$KIND" == "Normal" ]]; then
        RLO=$(awk -v m="$MEAN" -v s="$SD" 'BEGIN{ v=int(m-4*s); print (v<1?1:v) }')
        RHI=$(awk -v m="$MEAN" -v s="$SD" 'BEGIN{ print int(m+4*s) }')
    else
        RLO=$(head -1 "$WORK/model.tsv" | awk '{print $1}')
        RHI=$(tail -1 "$WORK/model.tsv" | awk '{print $1}')
    fi
    awk -v lo="$RLO" -v hi="$RHI" '$1>=lo && $1<=hi' "$WORK/truth.tsv" > "$WORK/truth_r.tsv"
    read -r RMEAN RSD RSKEW RP05 RP50 RP95 RP99 RBINS RGAPS < <(moments < "$WORK/truth_r.tsv")
    RPAIRS=$(awk '{s+=$2} END{print s+0}' "$WORK/truth_r.tsv")
    EXCL=$(awk -v a="$RPAIRS" -v b="$NPAIRS" 'BEGIN{printf "%.4f", (b-a)/b*100}')

    printf "%-24s %10.2f %9.2f %9.3f %6d %6d %6d %6d %7d %6d\n" \
           "REAL BAM in $KIND range" "$RMEAN" "$RSD" "$RSKEW" "$RP05" "$RP50" "$RP95" "$RP99" "$RBINS" "$RGAPS"
    printf "%-24s %10.2f %9.2f %9.3f %6d %6d %6d %6d %7d %6d\n" \
           "model: $KIND" "$MMEAN" "$MSD" "$MSKEW" "$MP05" "$MP50" "$MP95" "$MP99" "$MBINS" "$MGAPS"
    printf "%-24s  support %d-%d, excludes %s%% of pairs as outliers\n" "" "$RLO" "$RHI" "$EXCL"

    REPORT+=("$(awk -v k="$KIND" -v mm="$MMEAN" -v ms="$MSD" -v mk="$MSKEW" -v m9="$MP99" \
                    -v tm="$RMEAN" -v ts="$RSD" -v tk="$RSKEW" -v t9="$RP99" -v g="$MGAPS" '
        BEGIN {
            dm=(mm-tm)/tm*100; if(dm<0)dm=-dm
            ds=(ms-ts)/ts*100; if(ds<0)ds=-ds
            dk=mk-tk;          if(dk<0)dk=-dk
            d9=(m9-t9)/t9*100; if(d9<0)d9=-d9
            printf "%s|%.2f|%.2f|%.3f|%.2f|%d", k, dm, ds, dk, d9, g
        }')")
done

echo; echo "Agreement with the real BAM, over each model's own support:"
for r in "${REPORT[@]}"; do
    IFS='|' read -r k dm ds dk d9 g <<<"$r"
    printf "  %-10s mean %6s%%   sd %6s%%   skew off by %6s   p99 %6s%%\n" "$k" "$dm" "$ds" "$dk" "$d9"
    # Only the discrete model claims to reproduce the shape. A Normal is EXPECTED to miss
    # the skew -- asserting against it would just re-measure the thing we already know.
    [[ "$k" == "Discrete" ]] || continue
    awk -v v="$dm" 'BEGIN{exit !(v>2.0)}'  && { echo "      FAIL mean off by ${dm}% (>2%)";  FAIL=1; }
    awk -v v="$ds" 'BEGIN{exit !(v>5.0)}'  && { echo "      FAIL sd off by ${ds}% (>5%)";    FAIL=1; }
    awk -v v="$dk" 'BEGIN{exit !(v>0.15)}' && { echo "      FAIL skew off by ${dk} (>0.15)"; FAIL=1; }
    awk -v v="$d9" 'BEGIN{exit !(v>5.0)}'  && { echo "      FAIL p99 off by ${d9}% (>5%)";   FAIL=1; }
    [[ "$g" -gt 0 ]] && { echo "      FAIL the built model still has $g gaps"; FAIL=1; }
done

echo
if [[ "$FAIL" -eq 0 ]]; then
    echo "VERDICT: the discrete model reproduces the distribution it was built from"
else
    echo "VERDICT: DISAGREEMENT ABOVE -- the built model does not match its input"
fi
exit "$FAIL"
