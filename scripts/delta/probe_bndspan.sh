#!/usr/bin/env bash
# Rule 0 for the BND-span representation bridge.
#
# `selftest_truvari` compares the derived spans to THEMSELVES. That matches trivially and
# proves nothing about the actual question, which is cross-representation:
#
#     does truvari -t match a derived <BND_SPAN> to the <DUP:TANDEM>/<DEL>/<INV>
#     a caller emits for the same junction?
#
# That match is the entire premise of BNDspan scoring, it is the reason build_bnd_spans
# exists, and nothing in the pipeline tests it. If it does not hold, BNDspan recall is 0
# for reasons that are ours rather than the caller's — which is the mistake this whole
# line of work started from.
#
# Known-answer, both directions, because a configuration loose enough to match anything
# is as useless as one that matches nothing:
#   * a synthetic caller record at the span's own coordinates MUST match
#   * the same record shifted far away MUST NOT
#
# Standalone use against an existing run directory (seconds, no simulation). Either
# point it at derived spans, or at the BND truth and let it derive them with the real
# build_bnd_spans extracted from sv_pipeline.sbatch — which also checks that the ALT
# parser works on actual eidolon output rather than on a fixture:
#   scripts/delta/probe_bndspan.sh --spans "$OUTDIR/truth_sv_BNDspan.vcf.gz"
#   scripts/delta/probe_bndspan.sh --bnd   "$OUTDIR/truth_sv_BND.vcf.gz"
#
# Exit 0 = the bridge works as assumed. Exit 1 = it does not, and any BNDspan number
# produced with this configuration is meaningless.

set -euo pipefail

SPANS=""
BND_TRUTH=""
REFERENCE="${REFERENCE:-}"
# Offset applied to the "must match" case. Real callers do not hit a breakpoint exactly;
# job 20636298 had Manta 1 bp off. Non-zero on purpose — an exact-coordinate-only match
# would pass here and still fail on real data.
NEAR_OFFSET="${NEAR_OFFSET:-1}"
# Distance for the negative control. Far beyond any breakpoint imprecision.
FAR_OFFSET="${FAR_OFFSET:-50000}"
WORK="${WORK:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --spans)     SPANS="$2"; shift 2 ;;
        --bnd)       BND_TRUTH="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --work)      WORK="$2"; shift 2 ;;
        -h|--help)   sed -n '2,30p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$SPANS" && -z "$BND_TRUTH" ]]; then
    echo "ERROR: one of --spans <truth_sv_BNDspan.vcf.gz> or --bnd <truth_sv_BND.vcf.gz>" >&2
    exit 2
fi

if [[ -z "$WORK" ]]; then
    WORK="$(mktemp -d)"
    trap 'rm -rf "$WORK"' EXIT
fi
mkdir -p "$WORK"

# Derive the spans here when given the BND truth, using the REAL function rather than a
# reimplementation — so this also exercises the ALT parser against actual eidolon output.
if [[ -z "$SPANS" ]]; then
    [[ -f "$BND_TRUTH" ]] || { echo "ERROR: no such BND truth: $BND_TRUTH" >&2; exit 2; }
    SBATCH="$(cd "$(dirname "$0")" && pwd)/sv_pipeline.sbatch"
    [[ -f "$SBATCH" ]] || { echo "ERROR: cannot find $SBATCH to extract build_bnd_spans" >&2; exit 2; }
    OUTDIR="$WORK"          # build_bnd_spans writes its temporaries under $OUTDIR
    eval "$(awk '/^build_bnd_spans\(\) \{/,/^\}$/' "$SBATCH")"
    SPANS="$WORK/derived_BNDspan.vcf.gz"
    span_rc=0
    build_bnd_spans "$BND_TRUTH" "$SPANS" || span_rc=$?
    case "$span_rc" in
        0) ;;
        2) echo "No intra-contig junctions to span — the bridge is untestable here, and" >&2
           echo "  BNDspan scoring would legitimately not run." >&2
           exit 0 ;;
        *) echo "ERROR: span derivation failed (rc=$span_rc) — cannot probe the bridge." >&2
           exit 1 ;;
    esac
fi

[[ -f "$SPANS" ]] || { echo "ERROR: no such spans VCF: $SPANS" >&2; exit 2; }

n_spans=$(bcftools view -H "$SPANS" 2>/dev/null | wc -l)
if [[ "$n_spans" -eq 0 ]]; then
    echo "ERROR: $SPANS contains no records — nothing to probe." >&2
    exit 1
fi

# One truth record, so the arithmetic is unambiguous: TP=1/FN=0 or TP=0/FN=1.
TRUTH1="$WORK/probe_truth.vcf.gz"
{ bcftools view -h "$SPANS"; bcftools view -H "$SPANS" | head -1; } \
    | bcftools sort -O z -o "$TRUTH1" 2>/dev/null
bcftools index -f -t "$TRUTH1"

read -r P_CHROM P_POS P_END < <(bcftools query -f '%CHROM\t%POS\t%INFO/END\n' "$TRUTH1")
echo "probe junction: $P_CHROM:$P_POS-$P_END  (span $((P_END - P_POS)) bp, of $n_spans derived)"

# A caller's representation of that same junction: a tandem duplication over the span.
# This is what Manta emitted in job 20636298 for a truth junction — POS 1 bp off, END
# exactly on the mate — and what truvari refused to match as a BND.
mk_dup() {  # <out.vcf.gz> <offset>
    local out="$1" off="$2" hdr="$WORK/.hdr"
    bcftools view -h "$SPANS" | grep -vE '^#CHROM' > "$hdr"
    grep -q '^##ALT=<ID=DUP,' "$hdr" \
        || echo '##ALT=<ID=DUP,Description="Duplication">' >> "$hdr"
    { cat "$hdr"
      bcftools view -h "$SPANS" | grep -E '^#CHROM'
      awk -v c="$P_CHROM" -v p="$P_POS" -v e="$P_END" -v o="$off" 'BEGIN{OFS="\t"
          print c, p+o, ".", "N", "<DUP>", 60, "PASS",
                "SVTYPE=DUP;END=" (e+o) ";SVLEN=" (e-p), "GT", "0/1" }'
    } | bcftools sort -O z -o "$out" 2>/dev/null
    bcftools index -f -t "$out"
}

# -t (--typeignore) is the point: the caller says DUP, the truth says BND. --pctsize/
# --pctseq 0 because a span and a duplication are the same interval but not the same
# sequence claim. Matching the settings score_caller uses for BNDspan.
run_probe() {  # <comp.vcf.gz> <out-dir>
    truvari bench -b "$TRUTH1" -c "$1" -o "$2" \
        --passonly -t --pctsize 0 --pctseq 0 --sizemax 100000000 \
        ${REFERENCE:+-f "$REFERENCE"} >/dev/null 2>&1 || true
    [[ -f "$2/summary.json" ]] || { echo "-1"; return; }
    python3 -c "import json;d=json.load(open('$2/summary.json'));print(d.get('TP-base') or 0)"
}

NEAR="$WORK/comp_near.vcf.gz"; FAR="$WORK/comp_far.vcf.gz"
mk_dup "$NEAR" "$NEAR_OFFSET"
mk_dup "$FAR"  "$FAR_OFFSET"

rm -rf "$WORK/out_near" "$WORK/out_far"
tp_near=$(run_probe "$NEAR" "$WORK/out_near")
tp_far=$(run_probe "$FAR"  "$WORK/out_far")

rc=0
echo "  positive control: <DUP> at +${NEAR_OFFSET}bp      -> TP=$tp_near  (must be 1)"
echo "  negative control: <DUP> at +${FAR_OFFSET}bp  -> TP=$tp_far  (must be 0)"

if [[ "$tp_near" == "-1" || "$tp_far" == "-1" ]]; then
    echo "ERROR: truvari produced no summary.json — the probe could not run." >&2
    exit 1
fi
if [[ "$tp_near" -ne 1 ]]; then
    echo "FAIL: truvari did not match a caller's <DUP> to the derived <BND_SPAN> for the" >&2
    echo "  SAME junction, ${NEAR_OFFSET}bp apart. The representation bridge does not hold, so" >&2
    echo "  BNDspan recall measures our comparison rather than the caller. Do not report it." >&2
    rc=1
fi
if [[ "$tp_far" -ne 0 ]]; then
    echo "FAIL: truvari matched a <DUP> ${FAR_OFFSET}bp away from the junction. The" >&2
    echo "  configuration is matching on position far more loosely than intended, so a" >&2
    echo "  BNDspan 'match' does not mean the caller found this junction." >&2
    rc=1
fi
[[ "$rc" -eq 0 ]] && echo "  BNDspan bridge: PASS (matches the same junction, rejects a distant one)"
exit "$rc"
