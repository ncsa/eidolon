#!/usr/bin/env bash
# Draw one genome's worth of germline structural variants from a population sites VCF.
#
# WHY. `cand_per_mb` compares a real genome against a simulation that plants no structural
# variants at all (#684). An SV breakpoint is the textbook source of a candidate site --
# reads spanning it clip at one fixed reference position, so they agree by construction and
# the clip clears the 20 bp threshold. Neither sequencing-error indels (1-2 bp) nor variant
# indel placement can do that (#672).
#
# WHY DRAW RATHER THAN GENERATE. `sv_rate_scale` generates de novo SVs at uniformly random
# positions. Real SVs sit in segmental duplications and repeat-rich sequence, which is where
# mappability problems live -- random placement cannot reproduce `mapq0_pct`, real placement
# might. Drawing on allele frequency also produces a realistic per-genome complement without
# anyone having to convert a cohort-union site density into a per-genome rate, which is the
# calibration trap recorded in lib_realism_config.sh.
#
# The output goes to `input_vcf`, the path sv_model_defaults.rs already recommends for
# specific rearrangements: records are preserved verbatim and junction reads are generated.
#
# POPULATION AF IS NOT EMITTED AS INFO/AF. eidolon parses INFO/AF into a variant's
# `allele_fraction` -- the share of reads carrying the alt. For a germline heterozygote that
# is 0.5, not the population frequency; writing gnomAD's AF there would render a common
# variant at 2% VAF and make it invisible. Frequency is recorded as INFO/POP_AF, and dosage
# is left to GT.
#
# Usage:
#   tools/draw_gnomad_sv_vcf.sh <sites.vcf.gz> <out.vcf> [contig ...]
#
# Env:
#   SEED       integer, default 1        reproducible draw
#   TYPES      default DEL,DUP,INS,INV,CNV   SVTYPEs to keep
#   MIN_AF     default 0                 skip sites below this frequency
#   MAX_LEN    default 1000000           skip SVs longer than this
set -euo pipefail

VCF="${1:?usage: draw_gnomad_sv_vcf.sh <sites.vcf.gz> <out.vcf> [contig ...]}"
OUT="${2:?usage: draw_gnomad_sv_vcf.sh <sites.vcf.gz> <out.vcf> [contig ...]}"
shift 2
CONTIGS=("$@")

SEED="${SEED:-1}"
TYPES="${TYPES:-DEL,DUP,INS,INV,CNV}"
MIN_AF="${MIN_AF:-0}"
MAX_LEN="${MAX_LEN:-1000000}"

command -v bcftools >/dev/null 2>&1 || {
    echo "ERROR: bcftools not found. On Delta: module load bcftools, or conda activate bioinf" >&2
    exit 1
}
[[ -f "$VCF" ]] || { echo "ERROR: sites VCF not found: $VCF" >&2; exit 1; }

REGION_ARG=()
if [[ ${#CONTIGS[@]} -gt 0 ]]; then
    REGION_ARG=(-r "$(IFS=,; echo "${CONTIGS[*]}")")
    echo "Contigs:  ${CONTIGS[*]}"
else
    echo "Contigs:  all"
fi
echo "Source:   $VCF"
echo "Seed:     $SEED   types: $TYPES   min_af: $MIN_AF   max_len: $MAX_LEN"

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# BND/CPX/CTX are deliberately absent from the default TYPES. sv_model_defaults.rs records
# that germline BND in gnomAD is dominated by mapping artifacts in repetitive regions -- one
# record reaches AF 0.114, and a translocation carried by 11% of humans would be a landmark
# cytogenetic finding rather than a variant.
echo "[1/3] streaming sites ..."
bcftools query "${REGION_ARG[@]}" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%INFO/AF\n' \
    "$VCF" > "$TMP/sites.tsv"
NSITES=$(wc -l < "$TMP/sites.tsv")
echo "      $NSITES sites read"
[[ "$NSITES" -gt 0 ]] || { echo "ERROR: no sites read -- wrong contig names for this VCF?" >&2; exit 1; }

echo "[2/3] drawing genotypes ..."
awk -v OFS='\t' -v seed="$SEED" -v types="$TYPES" -v min_af="$MIN_AF" -v max_len="$MAX_LEN" '
# Park-Miller rather than awk`s rand(): srand() seeding differs between gawk and mawk, and a
# draw that is not reproducible across implementations is not a fixture anyone can re-run.
#
# The multiplier matters. A glibc-style LCG (a=1103515245, m=2^31) computes a*x up to
# ~2.4e18, past the 2^53 where awk`s doubles stop being exact, so it silently loses bits --
# measured as het at 463 against an HWE expectation of 500 on a 1000-site AF=0.5 fixture,
# with both homozygote classes enriched. 16807 * (2^31-2) is 3.6e13 and stays exact.
function urand(   ) { lcg = (16807 * lcg) % 2147483647; return lcg / 2147483647 }
BEGIN {
    lcg = (seed % 2147483646) + 1     # Park-Miller state must be in [1, m-1]
    n = split(types, t, ","); for (i = 1; i <= n; i++) keep[t[i]] = 1
    for (i = 0; i < 8; i++) urand()   # discard the first few; low seeds start cold
}
{
    chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; filt=$6; svtype=$7; end=$8; svlen=$9; af=$10
    if (filt != "PASS" && filt != ".")      { skipped["not PASS"]++; next }
    if (!(svtype in keep))                  { skipped["type " svtype]++; next }
    if (af == "." || af + 0 <= 0)           { skipped["no AF"]++; next }
    if (af + 0 < min_af + 0)                { skipped["below min_af"]++; next }
    # Length: SVLEN when present, else END-POS. gnomAD gives SVLEN as a magnitude for DEL.
    L = (svlen != "." ? (svlen < 0 ? -svlen : svlen) : (end != "." ? end - pos : 0))
    if (L > max_len + 0)                    { skipped["over max_len"]++; next }

    # Hardy-Weinberg on the population frequency: hom-ref (1-af)^2, het 2af(1-af), hom-alt
    # af^2. A site the genome does not carry is simply not emitted.
    a = af + 0; u = urand()
    p_homref = (1-a)*(1-a); p_het = 2*a*(1-a)
    if (u < p_homref)                        { nref++; next }
    gt = (u < p_homref + p_het) ? "0/1" : "1/1"
    if (gt == "0/1") nhet++; else nhom++

    # INFO carries what eidolon reads (SVTYPE/END/SVLEN) plus POP_AF for provenance.
    info = "SVTYPE=" svtype
    if (end   != ".") info = info ";END=" end
    if (svlen != ".") info = info ";SVLEN=" svlen
    info = info ";POP_AF=" af
    if (id == ".") id = svtype "_" chrom "_" pos
    print chrom, pos, id, ref, alt, ".", "PASS", info, "GT", gt > "/dev/stdout"
    kept[svtype]++
}
END {
    total = nhet + nhom
    printf("      drawn: %d (het %d, hom %d); %d sites not carried\n", total, nhet, nhom, nref) > "/dev/stderr"
    for (k in kept)    printf("        %-6s %d\n", k, kept[k]) > "/dev/stderr"
    for (k in skipped) printf("        skip %-18s %d\n", k, skipped[k]) > "/dev/stderr"
    if (total == 0) {
        print "ERROR: drew zero variants -- every site was filtered or not carried." > "/dev/stderr"
        exit 1
    }
}' "$TMP/sites.tsv" > "$TMP/body.tsv"

echo "[3/3] writing $OUT ..."
{
    echo '##fileformat=VCFv4.2'
    echo "##source=tools/draw_gnomad_sv_vcf.sh (seed=$SEED)"
    echo "##reference=drawn from $(basename "$VCF")"
    echo '##ALT=<ID=DEL,Description="Deletion">'
    echo '##ALT=<ID=DUP,Description="Duplication">'
    echo '##ALT=<ID=INS,Description="Insertion">'
    echo '##ALT=<ID=INV,Description="Inversion">'
    echo '##ALT=<ID=CNV,Description="Copy number variant">'
    echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
    echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">'
    echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT">'
    echo '##INFO=<ID=POP_AF,Number=1,Type=Float,Description="Population allele frequency the draw used. NOT INFO/AF: eidolon reads that as the per-read alt fraction, which for a germline call is set by GT.">'
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    LC_ALL=C sort -k1,1 -k2,2n "$TMP/body.tsv"
} > "$OUT"

N=$(grep -vc '^#' "$OUT")
echo
echo "Wrote $N variants to $OUT"
echo
echo "Use it with:"
echo "  INPUT_VCF=$OUT sbatch scripts/delta/realism_panel.sbatch    # once the panel exposes it"
echo "  # or directly:  input_vcf: $OUT   in a gen-reads config"
echo
echo "SANITY CHECK before trusting a run: a human genome carries SVs on the order of 1e4"
echo "genome-wide. Scale by the fraction of the genome you drew over. A count far from that"
echo "means MIN_AF or TYPES is doing something unintended, not that the genome is unusual."
