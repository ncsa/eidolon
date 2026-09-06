#!/usr/bin/env bash
# Tests for tools/draw_gnomad_sv_vcf.sh.
#
# This script generates the fixture a whole class of runs is measured against, so a silent
# bias in it becomes a silent bias in every conclusion drawn downstream. The genotype draw
# is checked against Hardy-Weinberg rather than for plausibility: an earlier version used a
# glibc-style LCG whose multiplier overflowed awk's 2^53 exact-integer range, and it read
# het 463 against an expectation of 500 while looking entirely reasonable.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRAW="${DRAW:-$HERE/../../../tools/draw_gnomad_sv_vcf.sh}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
eq()  { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }
has() { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }

# Assert the tool is where we think before anything else. Several comparisons below read
# generated files, and a missing tool makes them compare empty against empty -- which passes.
[[ -f "$DRAW" ]] || { echo "FATAL: tool not found at $DRAW" >&2; exit 1; }

if ! command -v bcftools >/dev/null 2>&1 || ! command -v bgzip >/dev/null 2>&1; then
    echo "SKIP: bcftools/bgzip not on PATH (conda activate bioinf, or module load bcftools)"
    exit 0
fi

# ── fixture ──────────────────────────────────────────────────────────────────
# 1000 DEL at AF=0.5 gives a sharp HWE expectation (250/500/250). The other rows exist to
# be REJECTED, each for a different reason, so a filter that stops working is visible.
make_sites() {
    local out="$1"
    {
        echo '##fileformat=VCFv4.2'
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">'
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="e">'
        echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="l">'
        echo '##INFO=<ID=AF,Number=A,Type=Float,Description="af">'
        echo '##contig=<ID=chr22,length=50818468>'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        for i in $(seq 0 999); do
            p=$((1000 + i * 1000))
            printf 'chr22\t%s\tDEL_%s\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=%s;SVLEN=-500;AF=0.5\n' \
                   "$p" "$i" "$((p + 500))"
        done
        printf 'chr22\t2000000\tBND_x\tN\t<BND>\t.\tPASS\tSVTYPE=BND;AF=0.9\n'
        printf 'chr22\t2001000\tFAIL_x\tN\t<DEL>\t.\tLOWQUAL\tSVTYPE=DEL;END=2001500;SVLEN=-500;AF=0.9\n'
        printf 'chr22\t2002000\tRARE_x\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2002500;SVLEN=-500;AF=0.001\n'
        printf 'chr22\t2003000\tHUGE_x\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=9003000;SVLEN=-7000000;AF=0.9\n'
        printf 'chr22\t2004000\tNOAF_x\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2004500;SVLEN=-500\n'
    } > "$out"
    bgzip -f "$out"
    tabix -f -p vcf "${out}.gz" 2>/dev/null || true
}
make_sites "$WORK/sites.vcf"
SITES="$WORK/sites.vcf.gz"

run() { SEED="$1" bash "$DRAW" "$SITES" "$2" chr22 2>"${2}.err"; }
counts() { grep -oE 'het [0-9]+, hom [0-9]+\); [0-9]+' "$1" | grep -oE '[0-9]+' | tr '\n' ' '; }

echo "=== the draw follows Hardy-Weinberg ==="
# Pooled over twelve seeds so the test is about the generator, not one lucky stream.
tot_het=0; tot_hom=0; tot_ref=0
for s in 1 2 3 5 8 13 21 42 99 500 1234 7777; do
    run "$s" "$WORK/d_$s.vcf" >/dev/null
    read -r h m r <<<"$(counts "$WORK/d_$s.vcf.err")"
    tot_het=$((tot_het + h)); tot_hom=$((tot_hom + m)); tot_ref=$((tot_ref + r))
done
n=$((tot_het + tot_hom + tot_ref))
# 1001 sites reach the genotype draw per seed, not 1000: the AF=0.001 row exists to exercise
# MIN_AF, and under the default MIN_AF=0 it is correctly drawn against (essentially always
# hom-ref). It contributes ~12 hom-refs to the pool below, which moves the chi-square by far
# less than its 5.99 threshold.
eq "every site is accounted for across 12 seeds" "$n" "12012"
# chi-square on 2 df; 5% critical value is 5.99. A generator losing bits reads ~13 here.
chi=$(awk -v h="$tot_het" -v m="$tot_hom" -v r="$tot_ref" -v n="$n" 'BEGIN{
    e1=n*0.25; e2=n*0.5; e3=n*0.25
    printf "%.2f", (r-e1)^2/e1 + (h-e2)^2/e2 + (m-e3)^2/e3 }')
if awk -v c="$chi" 'BEGIN{exit !(c < 5.99)}'; then
    ok "genotype proportions fit HWE (chi-square $chi < 5.99 on 2 df)"
else
    bad "genotype proportions fit HWE" "chi-square < 5.99" "chi-square $chi (het $tot_het hom $tot_hom ref $tot_ref of $n)"
fi
# Non-vacuity: the test must be able to fail. A 50/50 het/hom split would be chi-square ~4000.
chi_bad=$(awk -v n="$n" 'BEGIN{ e1=n*0.25; e2=n*0.5; e3=n*0.25
    printf "%.0f", (0-e1)^2/e1 + (n/2-e2)^2/e2 + (n/2-e3)^2/e3 }')
if awk -v c="$chi_bad" 'BEGIN{exit !(c >= 5.99)}'; then
    ok "the HWE check would reject a biased draw (chi-square $chi_bad)"
else
    bad "the HWE check would reject a biased draw" "a large chi-square" "$chi_bad"
fi

echo "=== the draw is reproducible, and seeds are independent ==="
run 42 "$WORK/a.vcf" >/dev/null; run 42 "$WORK/b.vcf" >/dev/null; run 43 "$WORK/c.vcf" >/dev/null
strip() { grep -v '^##source' "$1"; }
if diff -q <(strip "$WORK/a.vcf") <(strip "$WORK/b.vcf") >/dev/null; then
    ok "the same seed gives the same VCF"
else
    bad "the same seed gives the same VCF" "identical" "differ"
fi
if diff -q <(strip "$WORK/a.vcf") <(strip "$WORK/c.vcf") >/dev/null; then
    bad "a different seed gives a different VCF" "different output" "identical"
else
    ok "a different seed gives a different VCF"
fi

echo "=== population frequency never lands in INFO/AF ==="
# eidolon reads INFO/AF as the per-read alt fraction. Writing gnomAD's AF there would render
# a common heterozygote at its population frequency instead of 0.5, which is invisible in the
# output and wrong everywhere downstream.
eq "no bare INFO/AF key is emitted" \
   "$(grep -oE '(^|;|\t)AF=' "$WORK/a.vcf" | wc -l | tr -d ' ')" "0"
has "frequency is recorded as POP_AF instead"  "$(cat "$WORK/a.vcf")" "POP_AF="
has "the header explains why"                  "$(cat "$WORK/a.vcf")" "NOT INFO/AF"

echo "=== each filter rejects its own row, and nothing else ==="
err="$WORK/a.vcf.err"
has "BND is excluded by default"        "$(cat "$err")" "skip type BND"
has "non-PASS sites are excluded"       "$(cat "$err")" "skip not PASS"
has "sites without AF are excluded"     "$(cat "$err")" "skip no AF"
# Only the 1000 AF=0.5 DELs can survive; every other row is rejected for a distinct reason.
drawn_n="$(grep -oE 'drawn: [0-9]+' "$err" | grep -oE '[0-9]+')"
[[ -n "$drawn_n" && "$drawn_n" -gt 0 ]] \
  && ok "the run drew a non-zero number of variants" \
  || bad "the run drew a non-zero number of variants" "> 0" "${drawn_n:-nothing}"
eq "the file holds exactly what was reported drawn" \
   "$(grep -vc '^#' "$WORK/a.vcf")" "$drawn_n"
kept_ids=$(grep -v '^#' "$WORK/a.vcf" | cut -f3 | grep -cvE '^DEL_[0-9]+$' || true)
eq "no rejected row reached the output" "$kept_ids" "0"

echo "=== MIN_AF and MAX_LEN are honored ==="
SEED=42 MIN_AF=0.9 bash "$DRAW" "$SITES" "$WORK/minaf.vcf" chr22 >/dev/null 2>"$WORK/minaf.err" || true
has "MIN_AF drops sites below the threshold" "$(cat "$WORK/minaf.err")" "below min_af"
SEED=42 MAX_LEN=100 bash "$DRAW" "$SITES" "$WORK/maxlen.vcf" chr22 >/dev/null 2>"$WORK/maxlen.err" || true
has "MAX_LEN drops oversized SVs"            "$(cat "$WORK/maxlen.err")" "over max_len"

echo "=== a draw that yields nothing is a hard failure ==="
# Rule 4: an empty fixture is not a small fixture. Silently writing a header-only VCF would
# produce a run that measures nothing and reports it as a result.
if SEED=42 MIN_AF=0.999999 bash "$DRAW" "$SITES" "$WORK/empty.vcf" chr22 >/dev/null 2>"$WORK/empty.err"; then
    bad "an empty draw exits non-zero" "failure" "exit 0"
else
    ok "an empty draw exits non-zero"
fi
has "and says why" "$(cat "$WORK/empty.err")" "drew zero variants"

echo "=== the output is a VCF the toolchain accepts ==="
bcftools view -h "$WORK/a.vcf" >/dev/null 2>&1 \
  && ok "bcftools parses the header" || bad "bcftools parses the header" "valid VCF" "rejected"
eq "every genotype is het or hom-alt" \
   "$(bcftools query -f '[%GT]\n' "$WORK/a.vcf" 2>/dev/null | grep -cvE '^(0/1|1/1)$')" "0"
eq "SVTYPE survives the round trip" \
   "$(bcftools query -f '%INFO/SVTYPE\n' "$WORK/a.vcf" 2>/dev/null | sort -u | tr '\n' ' ')" "DEL "
eq "records are sorted by position" \
   "$(bcftools query -f '%POS\n' "$WORK/a.vcf" 2>/dev/null | sort -c -n 2>&1 | wc -l | tr -d ' ')" "0"

MIN_ASSERTIONS=21
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi

printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
