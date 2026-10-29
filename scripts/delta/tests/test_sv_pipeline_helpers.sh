#!/usr/bin/env bash
# Tests for three sv_pipeline.sbatch functions with no coverage until #825, plus the index
# guard index_and_align relies on:
#   build_bnd_spans           three-way status: 0 = spans written, 2 = every junction is
#                             inter-contig (legitimate), 1 = a defect. An earlier version
#                             reported case 1 as case 2.
#   convert_manta_inversions  checks Manta's convertInversion.py did what it claims instead
#                             of trusting its exit code; 2 = converter absent (older Manta).
#   index_and_align           skips an aligned sample, derives R2 from R1, and indexes
#                             through index_reference_locked (lib_report.sh), which rebuilds
#                             a 0-byte index: a build killed mid-write leaves one that `-f`
#                             treats as present, and bwa-mem2 then fails every later run.
#
# The functions are EXTRACTED VERBATIM from the pipeline and lib_report.sh. bcftools is real
# (spans and conversions are really built and sorted); bwa-mem2, samtools sort/index, conda
# and Manta's converter are stubbed, which is the only way to provoke a converter that
# no-ops or contradicts itself.
#
# Usage:
#   scripts/delta/tests/test_sv_pipeline_helpers.sh            # run the suite
#   scripts/delta/tests/test_sv_pipeline_helpers.sh --mutate   # prove it is non-vacuous
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PIPELINE:-$HERE/../sv_pipeline.sbatch}"   # overridden by --mutate
LIBREPORT="${LIBREPORT:-$HERE/../lib_report.sh}"      # overridden by --mutate
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Mutations: <label>@<P|L>@<from>@<to>; P = the pipeline, L = lib_report.sh.
if [[ "${1:-}" == "--mutate" ]]; then
    survived=0
    while IFS='@' read -r label which from to; do
        [[ -n "$label" ]] || continue
        if [[ "$which" == "L" ]]; then src="$LIBREPORT"; else src="$PIPELINE"; fi
        cp "$src" "$WORK/mutant"
        FROM="$(printf '%b' "$from")" TO="$(printf '%b' "$to")" \
            perl -0pi -e 's/\Q$ENV{FROM}\E/$ENV{TO}/' "$WORK/mutant"
        if cmp -s "$src" "$WORK/mutant"; then
            printf '  ERROR   %-46s mutation did not apply — anchor not found\n' "$label"
            survived=$((survived+1)); continue
        fi
        if [[ "$which" == "L" ]]; then
            LIBREPORT="$WORK/mutant" bash "$0" >/dev/null 2>&1; rc=$?
        else
            PIPELINE="$WORK/mutant" bash "$0" >/dev/null 2>&1; rc=$?
        fi
        if [[ "$rc" -eq 0 ]]; then
            printf '  SURVIVED %-45s <- nothing caught this\n' "$label"
            survived=$((survived+1))
        else
            printf '  caught   %s\n' "$label"
        fi
    done <<'MUTATIONS'
inter-contig only reads as a defect@P@            return 2\n        fi\n        echo "  ERROR: BND spans: 0 junctions@            return 1\n        fi\n        echo "  ERROR: BND spans: 0 junctions
an unparsable truth reads as legitimate@P@if [[ "$bad" -eq 0 && "$inter" -gt 0 ]]; then@if [[ "$inter" -gt 0 ]]; then
a junction pair yields two spans@P@if (key in seen) { dup++; next }@if (0) { dup++; next }
inter-contig junctions become spans@P@if (mc != $1) { inter++; next }@if (0) { inter++; next }
missing header lines are not backfilled@P@grep -q "^${kind}=<ID=${id}," "$hdr" || echo "$decl" >> "$hdr"@:
a missing BND truth is not an error@P@[[ -f "$src" ]] || { echo "  ERROR: BND spans: no BND truth at $src" >&2; return 1; }@:
an absent converter reads as converted@P@        conda_activate bioinf || true\n        return 2@        conda_activate bioinf || true\n        return 0
a crashed converter is trusted@P@    if [[ "$rc" -ne 0 ]]; then\n        echo "  ERROR: convertInversion.py failed@    if false; then\n        echo "  ERROR: convertInversion.py failed
a no-op conversion is trusted@P@    if [[ "$inv_out" -eq 0 ]]; then@    if false; then
an inconsistent conversion is trusted@P@    if [[ "$bnd_out" -ge "$bnd_in" ]]; then@    if false; then
an aligned sample is realigned@P@    if [[ -f "$out_bam" && -f "${out_bam}.bai" ]]; then@    if false; then
R2 is not derived from R1@P@local in_r2="${in_fq/_r1./_r2.}"@local in_r2="$in_fq"
paired reads aligned as single-end@P@    if [[ "$PAIRED" == "true" ]]; then\n        local in_r2@    if false; then\n        local in_r2
alignment skips the index step@P@    index_reference_locked "$REFERENCE"\n    echo "  Aligning@    echo "  Aligning
a 0-byte index reads as present@L@        if [[ ! -s "${ref}.bwt.2bit.64" ]] && [[ ! -s "${ref}.bwt" ]]; then\n            echo "  Indexing reference for BWA-MEM2 (one-time, lock held)..."@        if [[ ! -f "${ref}.bwt.2bit.64" ]] && [[ ! -f "${ref}.bwt" ]]; then\n            echo "  Indexing reference for BWA-MEM2 (one-time, lock held)..."
MUTATIONS
    printf '\n──────── %d mutation(s) survived ────────\n' "$survived"
    [[ "$survived" -eq 0 ]]
    exit $?
fi

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); printf '  ok    %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  FAIL  %s\n     expected: %s\n     actual:   %s\n' "$1" "$2" "$3"; }
is()   { [[ "$2" == "$3" ]] && ok "$1" || bad "$1" "$3" "$2"; }   # is <name> <actual> <expected>
has()  { case "$2" in *"$3"*) ok "$1";; *) bad "$1" "contains: $3" "$2";; esac; }
hasnt(){ case "$2" in *"$3"*) bad "$1" "must NOT contain: $3" "$2";; *) ok "$1";; esac; }

for tool in bcftools bgzip; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        if [[ -n "${CI:-}" ]]; then
            echo "FATAL: $tool is not on PATH and CI is set; refusing to report success" >&2
            exit 2
        fi
        echo "SKIP: $tool not on PATH (set CI=1 to make this fatal)"
        exit 0
    fi
done

extract() { awk "/^$1\(\)/,/^}\$/" "$2"; }
for fn in build_bnd_spans convert_manta_inversions index_and_align; do
    src="$(extract "$fn" "$PIPELINE")"
    [[ -n "$src" ]] || { echo "FATAL: could not extract $fn from $PIPELINE"; exit 2; }
    eval "$src"
done
src="$(extract index_reference_locked "$LIBREPORT")"
[[ -n "$src" ]] || { echo "FATAL: could not extract index_reference_locked from $LIBREPORT"; exit 2; }
eval "$src"

OUTDIR="$WORK/out"; mkdir -p "$OUTDIR"

# A VCF from "CHROM POS ID ALT INFO" lines. $1 = output, $2 = 1 to declare END/SVLEN in the
# header (0 leaves them out, to test the backfill).
mkvcf() {
    local out="$1" full="$2"; shift 2
    {
        printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n'
        printf '##contig=<ID=chr2,length=1000000>\n'
        printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        if [[ "$full" == "1" ]]; then
            printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
            printf '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">\n'
        fi
        printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
        local rec
        for rec in "$@"; do
            set -- $rec
            printf '%s\t%s\t%s\tN\t%s\t.\tPASS\t%s\tGT\t0/1\n' "$1" "$2" "$3" "$4" "$5"
        done
    } | bcftools sort -O z -o "$out" 2>/dev/null
}
nrec() { bcftools view -H "$1" 2>/dev/null | wc -l | tr -d ' '; }

echo "=== build_bnd_spans: an intra-contig junction pair becomes one span ==="
mkvcf "$WORK/intra.vcf.gz" 0 \
    "chr1 100 a N[chr1:5000[ SVTYPE=BND" \
    "chr1 5000 b ]chr1:100]N SVTYPE=BND"
out="$(build_bnd_spans "$WORK/intra.vcf.gz" "$WORK/spans.vcf.gz" 2>&1)"; rc=$?
is   "it returns 0"                        "$rc" "0"
is   "the mate pair yields ONE span"       "$(nrec "$WORK/spans.vcf.gz")" "1"
is   "spanning the two breakends"          "$(bcftools query -f '%POS-%INFO/END,%INFO/SVLEN' "$WORK/spans.vcf.gz" 2>/dev/null | tr -d '\n')" "100-5000,4900"
has  "with END backfilled into a header that lacked it" "$(bcftools view -h "$WORK/spans.vcf.gz")" "##INFO=<ID=END,"

echo "=== build_bnd_spans: only inter-contig junctions is legitimate (rc 2) ==="
mkvcf "$WORK/inter.vcf.gz" 1 \
    "chr1 100 a N[chr2:500[ SVTYPE=BND" \
    "chr2 500 b ]chr1:100]N SVTYPE=BND"
out="$(build_bnd_spans "$WORK/inter.vcf.gz" "$WORK/s2.vcf.gz" 2>&1)"; rc=$?
is   "it returns 2, not 1"                 "$rc" "2"
has  "and says why"                        "$out" "all 2 junction record(s) are inter-contig"

echo "=== build_bnd_spans: an unparsable truth is a defect (rc 1), even beside inter-contig records ==="
mkvcf "$WORK/broken.vcf.gz" 1 \
    "chr1 100 a <BND> SVTYPE=BND" \
    "chr1 300 c N[chr2:500[ SVTYPE=BND"
out="$(build_bnd_spans "$WORK/broken.vcf.gz" "$WORK/s3.vcf.gz" 2>&1)"; rc=$?
is   "it returns 1, not 2"                 "$rc" "1"
has  "and names the unparsable ALT"        "$out" "unparsable ALT=1"

echo "=== build_bnd_spans: mixed intra and inter writes only the intra spans ==="
mkvcf "$WORK/mixed.vcf.gz" 1 \
    "chr1 100 a N[chr1:5000[ SVTYPE=BND" \
    "chr1 5000 b ]chr1:100]N SVTYPE=BND" \
    "chr1 9000 c N[chr2:500[ SVTYPE=BND"
build_bnd_spans "$WORK/mixed.vcf.gz" "$WORK/s4.vcf.gz" >/dev/null 2>&1; rc=$?
is   "it returns 0"                        "$rc" "0"
is   "one span, the inter-contig record skipped" "$(nrec "$WORK/s4.vcf.gz")" "1"

echo "=== build_bnd_spans: a missing truth file is an error ==="
out="$(build_bnd_spans "$WORK/nope.vcf.gz" "$WORK/s5.vcf.gz" 2>&1)"; rc=$?
is   "it returns 1"                        "$rc" "1"
# The status alone is not enough: without the check, an absent file falls through to "0
# junctions derived ... the ALT parser or the truth is wrong", blaming the wrong thing.
has  "and names the missing file"          "$out" "no BND truth at $WORK/nope.vcf.gz"

# ── convert_manta_inversions: stub conda, Manta's layout and its converter ─────
conda_activate() { return 0; }
MANTA_ENV=manta; REFERENCE="$WORK/ref.fa"; : > "$REFERENCE"
mkdir -p "$WORK/manta/bin" "$WORK/manta/libexec"
printf '#!/bin/sh\n' > "$WORK/manta/bin/configManta.py"; chmod +x "$WORK/manta/bin/configManta.py"
PATH="$WORK/manta/bin:$PATH"
# The converter, as a `python` function: $1 = convertInversion.py, $4 = the plain VCF.
python() {
    local vcf="$4"
    case "$CONV_MODE" in
        ok)           grep -v 'SVTYPE=BND' "$vcf"
                      printf 'chr1\t100\tinv1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=5000\tGT\t0/1\n' ;;
        noop)         cat "$vcf" ;;
        inconsistent) cat "$vcf"
                      printf 'chr1\t100\tinv1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=5000\tGT\t0/1\n' ;;
        crash)        echo "converter blew up" >&2; return 3 ;;
    esac
}
mkvcf "$WORK/manta_raw.vcf.gz" 1 \
    "chr1 100 m1 N]chr1:5000] SVTYPE=BND" \
    "chr1 5000 m2 N]chr1:100] SVTYPE=BND" \
    "chr1 20000 d1 <DEL> SVTYPE=DEL;END=20500"

echo "=== convert_manta_inversions: a real conversion is accepted ==="
: > "$WORK/manta/libexec/convertInversion.py"
CONV_MODE=ok out="$(convert_manta_inversions "$WORK/manta_raw.vcf.gz" "$WORK/conv.vcf.gz" 2>&1)"; rc=$?
is   "it returns 0"                        "$rc" "0"
is   "the output holds the INV"            "$(bcftools view -H -i 'INFO/SVTYPE="INV"' "$WORK/conv.vcf.gz" 2>/dev/null | wc -l | tr -d ' ')" "1"
has  "and reports the counts"              "$out" "2 BND -> 0 BND + 1 INV"

echo "=== convert_manta_inversions: a no-op conversion is refused ==="
CONV_MODE=noop out="$(convert_manta_inversions "$WORK/manta_raw.vcf.gz" "$WORK/conv2.vcf.gz" 2>&1)"; rc=$?
is   "it returns 1"                        "$rc" "1"
has  "and says no INV appeared"            "$out" "produced NO <INV> records"

echo "=== convert_manta_inversions: an inconsistent conversion is refused ==="
CONV_MODE=inconsistent out="$(convert_manta_inversions "$WORK/manta_raw.vcf.gz" "$WORK/conv3.vcf.gz" 2>&1)"; rc=$?
is   "it returns 1"                        "$rc" "1"
has  "and says the BND count did not drop" "$out" "BND count did not drop"

echo "=== convert_manta_inversions: a crashed converter is refused ==="
CONV_MODE=crash out="$(convert_manta_inversions "$WORK/manta_raw.vcf.gz" "$WORK/conv4.vcf.gz" 2>&1)"; rc=$?
is   "it returns 1"                        "$rc" "1"
has  "with the converter's rc and message" "$out" "convertInversion.py failed (rc=3)"
has  "and its stderr"                      "$out" "converter blew up"

echo "=== convert_manta_inversions: an older Manta without the converter is rc 2 ==="
rm -f "$WORK/manta/libexec/convertInversion.py"
CONV_MODE=ok out="$(convert_manta_inversions "$WORK/manta_raw.vcf.gz" "$WORK/conv5.vcf.gz" 2>&1)"; rc=$?
is   "it returns 2"                        "$rc" "2"
[[ ! -e "$WORK/conv5.vcf.gz" ]] && ok "and writes no output" || bad "and writes no output" "no file" "conv5.vcf.gz exists"

# ── index_and_align and its index guard: stub the aligner and samtools ─────────
LOG="$WORK/tools.log"
bwa-mem2() {
    echo "bwa-mem2 $*" >> "$LOG"
    if [[ "$1" == "index" ]]; then printf 'idx' > "$2.bwt.2bit.64"; fi
    if [[ "$1" == "mem" ]]; then echo "@SQ aligned"; fi
}
samtools() {
    echo "samtools $*" >> "$LOG"
    case "$1" in
        sort)  local out="" prev=""; for a in "$@"; do [[ "$prev" == "-o" ]] && out="$a"; prev="$a"; done
               cat > "$out" ;;
        index) : > "$2.bai" ;;
        faidx) : > "$2.fai" ;;
    esac
}
THREADS=2

echo "=== index_reference_locked: a 0-byte index is rebuilt, a real one is kept ==="
REF2="$WORK/ref2.fa"; : > "$REF2"; : > "$REF2.bwt.2bit.64"   # 0 bytes: a killed build
: > "$LOG"; index_reference_locked "$REF2" >/dev/null 2>&1
has  "a 0-byte index triggers a rebuild"   "$(cat "$LOG")" "bwa-mem2 index $REF2"
: > "$LOG"; index_reference_locked "$REF2" >/dev/null 2>&1
hasnt "a real index is not rebuilt"        "$(cat "$LOG")" "bwa-mem2 index"

echo "=== index_and_align: paired reads, R2 derived from R1 ==="
printf 'idx' > "$REFERENCE.bwt.2bit.64"; : > "$REFERENCE.fai"
PAIRED=true; : > "$LOG"
index_and_align "$WORK/tumor_r1.fastq.gz" "$WORK/tumor.bam" tumor >/dev/null 2>&1; rc=$?
is   "it returns 0"                        "$rc" "0"
has  "it aligns R1 with the derived R2"    "$(grep '^bwa-mem2 mem' "$LOG")" "$WORK/tumor_r1.fastq.gz $WORK/tumor_r2.fastq.gz"
has  "with a read group"                   "$(grep '^bwa-mem2 mem' "$LOG")" "ID:tumor"
[[ -f "$WORK/tumor.bam" && -f "$WORK/tumor.bam.bai" ]] && ok "and writes the BAM and its index" \
    || bad "and writes the BAM and its index" "tumor.bam + .bai" "missing"

echo "=== index_and_align: single-end reads align R1 alone ==="
PAIRED=false; : > "$LOG"
index_and_align "$WORK/normal_r1.fastq.gz" "$WORK/normal.bam" normal >/dev/null 2>&1
hasnt "no R2 on a single-end run"          "$(grep '^bwa-mem2 mem' "$LOG")" "normal_r2"

echo "=== index_and_align: an aligned sample is skipped; a fresh one indexes first ==="
PAIRED=true; : > "$LOG"
out="$(index_and_align "$WORK/tumor_r1.fastq.gz" "$WORK/tumor.bam" tumor 2>&1)"
has  "it says the BAM is already present"  "$out" "tumor BAM already present"
is   "and runs no aligner"                 "$(grep -c '^bwa-mem2 mem' "$LOG")" "0"
rm -f "$REFERENCE.bwt.2bit.64"; : > "$LOG"
index_and_align "$WORK/x_r1.fastq.gz" "$WORK/x.bam" x >/dev/null 2>&1
is   "a missing index is built before aligning" "$(grep -E '^bwa-mem2 (index|mem)' "$LOG" | cut -d' ' -f2 | tr '\n' ' ')" "index mem "

MIN_ASSERTIONS=34
TOTAL=$((PASS + FAIL))
if [[ "$TOTAL" -lt "$MIN_ASSERTIONS" ]]; then
    printf '\n  FAIL  only %d assertions ran, expected at least %d\n' "$TOTAL" "$MIN_ASSERTIONS"
    FAIL=$((FAIL+1))
fi
printf '\n──────── %d passed, %d failed ────────\n' "$PASS" "$FAIL"
[[ "$FAIL" -eq 0 ]]
