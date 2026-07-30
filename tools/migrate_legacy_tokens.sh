#!/usr/bin/env bash
# migrate_legacy_tokens.sh — normalize a pre-v3.0.0 eidolon VCF to the EIDOLON_* tokens.
#
# v3.0.0 renamed eidolon's emitted output tokens (see CHANGELOG):
#
#   INFO tags      NEAT_ORIGIN / NEAT_PROVENANCE / NEAT_REASON / NEAT_CCF / NEAT_VAF
#                    -> EIDOLON_*
#   sample column  NEAT_simulated_sample -> EIDOLON_simulated_sample
#
# Anything that filters on the new names (the benchmark harnesses, downstream
# scoring) cannot read a v2.0.0 truth VCF: bcftools rejects an undefined tag in a
# filter expression at PARSE time, so even
#
#   -i 'INFO/EIDOLON_ORIGIN="somatic" || INFO/NEAT_ORIGIN="somatic"'
#
# fails outright rather than falling back. Converting the file up front is the
# reliable fix, and it's what this script does.
#
# The conversion is IDEMPOTENT and quiet: running it on an already-migrated file (or a
# non-eidolon VCF) leaves every record and declaration untouched and exits 0, so callers
# can invoke it unconditionally without probing first.
#
# "Untouched" is about the DATA, not the bytes — every pass goes through bcftools, which
# adds a `##bcftools_viewCommand` provenance line and normalizes float formatting
# (`1.0` -> `1`). Two other cosmetic effects of `--rename-annots`: the tag's
# `Description=` text is not rewritten, so a converted header can still read
# "...is NEAT_CCF x allele dosage"; and output is always BGZF even if the input was
# plain gzip. None of these change record semantics.
#
# NOT handled: FASTQ/BAM read-name prefixes (@RNEAT_generated_ / RNEAT_chimeric_).
# Rewriting those means a full pass over files that are routinely hundreds of GB,
# to change ~19 bytes per record. `eidolon filter-reads` accepts both prefixes
# natively instead, so legacy FASTQs need no conversion.
#
# Usage:
#   tools/migrate_legacy_tokens.sh IN.vcf[.gz] OUT.vcf.gz
#   tools/migrate_legacy_tokens.sh --check IN.vcf[.gz]   # report only, convert nothing
#
# Exit: 0 on success (converted or already current). With --check, 0 = already
# current, 10 = legacy tokens present.
set -euo pipefail

die() { echo "migrate_legacy_tokens: $*" >&2; exit 1; }

CHECK_ONLY=0
if [[ "${1:-}" == "--check" ]]; then CHECK_ONLY=1; shift; fi
IN="${1:-}"
[[ -n "$IN" ]] || die "usage: $0 IN.vcf[.gz] OUT.vcf.gz  (or --check IN.vcf[.gz])"
[[ -f "$IN" ]] || die "input not found: $IN"
command -v bcftools >/dev/null || die "bcftools not on PATH"

if [[ "$CHECK_ONLY" -eq 0 ]]; then
    OUT="${2:-}"
    [[ -n "$OUT" ]] || die "missing OUT argument (use --check for a report-only run)"
    [[ "$OUT" == *.gz ]] || die "OUT must end in .gz (output is bgzipped): $OUT"
    [[ "$(readlink -f "$IN")" != "$(readlink -f "$OUT")" ]] || die "IN and OUT must differ"
fi

# ── Which legacy tokens does this file actually declare? ──────────────────────
# Read the header once. Only tags DECLARED in the header can be renamed by
# bcftools annotate, and only declared tags are what break filter expressions,
# so the header is the right thing to key off.
header="$(bcftools view -h "$IN")"

LEGACY_INFO=(NEAT_ORIGIN NEAT_PROVENANCE NEAT_REASON NEAT_CCF NEAT_VAF)
rename_map=""
found=()
for tag in "${LEGACY_INFO[@]}"; do
    # Exact ID match, so NEAT_ORIGIN doesn't also match a NEAT_ORIGINAL.
    if grep -q "^##INFO=<ID=${tag}," <<<"$header"; then
        found+=("$tag")
        rename_map+="INFO/${tag} EIDOLON_${tag#NEAT_}"$'\n'
    fi
done

# Sample columns: rename only names we know are ours, never an arbitrary sample.
#
# Take the names from the #CHROM line of the header we already captured, NOT from a
# second `bcftools query -l`. A process substitution's failure is invisible to set -e:
# if the command died, the loop body would simply never run, `legacy_samples` would
# stay 0, and we'd rename the INFO tags while silently leaving the sample column as
# NEAT_simulated_sample — a half-converted file returned with exit 0.
chrom_line="$(grep -m1 '^#CHROM' <<<"$header")" \
    || die "no #CHROM header line in $IN — not a VCF?"
legacy_samples=0
new_samples=""
# Fields 1-9 are the fixed VCF columns; 10+ are samples (absent in a sites-only VCF).
n_fields="$(awk -F'\t' '{print NF; exit}' <<<"$chrom_line")"
if [[ "$n_fields" -gt 9 ]]; then
    while IFS= read -r s; do
        [[ -n "$s" ]] || continue
        case "$s" in
            NEAT_simulated_sample) new_samples+="EIDOLON_simulated_sample"$'\n'; legacy_samples=1 ;;
            *)                     new_samples+="$s"$'\n' ;;
        esac
    done < <(cut -f10- <<<"$chrom_line" | tr '\t' '\n')
fi

if [[ "$CHECK_ONLY" -eq 1 ]]; then
    if [[ ${#found[@]} -eq 0 && "$legacy_samples" -eq 0 ]]; then
        echo "already current: no pre-v3.0.0 tokens in $IN"
        exit 0
    fi
    echo "legacy tokens in $IN:"
    [[ ${#found[@]} -gt 0 ]] && echo "  INFO: ${found[*]}"
    [[ "$legacy_samples" -eq 1 ]] && echo "  sample column: NEAT_simulated_sample"
    exit 10
fi

# ── Convert ──────────────────────────────────────────────────────────────────
if [[ ${#found[@]} -eq 0 && "$legacy_samples" -eq 0 ]]; then
    # Already current (or not an eidolon VCF). Normalize to bgzip+index anyway so
    # the caller gets a usable OUT either way — that's what makes an unconditional
    # call site safe.
    bcftools view -O z -o "$OUT" "$IN"
    bcftools index -f -t "$OUT"
    exit 0
fi

echo "migrate_legacy_tokens: converting pre-v3.0.0 tokens in $(basename "$IN")" >&2
[[ ${#found[@]} -gt 0 ]] && echo "  INFO: ${found[*]} -> EIDOLON_*" >&2
[[ "$legacy_samples" -eq 1 ]] && echo "  sample: NEAT_simulated_sample -> EIDOLON_simulated_sample" >&2

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
stage="$IN"

if [[ ${#found[@]} -gt 0 ]]; then
    # Guard the case of a half-converted file: renaming NEAT_X to EIDOLON_X when
    # EIDOLON_X is already declared would put a duplicate ID in the header.
    for tag in "${found[@]}"; do
        if grep -q "^##INFO=<ID=EIDOLON_${tag#NEAT_}," <<<"$header"; then
            die "both ${tag} and EIDOLON_${tag#NEAT_} are declared in $IN — refusing to
  create a duplicate header ID. Drop one with 'bcftools annotate -x INFO/${tag}' first."
        fi
    done
    printf '%s' "$rename_map" > "$tmpdir/rename.txt"
    bcftools annotate --rename-annots "$tmpdir/rename.txt" -O z -o "$tmpdir/renamed.vcf.gz" "$stage"
    stage="$tmpdir/renamed.vcf.gz"
fi

if [[ "$legacy_samples" -eq 1 ]]; then
    printf '%s' "$new_samples" > "$tmpdir/samples.txt"
    bcftools reheader -s "$tmpdir/samples.txt" -o "$tmpdir/reheadered.vcf.gz" "$stage"
    stage="$tmpdir/reheadered.vcf.gz"
fi

bcftools view -O z -o "$OUT" "$stage"
bcftools index -f -t "$OUT"
