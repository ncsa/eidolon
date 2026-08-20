#!/usr/bin/env bash
# Validate the v1.10 CNV/SV model pipeline against real gnomAD-SV data.
#
# Two modes:
#   smoke   — chr22 slice only (~2 min, good for first runs)
#   full    — full gnomAD-SV + GRCh38 (~10-30 min, real validation)
#
# What it does:
#   1. Downloads gnomAD-SV v4 sites VCF (or uses a local copy if provided)
#   2. Downloads / requires GRCh38 reference FASTA
#   3. Pre-filters: split multi-allelic, drop INV/CPX/CTX/BND/<INS:ME:*>
#      (informational — the trainer drops these anyway, but pre-filtering
#      keeps the warn-log volume manageable)
#   4. Times `eidolon gen-mut-model` and reports wall-clock duration
#   5. Pretty-prints the resulting `sv_model` block from the trained JSON
#   6. Compares the fitted parameters against the bundled defaults in
#      eidolon-core/src/models/sv_model_defaults.rs and flags large deviations
#
# Dependencies (must be on PATH):
#   - eidolon (built via `cargo build --release`)
#   - bcftools
#   - jq
#   - samtools (for FASTA indexing / slicing)
#   - curl (or wget) for downloads
#
# Usage:
#   ./tools/validate_with_gnomad.sh smoke
#   ./tools/validate_with_gnomad.sh full
#   GNOMAD_VCF=/path/to/local.vcf.gz REF_FASTA=/path/to/GRCh38.fa \
#       ./tools/validate_with_gnomad.sh full
#
set -euo pipefail

MODE="${1:-smoke}"
WORK="${WORK_DIR:-./gnomad_validation_work}"
mkdir -p "$WORK"
cd "$WORK"

# ── Inputs ──────────────────────────────────────────────────────────────
# Allow callers to supply local files via env vars; otherwise we download
# the smoke-friendly small slices. Full gnomAD-SV is large (~3.5GB) and
# GRCh38 is ~900MB — script will refuse to download those without an
# explicit ALLOW_LARGE_DOWNLOADS=1.
# Verified live 2026-08-02 (1.74 GB). The previous path — release/4/sv/gnomad_v4_sv.sites.vcf.gz
# — now 404s; gnomAD reorganised to release/<ver>/genome_sv/ with a dotted filename.
# A rotted URL makes this script fail at step 1, so check it before assuming the
# validation "has not been run lately" for any other reason.
GNOMAD_VCF_URL="${GNOMAD_VCF_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz}"
GRCH38_URL="${GRCH38_URL:-https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"

ensure_vcf() {
    if [[ -n "${GNOMAD_VCF:-}" && -f "$GNOMAD_VCF" ]]; then
        echo "Using local gnomAD VCF: $GNOMAD_VCF"
        return
    fi
    if [[ "$MODE" == "full" && "${ALLOW_LARGE_DOWNLOADS:-0}" != "1" ]]; then
        echo "Full gnomAD-SV is ~3.5GB. Set ALLOW_LARGE_DOWNLOADS=1 to download, or" >&2
        echo "supply GNOMAD_VCF=/path/to/file.vcf.gz pointing at a local copy." >&2
        exit 2
    fi
    if [[ ! -f gnomad_sv.vcf.gz ]]; then
        echo "Downloading gnomAD-SV sites VCF..."
        curl -sSL -o gnomad_sv.vcf.gz "$GNOMAD_VCF_URL"
    fi
    GNOMAD_VCF="$(pwd)/gnomad_sv.vcf.gz"
}

ensure_ref() {
    if [[ -n "${REF_FASTA:-}" && -f "$REF_FASTA" ]]; then
        echo "Using local reference: $REF_FASTA"
        return
    fi
    if [[ "${ALLOW_LARGE_DOWNLOADS:-0}" != "1" ]]; then
        echo "GRCh38 is ~900MB. Set ALLOW_LARGE_DOWNLOADS=1 to download, or" >&2
        echo "supply REF_FASTA=/path/to/GRCh38.fa pointing at a local copy." >&2
        exit 2
    fi
    if [[ ! -f GRCh38.fa ]]; then
        echo "Downloading GRCh38..."
        curl -sSL -o GRCh38.fa.gz "$GRCH38_URL"
        gunzip GRCh38.fa.gz
        samtools faidx GRCh38.fa
    fi
    REF_FASTA="$(pwd)/GRCh38.fa"
}

# ── Slice & filter ──────────────────────────────────────────────────────
prep_input() {
    local out="filtered.vcf.gz"
    local region=""
    if [[ "$MODE" == "smoke" ]]; then
        # Detect contig naming via the tabix index (fast — `bcftools index
        # --stats` reads the .tbi rather than scanning the data). The
        # earlier "grep first 200 lines of bcftools view" heuristic
        # silently broke on gnomAD-SV because the header alone is hundreds
        # of lines long and there are no data records in the first 200.
        local first_contig
        first_contig="$(bcftools index --stats "$GNOMAD_VCF" 2>/dev/null | head -1 | cut -f1)"
        if [[ -z "$first_contig" ]]; then
            echo "ERROR: bcftools index --stats returned nothing. Index missing?" >&2
            echo "       Try: tabix -p vcf $GNOMAD_VCF" >&2
            exit 1
        fi
        if [[ "$first_contig" == chr* ]]; then
            region="chr22"
        else
            region="22"
        fi
        echo "Smoke mode: restricting to '$region' (detected naming from '$first_contig')"
    fi
    # Split multi-allelic records (so MCNV doesn't get silently dropped at
    # the reader) and exclude SV types v1.11 can't generate so the trainer
    # warn log stays readable. Pipe through bcftools view -r FIRST so the
    # region restriction is applied before any further filtering.
    echo "Filtering: splitting multi-allelic, dropping CPX/CTX"
    if [[ -n "$region" ]]; then
        bcftools view -r "$region" "$GNOMAD_VCF" \
            | bcftools norm -m - \
            | bcftools view -e 'ALT~"^<CPX>$" || ALT~"^<CTX>$"' \
                -Oz -o "$out"
    else
        bcftools norm -m - "$GNOMAD_VCF" \
            | bcftools view -e 'ALT~"^<CPX>$" || ALT~"^<CTX>$"' \
                -Oz -o "$out"
    fi
    bcftools index -t "$out"
    local n
    n="$(bcftools view -H "$out" | wc -l)"
    echo "Filtered records: $n"
    if [[ "$n" -eq 0 ]]; then
        echo "ERROR: zero records survived filtering. Likely causes:" >&2
        echo "  - Contig naming mismatch (script chose '$region')" >&2
        echo "  - Input VCF is empty or unindexed" >&2
        echo "Debug commands:" >&2
        echo "  bcftools index --stats $GNOMAD_VCF | head -5" >&2
        echo "  bcftools view -H $GNOMAD_VCF | head -1 | cut -f1" >&2
        exit 1
    fi
    INPUT_VCF="$(pwd)/$out"
}

# ── Train ───────────────────────────────────────────────────────────────
train() {
    local config="train.yml"
    local output="trained_model.json.gz"
    cat > "$config" <<EOF
reference: $REF_FASTA
vcf_file: $INPUT_VCF
output_file: $output
overwrite_output: true
EOF
    echo "Running gen-mut-model..."
    local start=$(date +%s)
    /usr/bin/time -v eidolon gen-mut-model -c "$config" 2> train.timing.log
    local end=$(date +%s)
    local elapsed=$((end - start))
    echo "Wall clock: ${elapsed}s"
    echo "Peak RSS: $(grep 'Maximum resident set size' train.timing.log | awk '{print $NF}') KB"
    TRAINED_MODEL="$(pwd)/$output"
}

# ── Inspect ─────────────────────────────────────────────────────────────
inspect() {
    echo "=========================================="
    echo "Trained sv_model contents:"
    echo "=========================================="
    zcat "$TRAINED_MODEL" | jq '.sv_model'
    echo
    echo "=========================================="
    echo "Sanity checks:"
    echo "=========================================="
    local sv_model
    sv_model="$(zcat "$TRAINED_MODEL" | jq -c '.sv_model')"
    if [[ "$sv_model" == "null" ]]; then
        echo "FAIL: sv_model is null — no SV types cleared the 2-observation fit bar."
        echo "      This shouldn't happen on real gnomAD data; investigate the filter chain."
        return 1
    fi
    # Type probabilities should sum to ~1
    local type_sum
    type_sum="$(echo "$sv_model" | jq '.type_probabilities | to_entries | map(.value) | add')"
    echo "type_probabilities sum: $type_sum (expect ~1.0)"
    # Compare against bundled defaults — fail soft on large deviations.
    local del_p dup_p cnv_p bnd_p inv_p ins_p
    del_p="$(echo "$sv_model" | jq '.type_probabilities.Del // 0')"
    dup_p="$(echo "$sv_model" | jq '.type_probabilities.Dup // 0')"
    cnv_p="$(echo "$sv_model" | jq '.type_probabilities.Cnv // 0')"
    bnd_p="$(echo "$sv_model" | jq '.type_probabilities.Bnd // 0')"
    inv_p="$(echo "$sv_model" | jq '.type_probabilities.Inv // 0')"
    ins_p="$(echo "$sv_model" | jq '.type_probabilities.Ins // 0')"
    echo "Type breakdown:  DEL=$del_p     DUP=$dup_p     CNV=$cnv_p     BND=$bnd_p     INV=$inv_p     INS=$ins_p"
    echo "Default expects: DEL=0.6433  DUP=0.1470  CNV=0.0004  BND=0.1943  INV=0.0050  INS=0.0100 (v1.11.1)"
    local per_base
    per_base="$(echo "$sv_model" | jq '.per_base_rate')"
    echo "per_base_rate: $per_base (default 6.009e-4 from full gnomAD-SV v4.1)"
    local hom
    hom="$(echo "$sv_model" | jq '.homozygous_frequency')"
    echo "homozygous_frequency: $hom (default 0.20 — literature; gnomAD-SV is sites-only)"
}

# ── Optional: smoke gen-reads against the trained model ─────────────────
smoke_gen_reads() {
    if [[ "$MODE" != "smoke" ]]; then return; fi
    echo "=========================================="
    echo "Running gen-reads against the trained model (sv_rate_scale=1.0)"
    echo "=========================================="
    # Slice the reference to chr22 so gen-reads runs in seconds, not
    # hours. Using the full hg38 here was a footgun — it works, but a
    # smoke test that takes ~1h with a 10GB log isn't a smoke test.
    # Cache the chr22-only FASTA in the work dir; awk is portable and
    # doesn't require samtools/bgzip/.fai indexing on the source file.
    local smoke_ref="$(pwd)/smoke_ref_chr22.fa.gz"
    if [[ ! -f "$smoke_ref" ]]; then
        echo "Extracting chr22 from $REF_FASTA → $smoke_ref"
        if [[ "$REF_FASTA" == *.gz ]]; then
            zcat "$REF_FASTA"
        else
            cat "$REF_FASTA"
        fi | awk -v target="chr22" '
            /^>/ {
                # Strip ">" and take the first whitespace-delimited token.
                split(substr($0, 2), parts, /[ \t]+/);
                keep = (parts[1] == target);
            }
            keep { print }
        ' | gzip > "$smoke_ref"
        if [[ ! -s "$smoke_ref" ]] || [[ "$(zcat "$smoke_ref" | wc -l)" -lt 2 ]]; then
            echo "ERROR: chr22 extraction produced an empty file. The reference may use" >&2
            echo "       a different contig name. Check with:" >&2
            echo "       zcat $REF_FASTA | grep -m5 '^>'" >&2
            return 1
        fi
    fi
    local config="gen.yml"
    cat > "$config" <<EOF
reference: $smoke_ref
read_len: 150
coverage: 10
ploidy: 2
paired_ended: false
mutation_model: $TRAINED_MODEL
sv_rate_scale: 1.0
output_dir: $(pwd)/reads_out
output_filename: validation
produce_fastq: true
produce_vcf: true
overwrite_output: true
rng_seed: gnomad-validation
EOF
    mkdir -p reads_out
    local start=$(date +%s)
    eidolon gen-reads -c "$config" 2>&1 | tail -20
    local end=$(date +%s)
    echo "gen-reads wall clock: $((end - start))s"
    local n_sv
    n_sv="$(zcat reads_out/validation.vcf.gz 2>/dev/null \
        | grep -v '^#' \
        | grep -E '<DEL>|<DUP>|<CNV>|<INV>|[[][^]]*:[^]]*[[|]][^]]*:[^]]*[]]' \
        | wc -l)"
    echo "SVs in output golden VCF: $n_sv"
    if [[ "$n_sv" -eq 0 ]]; then
        echo "WARN: zero de novo SVs in output — sampler may have rejected everything."
    fi
}

# ── Main ────────────────────────────────────────────────────────────────
case "$MODE" in
    smoke|full) ;;
    *) echo "usage: $0 {smoke|full}" >&2; exit 1 ;;
esac

ensure_vcf
ensure_ref
prep_input
train
inspect
smoke_gen_reads
echo
echo "Validation run complete. See $WORK/train.timing.log for full /usr/bin/time output."
