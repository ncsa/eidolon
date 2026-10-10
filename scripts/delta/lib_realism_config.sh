#!/usr/bin/env bash
# Simulation-config generation for the realism panel, shared by realism_panel.sbatch and
# its tests.
#
# In its own file for the same reason select_contigs is: so the test exercises the SAME
# code the job runs. Grepping the sbatch for a `fragment_model:` line tests that the string
# is present, not that it is emitted under the right conditions -- and the bug this exists
# to prevent is precisely a conditional one.

# frag_model_ceiling <model.json.gz>
#
# Largest fragment length the model can produce. Discrete models have a hard maximum
# (the builder trims outliers); a Normal is unbounded, so mean + 4 sd is used as the
# practical ceiling. Prints nothing and returns non-zero when it cannot tell.
frag_model_ceiling() {
    local model="$1"
    [[ -s "$model" ]] || return 1
    command -v jq >/dev/null 2>&1 || return 1
    local v
    v="$(zcat "$model" 2>/dev/null \
         | jq -r 'if has("Discrete") then (.Discrete.distribution.values | max)
                  else ((.Normal.mean + 4 * .Normal.st_dev) | floor) end' 2>/dev/null)"
    [[ -n "$v" && "$v" != "null" ]] || return 1
    printf '%s\n' "$v"
}

# write_sim_config <out.yml> <reference> <outdir> <seed> <threads> <depth> <read_len> \
#                  <frag_mean> <frag_sd>
#
# Optional model paths are read from the environment: GC_BIAS_MODEL, FRAGMENT_MODEL,
# SEQ_ERROR_MODEL, QUALITY_MODEL, MUTATION_MODEL, GC_NORMALIZE. Unset means "let eidolon
# use its built-in default", which is a real choice and not the same as a trained model.
#
# ADAPTERS (unset | truseq | nextera) is the same kind of choice and belongs in the same
# place. Adapter read-through is a DIRECT source of soft clips -- eidolon's own
# `adapter_readthrough_is_soft_clipped_in_the_bam` asserts exactly that -- and the panel
# measured clip_pct 19x below real data with it switched off. Enabling it also flips
# `keep_short` in gen_reads, so fragments shorter than a read stop being rejected by the
# retry loop and become clipped reads instead; with a trained fragment model that is 0.557%
# of the mass on HCC1395 normal. Whether that lands near the real 0.19% or overshoots is a
# measurement, which is the point of exposing the knob rather than picking a value here.
# SV_RATE_SCALE (default 0.0) multiplies the SV model's per_base_rate. It is exposed
# because `cand_per_mb` compares a real genome -- which carries structural variants --
# against a simulation that plants none, so the metric was never measuring the same thing
# on both sides. See #684.
#
# CALIBRATION WARNING. The bundled model's per_base_rate (6.009e-4) is fit as
# `n / total_reflen_seen` from the gnomAD-SV SITES VCF: 1,832,360 records over ~3.05 Gb.
# A sites VCF is the UNION across the cohort, not one genome's complement, so scale 1.0
# asks for ~1.83M SVs in a single simulated genome against a real per-genome count on the
# order of 1e4. Sweep this empirically; do not assume 1.0 is realistic.
write_sim_config() {
    local out="$1" reference="$2" outdir="$3" seed="$4" threads="$5" depth="$6" \
          read_len="$7" frag_mean="$8" frag_sd="$9"

    cat > "$out" <<YML
reference: $reference
output_dir: $outdir
output_filename: sim
rng_seed: "$seed"
num_threads: $threads
coverage: $depth
read_len: $read_len
paired_ended: true
produce_fastq: true
sv_rate_scale: ${SV_RATE_SCALE:-0.0}
YML

    # fragment_mean/st_dev and fragment_model are two sources for the same thing, and
    # gen_reads/utils/runner.rs takes the explicit mean/st_dev path when they are present.
    # That is how the panel silently OVERRODE eidolon's own shipped empirical fragment
    # distribution with Normal(400, 90) on every run it ever did -- and then reported the
    # resulting symmetric insert distribution as a realism gap. Emit one or the other.
    if [[ -n "${FRAGMENT_MODEL:-}" ]]; then
        printf 'fragment_model: %s\n' "$FRAGMENT_MODEL" >> "$out"
    else
        printf 'fragment_mean: %s\nfragment_st_dev: %s\n' "$frag_mean" "$frag_sd" >> "$out"
    fi

    local pair key val
    # INPUT_VCF supplies variants directly rather than through a model. Drawn germline SVs
    # (tools/draw_gnomad_sv_vcf.sh) go here: cand_per_mb needs breakpoints at real loci, and
    # supplied variants are ADDED to the de novo mutations rather than replacing them, so the
    # run still carries the model's SNPs and small indels.
    for pair in "input_vcf:${INPUT_VCF:-}" \
                "gc_bias_model:${GC_BIAS_MODEL:-}" \
                "sequence_error_model:${SEQ_ERROR_MODEL:-}" \
                "quality_score_model:${QUALITY_MODEL:-}" \
                "mutation_model:${MUTATION_MODEL:-}" \
                "gc_bias_normalize_coverage:${GC_NORMALIZE:-}"; do
        key="${pair%%:*}"; val="${pair#*:}"
        [[ -n "$val" ]] && printf '%s: %s\n' "$key" "$val" >> "$out"
    done

    # A nested map, not a scalar, so it cannot ride the loop above.
    if [[ -n "${ADAPTERS:-}" ]]; then
        printf 'adapters:\n  enabled: true\n  preset: %s\n' "${ADAPTERS:-}" >> "$out"
    fi
    return 0
}

# write_settings <out.tsv> <eidolon_version> <git> <real_bam> <reference> <real_depth> \
#                <sim_depth>
#
# The effective settings of a panel run, as key/value lines archived beside panel.tsv. Every
# number the panel reports is conditional on these, and the job log that used to be their
# only record does not survive: job 22583881 could not be reproduced once its log was gone,
# and a later run differed in real depth, fragment source and adapters without anything on
# file saying so. The real BAM's size is recorded because a restaged BAM keeps its path.
# Everything else is read from the same environment write_sim_config reads.
write_settings() {
    local out="$1" version="$2" git="$3" real_bam="$4" reference="$5"
    local real_depth="$6" sim_depth="$7" frag bytes
    frag="Normal(${FRAG_MEAN:-}, ${FRAG_SD:-})"
    [[ -n "${FRAGMENT_MODEL:-}" ]] && frag="model ${FRAGMENT_MODEL}"
    bytes="$(stat -c %s "$real_bam" 2>/dev/null || echo unknown)"
    {
        printf 'key\tvalue\n'
        printf 'eidolon\t%s\n' "$version"
        printf 'git\t%s\n' "$git"
        printf 'real_bam\t%s\n' "$real_bam"
        printf 'real_bam_bytes\t%s\n' "$bytes"
        printf 'reference\t%s\n' "$reference"
        printf 'align_reference\t%s\n' "${ALIGN_REFERENCE:-same as reference}"
        printf 'contig\t%s\n' "${CONTIG:-every contig with reads}"
        printf 'n_regions\t%s\n' "${N_REGIONS:-}"
        printf 'region_bp\t%s\n' "${REGION_BP:-}"
        printf 'seed\t%s\n' "${SEED:-}"
        printf 'read_len\t%s\n' "${READ_LEN:-}"
        printf 'match_read_len\t%s\n' "${MATCH_READ_LEN:-}"
        printf 'read_len_tol\t%s\n' "${READ_LEN_TOL:-}"
        printf 'real_depth\t%s\n' "$real_depth"
        printf 'sim_depth\t%s\n' "$sim_depth"
        printf 'max_sim_depth\t%s\n' "${MAX_SIM_DEPTH:-}"
        printf 'fragment\t%s\n' "$frag"
        printf 'gc_bias_model\t%s\n' "${GC_BIAS_MODEL:-eidolon default}"
        printf 'sequence_error_model\t%s\n' "${SEQ_ERROR_MODEL:-eidolon default}"
        printf 'quality_score_model\t%s\n' "${QUALITY_MODEL:-eidolon default}"
        printf 'mutation_model\t%s\n' "${MUTATION_MODEL:-eidolon default}"
        printf 'gc_normalize\t%s\n' "${GC_NORMALIZE:-eidolon default}"
        printf 'adapters\t%s\n' "${ADAPTERS:-off}"
        printf 'min_clip\t%s\n' "${MIN_CLIP:-}"
        printf 'min_support\t%s\n' "${MIN_SUPPORT:-}"
        printf 'max_tlen\t%s\n' "${MAX_TLEN:-}"
        printf 'depth_lag\t%s\n' "${DEPTH_LAG:-}"
    } > "$out"
}
