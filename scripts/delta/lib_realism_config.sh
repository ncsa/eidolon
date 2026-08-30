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
sv_rate_scale: 0.0
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
    for pair in "gc_bias_model:${GC_BIAS_MODEL:-}" \
                "sequence_error_model:${SEQ_ERROR_MODEL:-}" \
                "quality_score_model:${QUALITY_MODEL:-}" \
                "mutation_model:${MUTATION_MODEL:-}" \
                "gc_bias_normalize_coverage:${GC_NORMALIZE:-}"; do
        key="${pair%%:*}"; val="${pair#*:}"
        [[ -n "$val" ]] && printf '%s: %s\n' "$key" "$val" >> "$out"
    done
    return 0
}
