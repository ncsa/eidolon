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
# SEQ_ERROR_MODEL, QUALITY_MODEL, MUTATION_MODEL, GC_NORMALIZE. PLOIDY sets the simulated
# genome's ploidy (gen-reads defaults to 2); a haploid sample such as E. coli needs 1 (#844).
# MUTATION_RATE overrides the mutation model's rate. MUTATION_RATE=0 with INPUT_VCF gives a
# simulated genome carrying exactly the supplied variants: input_vcf otherwise ADDS to the
# de novo mutations (#837).
# Unset means "let eidolon use its built-in default", which is a real choice and not the
# same as a trained model.
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
                "gc_bias_normalize_coverage:${GC_NORMALIZE:-}" \
                "ploidy:${PLOIDY:-}" \
                "mutation_rate:${MUTATION_RATE:-}"; do
        key="${pair%%:*}"; val="${pair#*:}"
        [[ -n "$val" ]] && printf '%s: %s\n' "$key" "$val" >> "$out"
    done

    # A nested map, not a scalar, so it cannot ride the loop above.
    if [[ -n "${ADAPTERS:-}" ]]; then
        printf 'adapters:\n  enabled: true\n  preset: %s\n' "${ADAPTERS:-}" >> "$out"
    fi
    return 0
}

# align_sim_reads <reference> <r1> <r2> <out.bam> <bwa.log> <threads>
#
# Aligns the simulated reads and refuses to hand on an incomplete BAM. Under `set -euo pipefail`
# a failed aligner used to end the job with exit 1 and no explanation: on job 23251798 bwa-mem2
# hit a transient read failure on /scratch after 48 of ~97 read batches, `samtools sort` wrote
# the half it had received (122,880,000 of 249,098,168 reads), and the job died with nothing in
# its log but a sort message. Diagnosing it took four round trips. This names the step and
# shows the aligner's last lines instead.
#
# BWA_BIN and SAMTOOLS_BIN exist so tests can substitute a failing aligner.
align_sim_reads() {
    local ref="$1" r1="$2" r2="$3" bam="$4" log="$5" threads="$6"
    local bwa="${BWA_BIN:-bwa-mem2}" samtools="${SAMTOOLS_BIN:-samtools}" st was_e=0
    [[ $- == *e* ]] && was_e=1
    set +e
    "$bwa" mem -t "$threads" "$ref" "$r1" "$r2" 2> "$log" \
        | "$samtools" sort -@ 8 -o "$bam" -
    st=("${PIPESTATUS[@]}")
    [[ "$was_e" -eq 1 ]] && set -e
    if [[ "${st[0]}" -ne 0 || "${st[1]}" -ne 0 ]]; then
        echo "ERROR: aligning the simulated reads failed (bwa-mem2 exit ${st[0]}, samtools sort exit ${st[1]})." >&2
        echo "  $bam is incomplete and will not be measured. Last lines of $log:" >&2
        tail -5 "$log" | sed 's/^/    /' >&2
        echo "  '[gzread] Success' there means a read from a FASTQ failed mid-stream; on job" >&2
        echo "  23251798 the file itself was intact (gzip -t passed), so check the filesystem and rerun." >&2
        return 1
    fi
    return 0
}
