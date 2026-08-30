#!/usr/bin/env bash
# Contig selection for the realism panel, shared by realism_panel.sbatch and its tests.
#
# In its own file rather than inline in the sbatch so the test exercises the SAME code the
# job runs. A test that re-extracts a block with sed tests a copy, and when the pattern
# stops matching it tests an empty string and passes — which is how the summariser's first
# test suite passed every assertion against a zero-byte awk program.
#
# The check this makes: the contig must carry reads in the real BAM AND exist in the
# reference AT THE SAME LENGTH. Neither half is optional.
#
#   * Picked from the BAM alone, a whole-genome BAM yields chr1 while REFERENCE defaults to
#     chr22 — the job would simulate one contig and measure another, and would not say so
#     until the panel hit an unknown contig in the simulated BAM, hours later, after the
#     whole simulate-and-align spend.
#   * Same name and different length means different builds. Every coordinate in the region
#     BED would then name different sequence on each side, and the job would still print a
#     tidy comparison table. "Comparing across references compares the references" is the
#     panel's own premise; this is where that gets enforced instead of documented.

# select_contig <real_bam> <reference_fai> <idxstats_file> [requested_contig]
# Prints "<contig>\t<length>" on stdout. Diagnostics go to stderr. Non-zero on any failure.
select_contig() {
    local bam="$1" fai="$2" idxstats="$3" requested="${4:-}"
    local contig bam_len ref_len

    if [[ -n "$requested" ]]; then
        awk -v c="$requested" '$1==c && $3>0 {found=1} END{exit !found}' "$idxstats" || {
            echo "ERROR: CONTIG=$requested carries no reads in $bam" >&2
            return 1
        }
        contig="$requested"
    else
        contig="$(awk 'NR==FNR {ref[$1]=1; next} $3>0 && ($1 in ref) {print $1; exit}' \
                  "$fai" "$idxstats")"
        if [[ -z "$contig" ]]; then
            echo "ERROR: no contig carries reads in $bam and also exists in the reference." >&2
            echo "       No sequenced contig is common to both, which usually means different" >&2
            echo "       builds or different naming conventions (chr22 vs 22). Comparing across" >&2
            echo "       references compares the references, so this is fatal rather than" >&2
            echo "       something to work around." >&2
            echo "       BAM contigs with reads: $(awk '$3>0 {printf "%s ", $1}' "$idxstats" | head -c 200)" >&2
            echo "       Reference contigs:      $(awk '{printf "%s ", $1}' "$fai" | head -c 200)" >&2
            return 1
        fi
    fi

    bam_len="$(awk -v c="$contig" '$1==c {print $2; exit}' "$idxstats")"
    ref_len="$(awk -v c="$contig" '$1==c {print $2; exit}' "$fai")"

    if [[ -z "$ref_len" ]]; then
        echo "ERROR: $contig is not present in the reference." >&2
        return 1
    fi
    if [[ "$bam_len" != "$ref_len" ]]; then
        echo "ERROR: $contig is $bam_len bp in $bam but $ref_len bp in the reference." >&2
        echo "       Same name, different length: these are different builds. The regions" >&2
        echo "       would not name the same sequence on each side, and the resulting table" >&2
        echo "       would be a comparison of two references rather than of two read sets." >&2
        return 1
    fi

    printf '%s\t%s\n' "$contig" "$ref_len"
}
