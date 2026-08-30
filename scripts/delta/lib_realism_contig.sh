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

# select_contigs <real_bam> <reference_fai> <idxstats_file> [requested_contig]
# Prints one "<contig>\t<length>" line per usable contig, in reference order. Diagnostics
# go to stderr. Non-zero on any failure.
#
# ALL shared contigs, not the first one. The simulation covers the whole reference whatever
# happens, so measuring one contig of a three-chromosome reference pays for three and reads
# one. Spreading the loci across every shared contig is free and gives the baseline its
# spread across GC and repeat contexts rather than within a single chromosome.
select_contigs() {
    local bam="$1" fai="$2" idxstats="$3" requested="${4:-}"
    local contigs contig bam_len ref_len emitted=0

    if [[ -n "$requested" ]]; then
        awk -v c="$requested" '$1==c && $3>0 {found=1} END{exit !found}' "$idxstats" || {
            echo "ERROR: CONTIG=$requested carries no reads in $bam" >&2
            return 1
        }
        contigs="$requested"
    else
        # Reference order, so the output is stable and reads the way the reference does.
        contigs="$(awk 'NR==FNR { if ($3 > 0) seq[$1] = 1; next } ($1 in seq) { print $1 }' \
                   "$idxstats" "$fai")"
        if [[ -z "$contigs" ]]; then
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

    for contig in $contigs; do
        bam_len="$(awk -v c="$contig" '$1==c {print $2; exit}' "$idxstats")"
        ref_len="$(awk -v c="$contig" '$1==c {print $2; exit}' "$fai")"

        if [[ -z "$ref_len" ]]; then
            echo "ERROR: $contig is not present in the reference." >&2
            return 1
        fi
        # A shared NAME with a different LENGTH means different builds. Fatal rather than
        # skipped: skipping it would quietly measure a subset of what was asked for, and
        # the mismatch is evidence about the whole pairing, not about one contig.
        if [[ "$bam_len" != "$ref_len" ]]; then
            echo "ERROR: $contig is $bam_len bp in $bam but $ref_len bp in the reference." >&2
            echo "       Same name, different length: these are different builds. The regions" >&2
            echo "       would not name the same sequence on each side, and the resulting table" >&2
            echo "       would be a comparison of two references rather than of two read sets." >&2
            return 1
        fi

        printf '%s\t%s\n' "$contig" "$ref_len"
        emitted=$((emitted + 1))
    done

    [[ "$emitted" -gt 0 ]] || { echo "ERROR: no usable contig found." >&2; return 1; }
}

# place_regions <contigs.tsv> <n_regions> <region_bp> <margin> <out.bed>
# Writes a BED of evenly spaced regions distributed across the contigs. Prints the number
# placed and the number of eligible contigs as "<placed>\t<eligible>". Non-zero if it
# cannot place any.
#
# Distribution, not concentration: the simulation covers the whole reference regardless of
# what gets measured, so on a three-chromosome reference measuring one contig pays for
# three and reads one. The extra loci also come from different GC and repeat contexts,
# which is the spread the baseline exists to record.
place_regions() {
    local contigs="$1" n_regions="$2" region_bp="$3" margin="$4" out="$5"
    local c l eligible base extra idx want usable i start placed n_elig skipped=0

    eligible="$(mktemp)"
    while IFS=$'\t' read -r c l; do
        [[ -n "$c" ]] || continue
        if [[ $(( l - 2 * margin - region_bp )) -gt 0 ]]; then
            printf '%s\t%s\n' "$c" "$l" >> "$eligible"
        else
            # Listed, not dropped in silence. A reference of short contigs would otherwise
            # reduce to a handful of loci with the table still reading as if it had all of
            # them, which is the shape of every quiet failure in this repo so far.
            echo "  NOTE: $c ($l bp) is too short for a ${region_bp} bp region with ${margin} bp margins — skipped" >&2
            skipped=$((skipped + 1))
        fi
    done < "$contigs"

    n_elig="$(wc -l < "$eligible")"
    if [[ "$n_elig" -eq 0 ]]; then
        rm -f "$eligible"
        echo "ERROR: no contig is long enough for a ${region_bp} bp region with ${margin} bp margins." >&2
        echo "       Lower REGION_BP, or use a reference with longer contigs. Measuring" >&2
        echo "       nothing and reporting 0.0 artifacts are indistinguishable in a table," >&2
        echo "       so this is fatal rather than a warning." >&2
        return 1
    fi

    # Remainder to the earliest contigs, so the total is exactly n_regions whenever the
    # eligible contigs can hold them.
    base=$(( n_regions / n_elig ))
    extra=$(( n_regions % n_elig ))
    idx=0
    : > "$out"
    while IFS=$'\t' read -r c l; do
        want=$base
        # An `if`, not `[[ ... ]] && want=$((want+1))`: that list returns 1 when the test is
        # false, and under `set -e` a bare failing list exits the job. It would have killed
        # every run with more contigs than remainder.
        if [[ "$idx" -lt "$extra" ]]; then want=$(( want + 1 )); fi
        idx=$(( idx + 1 ))
        if [[ "$want" -le 0 ]]; then continue; fi
        usable=$(( l - 2 * margin - region_bp ))
        for ((i=0; i<want; i++)); do
            start=$(( margin + (usable / want) * i ))
            printf '%s\t%d\t%d\n' "$c" "$start" "$(( start + region_bp ))" >> "$out"
        done
    done < "$eligible"
    rm -f "$eligible"

    placed="$(wc -l < "$out")"
    if [[ "$placed" -eq 0 ]]; then
        echo "ERROR: placed no regions despite $n_elig eligible contig(s)." >&2
        return 1
    fi
    printf '%s\t%s\n' "$placed" "$n_elig"
}
