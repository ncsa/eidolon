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

# n_mask <reference.fa> <contigs.tsv> [min_run] > mask.bed
#
# Emits "<contig>\t<start>\t<end>" (0-based half-open) for every run of N at least
# <min_run> long, for the contigs listed in contigs.tsv. This is what place_regions needs
# in order to not put a window on sequence that does not exist.
#
# WHY THIS EXISTS. Job 21622644 died in the measure step on chr21:2000000-2400000, which is
# inside the acrocentric p-arm — several megabases of N in GRCh38 because that stretch of
# the chromosome is not in the reference. place_regions had no way to know: it places by
# coordinate and only ever saw contig LENGTHS. Every acrocentric chromosome (13, 14, 15,
# 21, 22) has the same hole at the same place, so the panel could not run on a reference
# containing any of them. The panel itself was right to refuse — reader.rs treats a
# read-less region as an error precisely because zero artifacts and zero reads are
# indistinguishable downstream. The bug was upstream, in what it was asked to measure.
#
# RESOLUTION IS ONE FASTA LINE (~60 bp), NOT ONE BASE. A whole-line test is one regex per
# line; going per-base means 162 million substr() calls for a three-chromosome reference,
# which took minutes rather than seconds when measured. Boundaries can therefore be off by
# up to a line width. That is deliberate and it is enough: the blocks being avoided are
# megabases wide, and min_run keeps small gaps out of the mask entirely so a stray N in an
# otherwise good window does not cost us the window.
n_mask() {
    local ref="$1" contigs="$2" min_run="${3:-10000}"
    local c

    [[ -s "$ref" ]] || { echo "ERROR: n_mask: reference not found: $ref" >&2; return 1; }
    command -v samtools >/dev/null 2>&1 \
        || { echo "ERROR: n_mask: samtools is required to read the reference." >&2; return 1; }

    while IFS=$'\t' read -r c _; do
        [[ -n "$c" ]] || continue
        samtools faidx "$ref" "$c" \
          | awk -v ctg="$c" -v minrun="$min_run" '
                BEGIN { pos = 0 }
                NR == 1 { next }                       # the > header
                {
                    alln = ($0 ~ /^[Nn]+$/)
                    if (alln && !inrun)        { inrun = 1; rs = pos }
                    else if (!alln && inrun)   { if (pos - rs >= minrun) print ctg "\t" rs "\t" pos; inrun = 0 }
                    pos += length($0)
                }
                END { if (inrun && pos - rs >= minrun) print ctg "\t" rs "\t" pos }'
    done < "$contigs"
}

# place_regions <contigs.tsv> <n_regions> <region_bp> <margin> <out.bed> [mask_bed]
#
# Writes a BED of evenly spaced regions distributed across the contigs, and prints
# "<placed>\t<eligible>". Non-zero if it cannot place any.
#
# Distribution, not concentration: the simulation covers the whole reference regardless of
# what gets measured, so on a three-chromosome reference measuring one contig pays for
# three and reads one. The extra loci also come from different GC and repeat contexts,
# which is the spread the baseline exists to record.
#
# With no mask_bed the arithmetic is exactly what it has
# always been -- evenly spaced starts from the margin -- so existing runs do not shift.
#
# With a mask_bed, candidates are oversampled across the usable span, any candidate
# OVERLAPPING a masked interval at all is discarded, and what survives is thinned back to
# the requested count. Rejecting on any overlap rather than on a fraction is deliberate:
# min_run already keeps trivial gaps out of the mask, so anything left is a real hole, and
# there are far more candidates than needed. <placed> can now be LESS than <n_regions> --
# callers must read it rather than assume they got what they asked for.
place_regions() {
    local contigs="$1" n_regions="$2" region_bp="$3" margin="$4" out="$5" mask="${6:-}"
    local c l eligible base extra idx want usable i start placed n_elig skipped=0
    local cand kept ncand step nkept dropped=0

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

        if [[ -z "$mask" ]]; then
            for ((i=0; i<want; i++)); do
                start=$(( margin + (usable / want) * i ))
                printf '%s\t%d\t%d\n' "$c" "$start" "$(( start + region_bp ))" >> "$out"
            done
            continue
        fi

        # Oversample, drop anything touching a hole, then thin back to `want`. 6x the
        # request (floor 24) so a contig can lose most of its span to N -- chr22 loses its
        # first 10 Mb -- and still yield the full count from what is left.
        ncand=$(( want * 6 )); [[ "$ncand" -lt 24 ]] && ncand=24
        [[ "$ncand" -gt "$usable" ]] && ncand="$usable"
        cand="$(mktemp)"; kept="$(mktemp)"
        step=$(( usable / ncand )); [[ "$step" -lt 1 ]] && step=1
        for ((i=0; i<ncand; i++)); do echo $(( margin + step * i )); done > "$cand"

        # Mask first, candidates second. Half-open overlap: s < mask_end && s+bp > mask_start.
        awk -F'\t' -v C="$c" -v BP="$region_bp" -v MASK="$mask" '
            FILENAME == MASK { if ($1 == C) { m++; ms[m] = $2; me[m] = $3 } next }
            {
                s = $1; e = s + BP; bad = 0
                for (i = 1; i <= m; i++) if (s < me[i] && e > ms[i]) { bad = 1; break }
                if (!bad) print s
            }' "$mask" "$cand" > "$kept"

        nkept=$(wc -l < "$kept")
        dropped=$(( dropped + ncand - nkept ))
        if [[ "$nkept" -eq 0 ]]; then
            # Rule 4: a contig contributing nothing is a fact about the measurement, and
            # silence here would shrink the sample while the table read as if it had not.
            echo "  NOTE: $c contributed no region — every candidate fell in an N block" >&2
            rm -f "$cand" "$kept"; continue
        fi
        [[ "$nkept" -lt "$want" ]] && \
            echo "  NOTE: $c yielded only $nkept of $want regions clear of N" >&2

        awk -v k="$want" -v c="$c" -v w="$region_bp" '
            { a[NR] = $1 }
            END {
                if (NR == 0) exit
                n = (NR < k) ? NR : k
                for (i = 0; i < n; i++) {
                    idx = (n == 1) ? 1 : int(1 + i * (NR - 1) / (n - 1))
                    printf "%s\t%d\t%d\n", c, a[idx], a[idx] + w
                }
            }' "$kept" >> "$out"
        rm -f "$cand" "$kept"
    done < "$eligible"

    if [[ -n "$mask" && "$dropped" -gt 0 ]]; then
        echo "  NOTE: discarded $dropped candidate window(s) overlapping an N block" >&2
    fi
    rm -f "$eligible"

    placed="$(wc -l < "$out")"
    if [[ "$placed" -eq 0 ]]; then
        echo "ERROR: placed no regions despite $n_elig eligible contig(s)." >&2
        return 1
    fi
    printf '%s\t%s\n' "$placed" "$n_elig"
}
