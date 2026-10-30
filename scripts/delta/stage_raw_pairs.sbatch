#!/usr/bin/env bash
# Stage an aligned, duplicate-marked BAM from a run's FASTQs on ENA, for
# fit_overlap_matrix.sbatch (#779).
#
# WHY: the overlap fit counts disagreements between mates, so anything that edits read bases
# before alignment changes the answer. SEQC2's pre-aligned HCC1395 normal BAM was trimmed
# upstream (read lengths 32-151 bp) and gave a disagreement rate that FALLS 4x along the read,
# 0.059% -> 0.014%, where HG002 rises 0.10% -> 0.98%. Its @PG lines stop at bwa, so they
# cannot say what ran before it. Aligning a run whose reads are untouched removes the question.
#
# CHECK THE READS ARE RAW. An SRA upload can be processed too: SEQC2's own run, SRR7890943,
# is trimmed (min 35 bp, 36.5% under 151) though its mean length is 150.4. A mean hides this,
# so the pull prints the length distribution; untrimmed reads are all one length.
#
# Streams the first N_PAIRS pairs of each mate (no full download), aligns to REFERENCE with
# bwa-mem2, then fixmate -m | sort | markdup so duplicates carry the flag the fit skips.
# The head of a FASTQ is the first tiles of the flowcell, not a random sample; say so when
# quoting a number from it.
#
# USAGE
#   SRR=SRR14724533 TAG=hg001_novaseq sbatch scripts/delta/stage_raw_pairs.sbatch
# then fit, queued behind it:
#   BAM=$SCRATCH/raw_pairs/<TAG>/<TAG>.bam TAG=<TAG> sbatch --dependency=afterok:<id> scripts/delta/fit_overlap_matrix.sbatch
# and check the exit state, not just the log:
#   sacct -j <id> --format=JobID,State,ExitCode,Elapsed
#
# LOCAL READS. If $WORK/r1.fastq.gz and r2.fastq.gz already exist, nothing is downloaded and
# SRR is not needed: this aligns them through the same bwa-mem2 | fixmate | sort | markdup path,
# which is how simulated reads get compared with real ones on equal terms (#844, #861). N_PAIRS
# then defaults to the pair count of r1, checked against r2, so it need not be known at submit.
#   TAG=chm13_sim REFERENCE=... sbatch scripts/delta/stage_raw_pairs.sbatch

#SBATCH --job-name=eidolon-rawpairs
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=rawpairs-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
source "$REPO/scripts/delta/lib_report.sh"

TAG="${TAG:?set TAG to a short library name}"
REFERENCE="${REFERENCE:-$SCRATCH/neat_data/GRCh38.fa}"
WORK="${WORK:-$SCRATCH/raw_pairs/$TAG}"
LOCAL_READS=0
[[ -s "$WORK/r1.fastq.gz" && -s "$WORK/r2.fastq.gz" ]] && LOCAL_READS=1
THREADS="${SLURM_CPUS_PER_TASK:-16}"

module load samtools/1.22-cce19.0.0 2>/dev/null || true
setup_conda
conda_activate "${CONDA_ENV:-bioinf}"   # bwa-mem2
command -v samtools >/dev/null || { echo "FATAL: samtools not on PATH" >&2; exit 1; }
command -v bwa-mem2 >/dev/null || { echo "FATAL: bwa-mem2 not on PATH" >&2; exit 1; }
[[ -s "$REFERENCE.bwt.2bit.64" ]] || { echo "FATAL: no bwa-mem2 index for $REFERENCE" >&2; exit 1; }
mkdir -p "$WORK"

if [[ "$LOCAL_READS" -eq 1 ]]; then
    echo "=== stage_raw_pairs: $TAG, reads already in $WORK ==="
    if [[ -z "${N_PAIRS:-}" ]]; then
        n1=$(zcat "$WORK/r1.fastq.gz" | wc -l)
        n2=$(zcat "$WORK/r2.fastq.gz" | wc -l)
        [[ "$n1" -eq "$n2" && $((n1 % 4)) -eq 0 && "$n1" -gt 0 ]] \
            || { echo "FATAL: r1 has $n1 lines and r2 has $n2; not a matched pair of FASTQs" >&2; exit 1; }
        N_PAIRS=$((n1 / 4))
    fi
    echo "pairs: $N_PAIRS"
else
    N_PAIRS="${N_PAIRS:-40000000}"
    SRR="${SRR:?set SRR to an ENA/SRA run accession, or put r1/r2.fastq.gz in $WORK}"
    # ENA's report gives the FASTQ paths, so the SRR's directory layout need not be derived.
    urls=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=fastq_ftp" | awk -F'\t' 'NR == 2 { print $2 }')
    IFS=';' read -r u1 u2 <<< "$urls"
    [[ -n "${u1:-}" && -n "${u2:-}" ]] || { echo "FATAL: ENA lists no paired FASTQs for $SRR: '$urls'" >&2; exit 1; }
    echo "=== stage_raw_pairs: $TAG ($SRR), first $N_PAIRS pairs ==="
    echo "r1: $u1"
    echo "r2: $u2"

    # `head` closes the pipe early, so curl and zcat take SIGPIPE; pipefail comes off for these
    # pipelines only. The line counts below are what establish the pull finished.
    lines=$((N_PAIRS * 4))
    for m in 1 2; do
        u="u$m"
        fq="$WORK/r$m.fastq.gz"
        [[ -s "$fq" ]] && continue
        set +o pipefail
        curl -sf "https://${!u}" | zcat | head -n "$lines" | gzip -1 -c > "$fq.part"
        set -o pipefail
        got=$(zcat "$fq.part" | wc -l)
        [[ "$got" -eq "$lines" ]] || { echo "FATAL: r$m has $got lines, expected $lines" >&2; exit 1; }
        mv "$fq.part" "$fq"
    done
fi
echo "read lengths in r1 (count, bp), most common first:"
zcat "$WORK/r1.fastq.gz" | awk 'NR % 4 == 2 { print length($0) }' | sort -n | uniq -c | sort -rn | head -5

BAM="$WORK/$TAG.bam"
if [[ ! -s "$BAM.bai" ]]; then
    rg="@RG\tID:$TAG\tSM:$TAG\tLB:$TAG\tPL:ILLUMINA"
    bwa-mem2 mem -t "$THREADS" -R "$rg" "$REFERENCE" "$WORK/r1.fastq.gz" "$WORK/r2.fastq.gz" 2> "$WORK/bwa.log" \
        | samtools fixmate -m -@ 4 - - \
        | samtools sort -@ 8 -m 2G -T "$WORK/sort" - \
        | samtools markdup -@ 8 -s -f "$WORK/markdup.txt" - "$BAM"
    samtools index -@ 8 "$BAM"
fi

# Every pair must have come through: primary records = 2 x N_PAIRS.
primary=$(samtools view -c -F 0x900 "$BAM")
[[ "$primary" -eq $((N_PAIRS * 2)) ]] || { echo "FATAL: $primary primary records, expected $((N_PAIRS * 2))" >&2; exit 1; }
samtools flagstat -@ 8 "$BAM" > "$WORK/flagstat.txt"
cat "$WORK/flagstat.txt"

archive_run "raw_pairs/$TAG" "$WORK" flagstat.txt markdup.txt bwa.log
echo "BAM: $BAM"
echo "DONE."
