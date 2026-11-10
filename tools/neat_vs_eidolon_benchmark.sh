#!/usr/bin/env bash
# Quantitative benchmark: NEAT 4 (Python) vs eidolon (Rust), same input + params.
# Measures wall-clock time and peak RSS via /usr/bin/time -v, REPS times per
# config, reporting the median.
#
# Fairness: identical reference, coverage, read length, paired-end fragment
# params, and FASTQ-only output. RNG seeds cannot match across tools (different
# generators) so this compares resource usage, not output identity.
#
# Usage: tools/neat_vs_eidolon_benchmark.sh [outdir]
#   REPS=3 GENOMES=... THREAD_MODES=... override defaults via env.
set -uo pipefail

# ---- configuration ----------------------------------------------------------
EIDOLON_BIN="${EIDOLON_BIN:-$(cd "$(dirname "$0")/.." && pwd)/target/release/eidolon}"
NEAT_BIN="${NEAT_BIN:-$HOME/anaconda3/envs/neat/bin/neat}"
DATA_DIR="${DATA_DIR:-$HOME/code/data}"
GTIME="${GTIME:-/usr/bin/time}"

COVERAGE=10
READ_LEN=151
FRAG_MEAN=350
FRAG_SD=50
SEED=42
REPS="${REPS:-3}"

# genomes: "label:fasta_path"
GENOMES=(
  "ecoli:$DATA_DIR/ecoli.fa"
  "yeast:$DATA_DIR/yeast.fa"
  "c_elegans:$DATA_DIR/c_elegans.fa"
)
THREAD_MODES=(1 all)

OUTDIR="${1:-/tmp/neat_vs_eidolon_run}"
mkdir -p "$OUTDIR"
PERREP="$OUTDIR/per_rep.tsv"
MEDIAN="$OUTDIR/median.tsv"
echo -e "genome\tsize_mb\tthreads\ttool\trep\twall_s\tpeak_rss_mb\texit" > "$PERREP"
echo -e "genome\tsize_mb\tthreads\ttool\twall_s_med\tpeak_rss_mb_med\tn_ok" > "$MEDIAN"

NCORES=$(nproc)

# ---- helpers ----------------------------------------------------------------
parse_time() {  # /usr/bin/time -v file -> "wall_seconds peak_rss_mb"
  awk '
    /Elapsed \(wall clock\)/ {
      n=split($NF,a,":");
      if (n==3) wall=a[1]*3600+a[2]*60+a[3];
      else if (n==2) wall=a[1]*60+a[2];
      else wall=a[1];
    }
    /Maximum resident set size/ { rss=$NF/1024 }
    END { printf "%.2f %.1f", wall, rss }
  ' "$1"
}

median() {  # stdin: numbers, one per line -> median
  sort -n | awk '{v[NR]=$1} END{ if(NR==0){print "NA"; exit}
    if(NR%2){print v[(NR+1)/2]} else {printf "%.2f", (v[NR/2]+v[NR/2+1])/2} }'
}

run_config() {
  local label="$1" fasta="$2" mode="$3" tool="$4" cfg="$5" cmd_outdir="$6"
  local nthreads; [ "$mode" = "all" ] && nthreads=$NCORES || nthreads=$mode
  local size_mb; size_mb=$(awk -v b="$(stat -c%s "$fasta")" 'BEGIN{printf "%.1f", b/1048576}')
  local walls=() rsss=() ok=0

  for rep in $(seq 1 "$REPS"); do
    local tf="$OUTDIR/${tool}_${label}_t${mode}_r${rep}.time"
    local log="$OUTDIR/${tool}_${label}_t${mode}_r${rep}.log"
    rm -rf "$cmd_outdir"; mkdir -p "$cmd_outdir"
    echo ">>> [$tool] $label threads=$mode rep=$rep/$REPS"
    if [ "$tool" = "neat" ]; then
      "$GTIME" -v -o "$tf" "$NEAT_BIN" read-simulator -c "$cfg" -o "$cmd_outdir" -p bench > "$log" 2>&1
    else
      "$GTIME" -v -o "$tf" "$EIDOLON_BIN" gen-reads -c "$cfg" > "$log" 2>&1
    fi
    local ec=$? parsed wall rss
    parsed=$(parse_time "$tf"); wall=${parsed% *}; rss=${parsed#* }
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$size_mb" "$mode" "$tool" "$rep" "$wall" "$rss" "$ec" >> "$PERREP"
    echo "    wall=${wall}s rss=${rss}MB exit=$ec"
    if [ "$ec" -eq 0 ]; then walls+=("$wall"); rsss+=("$rss"); ok=$((ok+1)); fi
  done

  local wmed rmed
  wmed=$(printf "%s\n" "${walls[@]}" | median)
  rmed=$(printf "%s\n" "${rsss[@]}" | median)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$size_mb" "$mode" "$tool" "$wmed" "$rmed" "$ok" >> "$MEDIAN"
}

# ---- main loop --------------------------------------------------------------
for mode in "${THREAD_MODES[@]}"; do
  [ "$mode" = "all" ] && nthreads=$NCORES || nthreads=$mode
  for entry in "${GENOMES[@]}"; do
    label="${entry%%:*}"; fasta="${entry#*:}"
    [ -f "$fasta" ] || { echo "MISSING: $fasta"; continue; }

    neat_cfg="$OUTDIR/neat_${label}_t${mode}.yml"
    cat > "$neat_cfg" <<EOF
reference: $fasta
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_bam: false
produce_vcf: false
rng_seed: $SEED
threads: $nthreads
overwrite_output: true
EOF

    eidolon_out="$OUTDIR/eidolon_${label}_t${mode}"
    eidolon_cfg="$OUTDIR/eidolon_${label}_t${mode}.yml"
    cat > "$eidolon_cfg" <<EOF
reference: $fasta
read_len: $READ_LEN
coverage: $COVERAGE
paired_ended: true
fragment_mean: $FRAG_MEAN
fragment_st_dev: $FRAG_SD
produce_fastq: true
produce_vcf: false
produce_bam: false
overwrite_output: true
rng_seed: benchmark seed forty two
output_dir: $eidolon_out
output_basename: bench
EOF
    [ "$mode" != "all" ] && echo "num_threads: $nthreads" >> "$eidolon_cfg"

    run_config "$label" "$fasta" "$mode" "neat"  "$neat_cfg"  "$OUTDIR/neat_${label}_t${mode}"
    run_config "$label" "$fasta" "$mode" "eidolon" "$eidolon_cfg" "$eidolon_out"
  done
done

echo; echo "=== MEDIAN RESULTS (REPS=$REPS, ${NCORES} cores) ==="
column -t -s$'\t' "$MEDIAN"
