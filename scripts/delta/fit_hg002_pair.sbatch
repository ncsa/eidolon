#!/usr/bin/env bash
# Fit HG002's two mates into ONE model and report what #695 has been waiting for.
#
# WHY THIS EXISTS. #695 measured error_rate 0.002265 against the shipped 0.006638 and said
# plainly: do not adopt it, because it is R1 only and because a library-wide scalar may be the
# wrong shape when 12.22% of reads carry a collapsed tail. It named the sequence:
#
#   1. fit R2, report both
#   2. measure the rate separately for healthy vs degraded reads
#   3. adopt after #694, not before
#
# #723 shipped `fastq_file_r2:` and #720 shipped `fit_quality_degradation:`, so one job now
# answers 1 and 2 together. It also does the thing #723 recorded as unmeasured: whether a
# FOUR-population model (R1/R2 x healthy/degraded) reproduces a real library's mates, rather
# than a synthetic fixture's.
#
# REBUILD THE BINARY FIRST, and this script now checks rather than trusts. $EIDOLON is a build
# artifact in $SCRATCH, not the checkout: `git pull` does not rebuild it, and config is read by
# key lookup, so a binary predating #723 IGNORES `fastq_file_r2` instead of rejecting it and
# quietly fits one model for both mates. A binary predating #735 emits R2's quality profile
# backwards. Either produces a complete, plausible, wrong result. Step 0 refuses both.
#
#   cd $REPO && git pull && CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc \
#     CARGO_TARGET_DIR=$SCRATCH/cargo-target/eidolon cargo build --release
#
# WHY IT SUBSAMPLES, WHICH IS INTERIM. `gen-seq-error-model`'s `max_reads` takes the FIRST N
# records, and Illumina FASTQs are in flowcell order, so capping the fit reads one corner of one
# lane (#721). The stride below spans the file instead. Both mates are strided IDENTICALLY, so
# the subsamples remain mates. This goes away when the fitter can sample for itself.
#
# A NOTE ON WHAT CAN FAIL LOUDLY HERE, so it is not misread as a crash. Every population must
# reach the model's read length or the fit is REFUSED by name -- "the R2 population covers N bp,
# but the model covers M bp". That is the #720/#723 guard working, not a bug. If it fires, the
# two mates disagree on length and the subsample needs trimming, not the fitter.
#
# USAGE
#   sbatch scripts/delta/fit_hg002_pair.sbatch
#   STRIDE=200 sbatch scripts/delta/fit_hg002_pair.sbatch      # faster, coarser
#
# Then check the exit state, not just the log -- a killed job leaves partial output that reads
# like a short successful run:
#   sacct -j <id> --format=JobID,State,ExitCode,Elapsed

#SBATCH --job-name=eidolon-hg002pair
#SBATCH --partition=cpu
#SBATCH --account=bhrd-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=hg002pair-%j.log

set -euo pipefail

REPO="${REPO:-/projects/bhrd/jallen17/eidolon}"
# lib_report.sh resolves $SCRATCH, and this runs under `set -u`, so it is sourced BEFORE any
# path that mentions it -- which is every path below.
source "$REPO/scripts/delta/lib_report.sh"

# Defaults match where `stage_hg002.sh` WRITES, so the two scripts compose with no wiring:
#   ALIGN=0 sbatch scripts/delta/stage_hg002.sh
# puts hg002_R1/R2.fastq.gz in $SCRATCH/neat_data/hg002. (#695's earlier R1 fit read from
# /work/nvme/bhrd/jallen17/hg002 instead; pass DATA_DIR= to point at a library staged there.)
DATA_DIR="${DATA_DIR:-$SCRATCH/neat_data/hg002}"
R1="${R1:-$DATA_DIR/hg002_R1.fastq.gz}"
R2="${R2:-$DATA_DIR/hg002_R2.fastq.gz}"

# ~2M reads from the full 130M-record library. One chunk staged by stage_hg002.sh is ~20M
# pairs, where STRIDE=10 gives the same ~2M -- the script reports the count it actually got.
STRIDE="${STRIDE:-65}"
READ_LEN="${READ_LEN:-250}"
COVERAGE="${COVERAGE:-4}"
FRAGMENT_MEAN="${FRAGMENT_MEAN:-600}"
TAIL_CUT="${TAIL_CUT:-25}"
TAIL_WINDOW="${TAIL_WINDOW:-50}"
TAIL_CUT_DEEP="${TAIL_CUT_DEEP:-$((TAIL_CUT - 5))}"
# Floor on the subsample, so a truncated or mis-pathed input cannot produce a confident rate
# over a handful of reads. Lowered only for a smoke pass on a small fixture.
MIN_READS="${MIN_READS:-1000}"
WORK="${WORK:-$SCRATCH/hg002_pair}"
EIDOLON="${EIDOLON:-$SCRATCH/cargo-target/eidolon/release/eidolon}"
REFERENCE="${REFERENCE:-$REPO/eidolon/test_data/references/ecoli.fa}"

[[ -x "$EIDOLON" ]] || { echo "FATAL: no eidolon binary at $EIDOLON. Build it first (see header)." >&2; exit 1; }
for f in "$R1" "$R2" "$REFERENCE"; do
    [[ -f "$f" ]] || { echo "FATAL: not a file: $f" >&2; exit 1; }
done

mkdir -p "$WORK"
echo "=== fit_hg002_pair ==="
echo "repo:      $REPO ($(cd "$REPO" && git rev-parse --short HEAD), $(cd "$REPO" && git log -1 --format=%ci))"
echo "binary:    $EIDOLON ($(date -r "$EIDOLON" '+%Y-%m-%d %H:%M'))"
echo "R1:        $R1"
echo "R2:        $R2   stride $STRIDE"
echo "cut:       Q<$TAIL_CUT over the last $TAIL_WINDOW bases (fit and measurement both)"
echo "work:      $WORK"
echo

# ── 0. the binary is current, proven rather than assumed ────────────────────────
# Two silent failures this catches, both of which produce a full set of plausible numbers:
# a pre-#723 binary ignoring `fastq_file_r2`, and a pre-#735 binary emitting R2 backwards.
echo "[0/6] preflight: proving the binary understands fastq_file_r2 and emits R2 forwards ..."
PRE="$WORK/preflight"
mkdir -p "$PRE"
for m in 1 2; do
    src=$R1; [[ $m == 2 ]] && src=$R2
    # `zcat | head` makes zcat take SIGPIPE (141) and `set -e` would abort a step that
    # succeeded, so pipefail comes off for exactly this pipeline.
    set +o pipefail
    zcat -f -- "$src" | head -n 8000 | gzip -c > "$PRE/mate$m.fastq.gz"
    set -o pipefail
done
printf 'fastq_file: %s\nfastq_file_r2: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n' \
    "$PRE/mate1.fastq.gz" "$PRE/mate2.fastq.gz" "$PRE/model.json.gz" > "$PRE/fit.yml"
if ! "$EIDOLON" gen-seq-error-model -c "$PRE/fit.yml" > "$PRE/fit.log" 2>&1; then
    echo "FATAL: preflight fit failed. Last 20 lines:" >&2; tail -20 "$PRE/fit.log" >&2; exit 1
fi
python3 - "$PRE/model.json.gz" <<'PY' || exit 1
import gzip, json, sys
m = json.load(gzip.open(sys.argv[1], "rt"))
if "quality_score_model_r2" not in m:
    sys.exit("FATAL: the fitted model carries NO quality_score_model_r2. This binary predates\n"
             "       #723 and ignored `fastq_file_r2` instead of rejecting it -- both mates\n"
             "       would come from one fit and every per-mate number below would be a\n"
             "       duplicate. REBUILD (see the header) and resubmit.")
print("      ok: model carries a separate R2 population")
PY
printf 'reference: %s\noutput_dir: %s\noutput_filename: pre\noverwrite_output: true\nread_len: %s\ncoverage: 1\npaired_ended: true\nfragment_mean: %s.0\nfragment_st_dev: 60.0\nproduce_fastq: true\nproduce_bam: false\nproduce_vcf: false\nmutation_rate: 0.0\nrng_seed: preflight\nsequence_error_model: %s\n' \
    "$REFERENCE" "$PRE" "$READ_LEN" "$FRAGMENT_MEAN" "$PRE/model.json.gz" > "$PRE/gen.yml"
"$EIDOLON" gen-reads -c "$PRE/gen.yml" > "$PRE/gen.log" 2>&1
python3 - "$PRE/pre_r2.fastq.gz" <<'PY' || exit 1
import gzip, sys
head = tail = n = 0
for i, line in enumerate(gzip.open(sys.argv[1], "rt")):
    if i % 4 != 3: continue
    q = [ord(c) - 33 for c in line.strip()]
    if len(q) < 100: continue
    n += 1; head += sum(q[:25]) / 25; tail += sum(q[-25:]) / 25
if n < 100:
    sys.exit(f"FATAL: only {n} preflight R2 reads; cannot judge orientation.")
h, t = head / n, tail / n
if h <= t:
    sys.exit(f"FATAL: preflight R2 mean quality RISES along the read ({h:.2f} -> {t:.2f}).\n"
             f"       This binary predates #735 and emits R2's quality profile backwards.\n"
             f"       Every per-mate number below would be measured through that. REBUILD.")
print(f"      ok: preflight R2 degrades along the read ({h:.2f} -> {t:.2f}), n={n}")
PY
echo

# ── 1. strided subsample, both mates identically ────────────────────────────────
for m in 1 2; do
    src=$R1; [[ $m == 2 ]] && src=$R2
    SUB="$WORK/sub_r$m.fastq.gz"
    if [[ ! -s "$SUB" ]]; then
        echo "[1/6] subsampling R$m: keeping 1 record in $STRIDE ..."
        # pipefail STAYS ON: awk and gzip both read to EOF, so there is no early consumer and
        # no SIGPIPE. With it off, a truncated input makes zcat fail while awk and gzip
        # succeed, and the partial sample gets cached and fitted.
        rc=0
        zcat -f -- "$src" \
          | awk -v s="$STRIDE" 'NR%4==1 { rec++ } ((rec - 1) % s) == 0' \
          | gzip -c > "$SUB.part" || rc=$?
        [[ "$rc" -eq 0 ]] || { echo "FATAL: R$m subsample failed (rc $rc)" >&2; rm -f "$SUB.part"; exit 1; }
        mv -f "$SUB.part" "$SUB"
    fi
done
n1=$(zcat -f -- "$WORK/sub_r1.fastq.gz" | awk 'END { print NR/4 }')
n2=$(zcat -f -- "$WORK/sub_r2.fastq.gz" | awk 'END { print NR/4 }')
echo "      R1 $n1 reads, R2 $n2 reads"
[[ "${n1%.*}" -gt "$MIN_READS" && "${n2%.*}" -gt "$MIN_READS" ]] || { echo "FATAL: subsample holds $n1/$n2 reads, below the $MIN_READS floor" >&2; exit 1; }
# Rule 4: the denominator is part of the result. Mates that disagree in count mean the two
# files are not the same library in the same order, which invalidates every per-mate comparison.
[[ "${n1%.*}" -eq "${n2%.*}" ]] || { echo "FATAL: mate subsamples differ in size ($n1 vs $n2); the two FASTQs are not in the same order" >&2; exit 1; }

# ── 2. per-mate single fits: #695 step 1 ────────────────────────────────────────
# error_rate is ONE field shared across the pair in a paired model (#723 keeps the score set,
# offset and error rate shared), so the per-mate numbers #695 asks for come from separate fits.
for m in 1 2; do
    cfg="$WORK/fit_r$m.yml"
    printf 'fastq_file: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n' \
        "$WORK/sub_r$m.fastq.gz" "$WORK/model_r$m.json.gz" > "$cfg"
    echo "[2/6] fitting R$m alone ..."
    if ! "$EIDOLON" gen-seq-error-model -c "$cfg" > "$WORK/fit_r$m.log" 2>&1; then
        echo "FATAL: R$m fit failed. Last 20 lines:" >&2; tail -20 "$WORK/fit_r$m.log" >&2; exit 1
    fi
done

# ── 3. the four-population fit: #695 step 2, and #723's real-data follow-up ──────
cfg="$WORK/fit_pair.yml"
printf 'fastq_file: %s\nfastq_file_r2: %s\noutput_file: %s\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\nfit_quality_degradation: true\ndegradation_tail_cut: %s\ndegradation_tail_window: %s\n' \
    "$WORK/sub_r1.fastq.gz" "$WORK/sub_r2.fastq.gz" "$WORK/model_pair.json.gz" "$TAIL_CUT" "$TAIL_WINDOW" > "$cfg"
echo "[3/6] fitting both mates with the degraded population (four tensors) ..."
if ! "$EIDOLON" gen-seq-error-model -c "$cfg" > "$WORK/fit_pair.log" 2>&1; then
    echo "FATAL: paired fit failed. Last 30 lines:" >&2; tail -30 "$WORK/fit_pair.log" >&2
    echo "NOTE: a refusal naming a population that 'covers N bp' is the #720/#723 ragged guard," >&2
    echo "      not a crash -- the two mates disagree on read length." >&2
    exit 1
fi
sed 's/\x1b\[[0-9;]*m//g' "$WORK/fit_pair.log" \
  | grep -E "[Dd]egraded population|Separability|error_rate|read length|R2 " || true
echo

# ── 4. generate paired reads from the four-population model ─────────────────────
echo "[4/6] generating paired reads from the four-population model ..."
printf 'reference: %s\noutput_dir: %s\noutput_filename: gen_pair\noverwrite_output: true\nread_len: %s\ncoverage: %s\npaired_ended: true\nfragment_mean: %s.0\nfragment_st_dev: 60.0\nproduce_fastq: true\nproduce_bam: false\nproduce_vcf: false\nmutation_rate: 0.0\nrng_seed: hg002pair\nsequence_error_model: %s\n' \
    "$REFERENCE" "$WORK" "$READ_LEN" "$COVERAGE" "$FRAGMENT_MEAN" "$WORK/model_pair.json.gz" > "$WORK/gen_pair.yml"
"$EIDOLON" gen-reads -c "$WORK/gen_pair.yml" > "$WORK/gen_pair.log" 2>&1

# ── 5. model shapes and the fitted rates ────────────────────────────────────────
echo "[5/6] model shapes ..."
bash "$REPO/scripts/delta/describe_quality_model.sh" \
    "$WORK/model_r1.json.gz" "$WORK/model_r2.json.gz" "$WORK/model_pair.json.gz" || true
echo

# ── 6. the comparison ───────────────────────────────────────────────────────────
measure() {  # <label> <fastq>
    FASTQ="$2" MAX_READS=0 OUT="$WORK/meas_$1.txt" \
        TAIL_CUT="$TAIL_CUT" TAIL_CUT_DEEP="$TAIL_CUT_DEEP" TAIL_WINDOW="$TAIL_WINDOW" \
        bash "$REPO/scripts/delta/measure_quality_degradation.sh" >/dev/null
    awk -v l="$1" -v cut="$TAIL_CUT" -v deep="$TAIL_CUT_DEEP" '
        $1 == "tail" && $2 == "Q<" cut  { a = $5 }
        $1 == "tail" && $2 == "Q<" deep { b = $5 }
        END {
            if (a == "" || b == "") {
                printf "FATAL: no Q<%s / Q<%s row in the measurement of %s\n", cut, deep, l > "/dev/stderr"
                exit 1
            }
            printf "  %-24s Q<%-3s %-9s Q<%-3s %s\n", l, cut, a, deep, b
        }' "$WORK/meas_$1.txt"
}

echo "[6/6] === RESULT ==="
echo
echo "-- fitted error_rate, per mate (#695 step 1). Shipped default is 0.006638."
for m in 1 2; do
    printf '  R%s  %s\n' "$m" "$(grep -oE 'Computed error_rate: [0-9.]+' "$WORK/fit_r$m.log" | tail -1 | awk '{print $3}')"
done
echo
echo "-- degraded fraction, per mate (#695 step 2), from the four-population fit."
# Read from the MODEL, not the log. The two mates' log lines are worded differently -- R1's is
# "Degraded population:" and R2's is "R2 degraded population:" -- and a grep written against one
# silently drops the other, which is how a smoke pass reported R2 alone. The model file carries
# both fractions under one schema and cannot drift from the wording.
python3 - "$WORK/model_pair.json.gz" <<'PY2'
import gzip, json, sys
m = json.load(gzip.open(sys.argv[1], "rt"))
for label, key in (("R1", "quality_score_model"), ("R2", "quality_score_model_r2")):
    q = m.get(key)
    if q is None:
        print(f"  {label}  MISSING -- the model carries no {key}")
        continue
    d = q.get("degradation")
    print(f"  {label}  read_fraction {d['read_fraction']:.4f}" if d else
          f"  {label}  no degraded population in the model")
PY2
echo "      denominators and separability, from the fit log:"
sed 's/\x1b\[[0-9;]*m//g' "$WORK/fit_pair.log" \
  | grep -E "[Dd]egraded population: [0-9]+ of|Separability" | sed 's/^/        /' \
  || echo "        (none reported)"
echo
echo "-- collapsed-tail rate: real against simulated, per mate"
measure "real R1 (subsample)"  "$WORK/sub_r1.fastq.gz"
measure "sim  R1"              "$WORK/gen_pair_r1.fastq.gz"
measure "real R2 (subsample)"  "$WORK/sub_r2.fastq.gz"
measure "sim  R2"              "$WORK/gen_pair_r2.fastq.gz"
echo
echo "READ IT THIS WAY. The two real rows are the target and their RATIO is what #723 rests on"
echo "(HG002 measured R2 at 2.34x R1). The two sim rows should reproduce each mate's own rate,"
echo "and #694 predicts they will fall short -- the question this job answers is by how much,"
echo "per mate, on real data rather than a fixture."
echo
echo "Outputs in $WORK. Paste the RESULT block back, with the repo SHA from the header."
