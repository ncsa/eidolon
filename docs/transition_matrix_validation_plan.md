# SNP transition matrix — validation status and Delta plan

Status: **local tier complete; Delta tier specified and NOT yet run** (2026-08-12).

Covers `gen-seq-error-model`'s `bam_file:` / `transition_matrix_file:` inputs — the
BAM-derived SNP substitution matrix, which is one of the few model-builder capabilities
with no counterpart in NEAT 4 (its builders never open a BAM for this). Because it is
cited as a differentiator in the README comparison table, the evidence behind it should be
stated rather than assumed.

## Why this document exists

The feature shipped with tests that asserted only that a file appeared and deserialized.
Three of them constructed a known answer — an all-A→C BAM, a 10/5/5 count matrix — and
then asserted nothing about it. One said so in its own comment: *"We can't easily inspect
the matrix after serialization, so just verify that the runner completes."* Under that
coverage the matrix could have been garbage, or the default, or transposed, and every test
would have passed.

That is rule 1 in `CLAUDE.md` (assert content, not existence), and it had gone unnoticed
because the code is in fact correct. A correct implementation with vacuous tests reads
exactly like a verified one until something changes.

## What is now proven locally

Six unit tests in `eidolon/src/gen_seq_error_model/utils/runner.rs` and one output-fidelity
test in `eidolon/tests/model_output_fidelity.rs`:

| Claim | Test |
|---|---|
| MD-tagged mismatches produce the matching row | `test_runner_with_bam_md_tags_puts_all_a_weight_on_c` |
| Multi-target, asymmetric patterns are fitted per row | `test_runner_bam_transitions_track_a_mixed_mismatch_pattern` |
| The ref/read axes are not swapped | `test_build_transition_matrix_from_counts_is_not_transposed` |
| Observed counts drive the fitted row | `test_build_transition_matrix_from_counts_uniform_fallback` |
| A TSV overrides a BAM | `test_tsv_takes_precedence_over_bam` |
| No MD tags → hard error, not a silent default (#529) | `test_runner_bam_no_md_tags_is_an_error` |
| No `bam_file:` → default matrix, not an invented one | `test_runner_no_bam_file_uses_the_default_matrix` |
| **The matrix decides the substituted base in output reads** | `built_seq_error_transition_matrix_decides_the_substituted_base` |

Non-vacuity was established by mutation, not asserted:

| Mutation | Result |
|---|---|
| BAM branch returns `None` (counts ignored) | 2 tests fail |
| `transition_matrix_file` / `bam_file` precedence inverted | precedence test fails |
| `counts[i][j]` → `counts[j][i]` (axes transposed) | 4 tests fail |
| `generate_snp_error` reads a default instead of the model's matrix | fidelity test fails |

Two findings came out of that exercise:

- **The first version of the mixed-pattern test was itself vacuous.** Its fixture
  (`ref=ACGT` / `read=CATG`) produced A→C and C→A in equal number — a count matrix equal to
  its own transpose, so the transpose mutation passed it. The fixture is now asymmetric
  (`ref=AAAACCCC` / `read=CCGTAAAG`, A row 2C/1G/1T against C row 3A/1G). Single-entry rows
  are useless for this too, since normalization puts one nonzero cell at probability 1.0
  wherever it sits.
- **A forced matrix does not drive output to 100%.** `insertion_bias` is uniform over ACGT,
  so indel errors put inserted bases into the output that never consult the transition
  matrix. An absolute `>98% T` assertion fails on correct code. The fidelity test is
  therefore differential against a control run: default matrix gives a T share of
  **0.193**, forcing the A row to T gives **0.892**.

  This bullet previously read "`indel_probability` is 0.4, so ~40% of sequencing errors are
  indels". #660 corrected that constant to 0.01. The floor on this fixture is still ~39%,
  but for a different reason: #661 scales the indel share by local homopolymer run length,
  and the fixture is a 20,000-base poly-A reference — the curve's most enriched case
  (0.01 x 39.20 = 0.392). On a homopolymer-free reference the floor would be ~0.6%.

## What the local tier cannot establish

Every fixture above is synthetic and degenerate by design — one reference base, one forced
target, a hand-written count matrix. Three things remain unmeasured:

1. **Fit against a real substitution spectrum.** Real data exercises all four rows with a
   realistic skew (transition/transversion structure, strand and context effects). Nothing
   yet shows that a matrix fitted from a real BAM *resembles that BAM*.
2. **Round-trip through an aligner.** A user perceives this feature as "my simulated reads
   carry my instrument's error spectrum." That is only observable after alignment, at
   coverage, against a real reference.
3. **Scale.** `read_bam_transitions` streams a whole BAM through `walk_bam`. Runtime and
   peak RSS on a chromosome- or genome-scale BAM are unrecorded. The observer itself is a
   4×4 matrix, so RSS should be flat, but "should be" is not a measurement.

## Delta tier — proposed job

Fits the existing pattern: a `scripts/delta/run_transition_matrix_validation.sh` sibling to
`run_subclonal_vaf_validation.sh`, archiving through `archive_run` from `lib_report.sh` rather
than to a literal path. `RESULTS_DIR` already resolves to the one archive root
(`/projects/$ACCESS_PROJECT/$USER/eidolon-access-results`); naming a path here is how a doc
comes to reference `rneat-access-results/`, which was retired in the v2.0.0 rename and no
longer exists.

### Inputs

- GRCh38 chr22 (already staged by `fetch_validation_data.sh`).
- A real aligned BAM over chr22 **carrying MD tags**. The HG002 / GIAB alignment used by
  `stage_hg002.sh` is the natural choice. If MD is absent:
  `samtools calmd -b in.bam ref.fa > with_md.bam` — the builder requires it. Since #529 the
  builder **errors out** when a BAM yields no mismatches rather than falling back silently, so
  this failure now announces itself instead of arriving disguised as a result.

### Steps

1. **Independent ground truth.** Tally the training BAM's 12 off-diagonal substitution
   counts **without eidolon** — parse `MD`/`CIGAR` with `samtools view` + awk, or use
   `bcftools mpileup`. This must not come from `read_bam_transitions`, or the comparison is
   circular and proves only that the reader agrees with itself.
2. **Build two models.** `eidolon gen-seq-error-model` with `bam_file:` → *fitted*; the same
   config without it → *control* (default matrix). Same training FASTQ for both, so the
   quality model and `error_rate` are identical and the matrix is the only difference.
3. **Check the fit directly.** Compare the fitted model's 12 cells against step 1's tally.
4. **Simulate and align.** `gen-reads` over chr22 with `mutation_rate: 0.0` (so every
   mismatch in the output is a sequencing error, not a planted variant), 30×, paired-end;
   align with bwa-mem2. Run both models.
5. **Tally the simulated spectrum** from the aligned BAM, using the same step-1 tool.
   Insertion errors appear as CIGAR `I` rather than mismatches, so the uniform-insertion
   floor that constrains the local test does **not** contaminate this comparison — worth
   confirming in the output rather than assuming.
6. **Compare** fitted-vs-training and control-vs-training across the 12 cells.

### Pass criteria

- **Fit:** each of the 12 cells within 2% relative of the independent tally. A cell-wise
  check, not a summary statistic — a cosine similarity near 1.0 can hide one badly wrong
  cell.
- **Round trip:** max absolute deviation across the 12 proportions between simulated and
  training spectrum below 0.02, **and** the fitted model beating the control by at least
  5×. The control is what makes the number mean anything.
- **Denominators reported** on both sides: total mismatches counted, and per-cell counts.
  Per rule 4 a metric over an unknown denominator is not a result, and a zero or
  suspiciously small denominator is a **hard failure, not a warning** — the same shape as
  #450 and the `signature_check.sh` minimum-SNV gate.
- **Scale:** record wall-clock and peak RSS for the MD walk. No pass/fail; this becomes the
  first baseline entry, alongside `docs/model_builder_baseline.md`.

### Failure modes to guard explicitly

- MD tags absent → since #529 the builder errors out, so the job fails at the build step. Keep
  the identical-matrix check anyway: it is the backstop for any *other* route to a default
  matrix, and it is what would catch the guard being regressed. The job must never report a
  perfect control match as a pass.
- `mutation_rate: 0.0` not honored → planted variants inflate the mismatch tally and skew
  the spectrum toward the mutation model. Assert the output VCF is empty.
- The step-1 tool and eidolon disagreeing on how to count a multi-base MD run, or on
  soft-clipped/secondary records. Pin the filter set (`walk_bam` uses
  `BamWalkFilter::for_transitions()`, min MAPQ 0, secondary dropped) and mirror it in the
  awk pass, or the 2% tolerance is measuring a filter mismatch.

### Resource estimate

Alignment dominates: comparable to `germline_e2e.sbatch` on chr22 at 30×. The MD walk is a
streaming pass over the training BAM, so allow generously for I/O and expect flat RSS.
Both figures are guesses until the job runs, and should be replaced with measurements
rather than left as prose.

## Not verified

- No part of the Delta tier has run. Everything above the "Delta tier" heading is local,
  synthetic, and asserted; everything below it is a design.
- The BAM-derived matrix has never been compared against real data in any form.
- `read_bam_transitions`' handling of soft clips, secondary alignments and multi-base MD
  runs is exercised only by the fixtures in `bam_reader.rs`, which assert record counts —
  and, for the transition observer, an all-zero matrix. The non-zero real-BAM case is
  covered only transitively, through the runner tests.
- Long-read BAMs are untested here and are not in scope; see `docs/longread_epic_scope.md`.
