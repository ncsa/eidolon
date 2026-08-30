# Model-builder resource baseline & output fidelity

The model builders (`gen-seq-error-model`, `gen-frag-length-model`, `gen-gc-bias-model`,
`gen-mut-model`, `gen-bam-models`) are what a first-contact user runs on *their own*
data. `eidolon/tests/model_parity.rs` pins the fit algorithms on tiny fixtures; this
document records how they behave on real, full-size input — the resource envelope
and the proof that a fitted model actually shapes gen-reads output.

Harness: `scripts/delta/model_builders.sbatch`, fed by `scripts/delta/stage_soy.sh`.

## Resource envelope (Delta) — soybean ~978 Mb and human GRCh38 ~3.1 Gb

soy uses two staging modes: chr scopes to the largest contig (fast correctness
check); `FULL_GENOME=1` stages the whole ~200-contig assembly. human is the full
GRCh38 with the GIAB HG002 v4.2.1 truth VCF (`stage_hg002.sh`). Numbers are
wall-clock and peak RSS from `/usr/bin/time -v`.

| builder      | soy chr        | soy full-genome    | human GRCh38        | notes |
|--------------|----------------|--------------------|---------------------|-------|
| seq_error    | 0:06.7 / 19 MB | 0:53.8 / 18 MB     | 1:46.4 / **35 MB**  | streams the FASTQ — memory-flat at any scale |
| frag_length  | 0:03.3 / 50 MB | 0:42.0 / 748 MB    | 1:00.0 / 738 MB     | BAM TLEN scan |
| gc_bias      | 0:04.5 / 323 MB| 1:03.5 / 3.6 GB    | 2:18.0 / **11.4 GB**| loads reference + per-GC windows across all contigs |
| mut_model    | 0:03.4 / 228 MB| 0:39.9 / 1.3 GB    | 1:55.0 / 1.9 GB     | reference + VCF (human = ~4M-variant GIAB truth set) |
| bam_models   | 0:04.5 / 352 MB| 1:04.3 / 4.2 GB    | 2:19.3 / **11.4 GB**| frag + GC in one BAM pass — heaviest |
| roundtrip    | 177,275 reads  | 183,583 reads      | 187,966 reads       | gen-reads consumes the built models, emits valid pairs |

Source runs: soy chr = job 19887967, soy full = 19888625, human = 19899988. All `overall: PASS`.

### Memory scales with reference size — plan accordingly

The heavy builders (`gc_bias`, `bam_models`) hold the reference plus per-window
state, so peak RSS tracks **genome size**, not read count. Measured:

- ~4 GB peak at ~1 Gb (soybean full genome).
- **11.4 GB peak at ~3.1 Gb (human GRCh38)** — the soy→human extrapolation
  (predicted ~11–13 GB) held; `gc_bias` landed exactly, `bam_models` mildly sublinear.

Trivial on an HPC node, but a user building `bam_models`/`gc_bias` against a human
reference on a 16 GB laptop is near the edge — budget ~12 GB. `seq_error` is the
exception — it streams and stays flat (18–35 MB) regardless of genome or input size.

Nothing approached the harness's 128 GB / 8 h request at either scale; no OOM, no
runtime cliff, clean on both a ~200-contig assembly and full GRCh38.

### A real bug the human run caught

The first human run (job 19899126) `FAIL`ed `mut_model`: `gen-mut-model` rejected the
GIAB truth VCF with `MalformedVcf("FORMAT list and sample list different lengths")`.
The VCF is spec-compliant — trailing per-sample FORMAT fields may be dropped (all but
GT) — so eidolon's reader was too strict. Fixed in #364 (`read_open_vcf` + `extract_gt_str`
pad dropped trailing fields; over-long samples still rejected); the re-run (19899988)
passes all five. This is exactly the first-contact-on-real-data failure the harness
exists to surface — the H1N1 fixture never exercised a dropped-FORMAT record.

## Output fidelity — do the fitted numbers reach the reads?

A model that parses is not the same as a model that shapes the output. Verified
links (all are real, runnable tests, not assertions on paper):

| model → output | evidence |
|----------------|----------|
| explicit `fragment_mean` → BAM insert size (TLEN) | `gen_reads … test_paired_ended_insert_size_matches_model` |
| **built `fragment_model` file → BAM insert size** | `tests/model_fragment_fidelity.rs` — fit on a 200 bp-mean BAM, simulate with the file, output insert mean = **199.7 bp** (fitted 200.0) |
| **built `fragment_model` SHAPE → BAM insert size** | `tests/model_fragment_fidelity.rs::a_built_discrete_model_carries_the_right_tail_all_the_way_to_the_reads` — right-skewed input, measured on ecoli: input skew **0.684** → model **0.684** → generated reads **0.681**. The mean-only assertion above passed for years while a Normal fit discarded the entire shape, which is why this one exists. |
| **built `sequence_error_model` file → read qualities** | `tests/model_output_fidelity.rs` — train on a uniform-Phred-35 FASTQ, simulate, output mean quality = **35.0** |
| **built `mutation_model` file → output variant count** | `tests/model_output_fidelity.rs` — a 139× rate difference between two fitted models yields 0 vs **405** output variants (405 ≈ rate 3.17e-2 × 13,133 bp) |
| **built `gc_bias_model` file → which regions get sequenced** | `tests/model_output_fidelity.rs` — train on a BAM covering a 20%-GC contig but not an 80%-GC one (fitted `w[20]=1.98`, `w[80]=0.02`), simulate: output is **506 low-GC reads / 0 high-GC** |
| input VCF variant → output reads | `gen_reads … test_input_vcf_snp_appears_in_bam_reads` |

The `model_builders.sbatch` round-trip is only a *usability* smoke test (reads exist
and are correctly-lengthed); the fidelity tests above prove the build→simulate path
is numerically faithful through the model-file load path — for **all four** builders.

Because H1N1's ~40% GC can't populate distinct GC bins, the gc_bias test builds a
synthetic 2-contig reference (a 20%-GC and an 80%-GC contig) so the weight table is
deterministic; the others reuse the bundled H1N1 fixture.
