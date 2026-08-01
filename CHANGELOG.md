7/31/2026
=========
## eidolon v3.0.1 — measurement integrity

**Fixes only; no feature work and no change to simulated output.** Every item here is a
defect in something that *measures* eidolon, not in eidolon's simulation. They are
grouped into one patch release deliberately, so the Delta validation campaign that
follows attaches to a tagged version whose only difference from v3.0.0 is that the
measurements are now correct.

The common shape, and the reason these went unnoticed: each reported a plausible number
rather than failing. A harness that says `VERDICT: PASS` over a subset it silently
excluded, or `recall=0.000` for a comparison it could never make, is harder to catch
than a crash.

### `compare-vcfs` emitted VCFs with undeclared INFO tags (#444)

`FN_with_reasons.vcf` and `FP.vcf` pass each record's INFO column through verbatim but
wrote a fixed minimal header, so the artifacts carried tags their own header never
declared. bcftools tolerates that streaming VCF→VCF but hard-fails on BCF translation —
which `bcftools concat | sort` does internally — so it surfaced far downstream as an
opaque failure.

A fixed declaration set cannot work, because the two artifacts have different
provenance: FN records come from the **golden** VCF, FP records from the **called** VCF,
whose INFO is whatever the caller emitted (Mutect2 `TLOD`/`MBQ`, Strelka `SomaticEVS`).
Each artifact now inherits `##INFO`/`##FILTER`/`##FORMAT`/`##ALT`/`##contig` from its own
source, backfills an ID-only `##contig` for any record contig the source itself failed to
declare, and adds `GT` only when absent.

### Observed VAF was measured through a genotype call, destroying the low-VAF sites (#450)

The subclonal-VAF harness excluded the sites it exists to validate — 160 of 567 planted
sites, including an entire CCF cluster. Because the exclusion is VAF-dependent, the
reported bias was optimistic in the direction that flattered the result.

`bcftools mpileup | bcftools call -m -C alleles` was the cause. `-C alleles` constrains
which alleles are considered, but `call -m` still makes a diploid ML **genotype** call,
and a hom-ref call discards the uncalled allele together with its read count. On a
synthetic 200× BAM carrying the three CCF clusters, `call -m` destroyed **two of three**:

| site | intended | `mpileup \| call -m` | `mpileup` alone |
|---|---|---|---|
| 331 | 0.350 | `ALT=T` `AD=130,70` | `AD=130,70,0` → 0.350 |
| 631 | 0.175 | `ALT=.` `AD=165` | `AD=165,35,0` → 0.175 |
| 931 | 0.070 | `ALT=.` `AD=186` | `AD=186,14,0` → 0.070 |

Coverage cannot rescue it — identical at 100×/300×/600× — because the decision is driven
by the allele *fraction* against a diploid model, and no diploid genotype predicts 7%.

Both harnesses now read `FORMAT/AD` straight from mpileup. This also fixes
`run_scn_af_validation.sh` (#398's pool-AF validation), which shared the defect and was
worse: with no `-C alleles` at all, sites the caller declined to call vanished entirely —
the low end of the very spectrum it validates. `scn_af_compare.py` gained multi-allelic
expansion (so the truth's ALT matches against mpileup's `T,<*>` list and selects that
allele's own AD element), arity-checked per-allele fields (it had been taking `[0]`,
reporting the *first* allele's fraction for every allele of a multi-allelic record), and
zero-fill: a truth ALT the observed side never listed at a covered position now scores
**0.0 observed** rather than being dropped, since otherwise the VAF-dependent exclusion
just moves to a lower threshold.

### SV scoring reported numbers that were not measurements (#457)

Three defects:

- **Per-type runs filtered the truth but not the query**, so every record of another type
  counted as a false positive — producing arithmetically impossible output, with
  `manta_BND FP=39` against `manta_overall FP=20`. A subset cannot have more false
  positives than the whole. Per-type precision/FP/f1 from this harness were meaningless;
  only per-type recall was ever sound.
- **BND was unscoreable.** `truvari bench` matches by position *and size overlap*, and a
  breakend has no `SVLEN` with `END == POS`, so truvari cannot match one however correct
  the truth is. Adds `scripts/delta/bnd_proximity.py`, which scores junctions (not
  records — a junction is one event described by two mate records, so record-counting
  weighted BND 2× against DEL/DUP/INV) and is representation-agnostic, since callers
  legitimately re-represent a junction as `<DUP:TANDEM>` at the correct coordinates.
- **The aggregate was dominated by unmatchable types.** BND + CNV were 23 of 47 truth
  records in one run, all guaranteed misses, capping the aggregate near 0.5 regardless of
  caller quality — and making v3.0.0's *correct* BND fix look like a regression
  (0.667 → 0.413). The aggregate now covers only the truvari-scoreable subset and states
  what it excluded and where each type is scored instead.

### Harnesses now fail instead of reporting an unusable measurement

`signature_check.sh` warned `few SNVs — fit will be noisy` below 50 and fitted anyway.
SBS-96 has 96 trinucleotide contexts; a spectrum from fewer mutations than contexts
cannot populate them, and SigProfiler still returns an assignment that was reported as
though it measured something. Now fatal below `MIN_SNVS` (default 50 — the historical
bar, enforced rather than moved; §3.9's own analysis used 6,225 SNVs). This was the last
`WARNING`-on-an-unusable-measurement site from the harness audit.

### Also

- `mpileup -d` is now explicit (`MPILEUP_MAX_DEPTH`, default 2000). Its default of 250
  would silently downsample a 200×+ run.
- New read-the-artifact guards: abort when no site carries `FORMAT/AD`, and when a
  per-type query subset is empty report `recall=0.000` rather than emitting a precision
  figure computed from nothing.
- `eidolon/tests/vcf_validation.rs` gained whole-artifact validation of the
  `compare-vcfs` outputs.

Every fix in this release was reproduced and verified locally; no Delta time was spent
diagnosing or confirming any of them.

7/28/2026
=========
## eidolon v3.0.0 — subclonal heterogeneity (#405); output tokens renamed NEAT_* / RNEAT_* → EIDOLON_*; adopting SemVer

**Breaking output-format change** — no behavior change, only the names of emitted tokens.
Completes decision 1 of `docs/rename_eidolon_scope.md` (the post-2.0.0 output-token
migration). Downstream code that parses the old names must update.

### Why 3.0.0, and a versioning policy going forward

This release is **3.0.0**, not 2.1.0, and from here on eidolon follows
[Semantic Versioning 2.0.0](https://semver.org).

Earlier releases were versioned by feel, and by that habit this would have been a minor
bump: it completes a rename announced back in v2.0.0, and three prior releases shipped
breaking output changes as minor bumps (v1.11.0 and v1.12.0 both changed the FASTQ
read-name format; v1.6.0 removed config fields under a bolded "Breaking change"). That
pattern was not doing users any favours — a minor bump tells you it is safe to upgrade
without reading anything, and for these releases that was not true.

SemVer requires a project to **declare what its public API is**, so eidolon's is now
stated explicitly (also in the README):

**Part of the public API** — changes here require a MAJOR bump:
- Names and semantics of emitted VCF INFO tags, the VCF sample-column name, and the
  `FILTER`/`FORMAT` fields eidolon writes
- FASTQ/BAM **read-name (QNAME) format**, including the `EIDOLON_generated_` /
  `EIDOLON_chimeric_` prefixes and the positional fields encoded after them
- CLI subcommand names, flags, and configuration-YAML keys
- Model-file format compatibility (a model built by version *N* being readable by *N+1*)

**Not part of the public API** — these may change in a MINOR or PATCH release:
- The Rust library (`eidolon-core`) API surface; this crate exists to serve the binary
- Log output, progress reporting, and human-readable messages
- Exact simulated *content* for a given seed — reads are a random draw, and any change
  to the sampler or models legitimately changes output while preserving the format

For a read simulator the emitted format is the interface downstream code actually binds
to, which is why it heads the list. That is also the reasoning that makes this release
major: renaming the tokens is a breaking change to that interface, notwithstanding that
it was pre-announced and that no behavior changed.

### New feature — subclonal heterogeneity for gen-cancer-reads (#405)

The tumor is no longer a single global purity split, which collapsed every somatic variant
to ~one effective VAF. A tumor is now a mixture of **subclones** at distinct cancer-cell
fractions (CCF), and each somatic variant's observed alt fraction composes as
`purity × dosage × CCF` — the CCF *composes* with the variant's dosage and purity rather
than replacing them, so a het somatic at CCF *f* lands at ~*f*/2 in the merged output.
This is the intra-tumor heterogeneity axis the predecessor could not represent at all.

Three ways to specify the architecture (mutually exclusive):

- **Inline** — `subclones:`, a list of `{ccf, weight}`; the somatic burden is split across
  populations in proportion to `weight` (default 1.0 = equal share), `ccf ∈ (0, 1]`.
- **Fitted from real data** — `subclones_file:` accepts a **PyClone-VI** per-mutation table
  (grouped by cluster, weight = mutation count) or a **PCAWG-11 / DPClust** cluster table
  (`cluster` + `ccf` + a size column: `n_ssms` / `size` / `weight`).
- **Replayed** — `somatic_vcf:` honors a real tumor's *observed* VAF per site: each input
  site's VAF is divided by purity so the merged reads reproduce it, and the record is
  tagged `somatic`. Observed VAF above purity is clamped to 1.0; sites with no usable VAF
  fall back to genotype dosage.

**Ground truth in the truth VCF.** Every somatic record carries `INFO/EIDOLON_CCF` (the
intended CCF) and `INFO/EIDOLON_VAF` (the intended *observed* VAF — what a caller should
measure on the merged sample). Note `FORMAT/AF` in the truth is measured per-pass
(tumor-only, = dosage × CCF); `EIDOLON_VAF` is that × purity, so `EIDOLON_VAF` is the
correct target when scoring a somatic caller against the merged output.

**Real-data validation (Delta).** A round-trip the unit tests cannot reach: simulate, align
the merged reads with bwa-mem2, and measure the observed VAF at each somatic site against
the intended `EIDOLON_VAF`. Two complementary runs:

| Metric | Generative — soy scaffold | Reproductive — HCC1395 chr1 |
|---|---|---|
| Input VAF spectrum | 3 synthetic CCF clusters (≈0.04 / 0.16 / 0.40) | real tumor's empirical VAF, all deciles 0.1–1.0 |
| Coverage / purity | 200× / 0.8 | 200× / 0.9 |
| Somatic sites (compared) | 567 planted (387 scored) | 2,948 (2,698 scored) |
| Bias, mean(observed − intended) | −0.003 | −0.004 |
| MAE (vs ~200× noise floor) | 0.026 | 0.024 |
| Pearson r | 0.95 | **0.99** |

The generative run confirms the `purity × dosage × CCF` composition is unbiased and
accurate to the sampling limit; per-site Pearson is a spread-dependent summary for tight
discrete clusters, so fidelity is gated on bias + MAE-vs-noise-floor rather than r. The
reproductive run replays SEQC2 **HCC1395** (a triple-negative breast cancer cell line) at
its observed VAF over GRCh38 chr1 — a real tumor's full *continuous* distribution, where
Pearson is fully meaningful: r = 0.99, unbiased (−0.004), per-decile MAE uniform at
0.009–0.027 across the whole 0.1–1.0 range.

Harnesses: `scripts/delta/run_subclonal_vaf_validation.sh` (generative) and
`scripts/delta/run_hcc1395_reproductive.sh` (reproductive). Config template:
`template_config/gen_cancer_reads_template.yml`; details in `docs/cancer_simulator.md` and
`docs/cancer_howto.md`.

### Output tokens renamed NEAT_* / RNEAT_* → EIDOLON_*

- **VCF INFO tags:** `NEAT_ORIGIN` → `EIDOLON_ORIGIN`, `NEAT_PROVENANCE` → `EIDOLON_PROVENANCE`,
  `NEAT_REASON` → `EIDOLON_REASON`, `NEAT_CCF` → `EIDOLON_CCF`, `NEAT_VAF` → `EIDOLON_VAF`.
- **VCF sample column:** `NEAT_simulated_sample` → `EIDOLON_simulated_sample`.
- **FASTQ/BAM read-name prefixes:** `RNEAT_generated_*` → `EIDOLON_generated_*`,
  `RNEAT_chimeric_*` → `EIDOLON_chimeric_*` (these still carried the older `rusty-neat`
  prefix; now aligned with the tool name).
- Bundled scripts moved in lockstep (`tools/cancer_simulate.sh`, `cancer_benchmark.sh`,
  `cancer_sv_benchmark.sh`, `scripts/delta/*`); the `filter_reads` read-name parser updated
  (now length-agnostic); parity/baseline test assertions re-blessed.
- **Migration path for pre-v2.1.0 files.** The rename would otherwise strand existing
  artifacts, so:
  - **`tools/migrate_legacy_tokens.sh`** converts a v2.0.0 VCF in place of a manual fix-up
    (`--rename-annots` for the INFO tags, `reheader -s` for the sample column). It is
    **idempotent** — a no-op on an already-current or non-eidolon VCF — so call sites can
    invoke it unconditionally. A `--check` mode reports without converting (exit 10 = legacy
    tokens present). Note a dual-name bcftools filter is *not* a workable alternative:
    bcftools rejects an undefined tag in a `-i` expression at parse time, so even
    `EIDOLON_ORIGIN=… || NEAT_ORIGIN=…` fails outright on a file carrying only one of them.
  - `scripts/delta/{cancer_pipeline,sv_pipeline}.sbatch` run that converter automatically when
    the truth VCF predates the rename. `tools/cancer_benchmark.sh` and `cancer_sv_benchmark.sh`
    run bcftools inside a container (the host needn't have it), so they instead pre-flight the
    truth header and print an actionable message naming the converter. That check also makes
    the pre-existing "pass `--truth-filter ''`" hint reachable — bcftools' parse-time tag
    error previously aborted the run before it could be printed.
  - **`eidolon filter-reads` accepts both `@EIDOLON_generated_` and the v2.0.0
    `@RNEAT_generated_` read-name prefix**, so legacy FASTQs need no conversion at all —
    rewriting hundreds of GB to change a 19-byte prefix is not a reasonable ask. (The two
    prefixes differ in length, so the parser slices by the *matched* prefix.)
  - `write_vcf` now **drops eidolon's own pre-v2.1.0 INFO tokens** from a passed-through INFO
    column instead of reproducing them. Feeding a v2.0.0 golden VCF back in via `input_vcf:`
    or `somatic_vcf:` previously preserved e.g. `NEAT_PROVENANCE` verbatim on symbolic records
    while the header declared only the `EIDOLON_*` names — an undeclared tag, which streams
    fine VCF→VCF but hard-fails BCF translation (`bcftools concat | sort` does that
    internally). Dropping rather than renaming is deliberate: the writer re-emits the current
    equivalent, so a stale tag would contradict it.
- **Not changed:** NEAT pedigree in README/`CITATION.cff`/comparison text; the `NEAT_DATA`
  staging env var; and the benchmark harnesses' env vars that point at the *predecessor*
  NEAT 4 tool (`NEAT_BIN`, `NEAT_ENV`, `NEAT_SEED`, `NEAT_MAX_GENOME_MB`) — those name NEAT,
  not eidolon.
- **Also not changed (deliberate):** the sample-column names in *intermediate corpus* VCFs
  synthesized by the model-building scripts — `NEAT_tumor_corpus` (`tools/fetch_tumor_corpus.sh`,
  `fetch_cosmic_corpus.sh`, `fetch_cosmic_per_tissue_corpus.sh`) and `NEAT_pcawg_corpus`
  (`tools/build_pcawg_sv_vcf.py`). These are emitted tokens, not env vars, but they appear only
  in throwaway inputs to `gen-mut-model` — never in simulation output — so renaming them would
  force a re-fetch of every corpus for no downstream benefit.

7/21/2026
=========
## eidolon v2.0.0 — renamed from rusty-neat / rneat

The project is renamed **rusty-neat → eidolon** (binary `rneat` → `eidolon`). This is a
breaking change to the command name and crate identity; there are **no behavior changes**.

- **New name.** The tool, binary, and bioconda package are now `eidolon`; the library crate
  is `eidolon-core`. Invoke `eidolon gen-reads …`, `eidolon gen-cancer-reads …`, etc.
- **`rneat` deprecation shim.** A `rneat` alias binary still ships: it prints a rename warning
  and forwards all arguments to `eidolon`, so existing scripts keep working. It will be removed
  in a future release — please migrate to `eidolon`.
- **Output formats unchanged.** Golden-VCF INFO tags (`NEAT_PROVENANCE`, `NEAT_ORIGIN`), the
  `NEAT_simulated_sample` sample column, and FASTQ read-name prefixes (`RNEAT_generated_*`,
  `RNEAT_chimeric_*`) are **retained** for output stability and to preserve the NEAT lineage in
  the format. A future release may migrate these to `EIDOLON_*`.
- **NEAT pedigree retained** in the README and `CITATION.cff` — eidolon remains a Rust port of
  NEAT.
- Repository moves to `github.com/ncsa/eidolon`; GitHub redirects the old `rusty-neat` URLs.

7/20/2026
=========
## rneat v1.21.1 — realistic default het/hom ratio

The bundled default mutation model shipped `homozygous_frequency = 0.01` — a faithful
port of the NEAT lineage's own default (NEAT2 0.010, NEAT4 0.001) — which forces each
variant site to be homozygous-alt with probability 0.01 and so drives a het/hom ratio
of ~99 in default-model gen-reads output. That is unrealistic for any diploid, where a
resequenced sample is ~1/3 homozygous-alt. Corrected to **0.333** (het/hom ≈ 2.0),
confirmed by re-running the germline pipeline on chr22 with the corrected default
(42,997 sites → het/hom 2.00, P(hom) 0.333) and locked by a regression test.

- The sampling logic was already correct: `homozygous_frequency` is fit as
  (homozygous sites)/(total sites) and consumed directly as P(homozygous), with no
  inversion — a properly *fit* model already produced a realistic ratio. Only the
  embedded default parameter was wrong.
- Runs that supply their own fit mutation model are unaffected. Default-model gen-reads
  output changes (genotype calls, and which reads carry the alt allele); model-parity
  and the full test suite stay green.

7/18/2026
=========
## rneat v1.21.0 — per-variant allele fraction

Adds an optional per-variant allele fraction so an input VCF can drive a continuous
allele-frequency spectrum in the simulated reads, instead of only the `{0.5, 1.0}` the
genotype-based fraction could produce. Backward-compatible: variants without an allele
fraction are unchanged (model-parity + full suite green).

### Read generation — allele-frequency realism (#398)
- gen-reads now honors an optional per-variant `allele_fraction` parsed from the input VCF —
  `INFO/AF` if present, else `FORMAT/AD` as `alt / (ref + alt)`. When set, the alt allele is
  emitted on that fraction of overlapping reads, and the golden VCF's measured `GT:AD:DP:AF`
  tracks it. This lets rneat replay pooled or somatic allele-frequency spectra from a real
  sample — the reproductive-replay case the model builders could not express before. Allele
  depth is used over genotype for pooled data, where `GT` is not diploid-meaningful.
- **Backward-compatible.** Variants with no allele fraction keep the genotype-based behavior
  (homozygous → all alt, heterozygous → ~0.5) with the identical RNG draw, so default runs are
  byte-identical to prior output.
- **Validated on real data.** On soybean-cyst-nematode Pool-seq, simulated-vs-real per-site
  allele frequency reached Pearson r = 0.98 at 150x coverage, unbiased, with the residual
  error tracking the binomial coverage-noise floor.

### Notes
- `Variant` no longer derives `Eq` / `Hash` — the new `Option<f64>` field is only `PartialEq`,
  and `Variant` was never used as a map key or set member (a separate `VariantKey` is hashed).


7/13/2026
=========
## rneat v1.20.1 — bug-fix

Patch release on top of v1.20.0. Fixes a core RNG bug where `NeatRng::random()` could
return values outside `[0,1)` for certain seeds. Output-preserving for all well-behaved
seeds (model-parity + full suite unchanged), so v1.20.0 outputs are unaffected except where
the RNG was previously producing out-of-range draws.

### Bug fixes
- **RNG: `NeatRng::random()` could return values outside `[0,1)`.** The `mash()` seed hasher
  dropped the u32 truncation the JS/Alea algorithm it ports applies to its return
  (`(n >>> 0) * 2^-32`) — `n` was an unmasked `u64`, so `mash()` could return `>= 1`. That let
  `new_from_seed` state escape `[0,1)` for certain seeds — notably RNGs from `derive_child(idx)`
  (one per contig/chunk in whole-genome gen-reads) — and `random()` then produced out-of-range
  draws (observed: a child whose first `random()` was `-8311.87`). Downstream this crashed
  `Normal::inverse_cdf` in fragment-length sampling on the SV-weighted path (`x must be in
  [0, 1]`) at whole-genome scale, and could have mis-fed other RNG consumers. Fixed at the
  source by masking `mash()`'s return to 32 bits, which is a no-op whenever `n < 2^32` (every
  seed that already behaved) — so existing seeds/outputs are unchanged (model-parity + full
  suite green); only previously-broken seeds change. Regression test added.


7/7/2026
========
## rneat v1.20.0 — realism update

Adds trinucleotide-context-aware SNP placement so simulated data reproduces
context-specific mutational signatures (e.g. APOBEC), not just the overall mutation
rate. On by default; indel placement and the model builders are unaffected. See the
output-compatibility note below.

### Read generation — mutational-signature realism (#372)
- gen-reads now weights *where SNPs land* by the local trinucleotide's fitted propensity
  `w(ctx)`, so simulated data reproduces context-specific mutational signatures (e.g.
  APOBEC SBS2/13), not just the overall mutation rate. SNP placement was previously
  context-independent, which flattened signatures. **On by default.** Validated on SEQC2
  HCC1395: the SBS-96 cosine of real-vs-simulated somatic SNVs rose from **0.72 to 0.99**.
- **SNPs only.** `w(ctx)` is a SNP substitution propensity, so it drives SNP placement
  alone; indels stay context-neutral (regional-rate placement, unchanged). Their positional
  biology — homopolymers, short-tandem-repeats, microhomology — is not the SNP trinucleotide
  signature, and applying the SNP weight to them concentrates indels and degrades indel
  recall. A separate indel-context placement model is tracked under the realism epic (#311).
- **Memory-bounded.** Context-weighted placement processes each contig in chunks
  (≤ 4M positions) rather than materializing a contig-length per-position weight array, so
  peak memory stays flat on large chromosomes.
- **Output-compatibility note:** for any run whose mutation model carries context bias —
  which includes the bundled default model and all bundled COSMIC models — simulated *SNP
  positions* differ from v1.19.1 and earlier (they are now context-realistic). Indel
  positions, the mutation *rate* / count, and `gen-mut-model` / the model builders are all
  unaffected. A workflow that needs byte-identical output to ≤1.19.1 should pin that
  version. (Only a genuinely context-flat model takes the unchanged fast path.)
- **Validation:** Delta Tier-1 regression gate passes — RSS bounded, all germline
  recall/precision/Ts-Tv hold, determinism/order-independence unchanged. A develop-vs-#372
  A/B on the cancer suite config (COSMIC model, chr22, 30x, purity 0.6; 3 seed reps/arm)
  showed **#372 has no measurable effect on somatic caller recall** (SNV 0.924±0.020 vs
  develop 0.927±0.013; indel within seed noise). The two somatic-recall regression baselines
  (`scripts/delta/baseline_metrics.tsv`) were re-frozen from that replication — their prior
  sd was seeded from a different config and underestimated true run-to-run variance.


7/6/2026
========
## rneat v1.19.1

Bug-fix + logging release on top of v1.19.0. The output-affecting changes are the VCF
dropped-FORMAT fix and the `fragment_model` input paths; the logging change is
output-preserving (verbosity/perf only), so fidelity is unchanged from v1.19.0.

### Read generation
- `gen-reads` accepts a `fragment_model` as a paired-end fragment source — no longer
  requires explicit `fragment_mean`/`fragment_st_dev` when a model is supplied; the model
  takes precedence at runtime, mean/st_dev are ignored if both are set (#355).
- `gen-cancer-reads` accepts a `fragment_model` for paired-end runs, at parity with
  gen-reads (#366).

### Model building / VCF input
- `gen-mut-model` (and every VCF input path) now accepts spec-compliant records that drop
  trailing per-sample FORMAT fields (all but GT). Previously these errored with "FORMAT
  list and sample list different lengths", silently failing the model build on standard
  GIAB / bcftools-mpileup VCFs. Samples with *more* fields than FORMAT are still rejected (#364).

### Logging / performance
- **Default `.neat.log` level is now `info`** (was effectively `trace`). The `--log-level`
  arg carried `default_value("trace")`, which silently overrode the intended `info`
  default, so every run wrote a trace-level log. On a human-sized reference the per-base
  trace/debug events produced a multi-GB `.neat.log` and burned most of the wall-clock on
  log I/O. Pass `--log-level debug`/`trace` to opt into verbosity; a bare `--log-level`
  now means `debug`. The on-screen log is unchanged (always `info`). This is the likely
  cause of the ~2× wall-clock seen in the v1.19.0 regression benchmark (#340).
- Removed the per-base sequencing-error debug logs (`"Creating sequencing error"`,
  `"Snp error"`, `"Deletion error"`, `"Insertion error"`, `"Generating basic SNP error"`) —
  constant strings emitted ~once per erroneous base, no diagnostic value.
- Delta harnesses (`scripts/delta/*`) now pass `--log-level warn` to `gen-reads` so a
  benchmark/validation run never writes a large log to the shared filesystem.

### Validation / tooling
- Adapter readthrough is now part of the regression baseline (adapter_* rows in
  `baseline_metrics.tsv`; `collect_adapter_validation.sh CANDIDATE_TSV=…` feeds
  `regression_gate.sh`). See #338.


7/2/2026
========
## rneat v1.19.0

Optional Illumina 3′ adapter readthrough (#125) — realistic short-insert data that
exercises adapter-trim QC — plus the correctness fix and Delta validation surfaced
while building it. Default off; output byte-identical to v1.18.0 when the `adapters`
block is omitted.

### Read generation
- Optional 3′ adapter readthrough (#125): `adapters: {enabled, preset: truseq|nextera|custom, r1, r2}`.
  **Default off** — output byte-identical when the block is omitted. When enabled, a fragment
  whose insert is shorter than `read_len` is kept and the read is padded to `read_len` at its
  3′ end with adapter sequence (R1→`r1`, R2→`r2`, appended after R2 is reverse-complemented),
  carrying model quality/error, soft-clipped (`S`) in the golden BAM, introducing no variants.
  Incompatible with `long_reads`. Amount of readthrough scales with the insert distribution
  relative to `read_len`.
- `keep_short_fragments`: keep short inserts as insert-length **genomic** reads *without*
  adapter padding — the adapter-free short-insert mode.

### Fixed
- Adapter path emitted a malformed FASTQ record — quality string one character longer than
  the sequence — for degenerate zero-length inserts (`start == end`). Invisible to `zcat`/`wc`
  but bwa-mem2's parser halts at the first such record, silently truncating alignment to a few
  thousand reads. Zero-length inserts are now skipped (both mates, streams stay in sync); no
  effect when adapters are off.

### Validation / tooling
- 4-arm adapter validation harness (off / short_ctrl / on_raw / on_trim) with a contrast-based,
  data-driven collector that isolates the adapter effect from the short-insert coverage effect
  (`scripts/delta/`). Result: adapter readthrough is detectable and **fully callable** — GATK
  SNP/indel fidelity is unchanged (within replication noise) whether adapters are fastp-trimmed
  or left for BWA-MEM2 to soft-clip; the only fidelity difference vs a long-insert baseline is
  the reduced effective coverage of short inserts, present with adapters entirely absent.
- `germline_e2e.sbatch`: node-local staging (`LOCAL_STAGE`, immune to shared-Lustre OST drops)
  and a fastp thread cap (avoids a livelock at high `--cpus-per-task`).

### Docs
- README: new "3′ Adapter Readthrough" section (config, presets, use cases, validated behavior).


6/29/2026
=========
## rneat v1.18.0

Releases work that accumulated on `develop` after v1.17.4 (June 17) but was never
cut into a release — performance, memory, and read-generation refinements — plus
the ACCESS validation suite. Supersedes v1.17.4; includes everything in it.

### Performance & memory
- Eliminated the per-base heap allocation in the read-generation hot loop and made
  the fetched fragment zero-copy (`get_subseq_slice`) — faster single-thread read
  generation, the mode that dominates rneat's throughput.
- `target_bed`-aware reference loading: only the targeted contigs are materialized,
  so a region-sharded whole-genome run no longer holds the full reference in RAM
  (per-shard footprint dropped from ~GBs to ~MBs). Output is byte-identical —
  per-contig RNG keys off the original file position regardless of which contigs
  load.

### Read generation
- Reverse (R2) reads are generated forward over the right-end window and then
  reverse-complemented as a whole record, with variant offsets indexed by
  `seq_index` for both strands — correct SNP/insertion/deletion handling on R2.
- PE fragments are padded by one read length so a deletion near R2's tail consumes
  reference instead of truncating and dropping the pair.

### CLI
- `rneat --version` / `-V`.

### Validation
- Extensive Delta (NCSA ACCESS) validation campaign documented under `docs/`
  (report + figures), with a reusable regression suite (`scripts/delta/`,
  `docs/regression_protocol.md`) that gates future changes on memory, speed,
  fidelity, and determinism against replication-derived tolerances.

> Note: the develop⇄main history had diverged before this release; the
> read-generation bullet describes the code delta vs v1.17.4 — review the
> attribution against prior CHANGELOG entries when finalizing.



### Packaging: bundle third-party Rust license texts (Bioconda)

Per conda-forge/Bioconda policy, the licenses of all (transitive) Rust
dependencies are now bundled with the package. The Conda build runs
`cargo-bundle-licenses` to produce `THIRDPARTY.yml`, which is shipped via
`about.license_file` alongside `LICENSE.md`. A curated `THIRDPARTY.yml` is
committed at the repo root and reused by the build with `--previous`, so the
nine crates that publish no LICENSE file in their crate (the `noodles` family,
`r-efi`, `vectorize`, `wasip2`) carry hand-filled canonical license texts;
`--check-previous` fails the build if a future dependency change leaves a crate
without a recorded license. No source or runtime behavior changes. Addresses
the Bioconda recipe review (PR #66122).

6/16/2026
=========
## rneat v1.17.3

### Fix: Windows release build (VS2026 runner)

The `windows-latest` GitHub runner was upgraded to the `windows-2025-vs2026`
image (Visual Studio 2026), and the v1.17.2 Windows build failed with
`couldn't determine visual studio generator` while compiling `libz-ng-sys` (the
`zlib-ng` backend for `flate2`). The transitive `cmake` build dependency was
pinned at 0.1.54, which predates VS2026 support. Bumped `cmake` to 0.1.58
(VS2026 support landed in 0.1.55) via `Cargo.lock`; no source or other platform
changes. Re-released as v1.17.3 to produce the missing Windows binary.

6/15/2026
=========
## rneat v1.17.2

### Fix: reference N bases are now treated as gaps, not filled with random sequence

`gen-reads` previously overwrote **every** reference `N` with a random base at
FASTA load time (`apply_n_substitution`). Because the region machinery treats
masked bases as ordinary sequence (only `N`/`X` mark gaps), the loaded reference
contained no `N` — so each contig became a single non-N region and the gap
exclusions were silently bypassed:

- read fragments were anchored across centromere/telomere N-tracts,
- de-novo mutations were sampled inside N regions, and
- the v1.13.1 SV alignability gate (#224), which counts only `N`, always passed
  — making that fix inert in the real pipeline.

N bases are now left as-is at load. The existing non-N machinery (`map_buffer` /
`get_non_n_regions`) governs read anchoring, mutation placement, and SV
anchoring, so assembly gaps remain coverage dropouts. Output FASTQ is still
ACGTN-only (bases are unmasked at write time); reads spanning short interior N
runs may emit `N`. The model-training subcommands (`gen-mut-model`,
`gen-gc-bias`) already read raw sequence and were unaffected. Added an
end-to-end regression test asserting no variants are placed inside an N tract.

6/9/2026
========
## rneat v1.17.1

### Faster compression across all subcommands

v1.17.0 made `gen-reads` FASTQ output fast; this extends the same libdeflate
speedup to the remaining output paths:

- **BAM and `.vcf.gz`** now use libdeflate via the `noodles-bgzf` `libdeflate`
  feature (the golden VCF and BAM writers previously used slower zlib).
- **`gen-cancer-reads`** (the tumor/normal FASTQ merge) and **`filter-reads`**
  FASTQ output now use the libdeflate `BlockGzWriter`.
- **BAM-only runs** no longer compress the generated FASTQ into a discarded
  buffer (it now goes to a null sink), removing wasted work.

## rneat v1.17.0

### Faster FASTQ compression (libdeflate + single-pass combine)

Gzip was ~25–27% of single-threaded wall time. Two changes cut that:

- **libdeflate block-gzip** — a new `BlockGzWriter` compresses FASTQ output with
  libdeflate (~2–3× faster than zlib for the same gzip) in independent ~256 KiB
  gzip members. The result is a standard multi-member gzip file (read
  transparently by `gunzip`, `zcat`, and downstream tools); memory stays bounded
  to one block.
- **Single-pass combine** — per-contig temp files are now byte-concatenated into
  the final FASTQ instead of being decompressed and re-compressed, so each read
  is compressed once (in the parallel per-contig write) instead of twice.

Net effect: ~22% faster single-threaded read generation (e.g. the c_elegans
benchmark dropped from ~365 s to ~283 s), with no multi-threaded regression and
unchanged memory.

### Deterministic, byte-identical FASTQ output (#251)

For a given seed, FASTQ output is now byte-identical regardless of `num_threads`
(previously the read *set* was stable but record order varied). Read names are
unique per fragment and the final files are assembled in contig order.

### Opt-in sub-contig chunking (#251)

A new `chunk_size` config key can split contigs into sub-contig chunks so a
single large chromosome can be worked by multiple cores. It is **disabled by
default**: benchmarking showed read generation is memory-bandwidth bound, so
chunking does not improve wall time on commodity hardware (and can mildly
regress it). It is kept for CPU-bound or many-memory-channel scenarios. See the
README "Parallel Processing" section, including guidance that **fewer threads
can be faster** on bandwidth-limited desktops.

### Packaging and docs

- **`CITATION.cff`** and a **Bioconda recipe** (`conda-recipe/`) for one-command
  installation (Bioconda submission in review).
- A **NEAT 2 vs NEAT 4 vs rneat** comparison table in the README, framed around
  rneat's durable advantages (native cancer workflow, low/flat memory,
  reproducibility, single-binary deployment).
- A **feedback form** for soliciting real-world `rneat` usage notes that don't
  require direct action. For bugs and features, continue to use the GitHub
  templates.

6/6/2026
========
## rneat v1.16.0

### Per-tissue tumor models completed — SNP/indel half (#202)

The per-tissue cancer models are now tissue-specific in **both** halves. #237
shipped the per-tissue `sv_model`; this adds the per-tissue SNP/indel spectrum,
completing #202. Three full models are bundled:
`tools/cosmic_per_tissue_{BRCA,skin,lung}.json.gz`.

- **`tools/fetch_cosmic_per_tissue_corpus.sh`** — builds a per-tissue COSMIC
  SNP/indel corpus. The COSMIC GenomeScreensMutant VCF is tissue-aggregated, so the
  split joins the sample-level **TSV** (`COSMIC_PHENOTYPE_ID`) to the Classification
  file (`PRIMARY_SITE`) to collect each tissue's COSV mutation IDs, then filters the
  VCF (which already carries VCF-anchored alleles — no indel re-anchoring) to those
  IDs. Reuses the pan-cancer adapter's filtering (drop complex/ambiguous, chr-prefix,
  dedup).
- **`tools/graft_sv_model.py`** — grafts a per-tissue `sv_model` onto a per-tissue
  SNP/indel base. **`tools/build_per_tissue_models.sh`** orchestrates the full chain
  (fetch → `gen-mut-model` train → graft).
- The fits reproduce known tissue biology: **skin**/melanoma is the most
  SNV-dominated (UV C>T; ~1.6% indels), **breast** carries the most indels (~5.2%),
  **lung** sits between (~3.6%).
- `cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz` (pan-cancer SNP/indel + per-tissue
  SV) are retained as the `sv_model` donor for the graft step.
- `docs/cancer_howto.md` reorganized around training a model from your own somatic
  VCF (`gen-mut-model` → `tumor_model:`), with the public-corpus adapters as a
  convenience path.

This completes the core cancer-modeling epic (#129); the remaining advanced-SV
items (#191/#192) are tracked separately as Expanded-realism enhancements.

6/6/2026
========
## rneat v1.15.1

### Cancer wrap-up: parity test, output-dir fix, how-to guide

Tying off the native `gen-cancer-reads` work (#239):

- **`gen-cancer-reads` now creates `output_dir` if it doesn't exist** (via
  `create_dir_all`), matching `gen-reads`. Previously the cancer path built its
  per-pass configs directly and skipped the YAML parser's dir-creation step, so a
  non-existent `output_dir` failed mid-run with a bare `NotFound`.
- **Parity test** (`rneat/tests/cancer_parity.rs`): runs the native subcommand and
  `tools/cancer_simulate.sh` over the same reference/seed/purity and asserts they
  produce identical merged-FASTQ record multisets, per-pass golden VCF bodies, and
  origin-tagged truth-VCF `(chrom,pos,ref,alt,NEAT_ORIGIN)` sets. This is the proof
  gating eventual retirement of the shell script.
- **`fastq_merge` unit test** covering read-name tagging, including the case where
  a quality line begins with `@` (must not be tagged).
- **New copy-paste guide** `docs/cancer_howto.md` with worked examples
  (WGS tumor/normal, per-tissue SV models, purity sweeps, low-purity samples),
  output reference, benchmarking, and model training. The README cancer section is
  trimmed to a short intro pointing at it.

6/6/2026
========
## rneat v1.15.0

### Native `rneat gen-cancer-reads` subcommand (#239)

Tumor/normal cancer simulation is now a first-class subcommand — a native port
of `tools/cancer_simulate.sh` with **no `bcftools`/`awk` runtime dependency**.
`rneat gen-cancer-reads -c <config.yml>` runs two `gen-reads` passes over the
same reference (normal at `(1-purity)*C`; tumor at `purity*C`, sharing the
normal's germline) and merges them:

- **Tagged FASTQs** — each pass's read names are prefixed `N_`/`T_` before
  concatenation, so normal/tumor reads at the same coordinate don't collide as
  QNAMEs (which MarkDuplicates silently drops).
- **Origin-tagged truth VCF** — the two golden VCFs are merged into
  `<prefix>_merged_truth.vcf.gz` with `INFO/NEAT_ORIGIN ∈ {germline, somatic,
  shared}`, resolved from `NEAT_PROVENANCE` via an exact `(contig,pos,ref,alt)`
  key (tumor `denovo` → somatic, tumor `input` → shared, normal-only → germline).

The tumor pass defaults `tumor_mutation_rate` to **1e-5** (typical somatic
burden; `model` defers to the model's fitted rate), and automatically inherits
the breakpoint double-count fix (it lives in the shared `run_neat` path). v1
models a single tumor/normal split; N-way subclonal mixtures are a future
additive extension. Config template: `template_config/gen_cancer_reads_template.yml`;
design + decisions: `docs/archive/cancer_simulator_native_plan.md`.

`tools/cancer_simulate.sh` is retained as the reference implementation and for
the Docker-based caller benchmarks, pending a same-seed parity test.

### Configurable de novo SV length cap (#229)

The de novo SV sampler rejected drawn lengths above a hardcoded `contig_len / 4`
— a cap that keeps its overlap-rejection search tractable on chromosome-scale
references. On small contigs (bacteria, viruses) that cap sits below the trained
DEL median, starving DEL/DUP/INV draws and routing structural signal into BND;
it also blocked any attempt to model large or aneuploid events.

The cap is now a `gen-reads` config knob, **`sv_max_length_fraction`** (default
`0.25`). The default reproduces the historical `contig_len / 4` behavior
bit-for-bit, so existing configs are unaffected. Raise it toward `1.0` for small
contigs or to admit large de novo events — the sampler scales its retry budget
proportionally and its starvation warning now reports the active cap and
fraction. Only de novo SVs are affected; `input_vcf` SVs are untouched. See the
new structural-variants section in `template_config/gen_reads_template.yml` and
the updated aneuploidy notes in `docs/cancer_simulator.md`.

6/6/2026
========
## rneat v1.14.2

### Per-tissue cancer SV models (#202)

Tissue-specific tumor SV models, stratifying the v1.14.0/v1.14.1 PCAWG SvModel by
primary site. Three are bundled alongside the pan-cancer model:
`tools/cosmic_pancancer_sv_{BRCA,skin,lung}.json.gz`. Each is the pan-cancer
COSMIC SNP/indel base with an `sv_model` fitted from that tissue's PCAWG donors.

- `tools/build_pcawg_sv_vcf.py` gains `--sample-sheet` / `--projects` /
  `--tissue-label`: filter the somatic BEDPE/CNA donors to a tissue via
  `aliquot_id → dcc_project_code` (from `pcawg_sample_sheet.tsv`). gnomAD INS is
  germline and never filtered.
- `tools/build_per_tissue_sv_models.sh` drives build → `gen-mut-model` →
  normalize per tissue.
- `tools/fetch_pcawg_sv_corpus.sh` now also fetches `pcawg_sample_sheet.tsv`.

Tissue → PCAWG projects (SV-donor counts): **BRCA** → `BRCA` (211); **skin** →
`SKCM-US,MELA-AU` (106); **lung** → `LUAD-US,LUSC-US` (85). PCAWG's project codes
diverge from TCGA's and SKCM/LUAD alone are underpowered (~37 each), so skin/lung
group related project codes. The fits show expected tissue signal: BRCA is
**DUP-dominant** (Dup 0.317, the tandem-duplicator phenotype) at ~2× SV burden;
skin is **BND-enriched** (0.306); lung is **DEL-dominant** (0.349).

Pick one at run time:
`cancer_simulate.sh --tumor-model tools/cosmic_pancancer_sv_BRCA.json.gz …`.

**Scope:** only the **SV component** is tissue-specific — the SNP/indel side stays
pan-cancer COSMIC (per-tissue SNP/indel via MC3/COSMIC-classification is the
remaining open half of #202). INS length is shared from the gnomAD-derived base
(germline, tissue-agnostic); INS and focal-CNV *rates* stay pan-cancer literature
values, so a tissue's CNV fraction shifts mainly because its SV burden differs.

6/6/2026
========
## rneat v1.14.1

### Patch release: cancer-simulation fixes following the v1.14.0 PCAWG refit

Three corrections to the v1.14.0 cancer-simulation path: a read-generation
performance fix, a realistic default somatic mutation rate, and the long-standing
breakpoint double-counting fix.

#### Read-generation performance: O(log V) variant lookup (#233)

`write_block_fastq` scanned the entire contig's `variant_map` for every fragment
to find the variants overlapping each read — O(fragments × variants). With dense
models this dominated read generation: the bundled COSMIC tumor model's
corpus-aggregated rate (~5.5e-3 → ~220k SNVs on chr22) made the
`cancer_simulate.sh` tumor pass run ~22–27 min where the normal pass took
minutes. `MutatedMap.flagged_positions` is now kept sorted and the two read
windows are located by `partition_point` binary search (and `is_flagged` uses
`binary_search`), turning the per-read lookup into O(log V + hits). Output is
byte-identical (same variant set per read, consumed by key). Measured on chr22
cov2 with the COSMIC model: read-gen **23–27 min → ~27 s (~55×)**. Benefits any
high-variant-density run (deep coverage, dense `input_vcf`), not just cancer.

#### Realistic default tumor somatic mutation rate (#235)

The de-novo mutations added in `cancer_simulate.sh`'s tumor pass are somatic
(germline is carried through via `input_vcf`), but most tumor models — including
the bundled COSMIC model (~5.5e-3) — carry a corpus-aggregated rate that
overstates single-tumor burden by ~50–500×. `--tumor-mutation-rate` now
**defaults to 1e-5** (typical solid tumor, ~10 SNVs/Mb) instead of the model's
fitted rate. On the chr22 fixture this drops somatic SNVs from ~278k to ~507.
Pass `--tumor-mutation-rate model` to defer to the model's own fitted rate. This
is the SNP/indel rate only — SV burden (`sv_model.per_base_rate × --sv-rate-scale`)
and the normal-pass germline rate are unchanged.

**Behavioral change:** runs that relied on the old model-fitted default will now
emit far fewer (realistic) somatic SNVs unless `--tumor-mutation-rate model` is
passed.

#### Breakpoint double-counting fix (#236)

The chimeric pass emits junction-spanning reads for every chimeric SV junction,
but the regular per-contig pass also covered those breakpoints from the unbroken
reference, so a homozygous junction sat at ~2× coverage (het ~1.5×) — the caveat
carried since v1.12.0. The regular pass now drops the broken-allele fraction of
read-pairs that cross a junction (a pair crosses when its R1 `[start,start+rl)`
or R2 `[end−rl,end)` window contains the junction; gap junctions are left alone),
so total junction depth ≈ coverage. Covers **BND, INV, DEL, and CNV-loss**;
**DUP and CNV-gain are intentionally not suppressed** — their chimeric junction
is the tandem boundary (END→POS), a novel adjacency the regular linear pass never
reproduces, so it isn't double-counted. No-op (and no RNG drawn) for contigs
without suppressible SVs, so non-SV runs are byte-for-byte unchanged.

**Behavioral change:** simulated depth at homozygous BND/INV/DEL/CNV-loss
junctions drops from ~2× to ~1× coverage. SV-caller recall is unaffected (the
chimeric junction reads — the actual detection signal — are retained; only the
redundant reference reads are removed).

6/5/2026
========
## rneat v1.14.0

### Data-derived cancer SvModel refit from the PCAWG consensus callsets

Closes #218. The bundled tumor model's `sv_model`
(`tools/cosmic_v104_pancancer_model.json.gz`) was, through v1.13.x, a set of
hand-rolled literature values injected by `tools/inject_cancer_sv_model.py`.
This release replaces them with parameters **counted from real data**: the
PCAWG consensus structural-variant and copy-number callsets (Li et al., *Nature*
578, 2020; n=2,748 donors), with gnomAD-SV supplying the insertion length
distribution.

#### New data pipeline (replaces the hand-injection)

| Tool | Role |
|------|------|
| `tools/fetch_pcawg_sv_corpus.sh` | Download the open-tier PCAWG consensus SV (BEDPE) + CNV (CNA) callsets from the ICGC-ARGO object store (`--endpoint-url https://object.genomeinformatics.org`) and verify checksums. |
| `tools/build_pcawg_sv_vcf.py` | Convert BEDPE (DEL/DUP/INV/BND) + CNA (focal CNV) + gnomAD-SV (INS length) into one symbolic-SV VCF `gen-mut-model` can fit, plus a per-donor counts sidecar. |
| `tools/normalize_pcawg_sv_model.py` | Apply per-tumor rate/type-mix corrections to the fit and splice the result into the bundled model. Data-derived successor to `inject_cancer_sv_model.py`. |

`tools/inject_cancer_sv_model.py` is **deprecated** (banner added); it is kept
only for historical reference / bootstrapping a null-`sv_model` model.

#### What changed in the model

Type probabilities, now counted rather than assumed:

| type | v1.13.x (heuristic) | v1.14.0 (PCAWG) |
|------|--------------------:|----------------:|
| DEL  | 0.350 | 0.292 |
| DUP  | 0.170 | 0.258 |
| INV  | 0.080 | 0.134 |
| BND  | 0.350 | 0.232 |
| CNV  | 0.010 | 0.056 |
| INS  | 0.040 | 0.028 |

The refit ordering (DEL > DUP > BND > INV) now matches the published PCAWG
ranking — deletions most common, tandem duplications second, with inversions an
uncommon class. The old heuristic's "~35% BND / 8% INV" was a misreading of the
PCAWG breakdown. The fitted per-base rate (3.46e-8) lands within 10% of the
prior literature estimate (3.8e-8) — independent confirmation that value was
well-calibrated, now grounded in counted observations. Length distributions and
the CNV copy-number distribution are likewise fit from data; they are
substantially larger than the old values (e.g. DEL median ~115 kb vs the prior
~5 kb), reflecting genuine somatic SV sizes.

##### Inversion double-count fix

A balanced inversion appears in PCAWG BEDPE as **two** junction rows (one
head-to-head `h2hINV`, one tail-to-tail `t2tINV`). The first adapter pass counted
each as a separate inversion, ~2x over-representing INV (0.235 vs the corrected
0.134). The per-donor `h2h ~= t2t` counts confirm the pairing; `build_pcawg_sv_vcf.py`
now counts one INV per `h2hINV` and skips its `t2tINV` partner.

#### chr22 caller-recall validation (no regression)

Same fixture family as v1.13.1 (chr22, PE 30x, purity 0.5, sv-rate-scale 50),
refit tumor model, Manta in tumor/normal mode. Pooled across 3 RNG seeds
(DEL/DUP/CNV via truvari with `--pctseq 0 --typeignore --sizemax 3e8`;
INV/BND via the positional cross-type pass, since truvari cannot score BND):

| SV type | pooled recall (3 seeds) | per-seed range | acceptance |
|---------|-------------------------|----------------|------------|
| DEL | 42/71 = 59% | 48-68% | >=50% PASS (pooled; one seed marginal at 48%, within binomial noise at n~21) |
| DUP | 46/70 = 66% | 56-76% | >=30% PASS |
| CNV | 17/22 = 77% | 67-80% | — |
| INV | 15/25 = 60% (positional) | 43-71% | detected PASS |
| BND | 34/87 = 39% (positional, all somatic) | 32-44% | alignable >=30% PASS |

#### Benchmark scoring fixes (`tools/cancer_sv_benchmark.sh`)

The realistic (large) SV sizes exposed scoring gaps in the stock truvari path
that silently reported 0% recall. Fixed:

- **Truth filter** now keeps only symbolic SV records (`INFO/SVTYPE`), excluding
  literal large indels that belong to the SNV/indel benchmark and otherwise
  count as false negatives for a structural caller.
- **truvari flags** `--pctseq 0` (symbolic ALTs carry no sequence), `--sizemax
  3e8` (default 50 kb excluded realistic SVs), `--typeignore`, `--refdist 500`.
- **New positional BND + INV recall pass** (step 5c): truvari does not score
  BND-type records, and Manta emits translocations/inversions as BND-form
  junctions. A ±500 bp cross-type matcher (reproducing v1.13.1's method)
  reports their recall.

#### Known limitations / caveats

- **INS rate is literature, not counted.** PCAWG retrotransposition (MEI) calls
  are controlled-access; gnomAD-SV supplies the INS *length* distribution
  (germline, but mobile-element sizes transfer), while the INS *fraction*
  (~0.028) remains a somatic literature estimate. INS content is still random
  novel sequence (#190).
- **CNV rate is literature, not counted.** Raw focal-CNA segments (~135/donor)
  are segmentation, not discrete events; the CNV `type_probability` uses a
  literature per-tumor count while the CN distribution (capped at 10) and length
  are data-derived.
- **Aneuploidy** (whole-arm/chromosome CNV) remains on the depth-modulation /
  `input_vcf:` path (#189), not the SvModel sampler.

6/2/2026
========
## rneat v1.13.1

### Patch release: BND chimeric reads now produce caller-detectable split-read signal

Closes #224. v1.13.0's validation found Manta got 0% recall on somatic BND truth records despite the BND chimeric path emitting reads to the FASTQ. Investigation traced this to two compounding issues, both fixed here.

#### Issue 1 — Phantom N-region truth records (`sample_variants` anchor check too narrow)

`SvModel::sample_variants` rejected positions where the literal anchor base was N (line 473), but **passed positions where the anchor was a non-N base sitting at the boundary of a long N-tract**. The chimeric pair generation then sampled a `read_len`-sized piece extending into the N-tract, producing reads BWA couldn't place uniquely. On the v1.13.0 chr22 fixture, ~37% of de novo BND truth records were unrecoverable for this reason (anchor passed the single-base check; chimeric piece landed in the N-tract).

**Fix:** new helper `anchor_window_alignable(sequence, pos, window=200)` in `common/src/structs/sv_model.rs` rejects positions where the ±200bp window is >20% N. Applied to both the anchor and (for non-BND types) the END position. For de novo BNDs, the mate-position picker now retries up to 32 times to find a clean mate before skipping the record.

#### Issue 2 — Imbalanced chimeric offset (`process_chimeric_variants`)

The chimeric junction offset was sampled uniformly from `[1, min(frag_len-1, read_len-1)]`, producing tiny pieces (offset=1 or 2) ~1-2% of the time. Small offsets yield a 1-2 bp local piece and ~150bp mate piece — BWA aligns the longer piece confidently but the supplementary alignment for the shorter piece is too short to anchor uniquely, so BWA emits it with `MAPQ=0`. Manta filters `MAPQ<20` supplementary alignments out of its candidate pool, so these "weak" split-reads contribute no signal.

**Fix:** new helper `balanced_chimeric_offset()` in `rneat/src/gen_reads/utils/runner.rs` samples from `[read_len/4, frag_len - read_len/4]`, guaranteeing both junction pieces are ≥ `read_len/4` (~38bp at read_len=151). That clears BWA-MEM's anchor threshold for confident split alignment and pushes the supplementary MAPQ past Manta's filter. Applied uniformly to all five chimeric paths (BND/INV/DEL/DUP/CNV).

#### Cancer SV benchmark: Manta BND recall now 77% on alignable subset

Same chr22 PE 30× / purity 0.5 / sv-rate-scale 50 fixture (different RNG seed than v1.13.0, so absolute truth counts differ slightly). Cross-type scoring (a truth event matches any Manta call within ±500bp at the same span endpoints; necessary because Manta classifies the v1.13.1 chimeric junctions as DEL/DUP based on the read-orientation pattern that BWA emits, the same scoring artifact as v1.13.0's INV-as-BND-quartet):

| SV type | v1.13.0 recall | v1.13.1 recall |
|---------|----------------|----------------|
| DEL     | 21/39 (54%)    | **24/37 (65%)** |
| DUP     | 12/22 (55%)    | 6/15 (40%) — different truth draw; see note |
| INV     | 4/4 (100%)     | 5/9 (56%) — larger sample, more realistic |
| BND (all)        | 0/35 (0%)  | **17/39 (44%)** |
| BND (alignable)  | 0/22 (0%)  | **17/22 (77%)** ✅ |

A note on the DUP "regression": v1.13.0's truth draw had 22 DUPs (some of which were small / easy); v1.13.1's draw at a different RNG seed had 15 DUPs of a different distribution. Without a side-by-side same-seed comparison, the 55% vs 40% is not strictly apples-to-apples. The chimeric DUP path itself is unchanged from v1.13.0 (Fix B applies symmetrically; Fix A only affects which positions get truth records).

Same caveat for INV: v1.13.0's 4/4=100% was a small-sample fluke; v1.13.1's 9-record sample at ~56% is more representative.

The BND alignable result (77%) is the headline: clears the v1.13.1 acceptance criterion (≥30%) by a wide margin.

#### Cross-type matching: what changed

Manta represents simulated junctions according to the read-orientation pattern BWA observes, not the truth-VCF's `SVTYPE` field. For a `<DEL>`-style chimeric junction Manta calls `SVTYPE=DEL`; for a tandem-DUP-style junction Manta calls `SVTYPE=DUP`. For inversion-like junctions (revcomp piece on the same chromosome) Manta sometimes calls `SVTYPE=BND` quartet, sometimes `SVTYPE=DUP:TANDEM` depending on the strand-flip pattern that lands at BWA's split-alignment site. The v1.13.1 BND chimeric pair (case `t]p]`: forward local + revcomp mate) produces a read pattern Manta interprets as DEL or DUP depending on the relative coordinates — that's what the "10 as DEL, 7 as DUP" breakdown reflects.

This is a known characteristic of Manta's somatic-mode classification and parallels the v1.13.0 INV-as-BND-quartet finding. Caller benchmarks should use cross-type matching for SV simulation truth.

#### Tests

- 328 lib + 223 binary unit tests pass.
- NEW `anchor_window_alignable_basics` (helper unit test) and `sample_variants_skips_n_adjacent_clean_anchors` (regression test for the wider N-window gate).
- All v1.12.x and v1.13.0 SV integration tests still passing (`bnd_fastq.rs`, `bnd_roundtrip.rs`, `inv_fastq.rs`, `del_chimeric.rs`, `dup_chimeric.rs`, `cnv_chimeric.rs`, `multi_sv_integration.rs`, `cosmic_bundle.rs`).

#### Known limitations carried into v1.13.1

- **Mid-chromosome N-tracts.** The ±200bp anchor-window check catches the common case where a position sits at the boundary of a long N-tract. It doesn't catch the rarer case where the chimeric piece (`read_len`-sized) reaches into a more distant N-tract while the immediate window is clean. A future refinement could check the actual chimeric piece extent based on SV type. For now, a small number of N-region phantom truth records may still appear (~17 of 39 BNDs in the chr22 fixture).
- **Breakpoint double-counting** remains the same caveat as v1.12.0 — regular per-contig reads at the breakpoint still co-exist with chimeric junction reads. Doesn't block Manta detection at this coverage but could in lower-coverage / lower-purity scenarios.

---

6/1/2026
========
## rneat v1.13.0

### Minor release: symbolic SVs now emit split-read + discordant-PE signal, literal long deletions land in CIGAR

Closes #220 and #221. v1.12.1's validation revealed that the v1.12.0 SV pipeline produced truth records but no detectable signal for split-read/PE-discordant SV callers like Manta — 0% recall on the bundled-model somatic SVs. This release closes that gap for DEL/DUP/CNV and fixes a separate pre-existing bug where literal long deletions from the standard indel_model appeared in the truth VCF but not in the FASTQ.

#### #221 — Literal long deletions now produce CIGAR D-ops
- **Bug.** In `generate_read`, the alt-allele branch wrote the 1-byte ALT for a literal Deletion but never advanced `seq_index` past the deleted REF bases. For a REF=77bp / ALT=1bp variant the read continued transcribing the deleted bases from the unbroken reference and emitted CIGAR `<read_length>M` with no D op — leaving the deletion completely invisible to downstream callers despite rneat's own AD counter reporting alt support.
- **Fix.** Apply the same skip+D-op pattern that `SequencingErrorType::DeletionError` already uses, capped at the remaining fragment buffer. Pure 5-line change in `common/src/file_tools/fastq_tools.rs`. Pinned by `test_literal_long_deletion_emits_d_ops` (asserts ≥49 D ops for a 50-bp homozygous deletion).
- **Impact.** Pre-existing bug since v1.10+. Affected every gen-reads run whose indel_model produced deletions long enough to be caller-relevant (typically ≥10bp). Most visible in cancer benchmarks where the COSMIC-trained indel distribution has a long tail.

#### #220 — Symbolic DEL/DUP/CNV gain chimeric junction reads
- **Gap.** Pre-v1.13.0 symbolic `<DEL>` / `<DUP>` / `<CNV>` records only modulated depth via `build_coverage_multipliers`. Manta, DELLY, GRIDSS and similar split-read/PE-discordant SV callers need junction-spanning reads to call SVs reliably; depth alone isn't enough. v1.12.0 already had chimeric paths for BND (#187) and INV (#188), but DEL/DUP/CNV were depth-only.
- **Implementation.** Extends `process_chimeric_variants` in `rneat/src/gen_reads/utils/runner.rs` with three new branches:
  - **DEL.** `generate_del_pair` + `get_del_pieces` stitch REF[..=POS] (anchor included) with REF[END..] (post-deletion). BWA aligns the result as discordant PE pairs (mate distance overshoots by END-POS) and split-read alignments (soft-clip at POS realigns to REF[END..]).
  - **DUP (tandem).** `generate_dup_pair` + `get_dup_pieces` stitch REF[end-k..end] with REF[location..location+k] — the tandem-copy boundary. BWA sees the inverted-coordinate signature (end before start) and emits split-read or wrong-orientation PE signal.
  - **CNV.** Dispatches to the DEL or DUP path based on whether INFO/CN is below or above the diploid baseline. CN=0 → DEL-like; CN=4 → DUP-like with the magnitude of CN deviation driving coverage scale. CNVs without INFO/CN are skipped with a debug log.
- **Pinned by** new integration tests `rneat/tests/del_chimeric.rs`, `rneat/tests/dup_chimeric.rs`, and `rneat/tests/cnv_chimeric.rs` (the last covers both CN<ploidy and CN>ploidy dispatch with cross-type negative assertions).

#### Bug fix: PE pair-sync regression exposed by #221
- **Bug.** In `write_block_fastq`, if `generate_read` returned `TruncatedRead` for r2, the function would `continue` the fragment loop — but r1 had already been written to buffer1. This desynced the r1/r2 streams and BWA-MEM aborted with `paired reads have different names`. Always reachable in principle; became common after #221's literal-DEL skip+D-op fix legitimately advances seq_index past the deleted bases and can exhaust the buffer for long deletions near a fragment edge.
- **Fix.** Generate r2 BEFORE writing r1; if r2 returns TruncatedRead, drop r1 alongside it. r1 and r2 now succeed together or fail together. Caught and resolved during chr22 end-to-end validation (BWA error on first pass, clean run after the fix).

#### Cancer SV benchmark: Manta recall now meets v1.13.0 acceptance criterion (≥50% on DEL/DUP)
Same chr22 PE 30× / purity 0.5 / sv-rate-scale 50 fixture used in v1.12.1 validation. Comparison of Manta somatic SV recall (±500bp position match, with INV→BND-quartet matching since Manta reports inversions as paired BND records, not `<INV>`):

| SV type | v1.12.1 recall | v1.13.0 recall |
|---------|----------------|----------------|
| DEL     | 0/32 (0%)      | **21/39 (54%)** |
| DUP     | 0/13 (0%)      | **12/22 (55%)** |
| CNV     | 0/2 (0%)       | **1/1 (100%)**  |
| INV     | 0/5 (0%)       | **4/4 (100%)** — detected as Manta BND quartets |
| BND     | 0/32 (0%)      | 0/35 (0%) — caller doesn't pick up the chimeric junctions |
| **TOTAL** | **0/84 (0%)** | **38/101 (38%)** |

Precision: 100% — every Manta PASS somatic call matched a real truth record at ±500bp.

A note on the INV scoring: Manta emits inversions as a pair of paired BND records (4 total per INV — start junction and end junction, each described from both sides) rather than as a single `<INV>` record. A type-strict matcher scores INV truth against `<INV>` calls and gets 0%, which is what we initially reported. The corrected matcher recognizes the BND-quartet representation and finds all 4 truth INVs were detected.

#### Known limitations carried into v1.13.0
- **BND still at 0% Manta recall.** Of the 35 somatic BND truth records, ~15 have at least one endpoint in chr22's N-rich centromere/telomere region (where coverage drops to zero because BWA can't align all-N reads — this is a reference-structure limitation, not an rneat bug). Of the remaining ~20 alignable BNDs, BND chimeric reads ARE present in the tumor BAM at the breakpoint positions, but Manta doesn't pick them up as somatic calls. Tracked at #224 — likely needs either stronger junction signal (higher num_frags multiplier for BND) or a caller that's less conservative for short-fragment BND signatures (GRIDSS).
- **Breakpoint double-counting.** The regular per-contig pass still generates unbroken-reference reads at SV breakpoints. For homozygous SVs the breakpoint locus ends up with regular reads PLUS junction reads — roughly 2× coverage at the boundary. Doesn't affect Manta's call decisions for DEL/DUP (verified by the recall numbers above) but inflates depth at the junction. Tracked as a v1.13.1 follow-up.

#### Tests
- 328 lib + 223 binary unit tests pass.
- NEW unit test: `test_literal_long_deletion_emits_d_ops`.
- NEW integration tests: `del_chimeric.rs` (1), `dup_chimeric.rs` (1), `cnv_chimeric.rs` (2).
- v1.12.0 SV integration tests still passing (`bnd_fastq.rs`, `bnd_roundtrip.rs`, `inv_fastq.rs`, `multi_sv_integration.rs`, `cosmic_bundle.rs`).

---

5/31/2026
=========
## rneat v1.12.1

### Patch release: cancer SV model bundled into the COSMIC tumor artifact

Closes #217. Before this patch the bundled `tools/cosmic_v104_pancancer_model.json.gz` had `sv_model: null`, so `cancer_simulate.sh`'s tumor pass produced zero symbolic SVs even with `--sv-rate-scale > 0`. The v1.12.0 BND / INV / INS pipeline was wired up but inert when driven by the bundled tumor model — the only way to exercise it was to train a custom model or supply BND / INV through `input_vcf:`.

#### What changed
- **NEW `tools/inject_cancer_sv_model.py`.** Injects a literature-derived `sv_model` component into any `MutationModel` JSON. Parameters are pan-cancer PCAWG defaults (Li et al., *Nature* 578, 112–121, 2020): `per_base_rate = 3.8e-8` (≈ 110 somatic SVs / tumor over 2.9 Gb), type distribution 35% DEL / 17% DUP / 35% BND / 8% INV / 4% INS / 1% CNV (BND over-represented vs gnomAD-derived germline default's 19%), length log-normals centered on cancer-typical medians (DEL ~5 kb, DUP ~30 kb, INV ~22 kb, CNV ~1 Mb, INS ~300 bp), and a bimodal copy-number distribution covering both homozygous loss (CN=0) and high-level amplification (CN=5+).
- **Re-packaged `tools/cosmic_v104_pancancer_model.json.gz`** (9.2 KB → 9.8 KB) with the injected `sv_model`. Existing SNV/indel parameters are unchanged — this is purely additive.
- **Regression test `rneat/tests/cosmic_bundle.rs`.** Two tests: one asserts the bundled file deserializes with a non-null `sv_model` carrying all six SV types and cancer-typical rates (catches accidental copy-from-germline-default); a second runs `gen-reads` end-to-end with the bundle wired in via `mutation_model:` + high `sv_rate_scale` and confirms symbolic / BND records reach the output VCF.
- **`tools/cancer_simulate.sh` help text** updated to reflect that the bundled tumor model now carries an SV component.

#### Why heuristic literature values, not a data refit
COSMIC GenomeScreensMutant is SNV/indel-only. A real cancer-SV refit needs a separate corpus — PCAWG consensus SV calls or COSMIC StructVar — and is tracked at **#218** for v1.13. Shipping literature-derived defaults now means users can exercise the v1.12.0 BND / INV / INS pathways with `cancer_simulate.sh --sv-rate-scale 1.0` today; the v1.13 refit will fine-tune the rates without changing any code paths.

#### Tests
- 327 lib + 223 binary unit tests pass (unchanged from v1.12.0).
- NEW `rneat/tests/cosmic_bundle.rs`: 2 tests, both passing.
- v1.12.0 SV integration tests still passing (`bnd_fastq.rs`, `bnd_roundtrip.rs`, `inv_fastq.rs`, `multi_sv_integration.rs`).

---

5/30/2026
=========
## rneat v1.12.0

### Minor release: stage-2 cancer-SV credibility — BND, INV, de novo INS

Closes the three foundational SV gaps that v1.11's cancer MVP couldn't yet generate. With these shipped, the rneat simulator covers the structural axis that most clinically-recognizable cancer SVs depend on: BCR-ABL / PML-RARA / EWSR1-FLI1 translocations (#187 BND), inv(16) AML / inv(3) MDS inversions (#188 INV), and L1 / Alu / SVA mobile-element insertions (#190 INS).

End-to-end verified on a synthetic 5 kb fixture with `sv_rate_scale=50.0` plus an input VCF carrying one BND pair and one INV: 125 BND records + 1 INV + 138 de novo literal-INS records in the output VCF, with ~2,245 BND-chimeric + 60 INV-chimeric reads in the FASTQ. The first de novo INS's 20-mer prefix appears 35× in the FASTQ — confirming the inserted bases actually reach the simulated reads, not just the truth VCF.

#### #187 — De novo `<BND>` translocations with junction reads (PR #198)
- **`SvModel` learns and samples `<BND>` records.** BND records bypass the 50 bp SV-length minimum (junction variants are conceptually 1 bp); `gen-mut-model` fits the frequency and genotype distribution from training VCFs; `gen-reads` samples and emits BND records when `sv_rate_scale > 0`.
- **Junction reads via a new chimeric-pair path.** `process_chimeric_variants` → `generate_chimeric_pair` → `get_bnd_pieces` + `get_stitched_sequence` synthesizes chimeric read pairs that span both sides of the breakend, correctly reverse-complementing per VCF 4.2 breakend orientation (the four `t[p[`, `t]p]`, `[p[t`, `]p]t` forms).
- **Paired-BND deduplication.** Canonical-ID HashSet keyed by `format!("BND_{:?}", (min_contig, min_pos, max_contig, max_pos))` ensures each junction is processed exactly once even though VCFs encode it from both sides.

#### #188 — `<INV>` inversions with read-strand flipping (PR #200)
- **Two-junction chimeric model.** Each de novo INV generates reads spanning BOTH breakpoints (start of inversion + end of inversion), each junction independently sampling `num_frags` fragments. `generate_inv_pair` + `get_inv_pieces` handle the orientation flipping: `REF[..POS-1] | RC(REF[POS..END])` at junction 1, `RC(REF[POS..END]) | REF[END+1..]` at junction 2.
- **Smarter fragment-length sampling.** Paired-end mode uses attempt-bounded retries for minimum length; single-end uses `read_len + 32 bp` padding so sequencing-error deletions don't truncate the read. Applied to BND too.
- **Graceful `TruncatedRead` handling.** Edge-case junctions that exhaust the available sequence are logged and skipped rather than aborting the whole simulation.

#### #190 — De novo `<INS>` with novel-sequence generation (PR #216)
- **Literal-Insertion representation, not symbolic.** De novo INS records are emitted as `REF="A"`, `ALT="A<novel-bases>"` — the same shape Manta / DELLY produce when they resolve the inserted sequence. Routes through gen-reads' existing literal-insertion machinery; the novel bases appear in fragment-loop reads at the locus through the existing per-base coin-flip.
- **Composition-weighted novel sequence.** New helper `sample_novel_insertion_bases` draws each base from the single-base composition of a ±250 bp window around the anchor. Falls back to uniform ACGT for all-N regions. Trinucleotide-context sampling is a deferred refinement.

#### Behavioral changes — read before upgrading
- **`gen-reads` now generates BND / INV chimeric reads when `sv_rate_scale > 0`** AND the SvModel has BND / INV in its type pool. The bundled gnomAD-SV-derived `default_sv_model()` does **not** include BND or INV (they were filtered out of the gnomAD-SV training corpus). Users wanting these in default runs must either train a model from a corpus that includes them (e.g. PCAWG consensus calls) or supply BND / INV records via `input_vcf:` for round-trip.
- **De novo `<INS>` is now a literal-ALT Insertion record.** Pre-v1.12.0 the model picked INS as a type but emitted a symbolic `<INS>` with no allele. Output VCFs from v1.11.x with de novo INS records will look different in v1.12.0 (resolved sequence in REF/ALT instead of empty symbolic).
- **Chimeric read QNAMEs carry a `junction` tag and a `frag_idx` 16-hex tag** so reads from the same BND / INV at the same fragment iteration don't collide. Anything pattern-matching the read-name format needs updating to handle `RNEAT_chimeric_<c1>_<pos>_<c2>_<mate>_<16-hex>/<mate-id>` (BND) and `RNEAT_chimeric_INV_<contig>_<pos>_<end>_<junction>_<16-hex>/<mate-id>` (INV).
- **`gen-mut-model` now fits and emits BND / INV components** in trained models. Existing v1.11.x trained models still load — the BND/INV fields are absent so `gen-reads` won't generate those types unless the model is retrained.

#### Known limitations
Restated for v1.12.0 (carried over from earlier releases unless marked):
- **BND / INV breakpoint double-counting.** The regular per-contig pass still generates reads covering BND / INV breakpoint positions (reading from the unbroken reference), so the breakpoint locus ends up with regular reads PLUS junction reads — roughly 2× coverage for homozygous variants, closer to correct for heterozygous. A proper fix would teach the regular pass to skip the broken-allele fraction at these positions; deferred to v2.
- **Symbolic SVs emit `.,.:.:.` for FORMAT/AD/DP/AF.** That set just got bigger — BND and INV are symbolic and use the placeholder shape. INS is now literal and gets real counts, so the FORMAT field is asymmetric across SV types.
- **No mobile-element library for de novo INS.** Insertions are random novel sequence weighted by local composition, not drawn from L1 / Alu / SVA / HERV consensus. The library-driven path is a #190 v2 extension.
- **Viral-integration simulation is out of scope.** HPV / HBV / EBV integration into the host reference needs its own design (a viral FASTA as an optional input, chunked-integration sampler).
- **Per-tissue tumor models still pending (#202).** Both training adapters (MC3, COSMIC) produce pan-cancer corpora.
- **Cancer-specific SvModel defaults.** Bundled `default_sv_model()` is still gnomAD-SV-derived (germline). PCAWG-derived cancer defaults would be the natural complement; tracked at the v2 follow-up.

#### Companion: somatic SV benchmark pipeline
- **NEW `tools/cancer_sv_benchmark.sh`.** Sibling to `tools/cancer_benchmark.sh` (which scores Mutect2 with som.py); this script runs **BWA-MEM → Manta tumor/normal → truvari bench** against the rneat truth VCF and produces precision/recall/F1 by SV type. Closes the stage-2 exit criterion ("a current somatic SV caller … with high enough recall to use rneat as a benchmark fixture") with a runnable artifact rather than a deferred claim. Pinned Docker images (Manta 1.6.0, truvari 4.3.1, BWA 0.7.18, samtools 1.21, bcftools 1.21); BAMs are reused from `cancer_benchmark.sh` if you've already run it on the same output dir. **Requires paired-end FASTQ inputs** — Manta, DELLY, and GRIDSS all use fragment-orientation signal for SV detection and refuse to run on single-end data. The script's `--help` and pre-flight validation point users at the right `cancer_simulate.sh --paired-ended` invocation if they need to regenerate.

#### Cancer-simulate merge fix for BND/INV truths
- **`tools/cancer_simulate.sh`'s merge step now handles per-pass VCFs that carry BND/INV records.** Earlier, `bcftools concat | bcftools sort` did an internal BCF translation that errored out on undeclared INFO fields (`MATEID` on BND records is the canonical trigger). The merge step now injects `##INFO=` declarations for the standard SV fields (MATEID, SVTYPE, END, SVLEN, CN, IMPRECISE, CIPOS, CIEND, BND_DEPTH, MATE_BND_DEPTH) via `bcftools annotate -h` before the concat. Verified end-to-end on a synthetic 5 kb fixture with a 2-record BND pair + 1 INV in the germline — all three records survive the merge tagged `NEAT_ORIGIN=shared`, with MATEID preserved.

#### Tests
- 327 lib + 223 binary unit tests pass.
- Per-feature integration tests: `rneat/tests/bnd_fastq.rs`, `rneat/tests/bnd_roundtrip.rs`, `rneat/tests/inv_fastq.rs`.
- NEW combined-SV integration test: `rneat/tests/multi_sv_integration.rs` — exercises BND + INV + de novo INS in a single gen-reads run on a synthetic fixture and verifies all three pathways flow into the FASTQ / VCF correctly.

5/30/2026
=========
## rneat v1.11.1

### Patch release: BGZF-framed golden VCFs + spec-conformant `##source=` header

Two small VCF-output cleanups that came out of the v1.11.0 cancer-MVP work — both let `tools/cancer_simulate.sh` drop its `bcftools view -O z` transcode workaround.

- **gen-reads' golden VCF is now BGZF-framed, not plain gzip (#206).** The writer switched from `flate2::write::GzEncoder` to `noodles::bgzf::io::Writer` — same gzip magic at the byte level, so any gzip-aware reader still works, but it's now tabix-indexable directly (`bcftools index -t` previously refused with "format that cannot be usefully indexed"). The cancer-simulator merge step's `bcftools view -O z` transcode pass is gone — it was a v1.11.0 workaround for exactly this. `noodles::bgzf` was already a transitive dep via the BAM writer; no new crates added.
- **Header line `##Generated by rusty-neat` replaced with `##source=rusty-neat-<version>` (#207).** The old line wasn't VCF 4.2 spec-conformant (`##KEY=VALUE` is the required shape per §1.4.1); bcftools warned `"Could not parse the header line"` on every rneat-produced VCF. The new line carries the workspace version too (`##source=rusty-neat-1.11.1`), so any downstream consumer that inspects it gets useful provenance for free.

#### Regression tests
- **NEW** `test_write_vcf_output_is_bgzf_framed`: reads the first 16 bytes of a write_vcf output and pins the BGZF-defining markers (gzip magic + `FLG.FEXTRA=1` + `"BC"` subfield identifier at bytes 12–13). A regression to plain gzip output would otherwise be invisible to existing reader tests because both formats decode identically through MultiGzDecoder.
- `test_write_vcf` updated to expect `##source=rusty-neat-` rather than `##Generated by rusty-neat`.
- `test_read_vcf_gzipped_input` keeps writing plain gzip (locally scoped to that test) to verify the reader still handles non-bgzf gzip input — important for users feeding external gzipped VCFs into rneat.

#### Behavioral notes
- **Anything pattern-matching the old `##Generated by rusty-neat` header line** in downstream tooling needs an update to `##source=rusty-neat-` (or `##source=` more loosely).
- **Output VCFs are byte-different** from v1.11.0 because BGZF segments the gzip stream into ~64 KB blocks. Decompressed content is identical. Cache-by-hash systems using whole-file hashes will see invalidations.

5/30/2026
=========
## rneat v1.11.0

### Minor release: cancer-MVP — origin-tagged truth VCFs, FORMAT/AD/DP/AF, COSMIC adapter

First release with end-to-end cancer-simulation capability. The bundled COSMIC pan-cancer model produces a benchmark-grade somatic truth VCF (origin tags + per-variant VAF) that `hap.py` / `som.py` / `vcfeval` consume without preprocessing. Validated on chr22 + GATK Mutect2 end-to-end: 67.5% SNV recall at 96.7% precision against the somatic-tagged truth set, at default settings.

#### New: tumor mutation training corpora (#183)
- **`tools/fetch_cosmic_corpus.sh` (#201).** Sibling to the existing `tools/fetch_tumor_corpus.sh` (TCGA MC3 MAF adapter). Converts COSMIC's GenomeScreensMutant VCF into a dedup, chr-prefixed, gen-mut-model-consumable corpus. v104 collapses ~50.6M raw records → ~16.7M unique (deduplicates ~11.8M multi-transcript annotation duplicates and remaps `MT → chrM`). `--train` chains directly into `rneat gen-mut-model`.
- **`tools/cosmic_v104_pancancer_model.json.gz`.** A 9.2 KB pre-trained pan-cancer mutation model bundled in `tools/` so the cancer simulator can be exercised end-to-end without a 1 GB COSMIC download. Interim placement — final hosting decision tracked at **#186**. **Important caveat:** the fitted `mutation_rate` of ~5.5e-3 is corpus-aggregated ("fraction of bp with any observed mutation across the entire COSMIC catalog"), not per-tumor burden. Use `cancer_simulate.sh --tumor-mutation-rate 1e-5` (or similar) to shape per-sample burden — see "Behavioral changes" below.

#### New: variant-origin annotation in golden VCFs (#185, PR #204)
- **`INFO/NEAT_PROVENANCE` on every gen-reads output record.** Two-layer split: `gen-reads` itself emits `denovo` / `input` per record (a generic "where in this pass did this variant come from"), staying cancer-agnostic. The cancer-orchestration script then resolves cross-pass provenance into `INFO/NEAT_ORIGIN = germline | somatic | shared` via `bcftools isec` + an awk annotation pass — variants present only in the normal pass are `germline`, only in the tumor pass (de novo) are `somatic`, in both passes are `shared` (germline that round-tripped). The merged `<prefix>_merged_truth.vcf.gz` carries both `NEAT_PROVENANCE` (audit trail) and `NEAT_ORIGIN` (benchmark label).
- **`tools/cancer_simulate.sh` merge step.** New post-processing block produces `<prefix>_merged_truth.vcf.gz` automatically when `bcftools` and `bgzip` are available; gracefully skipped otherwise with per-pass VCFs preserved. Includes a `bcftools view -O z` transcode step to convert gen-reads' plain-gzip output to bgzip framing (so tabix indexing works); the transcode becomes a no-op once **#206** lands and `write_vcf` emits bgzip directly.

#### New: FORMAT/AD, FORMAT/DP, FORMAT/AF per variant (#176, PR #205)
- **Golden VCF FORMAT column expanded from `GT` to `GT:AD:DP:AF`.** Populated from a per-variant allelic-depth counter incremented by the gen-reads fragment loop at the existing per-base coin-flip site — for each emitted read that overlaps a variant locus, the counter records whether the haplotype draw landed on ref or alt. The counter accumulates per-contig through the rayon parallel join and is consumed by `write_vcf` at output time.
- **Literal variants** carry computed values (e.g. `0/1:18,12:30:0.4000` for a het with binomial 50/50 noise on 30× coverage).
- **Symbolic SVs** emit `.,.:.:.` placeholders — span-based depth semantics don't fit a point-based counter; tracked as a v2 follow-up.
- This is what makes the truth VCF consumable by `hap.py` / `som.py` / `vcfeval` for VAF-stratified TP/FP/FN scoring.

#### New: cancer simulator orchestration knobs (#184 follow-ups)
- **`--normal-mutation-rate` and `--tumor-mutation-rate` on `cancer_simulate.sh`.** Each forwards to gen-reads' existing `mutation_rate:` config field, overriding the model's fitted value. Per-tumor burdens are typically 1e-5 to 1e-4 per bp — the bundled COSMIC model's 5.5e-3 reflects aggregate catalog coverage and should be overridden for single-tumor simulation.
- **Position-sorted golden VCF output.** Pre-fix, `write_vcf` iterated `block_map.variant_map` (a HashMap), emitting records in nondeterministic order — `bcftools index -t` (tabix) refused to index the output. The writer now collects all records per contig into a sorted Vec before emission.
- **Per-pass read-name prefixing in the FASTQ merge step.** The normal + tumor passes were producing FASTQs whose read names collided whenever both passes happened to sample the same start coordinate. The merge step now prefixes each pass's reads with `N_` / `T_` before `cat`-merging so the resulting FASTQ has no cross-pass name collisions.

#### New: per-fragment uniqueness tag in read names (#210, PR #211)
- **gen-reads QNAMEs are now globally unique within a pass.** Pre-fix, the read name was `RNEAT_generated_<contig>_<start>_<end>/<mate>` — two fragments at the same coordinate (which the birthday paradox guarantees at 30×+ coverage; we observed ~250k such collisions per pass on chr22) shared an identical QNAME. Picard MarkDuplicates would silently classify them as PCR duplicates and remove them, halving effective coverage at the affected loci. The new format appends a 16-char hex frag-index tag: `RNEAT_generated_<contig>_<start>_<end>_<uniq>/<mate>`. `rneat filter-reads`'s record-name parser updated to read the last three underscore-separated tokens as `(start, end, uniq)` and ignore `uniq`.

#### New: somatic-caller benchmark pipeline (PR #208)
- **`tools/cancer_benchmark.sh`.** Drives the standard cancer-MVP validation: align per-pass FASTQs → tumor + normal BAMs (BWA-MEM 0.7.18), call somatic variants (GATK Mutect2 4.5.0.0 in tumor/normal mode), filter (FilterMutectCalls), and score against the rneat truth VCF (Illumina som.py via the GA4GH `jmcdani20/hap.py:v0.3.12` image). Every tool runs in a pinned Docker (or Podman-as-Docker) container; per-step outputs are idempotent (re-running after a single-step change is cheap). A `--truth-filter` flag defaults to `INFO/NEAT_ORIGIN="somatic"` so the recall metric measures real somatic-caller performance rather than the (correct) refusal to call germline-carry-through records.

#### Behavioral changes — read before upgrading
- **Golden VCF FORMAT column changed from `GT` to `GT:AD:DP:AF`.** Any downstream pipeline that hardcoded the FORMAT field as just `GT` (parsing the SAMPLE column as a single genotype string) will need updating. The new shape is the standard VAF-aware truth-set format expected by hap.py / som.py / vcfeval.
- **Output VCF records now have an INFO column.** Pre-v1.11.0 emitted `.` (empty INFO); now emits `NEAT_PROVENANCE=denovo|input`. Symbolic SVs that came in via `input_vcf:` have `NEAT_PROVENANCE` appended to their preserved INFO (SVTYPE / END / SVLEN / CN survive verbatim).
- **Output VCF records are now position-sorted within each contig.** Was nondeterministic HashMap iteration order; broke tabix indexing. The new sorted-Vec-per-contig emission produces tabix-indexable output once #206 (bgzip framing) also lands.
- **gen-reads read names now carry a uniqueness tag.** Format went from `RNEAT_generated_<contig>_<start>_<end>/<mate>` to `RNEAT_generated_<contig>_<start>_<end>_<16-hex>/<mate>`. Anything that pattern-matched the old format (regex, awk field splits) needs updating. `rneat filter-reads` was updated in-tree; external consumers may not have been.
- **`tools/cancer_simulate.sh` produces a merged truth VCF.** New file `<prefix>_merged_truth.vcf.gz` appears in the output directory alongside the per-pass VCFs whenever bcftools is installed. Per-pass VCFs are still produced and preserved.

#### Known limitations
Restated for v1.11.0 (carried over unless marked):
- **Symbolic SVs (`<DEL>` / `<DUP>` / `<CNV>`) emit `.,.:.:.` on FORMAT/AD/DP/AF.** Span-based depth semantics require a different counter shape; tracked as a v2 follow-up.
- **The bundled COSMIC pan-cancer model's `mutation_rate` is corpus-aggregated, not per-tumor.** Always override with `--tumor-mutation-rate` for single-tumor simulation. The model's docstring and `docs/cancer_simulator.md` flag this explicitly.
- **gen-reads writes plain gzip, not bgzip.** `tools/cancer_simulate.sh` works around it by transcoding via `bcftools view -O z` before indexing. Once **#206** lands, the transcode step becomes a no-op and the script can drop it.
- **`##Generated by rusty-neat` header line is malformed** (should be `##source=rusty-neat`). Tracked at **#207**; non-fatal — bcftools warns but proceeds.
- **No per-tissue tumor models yet.** Both adapters produce pan-cancer corpora. Per-tissue split (BRCA / SKCM / LUAD) tracked at **#202**.
- **`<INV>` / `<BND>` SVs not yet generated de novo.** Tier-1 SV gaps for cancer simulation; tracked at **#187** / **#188**.

#### Validation
End-to-end smoke run on chr22 + bundled COSMIC model + `--tumor-mutation-rate 1e-5`:
- Pass 1 (normal, 15× chr22): ~48 min, 55,811 germline variants
- Pass 2 (tumor, 15× chr22 with input_vcf=normal_golden): ~50 min, 56,318 variants (55,811 germline carried through + 507 somatic de-novo)
- Merge step: 55,811 `NEAT_ORIGIN=shared` + 507 `NEAT_ORIGIN=somatic`, 0 `germline`-only (every germline variant round-tripped, as expected when `input_vcf` is the normal pass's golden VCF)
- GATK Mutect2 + FilterMutectCalls scored via som.py against the somatic-filtered truth: **67.5% SNV recall, 96.7% precision; 12.1% indel recall, 66.7% precision.** Single-end alignment defaults; PE + panel-of-normals + contamination model would substantially lift indel numbers.
=======
5/27/2026
=========
## rneat v1.11.2

### Structural variant (SV) modeling improvements
- **Full support for symbolic `<INS>` (Insertions) in `SvModel`.** The statistical model now learns the rate and length distribution of symbolic insertions from training VCFs (gnomAD-SV, etc.). Sampled `<INS>` records now carry appropriate `SVLEN` and `END` tags in the output golden VCF.
- **Improved VCF output for structural variants.** De novo sampled SVs now include `SVLEN` in their INFO field (negative for deletions, positive for others).
- **VCF header completeness.** The output VCF header now includes missing `ALT` definitions for `<DUP>` and `<CNV>`, plus standard `INFO` definitions for `SVTYPE`, `SVLEN`, `END`, and `CN`.
- **Refined SV length statistics.** Introduced `SvData::event_length` to correctly distinguish between reference span (bases removed/replaced) and event length (bases inserted), ensuring more accurate log-normal fitting for insertions.
- **Bundled default SV model updated.** Updated `default_sv_model` with heuristic parameters for `<INS>` (median ~300bp, matching Alu MEIs) and `<INV>` (median ~4.9kb), and adjusted type probabilities to better reflect human genome diversity.
- **Fix:** `SvModel::affected_range_for_existing` now correctly handles `<INS>` and breakends as point events (1-bp reference span) rather than using their `SVLEN` as the span. This prevents point insertions from blocking nearby variants during de novo sampling.

5/25/2026
=========
## rneat v1.10.3

### Patch release: fix silent zero-SV output at human-chromosome scale (#194)
- **`SvModel::sample_variants` now produces de-novo SVs at chromosome-scale `λ`.** The Knuth multiplicative Poisson sampler at the core of `sample_variants` errored out when `exp(-λ)` underflowed to 0 — which happens around λ ≈ 745 in f64. With v1.10.2's full-genome `per_base_rate = 4.841e-4`, every human chromosome at `sv_rate_scale = 1.0` lands well past that threshold (chr22: λ ≈ 24,205; chr1: λ ≈ 121,000). The caller logged a `warn!` and emitted **zero** de-novo SVs — a silent regression introduced by v1.10.2's gnomAD-SV refit that nobody had hit because no test fixture used a contig larger than ~1.7 kb, and no one had run `gen-reads` with `sv_rate_scale > 0` against a real chromosome since v1.10.1 landed.
- **Fix:** `sample_poisson` is now hybrid. For `λ < 30` it keeps Knuth's exact algorithm (preserves the RNG-draw sequence that small-λ pipeline tests are seeded against). For `λ ≥ 30` it switches to a Gaussian approximation `N(λ, √λ)` sampled via inverse-CDF on a uniform draw — the same pattern `sample_log_normal_usize` uses for log-normal. At λ = 30 the approximation error is below 1 part in 1000; by λ = 1000 it's indistinguishable from exact Poisson in practice. No new dependencies; ~25 LoC change.
- **Regression tests added:** `sample_poisson_{large,extreme}_lambda_*` exercise λ from 1,000 to 121,000 explicitly. `sample_variants_works_at_chromosome_scale` is the end-to-end pin: hand-rolls a 50-Mb synthetic contig with the v1.10.2 default `per_base_rate` and asserts the returned SV count is well above zero (used to be silently zero before this fix).

5/25/2026
=========
## rneat v1.10.2

### Patch release: documentation, CI maintenance, known-limitations cleanup
- **Cancer simulator design doc.** Consolidated the plan that emerged from the cancer-modeling investigation (umbrella #129; cancer-MVP sub-tickets #183–#186; SV-gap tickets #187–#192) into a single reference at `docs/cancer_simulator.md`. Covers the orchestration-first architecture, GDC TCGA MAFs as the v1 training corpus, the SV gap analysis tiered by cancer relevance, and a staged roadmap. Establishes `docs/` as the convention for design notes that should survive past their originating issue threads.
- **CI: bumped `actions/checkout@v4` → `@v5` (#175).** Node.js 20 actions are deprecated by GitHub Actions — forced upgrade June 2026, removal September 2026. `actions/checkout@v5` runs on Node.js 24. Two workflow files updated: `.github/workflows/rust_binaries.yml` (release-binary build, runs on `v*.*.*` tags) and `.github/workflows/rusty-neat-tests.yml` (PR test workflow).

### Known limitations
Restated for v1.10.2 (carried over from earlier releases unless marked otherwise):
- **`<INV>` and `<BND>` records round-trip from input but do not yet drive read-orientation flipping or junction-read modeling.** Reads from inversion / breakend spans still come from the forward-strand reference. This is the single biggest gap for cancer simulation — tracked at **#187** (BND with junction reads) and **#188** (INV with read-strand flipping).
- **Multi-allelic VCF records (`,`-separated ALTs) are skipped at the reader,** both literal and symbolic. Real somatic VCFs do contain multi-allelic indels; the cancer MVP works around this by training from MAF-derived bi-allelic VCFs.
- **The SV sampler's overlap-rejection retry budget caps at 10 × n.** Very high `sv_rate_scale` values can saturate retries and warn-and-stop with fewer SVs than the Poisson draw requested. The full-genome gnomAD-derived defaults (`per_base_rate ≈ 4.8e-4`) at `sv_rate_scale = 1.0` are well below this ceiling; documented for users who push the rate substantially higher.
- **The bundled `default_sv_model()` is germline-derived (gnomAD-SV v4.1).** A cancer-specific SvModel default trained against PCAWG SV consensus is a natural follow-up once **#187** / **#188** land. Until then, cancer simulations should override the default via an explicit `mutation_model:` config field.

### Notes
- **Phantom de-novo SNPs in `<DEL>` spans is no longer a limitation.** v1.10.0's depth-modulation work added a `mult == 0.0` override that zeroes the mutation rate over deleted spans before SNP sampling, and the fix is exercised by `pipeline_e2e::gen_reads_with_symbolic_del_modulates_depth_and_round_trips_to_vcf`. Mentioned here because the v1.9.0 known-limitations section flagged it explicitly and readers chasing that note will otherwise have to read code to confirm the resolution.

5/24/2026
=========
## rneat v1.10.1

### Patch release: SV-model refinements from real-data validation
- **`<CNV>` records sampled from a model with empty `cnv_copy_number_distribution` now carry a fallback CN.** Trained models that observed `<CNV>` records but had no `INFO/CN` (gnomAD-SV is the canonical case — per-sample CN lives in FORMAT) used to produce `<CNV>` outputs with `copy_number = None`, which `gen-reads` warn-and-skipped for coverage modulation, so every CNV in the output VCF was uninterpretable for downstream callers. `SvModel::sample_variants` now falls back to a hard-coded `FALLBACK_CN_DISTRIBUTION` (mirroring the bundled-default CN shape) when the trained distribution is empty but `Cnv` is in the type list. An `info!` log at sampler entry flags when the fallback is firing.
- **`gen-mut-model` now ingests whole-genome VCFs without OOM.** `read_vcf` loaded every `Variant` with its full `info`, `id`, `filter`, `genotype_str`, `format`, `sample`, and `quality_score` fields populated — none of which the trainer reads. On a full-genome gnomAD-SV ingest (~1.5M records after the standard INV/CPX/CTX/BND filter), the `info` column alone (population AFs, VEP annotations, etc.) carried ~2–3 GB of unused string data and OOM'd the trainer with the desktop session as collateral damage. Added `read_vcf_lean` / `Variant::from_file_lean` which leave the unused fields empty for callers that never touch them. Peak RSS on the full-genome run dropped from 3.6 GB (OOM'd) to 1.3 GB; the full-fat `read_vcf` stays unchanged for `gen-reads` and `compare-vcfs`, both of which still need the round-trip fields. Also demoted the per-record `Found genotype in vcf file` log from `debug` to `trace` — under `--log-level debug` it was emitting one line per VCF record (44 MB on the OOM'd attempt).
- **Bundled `sv_model_defaults` refit from the full-genome gnomAD-SV v4.1 validation run.** The lean reader unblocked the `full`-mode validation the v1.10.1 PR had originally deferred. The fitted values supersede both the v1.10.0 literature defaults and the interim chr22-only fits committed earlier in this release cycle. Three of the five parameters are now genome-wide MLEs: `per_base_rate = 4.841e-4`, type proportions DEL/DUP/CNV = 0.8171/0.1824/0.0005, and per-type log-normal lengths (DEL median 703 bp, DUP median 864 bp, CNV median 15.3 kb). The `per_base_rate` change from 8.9e-6 (chr22 interim) to 4.8e-4 (full genome) is **not** sampling drift — the trainer normalizes by `total_reflen_seen` (the full FASTA), so the chr22 interim default effectively reported "fraction of all-genome SVs that happen to fall on chr22," not a true per-base rate. The full-genome refit corrects that. `homozygous_frequency` and `cnv_copy_number_distribution` stay at literature-derived values: gnomAD-SV is sites-only with no GT, and per-sample CN lives in FORMAT — neither parameter is fittable from this corpus. `tools/validate_with_gnomad.sh full` reproduces the run that generated these numbers.

### Behavioral change — read before upgrading
- **Default SV generation rate increased ~54× relative to v1.10.0.** The refit above raises `MutationModel::default()`'s `per_base_rate` from `3.0e-6` (v1.10.0 literature default) to `4.841e-4` (v1.10.1 full-genome MLE). Anyone with a saved gen-reads YAML that sets `sv_rate_scale > 0` *and* doesn't supply an explicit `mutation_model:` will see roughly **54× more `<DEL>` / `<DUP>` / `<CNV>` records per base** after upgrading. The previous default was meaningfully below the true gnomAD-SV rate and produced unrealistically sparse SV output — the new default is correct for genome-wide simulation, but the change is large enough that downstream pipelines tuned against v1.10.0 numbers may need recalibration.
- **Migration paths:** (a) leave `sv_rate_scale` unset / `0.0` to disable de novo SV generation entirely (this remains the unchanged default and is unaffected); (b) scale your existing `sv_rate_scale` down by ~50× to approximate v1.10.0's output volume; or (c) pin the v1.10.0 defaults by passing an explicit `mutation_model:` file trained against your own corpus.
- `<CNV>` records also drop from ~10% of all generated SVs (v1.10.0 hand-rolled) to ~0.05% (full-genome MLE) — anyone exercising CNV-specific downstream paths will see ~200× fewer of them. If you need denser CNVs for testing, fit a model from a CNV-enriched VCF and override the default.

5/24/2026
=========
## rneat v1.10.0

### De novo CNV / SV generation from a learned model (#105)
- **New `SvModel` on `MutationModel`.** Optional field carrying five learned distributions: `per_base_rate` (Poisson λ per base), `type_probabilities` over `<DEL>` / `<DUP>` / `<CNV>`, per-type log-normal `length_log_normal` parameters, `cnv_copy_number_distribution`, and `homozygous_frequency`. Marked `#[serde(default)]` so v1.9 mutation-model JSON files load unchanged with `sv_model = None`.
- **`gen-mut-model` fits the SV component.** Symbolic ALTs in the training VCF feed a parallel accumulator; after the SNP/indel pass completes, `SvModel::fit_from_observations` returns `Some(SvModel)` if any type cleared the 2-observation per-type fit bar, else `None`. Filter chain (each stage logs a count): non-symbolic ALTs ignored → `<INV>` / breakends / unknown tags dropped (v1 only generates `<DEL>` / `<DUP>` / `<CNV>`) → records with no derivable span dropped → spans below 50bp dropped (they're indels, not SVs) → types with < 2 surviving samples dropped from the model.
- **`gen-reads` samples and emits SVs.** Per-contig Poisson draw → weighted type pick → log-normal length (rejected outside `[MIN_SV_LENGTH_BP, contig_len / 4]`) → uniform anchor with overlap and N-gap rejection → `INFO/CN` draw for `<CNV>` → Bernoulli for genotype using the trained homozygous fraction. De novo records merge with `input_vcf` SVs and flow through the existing depth-modulation path, so simulated coverage and the golden VCF round-trip behave identically for both sources.
- **`sv_rate_scale` config knob, default `0.0`.** De novo SV generation is **opt-in** — v1.9 users see no behavior change. `1.0` reproduces the rate of the trained corpus; larger values scale proportionally for stress testing.
- **Bundled default `SvModel`.** `MutationModel::default()` (the embedded "no model file supplied" path) now attaches a hand-rolled SV model with approximate gnomAD-SV v2.1 parameters (DEL/DUP/CNV proportions, per-type log-normal length, CN distribution, homozygous fraction). Parameter provenance is documented inline in `common/src/models/sv_model_defaults.rs`. Users doing serious benchmarking should retrain via `gen-mut-model` against the SV-rich VCF of their choice; the default exists for tire-kicking.
- **`gen-mut-model` now handles SV-only training corpora.** v1.9.1's `NoLiteralVariants` guard rejected VCFs with zero literal SNPs/indels (e.g. the gnomAD-SV sites-only release) — surfaced during the v1.10 release validation against real gnomAD-SV data. v1.10 instead produces an SV-only model with `mutation_rate = 0` and a populated `sv_model`. gen-reads loaded with this model generates zero SNPs/indels and relies entirely on de novo SVs (via `sv_rate_scale > 0`) or `input_vcf` records. The `NoLiteralVariants` error variant is retained in the enum but no longer fires from this path.
- **Bgzip-compressed VCFs and FASTAs now read correctly.** `read_gzip_lines` previously used `flate2::GzDecoder`, which only reads the first gzip stream and silently stops. Bgzip files — produced by every standard bioinformatics tool (`bgzip`, `tabix`, `bcftools view -Oz`, `samtools faidx`) — are a series of concatenated gzip streams (one per ~64 KB block), so all but the first block were dropped. Practical effect: 100% of records lost from every bgzipped VCF (including gnomAD-SV, dbSNP, anything tabix-indexed). Switched to `MultiGzDecoder`, which is a strict superset — it handles both single-stream plain gzip and multi-stream bgzip identically. Bundled test fixtures (`H1N1_read1.fq.gz`, `H1N1_read2.fq.gz`) were bgzipped, so the `seq_error.canonical.json.gz` baseline was re-blessed (the trainer now sees the full FASTQ instead of just the first block).
- **Sites-only VCFs (no FORMAT / SAMPLE columns) now parse instead of being silently dropped.** Real-world callsets like gnomAD-SV ship without per-sample genotype columns. The reader used to skip every record because GT was absent; it now defaults to heterozygous (`0/1`) for sites-only records and records with FORMAT lacking GT. The trainer's `homozygous_frequency` will report `0` for a sites-only corpus — an honest "we don't know" rather than a synthesized 50/50 split.
- **Default log level lowered from `trace` to `info`.** The `--log-level` default was `trace`, which meant the file logger captured per-base debug events from `gen-reads` (sequencing-error generation, read placement, etc.) — on the order of `coverage × read_len × reference_bp` log lines. A single `gen-reads` run against hg38 at coverage 10 produced a multi-GB `.neat.log` and burned most of its wall-clock on disk writes. Users who actually need verbose output now opt in with `--log-level debug` or `--log-level trace`.

### Known limitations
- The bundled default model is **literature-derived**, not refit on the actual gnomAD VCFs. Users who care about distribution fidelity should retrain.
- `<INV>` and breakends still aren't generated — neither the read-content nor the junction modeling for them has landed.
- The sampler's overlap rejection retries up to `10 × n` times; very high `sv_rate_scale` values may saturate retries and warn-and-stop with fewer SVs than the Poisson draw requested.

5/24/2026
=========
## rneat v1.9.1

### Patch release on top of v1.9.0
- **Cargo workspace version catchup.** v1.9.0 shipped with `Cargo.toml`'s `workspace.version` still at `1.8.1` — anything inspecting `cargo` metadata or the compiled binary's reported version read the old value. Bumped to `1.9.1` so the package metadata matches the release.
- **Hardened `.as_literal().unwrap()` sites against future symbolic-ALT leaks (#169).** Added `debug_assert!(v.alternate.is_literal(), ...)` guards at the four call sites in `gen_mut_model::runner`, `fastq_tools::generate_read`, `compare_vcfs::equivalence::apply_variants`, and `compare_vcfs::runner::VariantKey::from`. Also added an explicit `is_symbolic()` skip at the top of `gen_mut_model::runner`'s variant loop — symbolic records were already falling into the `_` arm of the `variant_type` match, but the loop still incremented `homozygous_count` for them, biasing the trained model.
- **Fixed `SvData::span` silently accepting `INFO/END < POS` (#170).** Used to compute `end.saturating_sub(location).saturating_add(1)`, so a malformed `END` strictly less than `POS` silently produced `Some(1)`. The `<DEL>` path masked this via a downstream empty-range check, but `<DUP>` / `<CNV>` / `<INV>` got a bogus 1-base coverage modulation instead of warn-and-skip. `span()` now returns `None` for `end < location`; `Variant::from_file` also clears the bad `END` at parse time with a `warn!` naming both `POS` and `END`.
- **Fixed `gen-mut-model` NaN-poisoning on training VCFs with no literal variants.** A training VCF containing only symbolic SVs (zero SNPs / insertions / deletions) NaN-divided the SNP-frequency math (`0/0`), which then serialized to JSON as `null` and failed to deserialize on re-read. The runner now returns `GenMutationModelError::NoLiteralVariants` with a clear message and writes no output file. Real training corpora always carry literal variants alongside SVs, so this is an edge-case hardening rather than a fix for an observed user incident.
- **Comment and test-file cleanup.**
    - Scrubbed all remaining iterative-development "Phase N" references in source and test comments (most of them landed alongside #169, with one straggler in `compare_vcfs::runner` cleared up in this release).
    - Renamed `compare_vcfs_phase{1,2,3,4}.rs` integration test files to descriptive names (`compare_vcfs_exact_match.rs`, `compare_vcfs_equivalence.rs`, `compare_vcfs_attribution.rs`, `compare_vcfs_reasons_vcf.rs`) matching the doc comments already inside them.

5/24/2026
=========
## rneat v1.9.0

### Symbolic / structural variant support: input → depth modulation → output round-trip (#105, #106, #107)
- **New variant data model.** `AlternateType` is now an enum — `Literal(Vec<Nucleotide>)` for SNPs / indels / multi-base substitutions and `Symbolic(SvData)` for structural ALTs. `SvData` carries the verbatim ALT string, an `SvType` tag (`Del`, `Dup`, `Cnv`, `Ins`, `Inv`, `Bnd`, `Unknown`), and optional `INFO/END`, `INFO/SVLEN`, `INFO/CN`. Helper methods `as_literal` / `as_symbolic` / `is_literal` / `is_symbolic` let callers dispatch without unwrapping. Both VCF writers (`common::file_tools::vcf_tools::write_vcf` and `compare_vcfs::utils::vcf_writer`) emit `SvData.raw_alt` verbatim so the original symbolic notation survives the round-trip.
- **Input VCF parser handles symbolic ALTs.** `parse_alternate` classifies `<DEL>` / `<DUP>` / `<CNV>` / `<INS>` / `<INS:ME:ALU>` / `<INV>`, breakend notation (`G]17:198982]`, `]13:123456]T`), and other `<TAG>` forms (recorded as `SvType::Unknown` with the raw ALT preserved). `Variant::from_file` populates `SvData` from `INFO/END` / `INFO/SVLEN` / `INFO/CN` via a new `parse_sv_info` helper; `INFO/SVTYPE` disagreement vs the ALT tag warns but trusts the ALT (it's the more canonical signal). Malformed INFO values, bare flags, and `.` / empty INFO all parse without error.
- **Per-region depth modulation in `gen-reads`.** Symbolic SVs from `input_vcf` are no longer dropped; they drive per-region coverage scaling instead. `<DEL>` zeroes out reads in its span (hom → ×0, het → ×(ploidy−1)/ploidy), `<DUP>` inflates (hom → ×2, het → ×(ploidy+1)/ploidy), and any SV with `INFO/CN` is scaled by CN / ploidy. `<CNV>` records without `INFO/CN` warn and pass through unmodulated. `<INS>`, `<INV>`, breakends, and unknown tags pass through at ×1 — round-tripped to the output VCF but not yet biologically simulated (orientation flip / junction reads are out of scope for this release). Overlapping SVs compose multiplicatively. The fragment loop splits each region of interest by multiplier boundaries and calls `generate_fragments` / `generate_weighted_fragments` once per sub-region with scaled coverage.
- **`MutatedMap` separates symbolic SVs from per-base mutation flags.** New `sv_records: Vec<Variant>` lives alongside `variant_map`; `from_interval` routes any `AlternateType::Symbolic` variant there so its position is never flagged for per-base mutation. `write_vcf` emits `sv_records` after `variant_map` and preserves the original `INFO` field verbatim so `SVLEN` / `END` / `CN` survive the round-trip. `mutate_position` falls back to the reference if a symbolic record somehow reaches `variant_map`, so per-base read generation can't panic on `as_literal().unwrap()`.
- **`compare-vcfs` handles symbolic ALTs cleanly.** Previously `filter_vcf` / `VariantKey::from` would panic on `as_literal().unwrap()` if the golden VCF contained a round-tripped `<DEL>` / `<DUP>` / `<CNV>`. They're now counted into a new `skipped.symbolic` bucket (golden / called) and excluded from classification — `compare-vcfs` is byte-wise REF/ALT comparison and has no notion of equivalent symbolic edits. Report **schema version bumped to 1.3.0**; consumers parsing the JSON should default missing `symbolic` to 0.
- **End-to-end coverage.** New `pipeline_e2e.rs` test runs `gen-reads` through the real binary with a homozygous `<DEL>` at `H1N1_HA:500-1500` and asserts (a) the symbolic record appears verbatim in the output golden VCF with `INFO/END` preserved, and (b) zero reads start inside the deleted span while the surrounding region is still covered. Unit tests cover the new helpers (`build_coverage_multipliers`, `coverage_multiplier_for`, `apply_coverage_factor`, `split_region_by_multipliers`, `scale_coverage`) plus the `MutatedMap` routing change and the `compare-vcfs` symbolic-skip path.

### Known limitations
- `<INV>` and breakend records round-trip but do not yet drive orientation flipping or junction-read modeling — reads from inversion / breakend spans still come from the forward-strand reference.
- The mutation-rate map excludes the anchor base of each input variant from de-novo SNP sampling but does not yet exclude full SV spans. A homozygous `<DEL>` produces zero coverage in its span, so generated SNPs inside the deletion are invisible in the FASTQ, but they will appear in the output golden VCF.
- Multi-allelic SV sites (`,`-separated ALTs) are still skipped at the reader, as they are for literal multi-allelic records.
- NEAT does not yet *generate* CNVs / SVs de novo from a model — the v1.9.0 work is round-trip + depth modulation for user-supplied symbolic variants only. The mutation-model integration called out in #105 remains open.

5/23/2026
=========
## rneat v1.8.1

### Internal refactors (no user-visible behavior change)
- **Split `gen-gc-bias-model` `RunConfiguration` into walker params + `GcBiasModelParams` (#148).** `run_from_coverage` now takes only the FASTA-sweep / model-write fields it actually uses (`reference`, `bed_table`, `output_file`, `overwrite_output`, `window_size`, `window_stride`, `min_windows_per_bin`). The standalone CLI keeps `bam_file` and `min_mapq` for its own walk and embeds the model params as `config.model`. The unified `gen-bam-models` runner now constructs `GcBiasModelParams` directly from its `GcBiasSection` — no more dummy `bam_file` clone or hard-coded `min_mapq: 0`.
- **Factor `check_overwrite` into `common::file_tools::file_io` (#149).** Replaces four duplicated "output already exists and `overwrite_output` is false" checks in `gen_bam_models` (frag_length + gc_bias sub-sections), `gen_frag_length_model`, `gen_seq_error_model`, and `gen_gc_bias_model`. Standalone error messages are uniform; each module wraps the returned `io::Error` into its own error type.

5/22/2026
=========
## rneat v1.8.0

### `compare-vcfs`: new subcommand for validating a downstream caller against a NEAT golden VCF
- New `rneat compare-vcfs` subcommand. Classifies every variant in the called VCF as TP, FN, or FP relative to a NEAT-simulated golden VCF; produces `comparison_summary.json`, `comparison_summary.txt`, `FN_with_reasons.vcf`, and (optional) `FP.vcf`. Pure-Rust port of NEAT 2.1's `vcf_compare.py` algorithm plus NEAT 4's FN-attribution layer — no `hap.py` / pysam / Python subprocess required.
- **Exact-match classification.** Per contig, variants are keyed by `(position, ref, alt)`. Intersection → TP; golden-only → FN; called-only → FP. Per-VCF skip counters surface in the report so multi-allelic / hom-ref / non-PASS / outside-target-bed / outside-simulated-contigs records aren't silently lost.
- **Equivalence sweep** (ported from NEAT 2.1, default on). For each FP, take a ±`equivalence_window` bp window of the reference and apply both the FP set and FN set within it. If the resulting byte sequences match (case-insensitive), the two sets are denotation-different spellings of the same edit — every consumed FN is promoted to TP, every consumed FP is dropped. Catches left/right-aligned indels, compound SNPs vs single multi-base substitution, and other denotation differences exact matching would mis-classify. Disable with `fast: true`.
- **NEAT-aware FN attribution** (ported from NEAT 4). Every surviving FN gets one or more reason tags: `outside_simulated_contigs` (reported alone; BED checks are skipped since they presuppose simulation), `outside_mutation_bed`, `outside_target_bed`, `unknown`. Counts roll up into `fn_attribution` in the report; `FN_with_reasons.vcf` annotates each FN with a `NEAT_REASON` INFO tag.
- **Chrom aliases.** Optional two-column TSV (`bed_chrom\treference_chrom`) remaps BED chrom names at load time so users can compare against a reference using one convention (e.g., `chr1`) when their BED uses another (e.g., `1`). Applied to both `mutation_bed` and `target_bed`.
- **Chrom-naming-mismatch warnings.** Surfaced in `comparison_summary.warnings` when any BED chrom (post-alias) is absent from the reference's FASTA contig list. Catches single-typo BEDs (`chr_MP` in an otherwise correct H1N1 BED) that would otherwise silently misattribute every variant on that contig. Distinct from NEAT 4's "zero overlap" semantics, which would swallow the single-typo case.
- Report shape is schema-versioned at `1.2.0`. Outputs section enumerates every artifact path; precision / recall / F1 derived from TP/FN/FP counts (with `null` returned for undefined ratios). New `equivalents_promoted` counter in totals so users can see how many TPs were rescued by the sweep.
- See README → "Comparing a downstream caller's VCF against the golden VCF" for the end-to-end runbook; `template_config/compare_vcfs_template.yml` is the canonical config.

### Model-parity regression suite + release shakeout
- New `rneat/tests/model_parity.rs` integration suite. Builds each of the four data-driven models (`gen-mut-model`, `gen-seq-error-model`, `gen-frag-length-model`, `gen-gc-bias-model`) from fixed inputs and compares the gzipped JSON output against checked-in canonical baselines under `rneat/test_data/baseline_models/`. Also covers `gen-bam-models` GC-bias byte-equality vs the standalone runner (complements the in-tree fragment-length parity test which only covered frag-length).
- Cross-process stability is achieved by canonicalizing the parsed JSON before comparison: object keys are sorted via `serde_json::Map`'s BTreeMap backing, and arrays of tuple-key pairs (the shape serde produces for `HashMap<(K1,K2), V>` fields like `MutationModel::statistical_models`) are sorted by serialized key. Neutralizes Rust's per-process HashMap hasher randomization.
- Baselines are gzipped (~120 KB total). Re-bless with `BLESS_BASELINES=1 cargo test --test model_parity` or via `scripts/regenerate_model_baselines.sh` — only when an intentional algorithmic change is being made.
- Added `scripts/release_shakeout.sh` — full-size hg38 chr22 manual shakeout runbook for pre-publish QA. Not run in CI. Trains each model on real chr22 data, asserts byte-parity between the unified `gen-bam-models` and the standalone runners (one config per output, since the unified `min_mapq` gates both observers), runs `gen-reads` end-to-end with the trained models, and reports peak RSS + elapsed time per step for cross-release comparison.
- Added `scripts/RELEASE_CHECKLIST.md` — tickable gate-list for the publish workflow (lint, full test suite, shakeout, cross-version determinism, tag + GitHub release).
- Cross-version determinism verified for v1.7.0 → v1.8.0 on H1N1 paired-end with a fixed seed: byte-identical FASTQ multiset hashes across binaries built from each branch. The new `compare-vcfs` work doesn't touch any `gen-reads` code paths, so this remains true.

### Minor
- `cargo fmt --check` drift in `gen-bam-models` + `main.rs` (carried in from a merge resolution) is fixed; workspace formatting is once again clean.
- `gen-bam-models` template config / README explicitly call out that the unified `min_mapq` gates *both* observers — surfaced during the v1.7.0 chr22 shakeout when the standalone `gen-frag-length-model` (which hard-codes MAPQ > 10) couldn't be matched with a unified `min_mapq: 20` config.
- Added HPC friendly binaries to the auto workflow for releases.
- 
5/22/2026
=========
## rneat v1.7.0

### `gen-bam-models`: unified single-pass model builder
- New `rneat gen-bam-models` subcommand that walks a single BAM once and builds whichever subset of BAM-derived models the config requests. Currently supports `frag_length` and `gc_bias` sections; both are optional, at least one must be present. Saves a BAM iteration per additional model relative to running the per-tool commands separately.
- Config is two optional sub-sections under one top-level `bam_file` / `min_mapq`:
  ```yaml
  bam_file: /path/to/aligned.bam
  min_mapq: 0
  frag_length:
    output_file: /path/to/frag.json.gz
    min_reads: 2
  gc_bias:
    reference: /path/to/reference.fa
    output_file: /path/to/gc.json.gz
    window_size: 100
    window_stride: 100
    min_windows_per_bin: 10
  ```
- `FragLengthObserver` is now self-filtering — it applies its own paired/first-in-pair/mate-mapped/same-ref/MAPQ checks internally regardless of the walker filter, so it remains correct when paired with `BamWalkFilter::for_coverage()` in the unified runner. The standalone `gen-frag-length-model` command continues to pass `for_frag_length()` to the walker; the observer's internal checks are then redundant no-ops, so behavior on that path is unchanged.
- Internal refactor: `gen_frag_length_model::utils::runner::run_from_tlens` and `gen_gc_bias_model::utils::runner::run_from_coverage` are now public so the unified runner can call them with observer outputs from the shared BAM walk. The standalone runners are thin wrappers over the same helpers.
- A parity test asserts the unified runner and the standalone `gen-frag-length-model` produce numerically identical fragment-length models against the same BAM.

### Maintenance / cleanup
- Workspace-wide `cargo fmt` pass. `cargo fmt --check` now exits 0.
- `cargo clippy --workspace --all-targets` now passes with no deny-level errors. Two pre-existing hard failures fixed: `absurd_extreme_comparisons` in `common::rng::range_i64` (`if max > i64::MAX` is always false — replaced with `if min > max` so the function actually validates its range) and `unused_io_amount` in `filter_reads::filter_lib` tests (switched `write(...)` to `write_all(...)`).
- Idiomatic cleanups via `cargo clippy --fix` and a final manual pass: dropped redundant `.clone()` on Copy types, redundant `.into()` on identity conversions, redundant `&` on already-borrowed args; collapsed nested `if X { if Y { ... } }` into combined `if X && Y { ... }`; rewrote `(a + b - 1) / b` → `a.div_ceil(b)`; modernized `repeat(x).take(n)` → `repeat_n(x, n)`, `x >= a && x < b` → `(a..b).contains(&x)`, `x % 2 == 0` → `x.is_multiple_of(2)`; replaced three `assert!(true)` placeholders with `.expect(...)` calls; added a `Default` impl for `NucleotideSelector`; converted `impl Into<usize> for Nucleotide` to the canonical `impl From<Nucleotide> for usize` (same for `TrinucFrame`).
- CLI dispatch arms in `main.rs` (filter-reads, gen-mut-model, gen-seq-error-model, gen-frag-length-model, gen-gc-bias-model, gen-bam-models + the previously-touched gen-reads) now use a single `if let Some((name, cmd)) = subcommand && cmd.contains_id("configuration_yaml")` form instead of two nested ifs. Drive-by: fixed a copy-paste log message in `gen-mut-model` that said "Running rneat filter-reads".
- `gen-reads::run_neat` now takes `&RunConfiguration` instead of `&Box<RunConfiguration>`. Call sites no longer need to `Box::new(config.clone())`.
- `gen-reads`'s test VCF reader uses `.map_while(|l| l.ok())` instead of `.filter_map(|l| l.ok())` so a misbehaving reader can't loop forever on repeated `Err`.
- `gen-reads` config parser accepts `rng_seed` as either a YAML string (`"42 hello"`) or a bare integer (`42`); previously bare integers panicked with `ConfigReadError("rng_seed", "String")`. Two new tests pin both forms.
- `BedRecord` now has an `is_empty()` companion to `len()`.
- `SnpTrinucModel::from_file` no longer takes a dummy `&self` (it never used it).
- `Nucleotide::is_masked` rewritten with `matches!` instead of a boolean-returning `match`.
- Added explicit `#[allow(...)]` annotations for three cases where clippy's suggestion doesn't fit: `should_implement_trait` on the seven `default()` methods that return `Result` (can't implement infallible `std::Default`), `field_reassign_with_default` in the `gen_reads/utils/config.rs` test module (clearer than struct-init for sparse field overrides on a 30-field config), and `same_item_push` on `fastq_tools::generate_read` (pushing `'D'` N times is the CIGAR encoding, not a mistake).

5/21/2026
=========
## rneat v1.6.0

### `gen-gc-bias-model`: BAM input replaces pre-computed coverage files
- `gen-gc-bias-model` now reads an aligned BAM directly. Coverage is accumulated in process per reference base using `samtools depth`-equivalent defaults (skip unmapped/secondary/supplementary, configurable `min_mapq`, default 0). The external `samtools depth` / `bedtools genomecov` preprocessing step is no longer required.
- **Breaking change.** Config fields `coverage_file`, `coverage_format`, and the `CoverageFormat` enum (`samtools-depth`, `bedtools-genomecov-d`, `bedtools-genomecov-dz`) have been removed. Replace `coverage_file` with `bam_file` pointing at the aligned BAM the depth file was originally generated from. `min_mapq` is the new optional knob (default 0; set to 20 to mirror typical `samtools depth -q 20` filtering).
- Memory profile changes: peak memory is now ~4 bytes × total reference bases that received at least one record (~13 GB for hg38 WGS at 30×). The plain-text coverage-file path's per-contig-on-demand behavior is gone; if memory is a constraint on whole-genome runs, split the BAM by chromosome.
- The legacy `CoverageReader`, `CoverageData`, and `coverage_reader.rs` module (459 lines) have been deleted along with the samtools-depth / bedtools-genomecov line parsers.

### Shared BAM-walking infrastructure
- Introduced `common::file_tools::bam_reader::walk_bam` — a single-pass BAM iterator with a configurable `BamWalkFilter` and a `RecordObserver` visitor trait. Multiple observers receive each kept record in one pass, avoiding redundant BAM opens across consumers.
- Concrete observers: `FragLengthObserver` (TLEN histogram, used by `gen-frag-length-model`), `TransitionObserver` (4×4 read-vs-reference mismatch matrix from MD tags, used by `gen-seq-error-model`), and `CoverageObserver` (per-base reference depth, used by `gen-gc-bias-model`).
- `read_fragment_lengths` and `read_bam_transitions` are now thin wrappers over `walk_bam`; their signatures are unchanged so existing call sites continue to compile.
- `BamWalkFilter::for_frag_length()` / `for_transitions()` / `for_coverage()` centralize the per-tool filter policy that was previously duplicated across modules.
- An `on_start(&Header)` hook lets observers cache reference_sequence_id → contig name/length lookup tables once per walk.

### Migration
- In `gen_gc_bias_model.yml`, replace:
  ```yaml
  coverage_file: /path/to/coverage.txt
  coverage_format: samtools-depth
  ```
  with:
  ```yaml
  bam_file: /path/to/aligned.bam
  min_mapq: 20  # optional; default 0 matches `samtools depth` defaults
  ```
  No other config fields changed.

5/20/2026
=========
## rneat v1.5.1
- `gen-reads` now correctly handles reference genomes that contain IUPAC ambiguity codes (R/Y/M/K/S/W/H/B/V/D). Previously these bases were silently mapped to N, producing excess-N reads indistinguishable from assembly gaps. Each IUPAC code is now stochastically resolved to one of its constituent bases at reference-load time using the per-contig simulation RNG, so results are fully reproducible given the same seed. N in the reference retains its existing gap semantics and is unaffected. A single `WARN` line is emitted per contig listing the count of resolved bases. `gen-mut-model` and `gen-gc-bias-model` continue to treat IUPAC codes as N (appropriate since they analyze existing VCF/coverage data where variant callers typically skip ambiguous positions).
- Added `resolve_iupac_bases(raw, rng)` to `common::file_tools::fasta_stream`, available to any future consumer that needs RNG-seeded IUPAC resolution.
- `FastaStream` now yields raw sequence strings rather than `Vec<Nucleotide>`, allowing each caller to choose its own conversion strategy.
- Unit tests: all-code round-trip (output contains only ACGTN), ACGTN passthrough, masked-base passthrough, R two-way uniform distribution (1 000 draws), H three-way uniform distribution (3 000 draws).
- Integration test: `gen_reads_with_iupac_reference_produces_no_iupac_in_output` — runs the full pipeline against a synthetic FASTA with every IUPAC code and asserts no IUPAC letter appears in the output FASTQ.
- Added a binned scoring feature. Activated by user input when forming the quality score model, this will force rneat to bin the quality scores of the input data (if they are already binned, this is trivial) to user-defined bins. When this model is used in gen_reads, it creates a binned output.
- Added a number of unit tests and integration tests.

## rneat v1.5.0

### Memory improvements in `gen-reads`
- Replaced the per-position `Vec<f64>` bias map (one `f64` per chromosome base) with a compact `Vec<(start, end, rate)>` segment list. Measured peak RSS for rice (373 Mbp, 12 chromosomes) dropped from 3.06 GB → 1.05 GB (~3× reduction); major page faults dropped from 109 K → 0. Genomes previously too large to run safely (miscanthus ~11.9 GB estimated) now fit comfortably (~5.2 GB estimated).
- `cover_dataset` (single-ended path): fragment pool reduced from `num_frags` entries (`O(coverage × contig_length)`) to a single-entry pool with a cyclic index. Eliminates up to ~23 MB of VecDeque allocation per large chromosome and millions of redundant RNG shuffle calls.

### `gen-mut-model` fixes for real-world VCF data
- VCF parser no longer panics on missing allele calls (`./1`, `./.`) in the GT field. Missing alleles are skipped during genotype determination; variants with no allele data are treated as homozygous alt.
- Multi-allelic ALT fields (comma-separated, e.g. `T,TAATGGAATGG`) are now skipped with a debug log rather than being stored as garbled Complex variants.
- Variants whose FORMAT or SAMPLE fields are unparseable now warn-and-skip rather than aborting the model build, making the parser robust to mixed-format real-world VCFs.
- Soft-masked FASTA bases (`Maskeda/c/g/t`) are now canonicalized to uppercase (`A/C/G/T`) before comparing against the VCF REF allele in the trinucleotide context pass. This previously caused a hard `BaseMismatch` error on any variant in a repeat-annotated region of a soft-masked reference (e.g. hg38). Genuine mismatches (different base, not just masking) now warn-and-skip instead of aborting.
- Optional `bed_file` key in `gen-mut-model` config now handled safely; previously a missing key triggered a HashMap index panic.
- Per-variant "Found genotype in vcf file" log message demoted from `INFO` to `DEBUG`, eliminating millions of log lines when building models from whole-genome VCFs (NA12878 WGS produced ~4 M `INFO` lines before this fix).
- Human NA12878 whole-genome VCF (hg38) now builds successfully in ~2:15 (wall clock), peak RSS 5.6 GB, yielding `mutation_rate: 0.001457`, `homozygous_frequency: 0.390`, variant distribution 86% SNP / 6.6% INS / 7.4% DEL.

5/9/2026
=========

## rneat v1.4.1
- `gen-reads` new feature: `long_reads` configuration flag (boolean, default `false`). When `true`, fragments sampled from the fragment length model that are shorter than `read_len` are accepted rather than discarded, and the resulting read is truncated to the actual fragment length. This reproduces the realistic read-length distribution seen in real long-read datasets (PacBio HiFi, Oxford Nanopore), where size selection is imperfect and a tail of shorter reads is always present. When `false` (the default, appropriate for Illumina short reads), fragments shorter than `read_len` continue to be rejected — for Illumina, short fragments produce adapter read-through rather than truncated reads, which is outside the scope of this simulator. The flag affects both the uniform-coverage fragment path and the GC-bias-weighted path, as well as the FASTQ and BAM writers, which compute an effective read length of `min(fragment_len, read_len)` per fragment when `long_reads: true`. Added to `template_config/gen_reads_template.yml` with inline documentation.

5/9/2026
=========

## rneat v1.4.0
- `gen-reads` bug fix: `generate_read` in `fastq_tools.rs` was rejecting fragments whose length equals `read_length` (`<=` guard); all single-ended reads were silently discarded. Changed to `<` so fragments of exactly `read_length` are accepted.
- `gen-reads` bug fix: `map_buffer` in `fasta_stream.rs` was treating soft-masked bases (`Maskeda/c/g/t`) as N-regions, excluding repeat-annotated but valid sequence from read generation. Only `N` and `X` now delimit true gap regions.
- `gen-reads` bug fix: `generate_fragments` was drawing fragment lengths from the fragment model for single-ended runs, causing severe under-coverage when the model mean exceeds `read_length`. Single-ended fragments are now fixed at exactly `read_length`; the fragment model applies only to paired-ended runs.
- `gen-reads` bug fix: `cover_dataset` accumulated `start` across loop iterations instead of resetting to `temp_end`, causing reads to pile up at the beginning of each region rather than spanning it evenly.
- `gen-reads` / `mutation_model` bug fix: soft-masked nucleotides (`Maskeda/c/g/t`) were leaking into VCF REF/ALT fields and FASTQ read sequences as lowercase letters (`a/c/g/t`). `generate_mutation` now applies `check_base()` to the extracted `ref_base` and to each base of the deletion reference slice, unmasking all bases before they are written to output.
- New integration tests in `gen-reads`: `test_paired_ended_bam_flags_correct` (SAM flags on every R1/R2 record), `test_paired_ended_insert_size_matches_model` (TLEN mean within ±20% of configured fragment mean), `test_input_vcf_snp_appears_in_bam_reads` (seeded input-VCF SNP visible in BAM reads at the correct position).
- `gen-seq-error-model`: FASTQ input is now streamed one record at a time instead of being fully loaded into memory. Peak memory is bounded by a single record rather than the entire file, which matters for large (multi-GB) FASTQ inputs. The `max_reads` early-exit now also fires without reading the rest of the file.
- `gen-mut-model`: Trinucleotide probability computation is now O(n) instead of O(n²). Transition counts are pre-grouped by reference frame once before the probability loop, eliminating a redundant linear scan per trinucleotide context.
- `filter-reads` (`filter_fastq` and `filter_vcf`): Replaced `format!("{}\n", line)` heap allocations on every written line with two `write_all` calls (bytes + newline). Eliminated double HashMap lookups (`contains_key` + index) in favour of a single `get` call. Contig name reconstruction changed from a manual push loop to `join("_")`. The read-name prefix check simplified from `find` + index comparison to `starts_with`.
- `filter-reads` config: Removed a redundant identity `match` on a bool.
- `gen-frag-length-model`: `filter_lengths` no longer clones its input `Vec` before sorting — the owned value is sorted in place, eliminating one allocation proportional to the fragment count.
- `gen-gc-bias-model`: `overall_mean` computation replaced two separate filtered iterations with a single `fold`, scanning the 101-bin array once instead of twice.
- Streaming FASTA reader (`FastaStream`) moved from `gen-gc-bias-model` to the `common` crate, making it available to all subcommands. `gen-reads` now streams the reference one contig at a time instead of writing temp JSON block files, eliminating up to ~570 MB of disk I/O for human-scale genomes and capping peak memory to a single contig rather than the full genome.
- Contig-level parallelism added to `gen-reads` via rayon `par_bridge`. All contigs are now processed concurrently using rayon's work-stealing thread pool regardless of whether BAM output is requested; the output ordering is preserved by sorting results by contig index before writing.
- BAM creation is now fully parallel. Each contig worker writes its alignment records to a private temp BAM body file (no header, plain BGZF blocks). After all contigs finish, a single concatenation pass strips the BGZF EOF block from each intermediate file and assembles them in reference order into the final coordinate-sorted BAM. No `samtools sort` step is required.
- New `num_threads` config key in `gen-reads`. Set to an integer to cap the rayon thread pool used for parallel contig processing. Omit (or set to `.`) to use all available cores (default). Set to `1` to disable parallelism entirely.
- `NeatRng::derive_child(idx)` added to the common RNG: each parallel contig task derives a deterministic per-contig seed from the parent RNG state and the contig index, so results are fully reproducible regardless of thread scheduling.
- Removed `fasta_reader` and `fasta_map` modules. `SequenceBlock`, `SequenceMap`, and `RegionType` now live in `common::structs::sequence_block`; `map_buffer` and `apply_n_substitution` moved to `common::file_tools::fasta_stream`. Eliminates the temp JSON block files that were written to disk for every contig.
- `gen-mut-model` migrated to a single-pass `FastaStream` loop: trinucleotide counting and variant processing now share one iteration over each contig's sequence instead of two separate passes.

5/3/2026
=========

## rneat v1.3.0
- Added gen_mut_model mimicking the functionality in python Neat. The purpose is to allow users to create custom models based on real data.
- Added `gen-seq-error-model` subcommand. Reads FASTQ (.fastq or .fastq.gz) input and outputs a serialized `SequencingErrorModel` + `QualityScoreModel` (JSON, gzip-compressed). Supports optional BAM file input to infer a SNP transition matrix from real aligned reads using MD tags, and an optional TSV file to supply a fully custom transition matrix. Priority order: TSV override > BAM inference > default matrix.
- BAM-based transition matrix inference reads MD tags from aligned reads (via noodles), walking CIGAR + MD tokens to identify reference/alt base pairs at each SNP position and accumulating per-context counts into the 4×4 substitution matrix.
- Added FASTQ read shuffle in gen-reads: output reads are now randomly shuffled so alignment-tool ordering assumptions are not baked in. Simplified read naming to reduce verbosity.
- Added `template_config/gen_seq_error_model_template.yml` documenting all config fields with inline comments.
- Updated README with a full "Generating a Sequencing Error Model" section covering config schema, BAM MD-tag requirement, TSV format, and priority order.
- Bug fixes: degenerate zero-row handling in transition matrix normalization; edge-position SNP handling in MD-tag walker.
- Added `gen-frag-length-model` subcommand. Reads a paired-end BAM or SAM file and outputs a `FragmentLengthModel::Normal` (gzipped JSON) fitted to the observed template lengths. Filters using a median + 10×MAD outlier ceiling and a configurable per-length `min_reads` count floor (default 2; set 0 to disable). Added `template_config/gen_frag_length_model_template.yml` and README section.
- Implemented `produce_bam` / `output_bam` in `gen-reads` config. Setting `produce_bam: true` now writes a golden BAM file whose reads exactly match the generated FASTQs (same sequence with variants and sequencing errors applied, same quality scores). CIGAR strings reflect both variant positions and sequencing error indels: deletion errors contribute `D` ops and insertion errors contribute `I` ops, so the BAM alignment accurately represents the simulated read relative to the reference. The BAM is written in coordinate-sorted order using a windowed flush strategy (inspired by NEAT 2.1): records are buffered per block, flushed in sorted order at each block boundary, and any cross-block paired-end R2 reads are carried forward to the correct position. No `samtools sort` is needed; run `samtools index` directly on the output. Note: heterozygous variants use the same probabilistic application as the FASTQ, so reads drawn to the reference allele will not show the variant in either file.
- FASTQ shuffle is now global across all contigs rather than per-contig, so no chromosome ordering is visible in the output. All per-block temp files from every contig are collected into memory together and shuffled in a single pass before writing. For large genomes (> 500 Mbp) a runtime warning is emitted recommending `seqkit shuffle` as a memory-bounded post-processing alternative.
- Added `target_bed` config key to `gen-reads`. When supplied, reads and variants are generated only over regions in the BED; contigs absent from the BED are skipped entirely (no blocks read, no temp files written). This is the recommended approach for exome or targeted-panel runs on large genomes, where post-hoc `filter-reads` would otherwise process the full genome first. Both `.bed` and `.bed.gz` inputs are accepted.
- Added `input_vcf` config key to `gen-reads`. A single-sample VCF (`.vcf` or `.vcf.gz`) can now be supplied to force specific variants into the simulation. Each record must carry a `GT` field; SNPs and indels are applied at the given position before random variant generation runs on the remaining mutable bases. Complex variants (multi-base REF and ALT) are skipped with a warning. Random variants still fire at `mutation_rate` for positions not covered by the input VCF; set `mutation_rate: 0` to suppress them entirely. Caveats: only the first ALT allele of multi-allelic records is used; REF alleles are not verified against the reference sequence.
- Added `gen-gc-bias-model` subcommand. Reads a reference FASTA and a per-base coverage file (samtools depth, bedtools genomecov -d, or bedtools genomecov -dz), tiles the reference in sliding windows, and accumulates mean coverage per integer GC% bin to produce a `GcBiasModel` (gzipped JSON). Bins with fewer than `min_windows_per_bin` observations receive a neutral weight of 1.0 so sparse GC% values do not distort the model. Designed to scale to large genomes: coverage is indexed once and loaded one contig at a time; the reference is streamed without temp files; both GC counting and coverage summing use O(1)-amortized sliding windows so peak memory is bounded by the longest contig rather than total genome size. Added `template_config/gen_gc_bias_model.yml` documenting all config fields.
- Added GC bias simulation to `gen-reads`. When a `gc_bias_model` is supplied, fragment start positions are sampled directly from a per-position probability distribution derived from GC content, replacing the previous rejection-sampling approach. The model's `window_size` is embedded in the model file so the application-time window always matches the training-time window. Set `gc_bias_normalize_coverage: true` (default) to inflate fragment count to compensate for low-weight regions; set `false` to preserve the raw requested coverage.

9/19/2025
=========

## rneat v1.2.0
- Executable name change to `rneat`. In order to differentiate rusty-neat from the main neat, we are changing the binary name to `rneat` and we have expanded subcommands and made them more robust.
- Added bed filtering functions. Now if you supply a "filter_output" bed file to your gen-reads run, it will produce the original files plus a filtered version showing only reads and variants that overlap with the regions in the bed.
- It's also possible to run this utility separately with `rneat filter-reads`

9/13/2025
=========

## rust-neat v1.1.2
- I missed an intermediate version. Basically, it was converting to the command `neat gen-reads` instead of `generate-reads`, which will allow us to add subcommands, similar to how python NEAT worked. I also fixed some bugs in the fastq generation and paired-ended creation.
- In this release, we are adding a bed reader that can limit the fasta regions of interest. Basically, it will take an input bed and only generate reads over the sections in the bed.

9/11/2025
=========

## rusty-neat v1.1.0
- We have refactored the fastq code. Now a temporary fastq is written to file per block read in by the fasta reader. This should allow us to parallelize the temporary fastqs. Then it's simply a matter of concatenating them together. Works great! I got through single-ended E coli depth 10 in 36 seconds. I got through paired-ended E coli in 78 seconds.
- Caveats: the fastqs are sorted by read start. Included in each name is the exact location where the sequence was drawn from. Now, when you run an alignment, the reads will be named in such a way that will have <chromosome name>_start-pos_length in the name, so you can shuffle the reads and run an alignment, and easily check to see if your alignment is in the correct contig and location.

9/9/2025
=========

## rusty-neat v1.0.0
- It's been a couple years since I had an update to rusty-neat. This new version is version 1.0.0! It is not yet fully functional, as it is missing the BAM creation, though it does work on some test data: H1N1.fa and ecoli.fa. It's very fast for H1N1, but that's not too surprising. It's slow right now on fastq creation. I already know what I need to do to fix that part, I put an [issue](https://github.com/ncsa/rusty-neat/issues/65) up here that I am working on in develop. But for now, rusty-neat functions and can be executed on some arbitrary data. Will update the readme as well. Once I get the fastq generation idea to work (will require changes to fastq_tools and runner scripts), next step will be parallelizing the main contig loop in runner. That should give us some speedup. Also, we went from using way too much memory to using almost no memory. I'm not sure which is better. It would be nice to find a balance. I still think that a database could work. I'm wondering, maybe we index on position. Just use like a json-based database: {"contig": "chr1", "sequence": {"1":A, "2": C, etc}} We'd want something that had fast retrieval for writing. Something to investigate, maybe.
- The code is significantly different now, I read in all the old models from the original NEAT2.0, and stored them as default models. I used a lot of the same algorithms that they used in NEAT 2.0, figuring it was closer to rust than NEAT 3.X.
- rusty-neat writes out a ton of files. It should clean them up afterward, but I'm not sure it does yet, may need to investigate if there is an added command required. But it stores them in temp and linux cleans that up. Basically, we are splitting up the input fasta into smaller chunks (tbd on if there is an optimal size) then writing those to file for later retrival. Need to test various chunk sizes. I am basically doing blocked gzip-type things, but with some twists. If someone really understood bgzip format, it could probably be used directly. I'm not there yet.
- My runtimes running on a WSL2 Ubuntu setup with the latest rust installed: 311 milliseconds for H1N1.fa, 13 minutes for ecoli. It's all the IO during that step. I knew it was ugly as I was writing it but it's a rough draft. Obviously all the repeat code in fastq_tools is not ideal.
- The output files, so far have looked good. I realized I haven't tested paired-ended yet.
- I decided to map out the N-regions then overwrite them with random valid bases. The main reason being that it kept dipping into N regions and trying to use them for snp creation and that isn't going to work. NEAT 2.0 had a data structure that kept track of how long N-regions were, then would skip them if they were longer than a read length or overwrite them if they were short. That might work. Potential room for improvement, but it's fast so far.
- Eliminated bam creation for now. The option is in there, but it warns you that it doesn't work yet. Once fastqs are sorted and maybe after parallelism, it'll be next up.

9/11/2023
=========

## rusty-neat v0.1.0

Moved the code from a branch under NEAT (https://github.com/ncsa/neat) to its own repo https://github.com/ncsa/rusty-neat. This project will be considered a separate, related project. We will use the same underlying design ideas as NEAT, but it will essentially be different and won't be guaranteed in any way to be able to reproduce what the original NEAT did (though we aim to get it as close as possible, the differences in Python and Rust would make byte-by-byte compatibility very difficult). I have implemented an idea I prototyped in a private repo (https://github.com/joshfactorial/basic_neat), of reading in the file in what I hope is a straightforward way. I was looking for stuff like running out of memory and speed. It can read in the `E. coli` fasta file in about 30 seconds on my laptop, which is pretty fast, compared to Python. What's more is it can store the whole thing in memory, due to the efficient way Rust uses memory. Instead of copying the entire data set every time in memory, quickly eating it up, Rust forces you to decide when to pass a read-only version, copy the dataset, or use up the dataset once and lose it forever (all valid choices depending on context). So I'm feeling hopeful that this will work. I am publishing an Alpha called v0.1.0 for people to test, if they feel so inclined. See the README for more details.

## rusty-neat v0.1

Goal: To release a version of NEAT in rust capable of reading an arbitrary fasta file,
and outputting a version of that file with 1% of the bases mutated as SNPs.

9/5/2023
=========
- Initialized Cargo package with hello world default script.
- Added function to select a random integer from a list (for selecting positions)

9/9/2023
=========
- Created a first pass program that mutates a hardcoded DNA string at one random non-N position.
- Added a logging function to print out time-stamped log items both to stdout and to file.
- Need to improve error handling all around, taking advantage of Rusts error handling system.
- Design decision: After toying with some options to create a standardized RNG to use throughout the program, it seems the best option is to use an on-the-fly rng per function for now. This eliminates reproducibility between runs, but should enable parallelism of processing.
- Redesign from the ground up. Trying to use KISS as much as possible, and only introduce complexity as it becomes necessary, or as I learn Rust tricks that may help this endeavor.
