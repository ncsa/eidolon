# Cancer simulation how-to

The following guide to `eidolon gen-cancer-reads` should give you some ideas on how to use eidolon to simula cancer reads 
and test your cancer pipelines. This is based on previous work done by NCSA and ICGC/ARGO and needs real testing to 
validate it works as expected. Design rationale and the mutation-rate / SV calibration details live in
[`cancer_simulator.md`](cancer_simulator.md).

## Pipeline

`gen-cancer-reads` runs two `gen-reads` passes over one reference and merges them:

```
  normal pass   cov = (1 − purity)·C    germline only
  tumor  pass   cov = purity·C          germline (shared) + somatic (de novo)
        │
        merge: tag reads N_/T_, concatenate; combine the two golden VCFs with origin tags
        ▼
  <prefix>_merged_r1.fastq.gz    → aligner + somatic caller
  <prefix>_merged_truth.vcf.gz   → INFO/NEAT_ORIGIN ∈ {germline, somatic, shared}
```

The tumor pass consumes the normal pass's golden VCF as its germline (`input_vcf`),
so both populations carry the same germline. No `bcftools`/`awk` runtime
dependency — the merge is native.

## Quick start

Smoke test on the bundled H1N1 reference (seconds; proves the build works):

```bash
cat > cancer.yml <<'YAML'
reference: eidolon/test_data/references/H1N1.fa
output_dir: ./cancer_out
output_prefix: smoketest
total_coverage: 30
purity: 0.6
read_len: 70
paired_ended: true
fragment_mean: 250
fragment_st_dev: 30
rng_seed: demo
overwrite_output: true
YAML

eidolon gen-cancer-reads -c cancer.yml
```

`template_config/gen_cancer_reads_template.yml` is the fully-commented config.

## Output files

| File | Contents |
|---|---|
| `<prefix>_merged_r1.fastq.gz` (+ `_r2`) | `N_`/`T_`-tagged, concatenated reads — the simulated biopsy. Feed to your aligner. |
| `<prefix>_merged_truth.vcf.gz` | Origin-tagged truth: `INFO/NEAT_ORIGIN ∈ {germline, somatic, shared}`. Score against this. |
| `<prefix>_normal.vcf.gz` | Germline-only truth. |
| `<prefix>_tumor.vcf.gz` | Germline + somatic truth. |
| `<prefix>_{normal,tumor}_r1.fastq.gz` | Per-pass reads (`keep_per_pass: false` deletes after merge). |

The `N_`/`T_` read-name tags prevent same-coordinate QNAME collisions between the
two passes (which MarkDuplicates would otherwise drop).

## Train your own model

`gen-cancer-reads` is model-driven: each pass takes a `.json.gz` mutation model
(`tumor_model:` / `normal_model:`). The bundled models are starting points — for
real work, train from your own somatic calls.

`eidolon gen-mut-model` fits a model from a reference plus a single-sample VCF:

```bash
cat > my_tumor.yml <<'YAML'
reference: /path/to/GRCh38.fa
vcf_file: /path/to/your_somatic_calls.vcf.gz
output_file: my_tumor_model.json.gz
overwrite_output: true
YAML
eidolon gen-mut-model -c my_tumor.yml
```

Input VCF expectations:

- single-sample, `GT` in `FORMAT`; contig names match the reference;
- SNPs and indels are fit by REF/ALT length class — multi-base REF **and** ALT
  (complex) records are skipped;
- symbolic SV records (`<DEL>` / `<DUP>` / `<CNV>` / `<INV>` / `<BND>` with
  `SVTYPE` / `END` / `SVLEN`) are fit into an `sv_model` when present, which is
  what `sv_rate_scale` draws from at simulation time;
- `bed_file:` restricts the fit to regions; `transition_matrix_file:` overrides the
  inferred SNP transition matrix. See `template_config/gen_mut_model_template.yml`.

The fitted `mutation_rate` is `variant_count / reference_length` — corpus-aggregated
if the VCF pools many tumors, so treat it as a spectrum descriptor, not a per-tumor
rate. Set per-tumor somatic burden at simulation time with `tumor_mutation_rate`
(see config knobs). Then:

```yaml
tumor_model: my_tumor_model.json.gz
# normal_model: my_germline_model.json.gz   # optional; default = built-in germline model
```

### Public-corpus adapters

If you don't have your own calls, these convert public corpora into a trainable VCF
and chain into `gen-mut-model` (`--train --reference`):

| Adapter | Corpus | Notes |
|---|---|---|
| `tools/fetch_cosmic_corpus.sh` | COSMIC GenomeScreensMutant | SNV + indel; academic login |
| `tools/fetch_tumor_corpus.sh` | TCGA MC3 PUBLIC | SNV only; open |
| `tools/fetch_cosmic_per_tissue_corpus.sh` + `tools/build_per_tissue_models.sh` | COSMIC, per `PRIMARY_SITE` | builds the per-tissue models below |

Pre-bundled, ready to use as `tumor_model:` without any download:
`tools/cosmic_per_tissue_{BRCA,skin,lung}.json.gz` (per-tissue SNP/indel + SV) and
`tools/cosmic_v104_pancancer_model.json.gz` (pan-cancer).

## Worked examples

Assume `eidolon` on `PATH` and a reference at `~/code/data/GRCh38.fa`.

### Tumor/normal with your own model

```bash
cat > run.yml <<'YAML'
reference: ~/code/data/GRCh38.fa
output_dir: ./out
output_prefix: tumor70
total_coverage: 60
purity: 0.7
read_len: 151
paired_ended: true
fragment_mean: 350
fragment_st_dev: 50
tumor_model: my_tumor_model.json.gz
tumor_mutation_rate: 1e-5
rng_seed: run1
overwrite_output: true
YAML
eidolon gen-cancer-reads -c run.yml
```

### With structural variants

Enable de novo SVs with `sv_rate_scale` (requires an `sv_model` in the tumor model
— present in your fit if the training VCF carried symbolic SVs, and in the bundled
per-tissue models):

```yaml
tumor_model: tools/cosmic_per_tissue_BRCA.json.gz
sv_rate_scale: 1.0       # 1.0 = the model's nominal rate; higher stress-tests SV callers
```

### Subclonal architecture (somatic VAF spectrum)

By default the somatic burden sits at ~one effective VAF, set by `purity` alone.
Real tumors are mixtures of subclones at distinct **cancer-cell fractions (CCF)**.
Add a `subclones:` list to spread de-novo somatic variants across those fractions —
each variant is assigned a subclone (share ∝ `weight`) and takes its `ccf`:

```yaml
purity: 0.8
subclones:
  - {ccf: 1.0, weight: 0.6}   # clonal / truncal — present in every tumor cell
  - {ccf: 0.4, weight: 0.3}   # major subclone
  - {ccf: 0.15, weight: 0.1}  # minor subclone
```

CCF is a **cellular-fraction factor** that composes with the variant's dosage and with
purity — it does not replace them:

```text
observed VAF = purity · dosage · CCF
```

`dosage` is the alt-copy fraction from the genotype (0.5 for a heterozygous SNV, 1.0 for
homozygous; `alt_copies/ploidy` in general). So a heterozygous somatic variant at CCF `f`
is observed at `~purity·f/2` — the value a subclonal-deconvolution tool (PyClone,
SciClone, …) inverts back to `f`. These are orthogonal axes (purity = normal
contamination, dosage = per-copy multiplicity, CCF = subclonal fraction). Germline
(shared) variants are unaffected. `ccf ∈ (0, 1]`; `weight` defaults to `1.0` (equal
share). Omit the block for the dosage-only default (output is byte-identical to
pre-subclone runs).

Each somatic record in the golden/merged-truth VCF carries two ground-truth INFO tags:

- **`NEAT_CCF`** — the intended cellular fraction (subclone CCF).
- **`NEAT_VAF`** — the intended **observed** VAF after tumor/normal mixing
  (`purity × dosage × CCF`; for a reproductive `somatic_vcf` replay, the original
  input VAF). This is the number a somatic caller measures on the merged reads.

Mind the difference from `FORMAT/AF`: that field is measured **per-pass (tumor-only)**,
so it reads `dosage × CCF` (≈ `NEAT_VAF ÷ purity`) and carries sampling noise.
For scoring a caller's VAF against ground truth, compare to `NEAT_VAF`; use `FORMAT/AF`
only if you want the tumor-cell fraction. Germline/shared records carry neither tag.

```bash
# planted observed VAF vs caller VAF: score against NEAT_VAF, not FORMAT/AF
bcftools view -H -i 'INFO/NEAT_ORIGIN="somatic"' merged_truth.vcf.gz \
  | awk -F'\t' '{match($8,/NEAT_VAF=[0-9.]+/); print $2, substr($8,RSTART+9)}'
```

### Building the architecture from real data

Instead of hand-authoring `subclones:`, point at a subclonal-deconvolution tool's
cluster table with `subclones_file:` (tab-separated, header required; mutually
exclusive with the inline list). eidolon folds the two shapes real tools emit into
`{ccf, weight}`:

| Tool | Shape | Columns used |
|---|---|---|
| PyClone / PyClone-VI | per-mutation | `cluster_id`, `cellular_prevalence` → grouped by cluster, weight = mutation count |
| PCAWG-11 / CSR / DPClust | cluster table | `cluster`, `ccf`, size (`n_ssms` / `size` / `weight`) → used directly |

```yaml
subclones_file: pyclone_clusters.tsv
```

Other tools convert with a one-liner — emit a header plus the two/three columns
above. CCF > 1.0 (noisy clonal clusters) is clamped to 1.0 and non-positive-CCF
rows are dropped, both warned with a count. Extra columns (`cellular_prevalence_std`,
`cluster_assignment_prob`, …) are ignored.

**Round-trip validation.** Ingest a real tumor's clusters → simulate → the golden
VCF's `NEAT_CCF` carries exactly those planted CCFs → run the deconvolution tool on
the *simulated* reads → confirm it recovers the architecture you fed in.
`scripts/delta/run_subclonal_vaf_validation.sh` runs the read-level half on real
data (Delta): it aligns the merged reads and checks the observed VAF at each somatic
site tracks `NEAT_VAF` (fidelity gated on bias + MAE-vs-noise-floor). On a soybean
scaffold at 200× with three subclones spanning 4–40% VAF, the aligned reads reproduce
the planted spectrum **unbiased** (mean err −0.003) and **to the sampling-noise floor**
(MAE 0.026, Pearson r 0.95) — see report §3.12.

### Reproductive replay (from a real somatic VCF)

The subclonal options above are *generative* — they invent somatic variants matching
a target architecture. To instead **replay a specific real tumor's somatic calls**,
point at a somatic VCF:

```yaml
purity: 0.6
tumor_mutation_rate: 0.0      # pure replay — no de-novo somatic on top
somatic_vcf: tumor_somatic.vcf.gz
```

Each variant is honored at its **observed VAF** (`INFO/AF`, else derived from
`FORMAT/AD`), divided by `purity` so the merged reads reproduce that VAF after
tumor/normal mixing (a raw Mutect2/Strelka VCF works directly). Replayed variants are
tagged `NEAT_ORIGIN=somatic` / `NEAT_PROVENANCE=somatic_input` in the truth — distinct
from germline even though both come from files. Notes:

- Composes with generation: leave `tumor_mutation_rate > 0` (and optional `subclones`)
  to layer de-novo somatic variants on top of the replayed set.
- A variant with no AF falls back to its genotype dosage (het ≈ 0.5, hom = 1.0).
- Observed VAF above `purity` is physically impossible for a somatic variant; it's
  clamped to a tumor-cell fraction of 1.0 with a warning.

### One germline, many tumor scenarios

Fix the germline once and sweep purity/depth by pointing each run at the same
`germline_vcf:` (any normal-pass golden, or a real germline VCF):

```bash
for p in 0.3 0.5 0.8; do
  sed "s/^purity:.*/purity: $p/; s/^output_prefix:.*/output_prefix: p$p/" run.yml \
    | sed "/^tumor_mutation_rate:/a germline_vcf: ./out/tumor70_normal.vcf.gz" \
    > run_$p.yml
  eidolon gen-cancer-reads -c run_$p.yml
done
```

## Benchmarking

Docker-based scoring pipelines, reads → scored caller in one command:

```bash
# SNV/indel: BWA-MEM → Mutect2 → som.py
tools/cancer_benchmark.sh \
    --reference ~/code/data/GRCh38.fa \
    --normal-fastq ./out/tumor70_normal_r1.fastq.gz \
    --tumor-fastq  ./out/tumor70_merged_r1.fastq.gz \
    --truth-vcf    ./out/tumor70_merged_truth.vcf.gz \
    --output-dir   ./bench

# SV: BWA-MEM → Manta → truvari
tools/cancer_sv_benchmark.sh \
    --reference ~/code/data/GRCh38.fa \
    --normal-fastq ./out/tumor70_normal_r1.fastq.gz \
    --tumor-fastq  ./out/tumor70_merged_r1.fastq.gz \
    --truth-vcf    ./out/tumor70_merged_truth.vcf.gz \
    --output-dir   ./sv_bench
```

`INFO/NEAT_ORIGIN` lets you filter the truth to somatic-only before scoring
(`--truth-filter`).

## Config knobs

| Key | Effect |
|---|---|
| `purity` | tumor cell fraction in (0,1); tumor pass = `purity·total_coverage`. Drives somatic VAF. |
| `total_coverage` | combined merged depth. Keep high enough that `purity·C` doesn't round to a useless depth. |
| `tumor_mutation_rate` | per-base somatic rate. Default `1e-5`. `model` = use the model's fitted rate. |
| `normal_mutation_rate` | per-base germline rate. Default = the model's fitted rate. |
| `sv_rate_scale` | de novo SV multiplier; `0` = off, `1.0` = the model's `sv_model` rate. |
| `subclones` | optional list of `{ccf, weight}` subclones; spreads somatic variants across CCFs → observed VAF = `purity · dosage · ccf`. Omit for the dosage-only default. |
| `germline_vcf` | fixed shared germline instead of de-novo generation. |
| `rng_seed` | seeds both passes (suffixed `-normal`/`-tumor`); printed to the log. |
| `keep_per_pass` | keep per-pass FASTQs (`false` = merged only). |

## Relationship to `tools/cancer_simulate.sh`

The native subcommand is a drop-in replacement for the original shell orchestrator,
verified equivalent by `eidolon/tests/cancer_parity.rs` (identical merged-FASTQ record
multisets, per-pass golden VCFs, and origin classifications). The script is retained
for the Docker benchmark wiring; prefer the subcommand for new work.
