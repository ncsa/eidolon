# Comparing a downstream caller's VCF against the golden VCF
After running `gen-reads`, you can validate a downstream variant caller against the simulated truth using `eidolon compare-vcfs`. The subcommand classifies every variant into TP / FN / FP, catches denotation-different alternates that exact matching would mis-classify as both an FN and an FP (e.g., a left-aligned indel in the called VCF versus the same indel right-aligned in the golden), and attributes the surviving FNs to specific reasons drawn from the simulator's configuration.

```bash
$ eidolon compare-vcfs -c compare_vcfs_config.yml
```

Copy `template_config/compare_vcfs_template.yml` and fill in:

```yaml
golden_vcf: /path/to/golden.vcf.gz
called_vcf: /path/to/called.vcf.gz
reference:  /path/to/reference.fa
output_dir: /path/to/report_dir
overwrite_output: false

# Optional: restrict comparison to these regions.
target_bed: .
# Optional: comma-separated list. "." treats every called contig as simulated.
contigs_simulated: .
# Optional: BED used by the simulator. Drives FN attribution.
mutation_bed: .
# Optional: TSV remapping BED chrom names (e.g., 1<TAB>chr1).
chrom_aliases: .

# Filters
include_homs: false       # ALT == REF rows
include_filtered: false   # FILTER != PASS / "."

# Equivalence sweep (NEAT 2.1 algorithm, ported)
equivalence_window: 50    # ±N bp around each FP
fast: false               # skip the sweep entirely

# Outputs
write_fp_vcf: false       # also produce FP.vcf
```

The run produces four files under `output_dir`:

- `comparison_summary.json` — schema-versioned (currently `1.3.0`), machine-readable. Includes TP/FN/FP totals + per-contig breakdown, precision / recall / F1, the FN attribution roll-up (`outside_simulated_contigs`, `outside_mutation_bed`, `outside_target_bed`, `unknown`), per-VCF skip counters (`multiallelic`, `homozygous_ref`, `filtered`, `outside_target_bed`, `outside_simulated_contigs`, `symbolic`), and any chrom-naming-mismatch warnings. Symbolic / structural ALTs (`<DEL>`, `<DUP>`, `<CNV>`, ...) are byte-incomparable, so they're counted into `skipped.symbolic` and excluded from TP/FN/FP classification.
- `comparison_summary.txt` — same content, human-readable.
- `FN_with_reasons.vcf` — every surviving FN, annotated with a `EIDOLON_REASON` INFO tag listing the attribution reasons.
- `FP.vcf` (optional) — every surviving FP, as-is. Off by default; enable with `write_fp_vcf: true`.

**Equivalence detection.** For each false-positive variant, `compare-vcfs` takes a ±`equivalence_window` bp window of the reference and applies both the FP set and the FN set within that window. If the resulting byte sequences are identical, the two sets are alternative spellings of the same edit and every consumed FN is promoted to TP. Set `fast: true` to skip this pass; the report's `totals.equivalents_promoted` counts how many TPs were rescued by it.

**FN attribution.** Each surviving FN is tagged with at least one reason. If the FN's contig wasn't in `contigs_simulated`, that reason is reported alone (the BED checks are skipped because they presuppose the contig was simulated). Otherwise the configured BEDs are checked and `unknown` is the fallback. Aliases configured via `chrom_aliases` are applied to BED chrom names at load time, so you can compare against a reference that uses one convention (e.g., `chr1`) when your BED uses another (e.g., `1`). When a BED has any chrom that doesn't appear in the reference's FASTA contig list (post-alias), `compare-vcfs` surfaces a warning in the report — useful for catching single-typo BEDs that would otherwise misattribute every variant on that contig.
