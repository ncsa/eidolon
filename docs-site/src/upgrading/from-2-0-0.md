# Upgrading from 2.0.0

v3.0.0 renamed the tokens eidolon **emits**. There is no behavior change — the same
reads and variants are produced — but anything that *parses* eidolon output needs
attention. In-tree consumers were all migrated; your own scripts were not.

| Surface | v2.0.0 and earlier | v3.0.0+ | Migration |
|---|---|---|---|
| VCF INFO tags | `NEAT_ORIGIN`, `NEAT_PROVENANCE`, `NEAT_REASON`, `NEAT_CCF`, `NEAT_VAF` | `EIDOLON_*` | **converter provided** |
| VCF sample column | `NEAT_simulated_sample` | `EIDOLON_simulated_sample` | **converter provided** |
| FASTQ read names | `@RNEAT_generated_*`, `RNEAT_chimeric_*` | `@EIDOLON_generated_*`, `EIDOLON_chimeric_*` | **not converted** — see below |
| BAM read names (QNAME) | `RNEAT_generated_*` | `EIDOLON_generated_*` | **not converted** — see below |

## VCFs — convert them

```bash
tools/migrate_legacy_tokens.sh old_truth.vcf.gz new_truth.vcf.gz
tools/migrate_legacy_tokens.sh --check old_truth.vcf.gz   # report only; exit 10 = legacy
```

Idempotent: safe to run on a file that's already current, or on a non-eidolon VCF, so
you can call it unconditionally in a pipeline. It rewrites the header declarations and
the record fields together, so the result is a valid VCF.

Note you **cannot** work around this with a dual-name bcftools filter — bcftools rejects
an undefined tag in a `-i` expression at parse time, so even
`INFO/EIDOLON_ORIGIN="somatic" || INFO/NEAT_ORIGIN="somatic"` fails outright on a file
that carries only one of the two. Convert the file instead.

## Read names — NOT converted, and this is the part that can bite you

Read names are **not** rewritten, in either FASTQ or BAM. Rewriting files that are
routinely hundreds of GB to change ~19 bytes per record isn't a reasonable ask, so
instead `eidolon filter-reads` accepts **both** prefixes natively — legacy FASTQs keep
working with no conversion step, and it warns once when it sees an old prefix.

**Nothing protects your own scripts.** A `grep '^@RNEAT_generated_'`, `awk` field split,
or regex over read names will match **zero records and exit 0** against v3.0.0 output —
a silent empty result, not an error. Audit for the old prefix before upgrading:

```bash
grep -rn 'RNEAT_generated_\|RNEAT_chimeric_' your_scripts/
```

### Best fix: make your parser prefix-agnostic

**Only the prefix changed.** Everything after `_generated_` —
`<contig>_<start>_<end>_<uniq>/<mate>` — is byte-for-byte identical between versions
(verified across every read of a golden BAM re-prefixed both ways). So rather than
rewriting files, strip the prefix and your script works against *either* version, and
against any future rename:

```bash
# version-proof: drop whatever prefix is present, keep the encoded fields
samtools view x.bam | cut -f1 | sed -E 's/^[A-Z]+_generated_//'
zcat x_r1.fastq.gz  | awk 'NR%4==1' | sed -E 's/^@[A-Z]+_generated_//'
```

Same trick for chimeric reads: `sed -E 's/^@?[A-Z]+_chimeric_//'`.

If you are joining an old artifact against a new one on read name, normalize both sides
this way. One unrelated gotcha for QNAME joins: eidolon's golden BAM **keeps** the `/1`
`/2` mate suffix, while bwa strips it from aligned BAMs — so strip that too
(`sed 's|/[12]$||'`) when joining a golden BAM to an aligner's output. That mismatch
predates this release.

### If you must rewrite the files instead

```bash
# FASTQ (expensive — a full re-compress; prefer the prefix-agnostic parse above)
zcat old_r1.fastq.gz | sed 's/^@RNEAT_generated_/@EIDOLON_generated_/' | gzip > new_r1.fastq.gz

# BAM QNAMEs (QNAME is field 1, and no header line contains the prefix, so ^ is safe)
samtools view -h old.bam | sed 's/^RNEAT_generated_/EIDOLON_generated_/' | samtools view -b -o new.bam
```

BAM QNAMEs are covered by no in-tree tooling — `filter-reads` reads FASTQ, and the
converter is bcftools-based so it cannot touch a BAM. Nothing eidolon ships parses BAM
read names either, so this only affects your own scripts.
