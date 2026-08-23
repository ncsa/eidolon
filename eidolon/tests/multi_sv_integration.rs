//! Combined-SV integration test for v1.12.0.
//!
//! The per-feature integration tests in this directory cover each SV gap in
//! isolation: `bnd_fastq.rs` and `bnd_roundtrip.rs` for BND (#187),
//! `inv_fastq.rs` for INV (#188). #190 lands as a unit test in `sv_model.rs`
//! because de novo INS reuses the existing literal-Insertion machinery
//! end-to-end (no new code paths in `gen_reads`).
//!
//! This test exercises all three SV gaps in a single `gen-reads` run to
//! catch cross-type interaction bugs that the per-type tests can't surface:
//! e.g. the `processed_ids: HashSet<String>` is shared across BND and INV
//! dedup; the chimeric-pair writer's FASTQ/BAM tail is shared; the
//! sample_variants overlap-rejection loop has to interleave INS (literal,
//! routed to `variant_map`) with BND/INV/DEL (symbolic, routed to
//! `sv_records`).

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

#[test]
fn gen_reads_emits_bnd_inv_and_de_novo_ins_in_one_run() {
    let (_dir, work) = fresh_workdir();

    // Input VCF carries one BND pair (intra-contig) + one INV, so we can pin
    // round-trip behavior. The de novo INS records come from sample_variants
    // amplified by sv_rate_scale.
    let input_vcf = work.join("input_multi_sv.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
        )
        .unwrap();
        // BND pair: H1N1_HA:500 ↔ H1N1_HA:1500 (intra-contig translocation).
        writeln!(
            f,
            "H1N1_HA\t500\tbnd_A\tG\tG]H1N1_HA:1500]\t60\tPASS\tSVTYPE=BND;MATEID=bnd_B\tGT\t0/1"
        )
        .unwrap();
        writeln!(
            f,
            "H1N1_HA\t1500\tbnd_B\tC\t[H1N1_HA:500[C\t60\tPASS\tSVTYPE=BND;MATEID=bnd_A\tGT\t0/1"
        )
        .unwrap();
        // INV: H1N1_HA:200-400 (small inversion well inside the segment).
        writeln!(
            f,
            "H1N1_HA\t200\t.\tA\t<INV>\t60\tPASS\tSVTYPE=INV;END=400\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "multi_sv");
    config.coverage = 30;
    config.produce_fastq = true;
    config.produce_vcf = true;
    config.input_vcf = Some(input_vcf);
    // High sv_rate_scale so de novo INS gets sampled at small-contig scale
    // (default rate × ~2 kb contig × 1.0 is borderline). The model still
    // mostly samples DEL/DUP/CNV at typical rates; we just want a guaranteed
    // de novo INS record to verify the literal-Insertion path.
    config.sv_rate_scale = Some(50.0);
    // #516 HISTORY: de novo insertions longer than `read_len - 1` were once REFUSED rather
    // than planted, because reads could carry at most that much of an insertion's novel
    // sequence. That cap is GONE — insertions are sampled in altered-haplotype coordinates
    // now — and the assertion further down pins its removal. Do not reintroduce the old
    // wording here; it contradicted that assertion in this very file for one commit.
    //
    // The read length stays at 150 for a different reason: it is the largest that still fits
    // the helper's hardcoded 200 bp fragment.
    config.read_len = 150;
    config.rng_seed = "v1.12.0-multi-sv".to_string();

    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_fastq = work.join("multi_sv_r1.fastq.gz");
    let out_vcf = work.join("multi_sv.vcf.gz");
    assert!(out_fastq.exists(), "FASTQ not produced at {:?}", out_fastq);
    assert!(out_vcf.exists(), "VCF not produced at {:?}", out_vcf);

    // ── VCF assertions ─────────────────────────────────────────────────────
    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(&out_vcf).unwrap()));
        r.lines().map(|l| l.unwrap()).collect()
    };
    let body: Vec<&String> = vcf_lines.iter().filter(|l| !l.starts_with('#')).collect();
    assert!(!body.is_empty(), "no records in output VCF");

    // The synthetic "chimeric" contig that `process_chimeric_variants`
    // uses as a control-flow tag for organizing its FASTQ/BAM outputs
    // must NOT leak into the VCF's contig declarations. A bogus
    // `##contig=<ID=chimeric,length=0>` line breaks strict downstream
    // parsers (truvari refuses to score against it).
    assert!(
        !vcf_lines
            .iter()
            .any(|l| l.contains("##contig=<ID=chimeric")),
        "##contig=<ID=chimeric,...> leaked into output VCF — \
         see collect_contig_result's pseudo-contig skip. Got: {:?}",
        vcf_lines
            .iter()
            .filter(|l| l.starts_with("##contig"))
            .collect::<Vec<_>>()
    );

    // Round-trip: both BND records from input should appear in output.
    assert!(
        body.iter().any(|l| l.contains("G]H1N1_HA:1500]")),
        "first BND record didn't round-trip; got: {body:?}"
    );
    assert!(
        body.iter().any(|l| l.contains("[H1N1_HA:500[C")),
        "second BND record didn't round-trip; got: {body:?}"
    );
    // Round-trip: the input INV record should appear in output.
    assert!(
        body.iter().any(|l| l.contains("<INV>")),
        "INV record didn't round-trip; got: {body:?}"
    );

    // De novo INS records — emitted as literal Insertion (REF=1 base,
    // ALT >1 base). The bundled SvModel's INS length distribution has
    // mu=5.7, sigma=1.0, so the inserted bases land in the ~50-1000 bp
    // range. Filter for ALT > 30 to exclude single-base SNPs/MNPs and
    // short sequencing errors (which won't appear here anyway, but
    // defense in depth).
    let de_novo_ins_count = body
        .iter()
        .filter(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            fields.len() >= 5
                && fields[3].len() == 1
                && fields[3]
                    .chars()
                    .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
                && fields[4].len() > 30
                && fields[4]
                    .chars()
                    .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
        })
        .count();
    assert!(
        de_novo_ins_count > 0,
        "expected ≥1 de novo literal Insertion record (per #190); got 0 in: {body:?}"
    );

    // #516 END-TO-END, now INVERTED. This assertion used to pin the opposite — that no
    // de novo insertion could declare more novel bases than `read_len - 1`, because
    // fragments were placed in reference offsets and so could never begin inside a
    // zero-reference-width event. The interim cap that enforced that has been removed:
    // long insertions are sampled in altered-haplotype coordinates now, where the
    // inserted sequence has width, so there is no length beyond which the reads cannot
    // support the declared SVLEN.
    //
    // Keeping an assertion here (rather than deleting one) matters for the same reason
    // the original gave: it pins that the de novo path actually reaches the new
    // behaviour. If the cap were reinstated — or the sampler stopped being wired in —
    // every other assertion in this test would stay green while long de novo insertions
    // silently vanished from the truth VCF again.
    //
    // This pins PLANTING only. Whether the inserted sequence reaches the reads is a
    // separate question, pinned in eidolon/tests/long_insertion_rework.rs.
    let max_novel = config.read_len - 1;
    let oversized: Vec<(usize, &&String)> = body
        .iter()
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 5
                || f[3].len() != 1
                || !f[4].chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
            {
                return None;
            }
            let novel = f[4].len().saturating_sub(f[3].len());
            (novel > max_novel).then_some((novel, l))
        })
        .collect();
    assert!(
        !oversized.is_empty(),
        "no de novo insertion exceeds read_len-1={max_novel}. At this SV rate the model \
         draws them well past a read length, so seeing none means they are being dropped \
         — the #516 cap is back, or the de novo path is not reaching haplotype sampling."
    );

    // ── FASTQ assertions ───────────────────────────────────────────────────
    let fastq_names: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_fastq).unwrap(),
        ));
        r.lines()
            .enumerate()
            .filter(|(i, _)| i % 4 == 0)
            .map(|(_, l)| l.unwrap())
            .collect()
    };

    // BND chimeric reads: QNAME format is
    //   EIDOLON_chimeric_<c1>_<pos>_<c2>_<mate>_<16-hex>/<mate-id>
    // (NOT prefixed with INV — that's INV's format).
    let bnd_chimeric = fastq_names
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_") && !l.contains("EIDOLON_chimeric_INV_"))
        .count();
    assert!(
        bnd_chimeric > 0,
        "expected ≥1 BND chimeric read in FASTQ; got 0 in {} total reads",
        fastq_names.len()
    );

    // INV chimeric reads: QNAME starts with EIDOLON_chimeric_INV_. Both
    // junctions (start `_1_` and end `_2_`) should be represented since the
    // INV is homozygous in input (1/1 → full coverage on each junction).
    let inv_chimeric = fastq_names
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_INV_"))
        .count();
    assert!(
        inv_chimeric > 0,
        "expected ≥1 INV chimeric read in FASTQ; got 0"
    );
    let has_junction_1 = fastq_names
        .iter()
        .any(|l| l.contains("EIDOLON_chimeric_INV_") && l.contains("_1_0"));
    let has_junction_2 = fastq_names
        .iter()
        .any(|l| l.contains("EIDOLON_chimeric_INV_") && l.contains("_2_0"));
    assert!(has_junction_1, "missing INV junction-1 reads");
    assert!(has_junction_2, "missing INV junction-2 reads");

    // De novo INS bases must actually reach the FASTQ — not just appear in
    // the truth VCF. Grab the first 20 bases of one de novo INS's inserted
    // sequence and confirm it shows up in at least one read. The literal-
    // Insertion path emits the novel bases via the existing fragment-loop
    // coin-flip; if that wiring breaks the bases would be missing.
    let first_ins_kmer: Option<String> = body
        .iter()
        .filter_map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 5
                && fields[3].len() == 1
                && fields[4].len() > 30
                && fields[4]
                    .chars()
                    .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
            {
                // Skip the anchor base, take the next 20.
                Some(fields[4][1..21].to_string())
            } else {
                None
            }
        })
        .next();
    if let Some(kmer) = first_ins_kmer {
        // Read sequence lines (line index % 4 == 1 within the FASTQ).
        let read_seqs: Vec<String> = {
            use flate2::read::MultiGzDecoder;
            use std::io::{BufRead, BufReader};
            let r = BufReader::new(MultiGzDecoder::new(
                std::fs::File::open(&out_fastq).unwrap(),
            ));
            r.lines()
                .enumerate()
                .filter(|(i, _)| i % 4 == 1)
                .map(|(_, l)| l.unwrap())
                .collect()
        };
        let hits = read_seqs.iter().filter(|s| s.contains(&kmer)).count();
        assert!(
            hits > 0,
            "first de novo INS's 20-mer prefix '{kmer}' didn't appear in any read — \
             the literal-INS bases aren't reaching the FASTQ"
        );
    }
}
