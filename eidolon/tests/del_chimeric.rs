//! Integration test for symbolic <DEL> chimeric reads (#220).
//!
//! Before v1.13.0, gen-reads only modulated DEPTH for symbolic deletions —
//! the regular per-contig pass still emitted unbroken-reference reads
//! across the deletion span. This meant Manta, DELLY, GRIDSS and similar
//! split-read/PE-discordant SV callers got near-zero recall on the
//! simulated deletions even though the truth VCF correctly reported them.
//!
//! v1.13.0 adds a chimeric-junction path mirroring the BND/INV pattern:
//! `process_chimeric_variants` emits read pairs that stitch REF[..=POS]
//! (anchor included) with REF[END..] (post-deletion). When BWA aligns
//! these against the unbroken reference they surface as discordant-PE
//! pairs and split-read alignments — the signals Manta consumes.
//!
//! This test exercises a single forced <DEL> on H1N1_HA (a 2 kb segment)
//! at high coverage and pins:
//!   1. The FASTQ contains `EIDOLON_chimeric_DEL_*` reads.
//!   2. Their QNAME encodes the POS and END the test injected (so a
//!      future refactor that drops the DEL tag or shifts coordinates
//!      fails loudly).
//!   3. The chimeric read count is roughly proportional to coverage —
//!      catches a regression where the new branch silently emits zero.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

#[test]
fn gen_reads_with_symbolic_del_emits_chimeric_junction_reads() {
    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_del.vcf");
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
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">"
        )
        .unwrap();
        writeln!(
            f,
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">"
        )
        .unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
        )
        .unwrap();
        // Homozygous symbolic DEL at H1N1_HA:500-799 (300 bp deletion).
        // Hom so every fragment carries the junction → maximizes the
        // chimeric-read signal.
        writeln!(
            f,
            "H1N1_HA\t500\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=799\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "del_chimeric");
    // High coverage so the chimeric count is unambiguous.
    config.coverage = 100;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_fastq = work.join("del_chimeric_r1.fastq.gz");
    assert!(out_fastq.exists(), "FASTQ not produced at {:?}", out_fastq);

    let fastq_qnames: Vec<String> = {
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

    let del_chimeric: Vec<&String> = fastq_qnames
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_DEL_"))
        .collect();

    assert!(
        !del_chimeric.is_empty(),
        "expected ≥1 EIDOLON_chimeric_DEL_ read in FASTQ for a homozygous \
         300bp DEL at 100× coverage; got 0 of {} total reads. \
         Symbolic-DEL chimeric path may be silently emitting zero reads.",
        fastq_qnames.len()
    );

    // QNAME shape: EIDOLON_chimeric_DEL_<contig>_<pos>_<end>_<16-hex>/<mate>.
    // The injected DEL is at H1N1_HA POS=500 END=799, so QNAMEs should
    // carry "_500_799_" as a substring. This pins the coordinate encoding
    // so a future refactor that shifts position or end by ±1 fails the
    // test directly instead of silently moving the junction.
    let qname_has_coords = del_chimeric
        .iter()
        .any(|q| q.contains("EIDOLON_chimeric_DEL_H1N1_HA_500_799_"));
    assert!(
        qname_has_coords,
        "DEL chimeric QNAMEs must carry the injected POS=500 END=799 \
         coordinates. Got: {:?}",
        del_chimeric.iter().take(3).collect::<Vec<_>>()
    );

    // For homozygous DEL at 100× coverage we expect ~100 junction reads
    // (chimeric_pair emits one r1 per fragment_idx in 0..coverage). Allow
    // some slack for fragment-length retries and truncated-read skips.
    assert!(
        del_chimeric.len() >= 30,
        "expected ≳ 30 DEL chimeric reads at 100× hom coverage; got {}",
        del_chimeric.len()
    );
}
