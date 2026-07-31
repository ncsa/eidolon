//! Integration test for symbolic <CNV> chimeric reads (#220).
//!
//! CNV chimeric reads dispatch to the DEL or DUP path based on whether
//! INFO/CN is below or above the diploid baseline. This test exercises
//! both directions:
//!   - CN=0 (homozygous loss) → DEL-like EIDOLON_chimeric_DEL_ reads
//!   - CN=4 (amplification, +2 vs diploid) → DUP-like EIDOLON_chimeric_DUP_ reads
//!
//! Note: the QNAMEs reuse the DEL/DUP tags rather than introducing a
//! separate `CNV` tag because the underlying junction shape is identical
//! to a DEL or DUP — only the bookkeeping (which truth record gets
//! matched) differs, and that lives at the VCF level, not the read level.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

fn run_cnv(test_name: &str, cn: u32) -> Vec<String> {
    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_cnv.vcf");
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
            "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">"
        )
        .unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
        )
        .unwrap();
        // 400 bp CNV at H1N1_HA:400-799 with the supplied CN.
        writeln!(
            f,
            "H1N1_HA\t400\t.\tG\t<CNV>\t60\tPASS\tSVTYPE=CNV;END=799;CN={cn}\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), test_name);
    config.coverage = 100;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_fastq = work.join(format!("{test_name}_r1.fastq.gz"));
    assert!(out_fastq.exists(), "FASTQ not produced at {:?}", out_fastq);

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
}

#[test]
fn cnv_with_cn_below_ploidy_emits_del_like_chimeric() {
    let qnames = run_cnv("cnv_loss", 0);
    let chimeric: Vec<&String> = qnames
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_DEL_H1N1_HA_400_799_"))
        .collect();
    assert!(
        !chimeric.is_empty(),
        "CN=0 CNV must dispatch to DEL-like chimeric path; got 0 such reads in {} total",
        qnames.len()
    );
    // Ensure NO DUP-tagged reads appear (CN < ploidy is a loss, not a gain).
    let wrong: Vec<&String> = qnames
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_DUP_H1N1_HA_400_799_"))
        .collect();
    assert!(
        wrong.is_empty(),
        "CN=0 must NOT produce DUP-tagged chimeric reads; got {} of those",
        wrong.len()
    );
}

#[test]
fn cnv_with_cn_above_ploidy_emits_dup_like_chimeric() {
    let qnames = run_cnv("cnv_gain", 4);
    let chimeric: Vec<&String> = qnames
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_DUP_H1N1_HA_400_799_"))
        .collect();
    assert!(
        !chimeric.is_empty(),
        "CN=4 CNV (gain, +2 vs diploid) must dispatch to DUP-like chimeric path; \
         got 0 such reads in {} total",
        qnames.len()
    );
    // No DEL-tagged reads either.
    let wrong: Vec<&String> = qnames
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_DEL_H1N1_HA_400_799_"))
        .collect();
    assert!(
        wrong.is_empty(),
        "CN=4 must NOT produce DEL-tagged chimeric reads; got {} of those",
        wrong.len()
    );
}
