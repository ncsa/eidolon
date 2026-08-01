//! Integration test for symbolic <DUP> tandem chimeric reads (#220).
//!
//! Mirrors del_chimeric.rs. The tandem-DUP junction creates a reversed-
//! coordinate split-read signature: bases from the END of the duplicated
//! region butt up against bases from the START of the next copy. BWA
//! aligns the two halves of a junction-spanning fragment to "end of dup"
//! and "start of dup" in inverted order — Manta consumes this as a
//! tandem-DUP call.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

#[test]
fn gen_reads_with_symbolic_dup_emits_chimeric_junction_reads() {
    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_dup.vcf");
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
        // Homozygous symbolic DUP at H1N1_HA:400-799 (400 bp duplication).
        writeln!(
            f,
            "H1N1_HA\t400\t.\tG\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=799\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "dup_chimeric");
    config.coverage = 100;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_fastq = work.join("dup_chimeric_r1.fastq.gz");
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

    let dup_chimeric: Vec<&String> = fastq_qnames
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric_DUP_"))
        .collect();

    assert!(
        !dup_chimeric.is_empty(),
        "expected ≥1 EIDOLON_chimeric_DUP_ read in FASTQ for a homozygous \
         400bp DUP at 100× coverage; got 0 of {} total reads",
        fastq_qnames.len()
    );

    // QNAME shape pins POS=400 END=799 — but ONLY the tag and the echoed inputs. The
    // QNAME is formatted from the same `location`/`end` the test supplied, never from
    // the emitted geometry, so shifting either by ±1 leaves this assertion satisfied.
    // The junction bases are asserted in chimeric_sequence_content.rs.
    let qname_has_coords = dup_chimeric
        .iter()
        .any(|q| q.contains("EIDOLON_chimeric_DUP_H1N1_HA_400_799_"));
    assert!(
        qname_has_coords,
        "DUP chimeric QNAMEs must carry POS=400 END=799. Got: {:?}",
        dup_chimeric.iter().take(3).collect::<Vec<_>>()
    );

    assert!(
        dup_chimeric.len() >= 30,
        "expected ≳ 30 DUP chimeric reads at 100× hom coverage; got {}",
        dup_chimeric.len()
    );
}
