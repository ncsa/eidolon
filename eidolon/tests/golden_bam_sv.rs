//! Focused Gate 2b checks for golden BAMs with planted SVs.
//!
//! The variant-free gate proves that the BAM and FASTQ writers agree, but it cannot prove that
//! an SV survives into the BAM geometry. These small, deterministic cases require the expected
//! CIGAR signature for a literal INS and require chimeric records for DEL/DUP/INV/BND. Structural
//! DEL/DUP/INV records are stitched junction reads (their CIGAR is aligned to each anchor), not a
//! single `D`/`N` operation against the unmodified reference.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, gate2::ref_base, h1n1_reference};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;

struct Case {
    tag: &'static str,
    record: String,
    signature: Signature,
}

enum Signature {
    Insertion(usize),
    Chimeric,
}

fn cases() -> Vec<Case> {
    let base = ref_base("H1N1_HA", 500);
    vec![
        Case {
            tag: "del",
            record: "H1N1_HA\t500\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=550\tGT\t1/1".into(),
            signature: Signature::Chimeric,
        },
        Case {
            tag: "ins",
            record: format!(
                "H1N1_HA\t900\t.\t{base}\t{base}ACGTACGTAC\t60\tPASS\tSVTYPE=INS\tGT\t1/1"
            ),
            signature: Signature::Insertion(10),
        },
        Case {
            tag: "dup",
            record: "H1N1_HA\t600\t.\tG\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=900\tGT\t1/1".into(),
            signature: Signature::Chimeric,
        },
        Case {
            tag: "inv",
            record: "H1N1_HA\t600\t.\tG\t<INV>\t60\tPASS\tSVTYPE=INV;END=900\tGT\t1/1".into(),
            signature: Signature::Chimeric,
        },
        Case {
            tag: "bnd",
            record: format!(
                "H1N1_HA\t500\tbnd\t{base}\t{base}[H1N1_HA:1500[\t60\tPASS\tSVTYPE=BND\tGT\t1/1"
            ),
            signature: Signature::Chimeric,
        },
    ]
}

fn write_vcf(path: &Path, record: &str) {
    let mut f = std::fs::File::create(path).unwrap();
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
    writeln!(f, "{record}").unwrap();
}

fn run_case(case: &Case) -> Vec<(String, Vec<(Kind, usize)>)> {
    let (_dir, work) = fresh_workdir();
    let input = work.join(format!("{}.vcf", case.tag));
    write_vcf(&input, &case.record);

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), case.tag);
    config.coverage = 30;
    config.read_len = 100;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.produce_bam = true;
    config.produce_vcf = false;
    config.input_vcf = Some(input);
    config.mutation_rate = Some(0.0);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let bam_path = work.join(format!("{}.bam", case.tag));
    assert!(
        bam_path.is_file(),
        "{}: golden BAM was not produced",
        case.tag
    );
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    reader.read_header().unwrap();
    let mut records = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        let qname = String::from_utf8_lossy(record.name().unwrap().as_ref()).to_string();
        let cigar = record
            .cigar()
            .iter()
            .map(|op| {
                let op = op.unwrap();
                (op.kind(), op.len())
            })
            .collect();
        records.push((qname, cigar));
    }
    records
}

#[test]
fn planted_sv_signatures_reach_golden_bam() {
    for case in cases() {
        let records = run_case(&case);
        assert!(!records.is_empty(), "{}: BAM had no records", case.tag);

        let mut mates: HashMap<&str, usize> = HashMap::new();
        for (qname, cigar) in &records {
            assert!(
                !qname.ends_with("/1") && !qname.ends_with("/2"),
                "{}: QNAME suffix leaked",
                qname
            );
            assert!(
                !cigar.is_empty(),
                "{}: {qname} has an empty CIGAR",
                case.tag
            );
            *mates.entry(qname).or_default() += 1;
        }
        assert!(
            mates.values().all(|&n| n == 2),
            "{}: every paired QNAME must have two BAM records",
            case.tag
        );

        match case.signature {
            Signature::Insertion(expected) => assert!(
                records.iter().any(|(_, cigar)| cigar
                    .iter()
                    .any(|(k, n)| *k == Kind::Insertion && *n >= expected)),
                "{}: no CIGAR insertion of at least {}bp reached the BAM",
                case.tag,
                expected
            ),
            Signature::Chimeric => {
                let chimeric = mates
                    .keys()
                    .filter(|q| q.contains("EIDOLON_chimeric"))
                    .count();
                assert!(
                    chimeric >= 3,
                    "{}: only {chimeric} chimeric QNAMEs reached the BAM",
                    case.tag
                );
            }
        }
    }
}
