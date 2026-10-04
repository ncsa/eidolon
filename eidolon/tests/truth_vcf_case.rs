//! #790 — the truth VCF's bases are uppercase whatever the reference's case.
//!
//! CORRECTNESS CRITERION. A soft-masked reference marks repeats in lowercase; the bases are
//! the same bases. So the same run on a reference and on its lowercased copy writes the same
//! truth VCF, SVs and BND included. FALSIFIED by any record that differs between the two,
//! which before the fix was every SV's REF and every BND's anchor base.
//!
//! The comparison skips `##reference=`, which names the FASTA path and so differs by
//! construction. Everything else must match line for line.
//!
//! FIXTURE: H1N1, because de novo BND is inter-contig and needs at least two contigs. The
//! assertion is equality between two runs, not a number, so H1N1's size does not matter here.
//! The run must actually contain BND and other SV records, or equality would be vacuous.

mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use std::path::{Path, PathBuf};

/// Truth VCF lines from one gen-reads run, minus the `##reference=` header.
fn truth_vcf(work: &Path, name: &str, reference: PathBuf) -> Vec<String> {
    let mut config = GenReadsConfig::new(reference, work.to_path_buf(), name);
    config.produce_fastq = true;
    config.produce_vcf = true;
    config.coverage = 1;
    config.rng_seed = "truth vcf case".to_string();
    config.sv_rate_scale = Some(50.0);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    read_gzip_fastq_lines(&work.join(format!("{name}.vcf.gz")))
        .into_iter()
        .filter(|l| !l.starts_with("##reference="))
        .collect()
}

fn lowercased_h1n1(dir: &Path) -> PathBuf {
    let text: String = std::fs::read_to_string(h1n1_reference())
        .unwrap()
        .lines()
        .map(|l| {
            if l.starts_with('>') {
                format!("{l}\n")
            } else {
                format!("{}\n", l.to_ascii_lowercase())
            }
        })
        .collect();
    let path = dir.join("H1N1_lower.fa");
    std::fs::write(&path, text).unwrap();
    path
}

#[test]
fn a_lowercased_reference_writes_the_same_truth_vcf() {
    let (_g, work) = fresh_workdir();
    let upper = truth_vcf(&work, "upper", h1n1_reference());
    let lower = truth_vcf(&work, "lower", lowercased_h1n1(&work));

    // Denominator: the comparison means nothing unless SVs, BND included, were planted.
    let body: Vec<&String> = upper.iter().filter(|l| !l.starts_with('#')).collect();
    let bnd = body.iter().filter(|l| l.contains("SVTYPE=BND")).count();
    let other_sv = body
        .iter()
        .filter(|l| l.contains("SVTYPE=") && !l.contains("SVTYPE=BND"))
        .count();
    assert!(
        bnd > 0 && other_sv > 0,
        "the fixture planted {bnd} BND and {other_sv} other SV records; it needs both"
    );

    assert_eq!(upper.len(), lower.len(), "record counts differ");
    let differing: Vec<(&String, &String)> = upper
        .iter()
        .zip(lower.iter())
        .filter(|(u, l)| u != l)
        .collect();
    assert!(
        differing.is_empty(),
        "{} of {} truth VCF records differ only because the reference was lowercase (#790). \
         First:\n  upper: {}\n  lower: {}",
        differing.len(),
        body.len(),
        differing[0].0,
        differing[0].1
    );
}
