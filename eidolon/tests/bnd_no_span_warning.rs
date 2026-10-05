//! A well-formed breakend has no END or SVLEN (VCF 4.2 §5.4: it is a point), so the
//! coverage pass must not warn that one is missing. A symbolic SV that does need a span and
//! lacks it must still warn (#497).
//!
//! The breakend case asserts the DEBUG line as well as the missing WARN: without it, a
//! breakend that never reached the coverage pass would pass the test for the wrong reason.
//! Both tests read the written log: the terminal logger sends WARN to stdout, not stderr.

mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

const MISSING_SPAN_WARN: &str = "has no END/SVLEN";
const BND_DEBUG: &str = "Breakend at 1-based POS 600: a point event";

fn write_vcf(path: &std::path::Path, records: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
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
    for r in records {
        writeln!(f, "{r}").unwrap();
    }
}

/// Runs gen-reads on `records` at debug log level and returns the written log.
fn run(records: &[&str]) -> String {
    let (_dir, work) = fresh_workdir();
    let vcf = work.join("in.vcf");
    write_vcf(&vcf, records);
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "run");
    config.coverage = 5;
    config.produce_fastq = true;
    config.input_vcf = Some(vcf);
    config.sv_rate_scale = Some(0.0); // only the supplied records, no de novo SVs
    let yaml = config.write_yaml();
    let log = work.join("run.log");
    let out = eidolon()
        .args(["--log-level", "debug", "--log-dest"])
        .arg(&log)
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "gen-reads failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::read_to_string(&log).unwrap()
}

#[test]
fn a_breakend_without_end_or_svlen_does_not_warn() {
    // REF bases are the actual reference at those positions.
    let log = run(&[
        "H1N1_HA\t600\t.\tT\tT[H1N1_PB2:900[\t60\tPASS\tSVTYPE=BND\tGT\t1/1",
        "H1N1_PB2\t900\t.\tG\t]H1N1_HA:600]G\t60\tPASS\tSVTYPE=BND\tGT\t1/1",
    ]);
    assert!(
        log.contains(BND_DEBUG),
        "the breakend must reach the coverage pass and be skipped there:\n{log}"
    );
    assert!(
        !log.contains(MISSING_SPAN_WARN),
        "a breakend is not missing a span:\n{log}"
    );
}

#[test]
fn a_deletion_without_end_or_svlen_still_warns() {
    let log = run(&["H1N1_HA\t600\t.\tT\t<DEL>\t60\tPASS\tSVTYPE=DEL\tGT\t1/1"]);
    assert!(
        log.contains(&format!(
            "Symbolic SV at 1-based POS 600 {MISSING_SPAN_WARN}"
        )),
        "a <DEL> with no span is under-specified and must say so:\n{log}"
    );
}
