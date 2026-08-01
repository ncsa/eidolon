//! Tests for `eidolon compare-af`, closing #466 for this helper.
//!
//! Two kinds, deliberately:
//!
//!   * **Golden** — the Rust must reproduce, byte for byte, what the Python it replaced
//!     produced on the same inputs. That is the migration's evidence, preserved so it
//!     survives deleting the Python.
//!   * **Known-answer** — inputs whose correct output is computable by hand, so the
//!     tests are not merely "whatever the implementation does". 7 alt reads in 100 is a
//!     fraction of exactly 0.070 whatever any implementation thinks.

mod common;

use common::eidolon;
use std::io::Write as _;
use std::path::PathBuf;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_data/compare_af")
}

/// Run the subcommand, returning `(exit_ok, stdout, stderr)` with the logger's own
/// lines stripped — those carry timestamps and are not part of the measurement.
fn run(args: &[&str]) -> (bool, String, String) {
    let out = eidolon().arg("compare-af").args(args).output().unwrap();
    let strip = |s: &[u8]| -> String {
        String::from_utf8_lossy(s)
            .lines()
            .filter(|l| {
                !l.contains("[INFO]")
                    && !l.contains("Welcome to eidolon")
                    && !l.contains("Processing finished")
            })
            .map(|l| format!("{l}\n"))
            .collect()
    };
    (out.status.success(), strip(&out.stdout), strip(&out.stderr))
}

#[test]
fn reproduces_the_python_output_byte_for_byte() {
    let f = fixtures();
    for (sim, stem, should_pass) in [
        ("sim.vcf.gz", "py_reference_sim", true),
        ("sim_450.vcf.gz", "py_reference_sim450", false),
    ] {
        let (ok, stdout, stderr) = run(&[
            "--truth",
            f.join("truth.vcf.gz").to_str().unwrap(),
            "--sim",
            f.join(sim).to_str().unwrap(),
            "--min-depth",
            "25",
        ]);
        let want_out = std::fs::read_to_string(f.join(format!("{stem}.stdout"))).unwrap();
        let want_err = std::fs::read_to_string(f.join(format!("{stem}.stderr"))).unwrap();
        assert_eq!(
            stdout, want_out,
            "{sim}: stdout diverged from the Python reference"
        );
        assert_eq!(
            stderr, want_err,
            "{sim}: stderr diverged from the Python reference"
        );
        assert_eq!(
            ok, should_pass,
            "{sim}: exit status diverged from the Python reference"
        );
    }
}

/// Write a minimal VCF with FORMAT/AD, and return its path.
fn write_vcf(dir: &std::path::Path, name: &str, records: &[&str]) -> PathBuf {
    let p = dir.join(name);
    let mut fh = std::fs::File::create(&p).unwrap();
    writeln!(fh, "##fileformat=VCFv4.2").unwrap();
    writeln!(fh, "##contig=<ID=chr1,length=100000>").unwrap();
    writeln!(
        fh,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
    )
    .unwrap();
    writeln!(
        fh,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
    )
    .unwrap();
    for r in records {
        writeln!(fh, "{r}").unwrap();
    }
    p
}

#[test]
fn seven_alt_reads_in_one_hundred_is_a_fraction_of_exactly_seven_percent() {
    // Known answer, independent of the implementation: AD=93,7 over a total of 100 is
    // 0.070. Truth declares 0.070; a correct reader reports MAE 0.0000.
    let tmp = tempfile::tempdir().unwrap();
    let truth = write_vcf(
        tmp.path(),
        "t.vcf",
        &[
            "chr1\t100\t.\tA\tG\t60\tPASS\t.\tAD\t930,70",
            "chr1\t200\t.\tA\tG\t60\tPASS\t.\tAD\t825,175",
            "chr1\t300\t.\tA\tG\t60\tPASS\t.\tAD\t650,350",
        ],
    );
    let sim = write_vcf(
        tmp.path(),
        "s.vcf",
        &[
            "chr1\t100\t.\tA\tG\t60\tPASS\t.\tAD\t93,7",
            "chr1\t200\t.\tA\tG\t60\tPASS\t.\tAD\t165,35",
            "chr1\t300\t.\tA\tG\t60\tPASS\t.\tAD\t130,70",
        ],
    );
    let (ok, stdout, _) = run(&[
        "--truth",
        truth.to_str().unwrap(),
        "--sim",
        sim.to_str().unwrap(),
    ]);
    assert!(ok, "should pass:\n{stdout}");
    assert!(
        stdout.contains("MAE=0.0000"),
        "0.070/0.175/0.350 on both sides must give MAE exactly 0:\n{stdout}"
    );
    assert!(
        stdout.contains("n=3"),
        "all three sites must be compared:\n{stdout}"
    );
}

#[test]
fn multiallelic_selects_the_matching_allele_not_the_first() {
    // mpileup emits `ALT=T,<*>` with `AD=93,7,0`. The truth's single ALT `T` must be
    // matched against the ALT LIST and ITS OWN AD element selected. Taking `[0]`
    // unconditionally reported the first allele's fraction for every allele — the #450
    // multi-allelic defect. Here ALT=G,T carries AD=900,30,70, so T's own fraction is
    // 70/1000 = 0.070 while the first allele's is 0.030.
    let tmp = tempfile::tempdir().unwrap();
    // Two sites, because a correlation needs at least two points — the same floor the
    // Python enforced. The second is plain; the first is the multi-allelic case.
    let truth = write_vcf(
        tmp.path(),
        "t.vcf",
        &[
            "chr1\t100\t.\tA\tT\t60\tPASS\t.\tAD\t9300,700",
            "chr1\t200\t.\tA\tG\t60\tPASS\t.\tAD\t800,200",
        ],
    );
    // The matching allele is deliberately SECOND. With it first, `counts[alt_idx + 1]`
    // and a hardcoded `counts[1]` are the same expression and the test proves nothing —
    // mutation testing caught exactly that.
    let sim = write_vcf(
        tmp.path(),
        "s.vcf",
        &[
            "chr1\t100\t.\tA\tG,T\t60\tPASS\t.\tAD\t900,30,70",
            "chr1\t200\t.\tA\tG\t60\tPASS\t.\tAD\t80,20",
        ],
    );
    let (ok, stdout, _) = run(&[
        "--truth",
        truth.to_str().unwrap(),
        "--sim",
        sim.to_str().unwrap(),
    ]);
    assert!(
        ok,
        "the `<*>` placeholder must not disqualify the real base:\n{stdout}"
    );
    assert!(
        stdout.contains("shared=2") && stdout.contains("MAE=0.0000"),
        "T should match its own AD element (70/1000 = 0.070), not the first allele's:\n{stdout}"
    );
}

#[test]
fn a_covered_site_with_no_observed_alt_scores_zero_rather_than_vanishing() {
    // This is the #450 failure mode in miniature. The sites most likely to have no alt
    // reads at all are the LOWEST-VAF ones, so dropping them biases the result
    // optimistically — exactly the direction that made a broken harness look clean.
    let tmp = tempfile::tempdir().unwrap();
    let truth = write_vcf(
        tmp.path(),
        "t.vcf",
        &[
            "chr1\t100\t.\tA\tG\t60\tPASS\t.\tAD\t900,100",
            "chr1\t200\t.\tA\tG\t60\tPASS\t.\tAD\t950,50",
        ],
    );
    // Position 200 is COVERED but lists no alt allele: zero reads carried it.
    let sim = write_vcf(
        tmp.path(),
        "s.vcf",
        &[
            "chr1\t100\t.\tA\tG\t60\tPASS\t.\tAD\t90,10",
            "chr1\t200\t.\tA\t<*>\t60\tPASS\t.\tAD\t100,0",
        ],
    );
    let (_, stdout, _) = run(&[
        "--truth",
        truth.to_str().unwrap(),
        "--sim",
        sim.to_str().unwrap(),
    ]);
    assert!(
        stdout.contains("1 truth allele(s) had coverage but zero observed reads"),
        "the covered-but-unobserved site must be zero-filled, not dropped:\n{stdout}"
    );
    assert!(
        stdout.contains("shared=2"),
        "both sites must be compared:\n{stdout}"
    );
}

#[test]
fn the_coverage_gate_fires_and_the_metrics_alone_would_not_have() {
    // The #450 reproduction has BETTER headline numbers than the complete set — bias
    // +0.0011 vs +0.0008, MAE 0.0227 vs 0.0240 — because dropping the lowest-VAF
    // stratum flatters the result. Only the denominator distinguishes them.
    let f = fixtures();
    let (ok, stdout, stderr) = run(&[
        "--truth",
        f.join("truth.vcf.gz").to_str().unwrap(),
        "--sim",
        f.join("sim_450.vcf.gz").to_str().unwrap(),
        "--min-depth",
        "25",
    ]);
    assert!(!ok, "18.3% unscored must fail the gate");
    assert!(
        stderr.contains("FAIL: 18.3% of planted truth alleles went unscored"),
        "stderr must name the shortfall:\n{stderr}"
    );
    assert!(
        stdout.contains("NOTHING SCORED"),
        "the excluded stratum must appear as an explicit row, not vanish:\n{stdout}"
    );
    // ...and the metrics it reports look FINE, which is the point.
    assert!(
        stdout.contains("MAE=0.0227"),
        "the misleading-but-good MAE should still be shown:\n{stdout}"
    );
}

#[test]
fn a_wrong_arity_field_is_refused_rather_than_guessed() {
    // A per-allele field whose length disagrees with the ALT count cannot be indexed
    // safely. Guessing would silently attribute one allele's number to another, so the
    // allele is skipped instead — here that leaves nothing comparable.
    let tmp = tempfile::tempdir().unwrap();
    let truth = write_vcf(
        tmp.path(),
        "t.vcf",
        &["chr1\t100\t.\tA\tG\t60\tPASS\t.\tAD\t900,100"],
    );
    // Two ALTs declared, but AD carries only two values where Number=R wants three.
    let sim = write_vcf(
        tmp.path(),
        "s.vcf",
        &["chr1\t100\t.\tA\tG,T\t60\tPASS\t.\tAD\t90,10"],
    );
    let (ok, stdout, stderr) = run(&[
        "--truth",
        truth.to_str().unwrap(),
        "--sim",
        sim.to_str().unwrap(),
    ]);
    assert!(
        !ok,
        "with the malformed allele refused there is nothing to compare, which must fail \
         loudly rather than report a number over one site:\n{stdout}{stderr}"
    );
    assert!(
        stderr.contains("fewer than 2 comparable sites"),
        "should say why:\n{stderr}"
    );
}

#[test]
fn multiallelic_af_field_selects_the_matching_allele() {
    // The AD path and the FORMAT/AF path index differently (`Number=R` vs `Number=A`)
    // and are served by different code. A multi-allelic test that only exercises AD
    // leaves `pick_per_allele` — the function that does the AF selection — completely
    // untested, which is exactly what mutation testing caught here.
    //
    // ALT=G,T with AF=0.10,0.40: the truth's ALT is T, so the answer is 0.40. Taking
    // the first element would report 0.10.
    let tmp = tempfile::tempdir().unwrap();
    let truth = write_vcf(
        tmp.path(),
        "t.vcf",
        &[
            "chr1\t100\t.\tA\tT\t60\tPASS\tAF=0.40\tAD\t600,400",
            "chr1\t200\t.\tA\tG\t60\tPASS\tAF=0.20\tAD\t800,200",
        ],
    );
    let sim = write_vcf(
        tmp.path(),
        "s.vcf",
        &[
            "chr1\t100\t.\tA\tG,T\t60\tPASS\tAF=0.10,0.40\tAD\t50,10,40",
            "chr1\t200\t.\tA\tG\t60\tPASS\tAF=0.20\tAD\t80,20",
        ],
    );
    let (ok, stdout, _) = run(&[
        "--truth",
        truth.to_str().unwrap(),
        "--sim",
        sim.to_str().unwrap(),
    ]);
    assert!(ok, "should compare cleanly:\n{stdout}");
    assert!(
        stdout.contains("MAE=0.0000"),
        "ALT T must take AF=0.40, its OWN element, not the first allele's 0.10:\n{stdout}"
    );
}
