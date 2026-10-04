//! #536 — every mapped golden BAM record carries an `NM` that is right.
//!
//! CORRECTNESS CRITERION. `NM` is the edit distance to the reference: mismatched aligned
//! bases plus inserted plus deleted bases, soft clips excluded. `samtools calmd` recomputes
//! it from the BAM and the FASTA independently of eidolon, so the two must agree on every
//! record. FALSIFIED by any mapped record without `NM`, or with an `NM` calmd disagrees with.
//!
//! The fixture must exercise what can go wrong, or agreement proves little: mismatches,
//! inserted and deleted bases, soft-clipped adapter (a clipped base is not an edit), reverse
//! reads, and SV junction reads, which are staged by a separate writer. Each is counted and
//! required to be present.
//!
//! `#[ignore]`: needs samtools, which the release-gates workflow installs and the PR
//! workflow's Rust job does not.

mod common;
use common::{eidolon, fresh_workdir, h1n1_reference};
use std::process::Command;

fn samtools(args: &[&str]) -> String {
    let out = Command::new("samtools")
        .args(args)
        .output()
        .expect("samtools not on PATH — this test needs it");
    assert!(
        out.status.success(),
        "samtools {args:?} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn nm_field(line: &str) -> Option<i64> {
    line.split('\t')
        .skip(11)
        .find_map(|f| f.strip_prefix("NM:i:"))
        .map(|v| v.parse().unwrap())
}

#[test]
#[ignore]
fn every_mapped_record_carries_the_nm_samtools_calmd_computes() {
    let (_g, work) = fresh_workdir();
    let reference = h1n1_reference();
    let yaml = work.join("config.yml");
    std::fs::write(
        &yaml,
        format!(
            "reference: {}\noutput_dir: {}\noutput_filename: nm\nread_len: 100\n\
             coverage: 20\npaired_ended: true\nfragment_mean: 150\nfragment_st_dev: 50\n\
             produce_fastq: true\nproduce_bam: true\nrng_seed: \"golden bam nm\"\n\
             sv_rate_scale: 30.0\nadapters:\n  enabled: true\n  preset: truseq\n",
            reference.display(),
            work.display()
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    let bam = work.join("nm.bam");
    let bam = bam.to_str().unwrap();
    let ours = samtools(&["view", "-F", "4", bam]);
    let theirs = samtools(&["calmd", bam, reference.to_str().unwrap()]);
    let theirs: Vec<&str> = theirs
        .lines()
        .filter(|l| !l.starts_with('@'))
        .filter(|l| {
            let flag: u32 = l.split('\t').nth(1).unwrap().parse().unwrap();
            flag & 4 == 0
        })
        .collect();
    let ours: Vec<&str> = ours.lines().collect();
    assert_eq!(
        ours.len(),
        theirs.len(),
        "calmd and view disagree on record count"
    );

    // Denominators: the cases that can go wrong must each be present.
    let cigar = |l: &str| l.split('\t').nth(5).unwrap().to_string();
    let flag = |l: &str| -> u32 { l.split('\t').nth(1).unwrap().parse().unwrap() };
    let n_with = |c: char| ours.iter().filter(|l| cigar(l).contains(c)).count();
    let n_mismatch = ours.iter().filter(|l| nm_field(l).unwrap_or(0) > 0).count();
    let n_reverse = ours.iter().filter(|l| flag(l) & 16 != 0).count();
    let (n_ins, n_del, n_clip) = (n_with('I'), n_with('D'), n_with('S'));
    let n_junction = ours
        .iter()
        .filter(|l| l.split('\t').next().unwrap().contains("_chimeric_"))
        .count();
    println!(
        "{} mapped records: {n_mismatch} with NM>0, {n_ins} with I, {n_del} with D, \
         {n_clip} with S, {n_reverse} reverse, {n_junction} junction",
        ours.len()
    );
    assert!(
        n_mismatch > 0 && n_ins > 0 && n_del > 0 && n_clip > 0 && n_reverse > 0 && n_junction > 0,
        "the fixture is missing a case it exists to cover"
    );

    let missing = ours.iter().filter(|l| nm_field(l).is_none()).count();
    assert_eq!(
        missing,
        0,
        "{missing} of {} mapped records carry no NM",
        ours.len()
    );

    let disagree: Vec<(i64, i64, &str)> = ours
        .iter()
        .zip(&theirs)
        .filter_map(|(o, t)| {
            let (a, b) = (nm_field(o).unwrap(), nm_field(t).unwrap());
            (a != b).then_some((a, b, *o))
        })
        .collect();
    assert!(
        disagree.is_empty(),
        "{} of {} records carry an NM samtools calmd disagrees with. First: ours {} calmd {}\n{}",
        disagree.len(),
        ours.len(),
        disagree[0].0,
        disagree[0].1,
        disagree[0].2
    );
}
