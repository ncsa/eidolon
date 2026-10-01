//! #780 — both mates of a fragment carry the same allele.
//!
//! CORRECTNESS CRITERION. A fragment is one molecule from one chromosome copy, so at a variant
//! both of its mates read the same allele. If this works, every read pair that covers a
//! heterozygous SNP with both mates shows the same base there, and about half of those pairs
//! carry the alternate. It is FALSIFIED by pairs whose mates disagree: a per-read allele draw
//! makes about half of het pairs disagree.
//!
//! KNOWN ANSWER, read straight off the FASTQ. Each SNP sits at the center of a 21-mer that
//! occurs once in the reference, and whose alternate form occurs nowhere, so a read's allele is
//! whichever 21-mer it contains (either orientation). Fragments are 150 bp and reads 100 bp, so
//! a pair's mates overlap on the middle 50 bp and many pairs cover each SNP with both mates.
//! Reads come from a Q40 model (error ~1e-4), so the tolerance for discordance is a few pairs.
//!
//!   600  T>C  het (0/1)             concordant; alt in ~half of pairs
//!   850  A>G  hom (1/1)             MUST NOT FIRE: every pair alt in both mates
//!   1100 G>A  het, INFO AF=0.2      the allele_fraction path: concordant; alt in ~a fifth
//!
//! The denominator (pairs with both mates covering the 21-mer) is asserted per site.
//!
//! Fixture: H1N1_HA, because nothing here depends on genome size.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use std::path::Path;

const READ_LEN: usize = 100;

/// (position, ref 21-mer, alt 21-mer)
const SITES: [(usize, &str, &str); 3] = [
    (600, "ACCATCCATCTACTAGTGCTG", "ACCATCCATCCACTAGTGCTG"),
    (850, "TGGTATTATCATTTCAGATAC", "TGGTATTATCGTTTCAGATAC"),
    (1100, "GGATGGTACGGTTATCACCAT", "GGATGGTACGATTATCACCAT"),
];

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// Which allele a read carries at a site: Some(true) alt, Some(false) ref, None not covered.
fn allele(read: &str, reference: &str, alt: &str) -> Option<bool> {
    let has = |k: &str| read.contains(k) || read.contains(&revcomp(k));
    match (has(reference), has(alt)) {
        (true, false) => Some(false),
        (false, true) => Some(true),
        _ => None,
    }
}

fn write_inputs(work: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let vcf = work.join("sites.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=H1N1_HA,length=1701>\n\
         ##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         H1N1_HA\t600\thet\tT\tC\t60\tPASS\t.\tGT\t0/1\n\
         H1N1_HA\t850\thom\tA\tG\t60\tPASS\t.\tGT\t1/1\n\
         H1N1_HA\t1100\taf\tG\tA\t60\tPASS\tAF=0.2\tGT\t0/1\n",
    )
    .unwrap();

    let fq = work.join("q40.fastq");
    let seq: String = "ACGT".chars().cycle().take(READ_LEN).collect();
    let body: String = (0..200)
        .map(|i| format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(READ_LEN)))
        .collect();
    std::fs::write(&fq, body).unwrap();
    let model = work.join("q40.json.gz");
    let cfg = work.join("fit.yml");
    std::fs::write(
        &cfg,
        format!(
            "fastq_file: {}\noutput_file: {}\noverwrite_output: true\nmax_reads: 0\nqual_offset: 33\n",
            fq.display(),
            model.display()
        ),
    )
    .unwrap();
    let out = eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(&cfg)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    (vcf, model)
}

#[test]
fn both_mates_of_a_fragment_carry_the_same_allele() {
    let (_g, work) = fresh_workdir();
    let (vcf, model) = write_inputs(&work);

    let mut cfg = GenReadsConfig::new(h1n1_reference(), work.clone(), "conc");
    cfg.paired_ended = true;
    cfg.read_len = READ_LEN;
    cfg.coverage = 1000;
    cfg.fragment_mean = Some(150.0);
    cfg.fragment_st_dev = Some(5.0);
    cfg.mutation_rate = Some(0.0);
    cfg.rng_seed = "780".to_string();
    cfg.input_vcf = Some(vcf);
    cfg.sequence_error_model = Some(model);
    let yaml = cfg.write_yaml();
    let out = eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let r1 = read_gzip_fastq_lines(&work.join("conc_r1.fastq.gz"));
    let r2 = read_gzip_fastq_lines(&work.join("conc_r2.fastq.gz"));
    assert_eq!(
        r1.len(),
        r2.len(),
        "R1 and R2 must pair up record for record"
    );

    let mut failures = Vec::new();
    for (pos, reference, alt) in SITES {
        let (mut alt_pairs, mut ref_pairs, mut discordant) = (0usize, 0usize, 0usize);
        for (a, b) in r1.chunks(4).zip(r2.chunks(4)) {
            let (ra, rb) = (allele(&a[1], reference, alt), allele(&b[1], reference, alt));
            match (ra, rb) {
                (Some(true), Some(true)) => alt_pairs += 1,
                (Some(false), Some(false)) => ref_pairs += 1,
                (Some(_), Some(_)) => discordant += 1,
                _ => {}
            }
        }
        let n = alt_pairs + ref_pairs + discordant;
        let alt_share = alt_pairs as f64 / (alt_pairs + ref_pairs).max(1) as f64;
        eprintln!(
            "pos {pos}: {n} pairs cover it with both mates; {alt_pairs} alt, {ref_pairs} ref, \
             {discordant} discordant; alt share {alt_share:.3}"
        );
        if n < 100 {
            failures.push(format!(
                "pos {pos}: only {n} pairs cover the site with both mates"
            ));
            continue;
        }
        if discordant as f64 > 0.02 * n as f64 {
            failures.push(format!(
                "pos {pos}: {discordant} of {n} pairs have mates carrying different alleles"
            ));
        }
        let (lo, hi) = match pos {
            600 => (0.35, 0.65),
            850 => (0.98, 1.0),
            _ => (0.10, 0.30),
        };
        if !(lo..=hi).contains(&alt_share) {
            failures.push(format!(
                "pos {pos}: alt share {alt_share:.3} outside [{lo}, {hi}]"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
