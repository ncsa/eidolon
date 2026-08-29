//! A large deletion must not remove coverage anywhere except where it is (#625).
//!
//! WHY THIS FILE EXISTS: `generate_fragments` refuses a region shorter than
//! `read_len + 2 * max_del_len`. That padding used to come from a CONTIG-WIDE maximum applied
//! to every sub-region, so one 5 kb deletion raised the floor to 10,151 bp everywhere and
//! silenced every sub-region below it — whether or not it contained a deletion. Measured on
//! ecoli 4.64 Mb: 205 of 356 sub-regions dropped, **22.3% of the contig producing zero reads**,
//! against 0.0% with no deletions. A caller sees enormous deletions no record describes, and
//! nothing warns.
//!
//! FIXTURE CHOICE: this uses ecoli, not a synthetic contig, because the defect needs a contig
//! that splits into many small sub-regions. A 150 kb synthetic contig was tried first and does
//! NOT reproduce it even with a 40 kb deletion — it has too few sub-regions. That is exactly
//! the H1N1 trap in reverse, and it is why the slower fixture is worth its runtime here.
//!
//! Aligner-free: read counts come straight from the FASTQ.

mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir};
use std::io::Write as _;
use std::path::{Path, PathBuf};

fn ecoli_reference() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
}

fn contig_sequence(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect()
}

/// Total R1 reads produced, optionally with a single homozygous deletion of `del_len`.
fn reads_generated(name: &str, del_len: Option<usize>) -> usize {
    let (_g, work) = fresh_workdir();
    let reference = ecoli_reference();

    let mut config = GenReadsConfig::new(reference.clone(), work.clone(), name);
    config.coverage = 6;
    config.produce_fastq = true;
    config.rng_seed = "deletion region padding".to_string();
    config.sv_rate_scale = Some(0.0);

    if let Some(len) = del_len {
        let seq = contig_sequence(&reference);
        let pos = 200_000usize;
        let vcf = work.join("in.vcf");
        let mut f = std::fs::File::create(&vcf).unwrap();
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
        writeln!(
            f,
            "Chromosome\t{pos}\t.\t{}\t{}\t60\tPASS\tSVTYPE=DEL\tGT\t1/1",
            &seq[pos - 1..pos + len],
            &seq[pos - 1..pos]
        )
        .unwrap();
        config.input_vcf = Some(vcf);
    }

    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    use flate2::read::MultiGzDecoder;
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(work.join(format!("{name}_r1.fastq.gz"))).unwrap();
    BufReader::new(MultiGzDecoder::new(f)).lines().count() / 4
}

/// One 5 kb deletion must cost ~5 kb of sequence, not a fifth of the contig.
///
/// Measured separation: the contig-wide padding gives a ratio of **0.776**, the per-region
/// padding **0.999**. The threshold sits between them with room on both sides, and the
/// failure mode being guarded is a collapse, not a wobble.
#[test]
fn a_large_deletion_does_not_silence_the_rest_of_the_contig() {
    let baseline = reads_generated("pad_none", None);
    let with_del = reads_generated("pad_del", Some(5_000));

    assert!(baseline > 0, "control produced no reads at all");

    let ratio = with_del as f64 / baseline as f64;
    assert!(
        ratio >= 0.97,
        "a single 5 kb deletion cost {:.1}% of all reads on a 4.6 Mb contig \
         ({baseline} -> {with_del}). It should cost ~0.1%. A contig-wide deletion padding \
         is silencing sub-regions that contain no deletion (#625).",
        (1.0 - ratio) * 100.0
    );
}

/// MUST NOT FIRE: the deletion itself still has to reach the molecule. A "fix" that merely
/// stopped padding everywhere would sail through the test above while breaking deletions.
///
/// Uses a SYNTHETIC single-sub-region contig, deliberately, because the assertion is about
/// the molecule getting shorter and on ecoli that shortening is confined to the sub-region
/// holding the deletion — the furthest read start is set by a later sub-region and does not
/// move. Measured: 4641584 -> 4641585, i.e. nothing, against correct code.
///
/// Asserted via haplotype LENGTH, not a coordinate window: read names carry positions on the
/// molecule being sequenced, so after a homozygous deletion every downstream read is named
/// lower. "No read starts inside the deleted span" therefore FAILS against correct code — an
/// earlier version of this test asserted exactly that and had to be withdrawn.
#[test]
fn a_homozygous_deletion_shortens_the_sequenced_molecule() {
    const LEN: usize = 150_000;
    const DEL: usize = 5_000;

    let synthetic = |dir: &Path| -> PathBuf {
        let bases = b"ACGT";
        let mut seq = Vec::with_capacity(LEN);
        let mut x: usize = 7;
        for _ in 0..LEN {
            x = (x * 1103515245 + 12345) % 2147483648;
            seq.push(bases[(x >> 16) % 4]);
        }
        let p = dir.join("ref.fa");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, ">c1").unwrap();
        for chunk in seq.chunks(70) {
            f.write_all(chunk).unwrap();
            writeln!(f).unwrap();
        }
        p
    };

    let furthest = |name: &str, del: bool| -> usize {
        let (_g, work) = fresh_workdir();
        let reference = synthetic(&work);
        let mut config = GenReadsConfig::new(reference.clone(), work.clone(), name);
        config.coverage = 20;
        config.produce_fastq = true;
        config.rng_seed = "deletion shortening".to_string();
        config.sv_rate_scale = Some(0.0);
        if del {
            let seq = contig_sequence(&reference);
            let pos = 70_000usize;
            let vcf = work.join("in.vcf");
            let mut f = std::fs::File::create(&vcf).unwrap();
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
            writeln!(
                f,
                "c1\t{pos}\t.\t{}\t{}\t60\tPASS\tSVTYPE=DEL\tGT\t1/1",
                &seq[pos - 1..pos + DEL],
                &seq[pos - 1..pos]
            )
            .unwrap();
            config.input_vcf = Some(vcf);
        }
        let yaml = config.write_yaml();
        eidolon()
            .args(["gen-reads", "-c"])
            .arg(yaml.path())
            .assert()
            .success();

        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let f = std::fs::File::open(work.join(format!("{name}_r1.fastq.gz"))).unwrap();
        BufReader::new(MultiGzDecoder::new(f))
            .lines()
            .map(|l| l.unwrap())
            .filter(|l| l.starts_with("@EIDOLON_generated_c1_"))
            .filter_map(|l| {
                l.strip_prefix("@EIDOLON_generated_c1_")?
                    .split('_')
                    .next()?
                    .parse::<usize>()
                    .ok()
            })
            .max()
            .unwrap_or(0)
    };

    let plain = furthest("shorten_none", false);
    let deleted = furthest("shorten_del", true);
    assert!(plain > 0 && deleted > 0, "no reads generated");

    let shrink = plain.saturating_sub(deleted);
    assert!(
        (3_500..=6_500).contains(&shrink),
        "a homozygous {DEL} bp deletion shortened the molecule by {shrink} bp (furthest \
         read start {plain} -> {deleted}); expected ~{DEL}. The deletion is not reaching \
         the molecule being sequenced."
    );
}

/// The OTHER half of #625: a deletion must not silence its own neighbourhood either.
///
/// The alt-haplotype span is in HAPLOTYPE coordinates, where the deletions on that molecule
/// have already been removed — there is no gap left to span, so padding its minimum region
/// for one is wrong twice over. Padding it anyway drops the sub-region holding the deletion:
/// measured, a 100 bp deletion sitting beside a 5 kb one produced a **3,628 bp** zero-depth
/// hole (497913-501540) instead of 100 bp.
///
/// Probes just BEFORE the small deletion, inside its own sub-region, where read names are
/// still ordinary coordinates. Aligner-free.
#[test]
fn a_small_deletion_beside_a_large_one_keeps_its_neighbourhood() {
    let (_g, work) = fresh_workdir();
    let reference = ecoli_reference();
    let seq = contig_sequence(&reference);

    let vcf = work.join("in.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
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
    // A 5 kb deletion far away sets the contig-wide maximum; a 100 bp one is the victim.
    for (pos, len) in [(200_000usize, 5_000usize), (500_000, 100)] {
        writeln!(
            f,
            "Chromosome\t{pos}\t.\t{}\t{}\t60\tPASS\tSVTYPE=DEL\tGT\t1/1",
            &seq[pos - 1..pos + len],
            &seq[pos - 1..pos]
        )
        .unwrap();
    }
    drop(f);

    let mut config = GenReadsConfig::new(reference, work.clone(), "pad_neighbour");
    config.coverage = 20;
    config.produce_fastq = true;
    config.rng_seed = "deletion region padding".to_string();
    config.sv_rate_scale = Some(0.0);
    config.input_vcf = Some(vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    use flate2::read::MultiGzDecoder;
    use std::io::{BufRead, BufReader};
    let fh = std::fs::File::open(work.join("pad_neighbour_r1.fastq.gz")).unwrap();
    let near = BufReader::new(MultiGzDecoder::new(fh))
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| l.starts_with("@EIDOLON_generated_Chromosome_"))
        .filter_map(|l| {
            l.strip_prefix("@EIDOLON_generated_Chromosome_")?
                .split('_')
                .next()?
                .parse::<usize>()
                .ok()
        })
        .filter(|p| (499_000..499_900).contains(p))
        .count();

    assert!(
        near > 0,
        "no reads at all in the 900 bp immediately before a 100 bp deletion. Its whole \
         sub-region has been dropped because a 5 kb deletion 300 kb away padded the \
         haplotype span's minimum region (#625)."
    );
}
