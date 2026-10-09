//! #691 — a literal SV deletion that runs past its sub-region must not widen the alt haplotype.
//!
//! CORRECTNESS CRITERION. A sub-region's alt haplotype is sampled only over that sub-region's
//! own projection, plus the haplotype that genuinely follows it. So a heterozygous 500 bp
//! deletion costs about the reads its 500 bp would have carried at half depth, wherever it
//! sits. FALSIFIED by a read count that moves by more than that when the deletion straddles a
//! sub-region boundary.
//!
//! The defect: the alt window's end was projected from the sub-region's last base, and when
//! that base is deleted the projection is `None`, which fell back to `haplotype_len()`, the
//! end of the CONTIG. The alt haplotype was then sampled from the deletion to the contig end
//! at its half share of depth, on top of the reference coverage that region already gets.
//! Measured on ecoli at 10x: +73,518 R1 reads against a predicted +73,540 for exactly that,
//! and fragments planned past the haplotype end, which is #691's `span 561 materialized 356`.
//!
//! FIXTURE: ecoli (4.6 Mb, one contig). A `<DUP>` makes the sub-region boundary; the
//! deletion is placed to straddle its start. The DUP-only run is the control, so the DUP's
//! own extra reads cancel. Aligner-free: counts come straight from the FASTQ.

mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir};
use std::io::Write as _;
use std::path::PathBuf;

const DUP_START: usize = 200_000; // 1-based POS of the <DUP>; its sub-region starts here
const DUP_END: usize = 210_000;
const DEL_LEN: usize = 500;
const COVERAGE: usize = 4;

fn ecoli_reference() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
}

fn contig_sequence() -> String {
    std::fs::read_to_string(ecoli_reference())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect()
}

/// R1 reads from a run carrying the het `<DUP>` and, optionally, a het literal deletion of
/// `DEL_LEN` bases anchored at 1-based `del_pos`.
fn reads_generated(name: &str, del_pos: Option<usize>) -> usize {
    let (_g, work) = fresh_workdir();
    let seq = contig_sequence();

    let vcf = work.join("in.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "##contig=<ID=Chromosome,length={}>", seq.len()).unwrap();
    for (id, ty) in [
        ("SVTYPE", "String"),
        ("END", "Integer"),
        ("SVLEN", "Integer"),
    ] {
        writeln!(
            f,
            "##INFO=<ID={id},Number=1,Type={ty},Description=\"{id}\">"
        )
        .unwrap();
    }
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
    if let Some(pos) = del_pos {
        writeln!(
            f,
            "Chromosome\t{pos}\t.\t{}\t{}\t60\tPASS\t.\tGT\t0/1",
            &seq[pos - 1..pos + DEL_LEN],
            &seq[pos - 1..pos]
        )
        .unwrap();
    }
    writeln!(
        f,
        "Chromosome\t{DUP_START}\t.\t{}\t<DUP>\t60\tPASS\tSVTYPE=DUP;END={DUP_END};SVLEN={}\tGT\t0/1",
        &seq[DUP_START - 1..DUP_START],
        DUP_END - DUP_START
    )
    .unwrap();
    drop(f);

    let mut config = GenReadsConfig::new(ecoli_reference(), work.clone(), name);
    config.coverage = COVERAGE;
    config.read_len = 151;
    config.paired_ended = true;
    config.fragment_mean = Some(400.0);
    config.fragment_st_dev = Some(90.0);
    config.rng_seed = "straddling deletion window".to_string();
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
    let f = std::fs::File::open(work.join(format!("{name}_r1.fastq.gz"))).unwrap();
    BufReader::new(MultiGzDecoder::new(f)).lines().count() / 4
}

/// The read count a het `DEL_LEN` deletion may move by: its span at half depth is ~3 R1
/// reads at 4x. 0.1% of the control (~60 reads) leaves room for that and for placement noise,
/// and sits far below the defect, which adds ~48% of the contig's reads.
fn assert_costs_only_its_own_span(label: &str, control: usize, with_del: usize) {
    assert!(control > 0, "control produced no reads at all");
    let diff = with_del.abs_diff(control);
    let tol = control / 1000;
    assert!(
        diff <= tol,
        "{label}: a het {DEL_LEN} bp deletion moved the read count {control} -> {with_del} \
         ({diff} reads, tolerance {tol}). It should cost ~3 reads. An alt haplotype is being \
         sampled outside its own sub-region (#691)."
    );
}

#[test]
fn a_deletion_straddling_a_subregion_boundary_costs_only_its_own_span() {
    let control = reads_generated("dup_only", None);
    // Anchored 200 bp before the DUP: the deletion removes the last 199 bases of the
    // preceding sub-region and runs 301 bases into the DUP's.
    let straddle = reads_generated("straddle", Some(DUP_START - 200));
    assert_costs_only_its_own_span("straddling", control, straddle);
}

/// MUST NOT FIRE: the same deletion wholly inside a sub-region was always handled. A fix that
/// broke the ordinary case would show up here first.
#[test]
fn a_deletion_inside_a_subregion_costs_only_its_own_span() {
    let control = reads_generated("dup_only_inner", None);
    let inner = reads_generated("inner", Some(100_000));
    assert_costs_only_its_own_span("inside", control, inner);
}
