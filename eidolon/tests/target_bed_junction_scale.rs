//! Under a `target_bed`, `coverage` means mean on-target depth, so every fragment budget is
//! multiplied by the budget scale. An SV event's junction fragments must be too: otherwise a
//! heterozygous event's junction support stays at its unscaled count while local depth rises by
//! the scale, and reads at `1/scale` of what the depth implies (#889).
//!
//! Falsified when the junction-read count of a run whose budget scale is `s` is not `s` times
//! that of the same events in a run whose scale is 1: measured per SV type (DEL, INV, BND,
//! `<CNV>` loss), because each type has its own call site. The junction reads are counted from
//! the read names (`EIDOLON_chimeric_<TYPE>_`).

mod common;

use common::{eidolon, read_gzip_fastq_lines};
use std::path::{Path, PathBuf};

const CONTIG: &str = "synth1";
const CONTIG_LEN: usize = 1_000_000;
const WIDE: (usize, usize) = (100_000, 180_000);

fn synthetic_sequence(n: usize) -> String {
    let mut x: u64 = 0x0FED_CBA9_8765_4321;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ['A', 'C', 'G', 'T'][((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Six events of one type, spread through the wide target, as VCF records.
fn records(kind: &str) -> Vec<String> {
    let mut rows = Vec::new();
    for k in 0..6usize {
        let pos = 110_000 + k * 12_000;
        let end = pos + 1_000;
        match kind {
            "DEL" | "INV" => rows.push(format!(
                "{CONTIG}\t{pos}\te{k}\tN\t<{kind}>\t60\tPASS\tSVTYPE={kind};END={end};SVLEN={}\tGT\t0|1",
                if kind == "DEL" { -1000 } else { 1000 }
            )),
            "CNV" => rows.push(format!(
                "{CONTIG}\t{pos}\te{k}\tN\t<CNV>\t60\tPASS\tSVTYPE=CNV;END={end};CN=1\tGT\t0|1"
            )),
            "BND" => {
                rows.push(format!(
                    "{CONTIG}\t{pos}\tb{k}a\tG\tG]{CONTIG}:{end}]\t60\tPASS\tSVTYPE=BND;MATEID=b{k}b\tGT\t0|1"
                ));
                rows.push(format!(
                    "{CONTIG}\t{end}\tb{k}b\tC\tC]{CONTIG}:{pos}]\t60\tPASS\tSVTYPE=BND;MATEID=b{k}a\tGT\t0|1"
                ));
            }
            _ => unreachable!(),
        }
    }
    rows
}

/// Runs gen-reads and returns (junction reads of this SV type, the logged budget scale).
fn run(dir: &Path, name: &str, kind: &str, targets: &[(usize, usize)]) -> (usize, f64) {
    let out: PathBuf = dir.join(name);
    std::fs::create_dir_all(&out).unwrap();
    let mut fa = format!(">{CONTIG}\n");
    for chunk in synthetic_sequence(CONTIG_LEN).as_bytes().chunks(60) {
        fa.push_str(std::str::from_utf8(chunk).unwrap());
        fa.push('\n');
    }
    std::fs::write(out.join("ref.fa"), fa).unwrap();
    let bed: String = targets
        .iter()
        .map(|(s, e)| format!("{CONTIG}\t{s}\t{e}\n"))
        .collect();
    std::fs::write(out.join("t.bed"), bed).unwrap();
    let mut vcf = format!(
        "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={CONTIG_LEN}>\n\
         ##ALT=<ID={kind},Description=\"x\">\n\
         ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"x\">\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n\
         ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"x\">\n\
         ##INFO=<ID=CN,Number=1,Type=Integer,Description=\"x\">\n\
         ##INFO=<ID=MATEID,Number=.,Type=String,Description=\"x\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    );
    for r in records(kind) {
        vcf.push_str(&r);
        vcf.push('\n');
    }
    std::fs::write(out.join("in.vcf"), vcf).unwrap();
    std::fs::write(
        out.join("c.yml"),
        format!(
            "reference: {d}/ref.fa\ntarget_bed: {d}/t.bed\ninput_vcf: {d}/in.vcf\nread_len: 100\n\
             coverage: 30\nploidy: 2\npaired_ended: true\nfragment_mean: 250\nfragment_st_dev: 25\n\
             produce_vcf: false\nproduce_fastq: true\nproduce_bam: false\nshuffle_fastq: false\n\
             sv_rate_scale: 0.0\noverwrite_output: true\noutput_dir: {d}\noutput_basename: o\n\
             rng_seed: junction scale\nnum_threads: 4\n",
            d = out.display()
        ),
    )
    .unwrap();
    let output = eidolon()
        .args(["gen-reads", "-c"])
        .arg(out.join("c.yml"))
        .output()
        .unwrap();
    assert!(output.status.success(), "gen-reads failed");
    let log = String::from_utf8_lossy(&output.stdout).into_owned();
    let scale: f64 = log
        .lines()
        .find(|l| l.contains("target_bed:") && l.contains("scaled "))
        .and_then(|l| l.split("scaled ").nth(1))
        .and_then(|t| t.split('x').next())
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| panic!("no budget-scale line in the log:\n{log}"));
    let lines = read_gzip_fastq_lines(&out.join("o_r1.fastq.gz"));
    // A CNV loss is emitted by the DEL-style path, a breakend as `<contig>_<pos>_<contig>`.
    let tag = match kind {
        "CNV" => "EIDOLON_chimeric_DEL_".to_string(),
        "BND" => format!("EIDOLON_chimeric_{CONTIG}_"),
        _ => format!("EIDOLON_chimeric_{kind}_"),
    };
    let n = lines.iter().step_by(4).filter(|l| l.contains(&tag)).count();
    (n, scale)
}

fn junction_reads_follow_the_budget_scale(kind: &str) {
    let tmp = tempfile::tempdir().unwrap();
    let (one, s_one) = run(tmp.path(), "one", kind, &[WIDE]);
    let mut many = vec![WIDE];
    many.extend((0..800).map(|i| (300_000 + i * 600, 300_000 + i * 600 + 180)));
    let (scaled, s) = run(tmp.path(), "scaled", kind, &many);
    eprintln!("[junction-scale] {kind}: {one} reads at scale {s_one:.3}, {scaled} at scale {s:.3}");
    assert!(
        one > 20,
        "{kind}: only {one} junction reads, so this tests nothing"
    );
    assert!(
        s > 1.3 * s_one,
        "scale {s:.3} against {s_one:.3}: the BED did not scale the budget"
    );
    let ratio = scaled as f64 / one as f64;
    let expect = s / s_one;
    assert!(
        (ratio - expect).abs() < 0.08 * expect,
        "{kind}: junction reads x{ratio:.3} across a x{expect:.3} budget scale ({one} -> {scaled})"
    );
}

#[test]
fn deletion_junction_reads_follow_the_budget_scale() {
    junction_reads_follow_the_budget_scale("DEL");
}

#[test]
fn inversion_junction_reads_follow_the_budget_scale() {
    junction_reads_follow_the_budget_scale("INV");
}

#[test]
fn breakend_junction_reads_follow_the_budget_scale() {
    junction_reads_follow_the_budget_scale("BND");
}

#[test]
fn cnv_loss_junction_reads_follow_the_budget_scale() {
    junction_reads_follow_the_budget_scale("CNV");
}
