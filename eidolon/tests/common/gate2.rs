//! Shared machinery for **Gate 2** — realigning eidolon's FASTQ with a real aligner and asking
//! whether the result carries the evidence a structural-variant caller keys on.
//!
//! See `docs/sv_polish_roadmap.md` for what the gates are. The short version: every other SV
//! test inspects eidolon describing its own work (its FASTQ, its golden BAM, its truth VCF). A
//! caller sees none of those — it sees reads somebody else aligned. This module is the seam.
//!
//! Analysis is deliberately done in Rust over the SAM with `noodles` rather than by piping
//! `samtools` through `awk`: the arithmetic *is* the assertion, so it should be readable and
//! debuggable. That choice caught a bug on its first run — depth was being accumulated across
//! all eight H1N1 contigs, which backfilled a deleted window and turned a homozygous deletion
//! into an apparent 1.2x *enrichment*.

#![allow(dead_code)]

use super::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::io::{BufRead, Write as _};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Reference window used as the "unaffected" baseline in every gate. Well clear of the event
/// loci the gates plant, and on the same contig so it shares that contig's coverage.
pub const FLANK_WINDOW: (usize, usize) = (1_000, 1_400);

/// A symbolic SV to plant via `input_vcf`.
pub struct SvSpec {
    pub svtype: &'static str,
    pub contig: &'static str,
    pub pos: usize,
    pub end: usize,
    /// `"1/1"` for homozygous. Gates use hom so a signature is present or absent rather than
    /// halved — a het event turns every assertion below into a ratio judgement.
    pub gt: &'static str,
}

/// Locate bwa-mem2, or fail with something actionable. **Never returns a "skip"** — a gate that
/// silently passes when its aligner is missing is worth less than no gate at all.
pub fn bwa_mem2() -> String {
    if let Ok(p) = std::env::var("BWA_MEM2") {
        assert!(
            Path::new(&p).is_file(),
            "BWA_MEM2={p} does not exist. Gate 2 cannot run without an aligner."
        );
        return p;
    }
    let found = Command::new("bwa-mem2")
        .arg("version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    assert!(
        found,
        "bwa-mem2 not found on PATH.\n\
         Gate 2 realigns eidolon's FASTQ with a real aligner; there is nothing to assert \
         without one.\n\
         On this workstation it lives in the `aln` conda environment:\n\
         \n    conda activate aln\n\
         \nor point at it directly:\n\
         \n    BWA_MEM2=/path/to/bwa-mem2 cargo test --test <gate> -- --ignored\n"
    );
    "bwa-mem2".to_string()
}

/// Generate paired reads over H1N1, optionally planting one symbolic SV. The control run uses
/// the same reference, seed and coverage, so the two differ **only** by the variant.
pub fn generate_reads(work: &Path, tag: &str, sv: Option<&SvSpec>) -> (PathBuf, PathBuf) {
    let mut config = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), tag);
    config.coverage = 60;
    config.read_len = 100;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.produce_bam = false;
    config.produce_vcf = true;
    // No de novo variants: the only difference between run and control must be the planted SV.
    config.mutation_rate = Some(0.0);
    config.sv_rate_scale = Some(0.0);

    if let Some(sv) = sv {
        let input_vcf = work.join(format!("{tag}.vcf"));
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
        writeln!(
            f,
            "{}\t{}\t.\tG\t<{}>\t60\tPASS\tSVTYPE={};END={}\tGT\t{}",
            sv.contig, sv.pos, sv.svtype, sv.svtype, sv.end, sv.gt
        )
        .unwrap();
        config.input_vcf = Some(input_vcf);
    }

    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let r1 = work.join(format!("{tag}_r1.fastq.gz"));
    let r2 = work.join(format!("{tag}_r2.fastq.gz"));
    assert!(r1.is_file() && r2.is_file(), "expected {r1:?} and {r2:?}");
    (r1, r2)
}

/// Align a FASTQ pair with bwa-mem2 and return the path to the SAM.
pub fn align(bwa: &str, work: &Path, tag: &str, r1: &Path, r2: &Path) -> PathBuf {
    // Index into the work dir so the repo's test_data is never written to.
    let local_ref = work.join("ref.fa");
    if !local_ref.exists() {
        std::fs::copy(h1n1_reference(), &local_ref).unwrap();
        let out = Command::new(bwa)
            .arg("index")
            .arg(&local_ref)
            .output()
            .expect("bwa-mem2 index failed to spawn");
        assert!(
            out.status.success(),
            "bwa-mem2 index failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }

    let sam_path = work.join(format!("{tag}.sam"));
    let out = Command::new(bwa)
        .args(["mem", "-t", "2"])
        .arg(&local_ref)
        .arg(r1)
        .arg(r2)
        .output()
        .expect("bwa-mem2 mem failed to spawn");
    assert!(
        out.status.success(),
        "bwa-mem2 mem failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::write(&sam_path, &out.stdout).unwrap();

    // A SAM with a header and no alignments would satisfy every "no signal" branch below, so
    // establish that reads were actually placed before anything is measured.
    let mapped = std::io::BufReader::new(std::fs::File::open(&sam_path).unwrap())
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('@'))
        .count();
    assert!(
        mapped > 100,
        "{tag}: bwa-mem2 emitted only {mapped} alignment record(s) — the alignment step failed, \
         so nothing measured from it would mean anything"
    );
    sam_path
}

/// What a caller can see, extracted from a realigned SAM.
#[derive(Debug, Default)]
pub struct Signatures {
    /// Mean depth strictly inside the event (20 bp in from each breakpoint, so a read clipped
    /// at a junction cannot contribute and the window measures the event, not its edges).
    pub depth_inside: f64,
    /// Mean depth over `FLANK_WINDOW`, well clear of the event.
    pub depth_outside: f64,
    /// Soft clips whose clip point is within ±25 bp of either breakpoint.
    pub clips_at_breakpoints: usize,
    /// Soft clips anywhere else — the background the above must stand out from.
    pub clips_elsewhere: usize,
    /// Leftmost-read-of-pair with |TLEN| inflated well past the fragment mean. The
    /// **deletion** signature: a pair spanning a deletion aligns further apart than it was.
    pub long_pairs: usize,
    /// Leftmost read of the pair is on the reverse strand — "everted" / RF orientation. The
    /// **tandem-duplication** signature: a pair spanning the duplication junction reads out of
    /// the second copy into the first, so the mates appear swapped.
    pub everted_pairs: usize,
    pub reads: usize,
}

/// Accumulate the signatures over `contig` for an event spanning `pos..=end`.
pub fn analyse(sam_path: &Path, contig: &str, pos: usize, end: usize) -> Signatures {
    let mut reader = std::fs::File::open(sam_path)
        .map(std::io::BufReader::new)
        .map(sam::io::Reader::new)
        .unwrap();
    let _header = reader.read_header().unwrap();

    let contig_len = 2_000usize;
    let mut depth = vec![0usize; contig_len];
    let mut sig = Signatures::default();

    for result in reader.records() {
        let record = result.unwrap();

        // H1N1 has EIGHT contigs. Without this filter every contig's reads land in one depth
        // array, the event window gets backfilled by unrelated contigs, and the signal vanishes.
        let on_target = matches!(
            record.reference_sequence_name(),
            Some(n) if n == contig.as_bytes()
        );
        if !on_target {
            continue;
        }
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        sig.reads += 1;
        let start = usize::from(start); // 1-based

        let mut ref_pos = start;
        let ops: Vec<_> = record.cigar().iter().map(|o| o.unwrap()).collect();
        for (i, op) in ops.iter().enumerate() {
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for p in ref_pos..(ref_pos + len).min(contig_len) {
                        depth[p] += 1;
                    }
                    ref_pos += len;
                }
                Kind::Deletion | Kind::Skip => ref_pos += len,
                Kind::SoftClip => {
                    // A leading clip points at the read's start; a trailing clip at ref_pos.
                    let clip_at = if i == 0 { start } else { ref_pos };
                    if clip_at.abs_diff(pos) <= 25 || clip_at.abs_diff(end) <= 25 {
                        sig.clips_at_breakpoints += 1;
                    } else {
                        sig.clips_elsewhere += 1;
                    }
                }
                _ => {}
            }
        }

        // Pair geometry, counted once per pair from the leftmost mate (TLEN > 0).
        let tlen = record.template_length().map(|t| t as i64).unwrap_or(0);
        if tlen > 0 {
            let event_len = end.saturating_sub(pos);
            if tlen as usize > 250 + event_len / 2 {
                sig.long_pairs += 1;
            }
            if record.flags().is_ok_and(|f| f.is_reverse_complemented()) {
                sig.everted_pairs += 1;
            }
        }
    }

    let mean = |lo: usize, hi: usize| -> f64 {
        let slice = &depth[lo.min(contig_len)..hi.min(contig_len)];
        if slice.is_empty() {
            return 0.0;
        }
        slice.iter().sum::<usize>() as f64 / slice.len() as f64
    };
    sig.depth_inside = mean(pos + 20, end.saturating_sub(20));
    sig.depth_outside = mean(FLANK_WINDOW.0, FLANK_WINDOW.1);
    sig
}

/// Run a gate: generate with and without the SV, align both, and return `(with_sv, control)`.
pub fn run_gate(sv: &SvSpec) -> (Signatures, Signatures, tempfile::TempDir) {
    let bwa = bwa_mem2();
    let (dir, work) = fresh_workdir();

    let (vr1, vr2) = generate_reads(&work, "withsv", Some(sv));
    let (cr1, cr2) = generate_reads(&work, "control", None);

    let with_sv = analyse(
        &align(&bwa, &work, "withsv", &vr1, &vr2),
        sv.contig,
        sv.pos,
        sv.end,
    );
    let control = analyse(
        &align(&bwa, &work, "control", &cr1, &cr2),
        sv.contig,
        sv.pos,
        sv.end,
    );
    (with_sv, control, dir)
}
