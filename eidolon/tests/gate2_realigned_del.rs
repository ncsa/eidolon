//! **Gate 2 for DEL** — does a real aligner, given eidolon's FASTQ, produce the evidence a
//! structural-variant caller actually looks for?
//!
//! Every other SV test in this repo inspects eidolon's own output: the FASTQ bases, the golden
//! BAM, the truth VCF. All three are eidolon describing its own work. A caller never sees any of
//! them — it sees reads that someone else aligned. **Nothing has ever checked that step.** The
//! whole-genome campaigns jump straight from FASTQ to "did Manta find it", which conflates
//! "the evidence is correct" with "the evidence was sufficient", and cannot tell you which
//! failed when the answer is no.
//!
//! So this aligns eidolon's reads with **bwa-mem2** and asserts, on the resulting alignment, the
//! three signatures every read-based DEL caller keys on:
//!
//! | signature | why a caller needs it |
//! |---|---|
//! | depth collapses over the deleted interval | depth-based callers (CNVnator, GATK gCNV) |
//! | soft clips pile up at both breakpoints | split-read callers (Manta, Delly, LUMPY) |
//! | spanning pairs have TLEN inflated by ~SVLEN | discordant-pair callers (all of the above) |
//!
//! Each is asserted against a **no-variant control** built from the same reference, seed and
//! coverage, so "the signature is present" cannot be satisfied by background noise — which is
//! how a 13% artifact once read as a real DUP defect (see `docs/sv_support_matrix.md`).
//!
//! ## Why this test is `#[ignore]`d
//!
//! CI has no aligner. Rather than skip silently when `bwa-mem2` is missing — a pass that means
//! nothing, the exact shape this project keeps re-earning — the test **fails** with instructions
//! if the binary is absent, and is `#[ignore]`d so CI never reaches it. Run it deliberately:
//!
//! ```text
//! conda activate aln          # or: export BWA_MEM2=/path/to/bwa-mem2
//! cargo test --test gate2_realigned_del -- --ignored --nocapture
//! ```
//!
//! The analysis is done in Rust over the SAM with `noodles`, not by piping samtools through awk,
//! so it is debuggable and its arithmetic is visible.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::io::{BufRead, Write as _};
use std::path::{Path, PathBuf};
use std::process::Command;

/// The planted deletion. Homozygous so every fragment over the locus carries it — a het event
/// halves every signal below and would make the control comparison a ratio rather than a
/// presence/absence question.
const CONTIG: &str = "H1N1_HA";
const DEL_POS: usize = 500; // 1-based anchor; deleted bases are POS+1..=END
const DEL_END: usize = 799;
const DEL_LEN: usize = DEL_END - DEL_POS; // 299 deleted reference bases

/// Locate bwa-mem2, or fail with something actionable. Never returns a "skip".
fn bwa_mem2() -> String {
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
         \n    BWA_MEM2=/path/to/bwa-mem2 cargo test --test gate2_realigned_del -- --ignored\n"
    );
    "bwa-mem2".to_string()
}

fn run(cmd: &mut Command, what: &str) {
    let out = cmd
        .output()
        .unwrap_or_else(|e| panic!("{what}: failed to spawn: {e}"));
    assert!(
        out.status.success(),
        "{what} failed ({}):\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
}

/// Generate paired reads over H1N1, optionally with a homozygous symbolic DEL.
/// Same seed and coverage either way, so the control differs only by the variant.
fn generate_reads(work: &Path, tag: &str, with_del: bool) -> (PathBuf, PathBuf) {
    let mut config = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), tag);
    config.coverage = 60;
    config.read_len = 100;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.produce_bam = false;
    config.produce_vcf = true;
    // No de novo variants: the only difference between the two runs must be the planted DEL.
    config.mutation_rate = Some(0.0);
    config.sv_rate_scale = Some(0.0);

    if with_del {
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
            "{CONTIG}\t{DEL_POS}\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;END={DEL_END}\tGT\t1/1"
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

/// What a caller can see, extracted from a realigned SAM.
#[derive(Debug, Default)]
struct Signatures {
    /// Mean depth strictly inside the deleted interval.
    depth_inside: f64,
    /// Mean depth in a control window well clear of the event.
    depth_outside: f64,
    /// Soft-clipped bases whose clip point falls within ±25 bp of either breakpoint.
    clips_at_breakpoints: usize,
    /// Soft clips anywhere else — the background this must stand out from.
    clips_elsewhere: usize,
    /// Properly-paired reads whose |TLEN| exceeds the no-variant mode by ~DEL_LEN.
    discordant_pairs: usize,
    reads: usize,
}

/// Parse the SAM and accumulate the three signatures. Deliberately in Rust rather than
/// `samtools | awk`: the arithmetic below is the assertion, and it should be readable.
fn analyse(sam_path: &Path) -> Signatures {
    let mut reader = std::fs::File::open(sam_path)
        .map(std::io::BufReader::new)
        .map(sam::io::Reader::new)
        .unwrap();
    let header = reader.read_header().unwrap();

    // Depth is accumulated over the whole contig, then summarised over two windows.
    let contig_len = 1_800usize;
    let mut depth = vec![0usize; contig_len];
    let mut sig = Signatures::default();

    let _ = &header;
    for result in reader.records() {
        let record = result.unwrap();
        // H1N1 has EIGHT contigs. Without this filter every contig's reads land in one depth
        // array, the deleted window gets backfilled by unrelated contigs, and the deletion
        // vanishes into an apparent 1.2x ENRICHMENT. The control comparison is what caught it.
        let on_target = matches!(
            record.reference_sequence_name(),
            Some(n) if n == CONTIG.as_bytes()
        );
        if !on_target {
            continue;
        }
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        sig.reads += 1;
        let start = usize::from(start); // 1-based

        // Walk the CIGAR: M/=/X/D/N advance the reference, S records a clip point.
        let mut ref_pos = start;
        let mut first_op = true;
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
                    // A leading clip points at this read's start; a trailing clip at ref_pos.
                    let clip_at = if first_op && i == 0 { start } else { ref_pos };
                    let near = |bp: usize| clip_at.abs_diff(bp) <= 25;
                    if near(DEL_POS) || near(DEL_END) {
                        sig.clips_at_breakpoints += 1;
                    } else {
                        sig.clips_elsewhere += 1;
                    }
                }
                _ => {}
            }
            first_op = false;
        }

        // A fragment spanning the deletion aligns with its mate ~DEL_LEN further away than a
        // fragment that does not. Count only clearly-inflated pairs, once per pair.
        let tlen = record
            .template_length()
            .map(|t| t.unsigned_abs() as usize)
            .unwrap_or(0);
        if record.flags().is_ok_and(|f| f.is_first_segment()) && tlen > 250 + DEL_LEN / 2 {
            sig.discordant_pairs += 1;
        }
    }

    let mean = |lo: usize, hi: usize| -> f64 {
        let slice = &depth[lo.min(contig_len)..hi.min(contig_len)];
        if slice.is_empty() {
            return 0.0;
        }
        slice.iter().sum::<usize>() as f64 / slice.len() as f64
    };
    // Strictly interior: 20 bp in from each breakpoint, so a read clipped at the junction
    // cannot contribute and the window measures the deletion rather than its edges.
    sig.depth_inside = mean(DEL_POS + 20, DEL_END - 20);
    sig.depth_outside = mean(1_000, 1_400);
    sig
}

fn align(bwa: &str, work: &Path, tag: &str, r1: &Path, r2: &Path) -> PathBuf {
    // Index into the work dir so the repo's test_data is never written to.
    let local_ref = work.join("ref.fa");
    std::fs::copy(h1n1_reference(), &local_ref).unwrap();
    run(
        Command::new(bwa).arg("index").arg(&local_ref),
        "bwa-mem2 index",
    );

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

    // A SAM with a header but no alignments would sail through every assertion below as
    // "no signal", so establish that reads were actually placed.
    let mapped = std::io::BufReader::new(std::fs::File::open(&sam_path).unwrap())
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('@'))
        .count();
    assert!(
        mapped > 100,
        "{tag}: bwa-mem2 emitted only {mapped} alignment record(s) — the alignment step \
         failed, so nothing below would be measuring the deletion"
    );
    sam_path
}

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_del_produces_the_three_signatures_a_caller_needs() {
    let bwa = bwa_mem2();
    let (_dir, work) = fresh_workdir();

    let (dr1, dr2) = generate_reads(&work, "withdel", true);
    let (cr1, cr2) = generate_reads(&work, "control", false);

    let del = analyse(&align(&bwa, &work, "withdel", &dr1, &dr2));
    let ctl = analyse(&align(&bwa, &work, "control", &cr1, &cr2));

    println!("  DEL run: {del:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: depth collapses over the deletion ────────────────────────────────
    // Homozygous, so the interior should be ~empty rather than merely reduced. Asserted as a
    // ratio against the same run's own flanking depth, which cancels any coverage difference
    // between the two runs.
    let del_ratio = del.depth_inside / del.depth_outside.max(1e-9);
    let ctl_ratio = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    assert!(
        del_ratio < 0.10,
        "depth over a HOMOZYGOUS deletion should collapse; interior/flank was {del_ratio:.3} \
         (inside {:.1}x, outside {:.1}x)",
        del.depth_inside,
        del.depth_outside
    );
    // MUST-NOT-FIRE: the control has no deletion, so the same window must be normally covered.
    // Without this, a run that produced no reads over the region for an unrelated reason would
    // satisfy the assertion above.
    assert!(
        ctl_ratio > 0.80,
        "control has no deletion, yet the same window is depleted (ratio {ctl_ratio:.3}) — \
         the window itself is unreliable, so signature 1 proves nothing"
    );

    // ── Signature 2: split reads pile up at the breakpoints ───────────────────────────
    assert!(
        del.clips_at_breakpoints >= 5,
        "a split-read caller needs clipped reads at the junction; found {} within ±25bp of \
         {DEL_POS}/{DEL_END}",
        del.clips_at_breakpoints
    );
    // MUST-NOT-FIRE: clipping at the breakpoints must be specific to the deletion, not the
    // aligner's background rate on this fixture.
    assert!(
        del.clips_at_breakpoints > ctl.clips_at_breakpoints * 3,
        "breakpoint clipping ({}) is not clearly above the control's background ({}) — \
         a caller could not distinguish it either",
        del.clips_at_breakpoints,
        ctl.clips_at_breakpoints
    );

    // ── Signature 3: spanning pairs are discordant by ~SVLEN ──────────────────────────
    assert!(
        del.discordant_pairs >= 5,
        "a paired-end caller needs pairs whose insert size is inflated by the deletion; \
         found {}",
        del.discordant_pairs
    );
    assert!(
        del.discordant_pairs > ctl.discordant_pairs * 3,
        "discordant pairs ({}) are not clearly above the control's background ({})",
        del.discordant_pairs,
        ctl.discordant_pairs
    );
}
