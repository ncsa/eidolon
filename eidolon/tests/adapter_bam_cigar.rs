//! Adapter read-through must reach the BAM as a SOFT CLIP, not as an aligned match.
//!
//! Two gaps met here, both found by audit:
//!
//! 1. **No test in this repo read a CIGAR out of a produced BAM.** `produce_bam: true`
//!    appeared in exactly one integration test, which inspected only `template_length`.
//!    So the BAM's alignment geometry — the thing every downstream aligner and caller
//!    consumes — was asserted nowhere.
//! 2. **"adapter" appeared in zero files under `eidolon/tests/`**, despite #125 backing
//!    the largest single number in the ACCESS report's fix table (SNP recall
//!    0.0004 → 0.944).
//!
//! Between them they hid a real defect: `fastq_tools.rs` tagged adapter bases `'S'` and
//! `bam_writer.rs::char_to_cigar_kind` mapped everything except `I`/`D` to `Match`, so
//! adapter was written as aligned sequence. Each side was self-consistent and unit
//! tested; nothing asserted they agreed. A soft clip does not consume reference, so the
//! defect also overstated each read's reference span.
//!
//! The negative control matters as much as the positive one: with adapters disabled the
//! output must contain no soft clips at all, or a validator that stamped `S` on
//! everything would satisfy the positive assertion.

mod common;

use common::{eidolon, h1n1_reference};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::path::Path;

/// Per-record (soft-clipped bases, query-consuming bases) from a BAM.
fn clip_and_query_lengths(path: &Path) -> Vec<(usize, usize)> {
    let file = std::fs::File::open(path).expect("bam not produced");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().unwrap();
    reader
        .records()
        .map(|r| {
            let rec = r.unwrap();
            let (mut clipped, mut query) = (0usize, 0usize);
            for op in rec.cigar().iter() {
                let op = op.unwrap();
                match op.kind() {
                    Kind::SoftClip => {
                        clipped += op.len();
                        query += op.len();
                    }
                    Kind::Match
                    | Kind::Insertion
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch => query += op.len(),
                    _ => {}
                }
            }
            (clipped, query)
        })
        .collect()
}

/// Run gen-reads on the H1N1 fixture with a short-insert library.
/// `adapters` toggles 3' read-through; everything else is held constant.
fn run(dir: &Path, name: &str, adapters: bool) -> Vec<(usize, usize)> {
    let cfg_text = format!(
        "reference: {ref}\nread_len: 100\ncoverage: 20\nploidy: 1\npaired_ended: true\n\
         fragment_mean: 60\nfragment_st_dev: 5\n\
         produce_bam: true\nproduce_fastq: false\nproduce_vcf: false\n\
         overwrite_output: true\noutput_dir: {out}\noutput_filename: {name}\n\
         rng_seed: adapter cigar\nnum_threads: 1\n{ad}",
        ref = h1n1_reference().display(),
        out = dir.display(),
        ad = if adapters {
            "adapters:\n  enabled: true\n  preset: truseq\n"
        } else {
            "keep_short_fragments: true\n"
        },
    );
    let cfg = dir.join(format!("{name}.yml"));
    std::fs::write(&cfg, cfg_text).unwrap();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    clip_and_query_lengths(&dir.join(format!("{name}.bam")))
}

#[test]
fn adapter_readthrough_is_soft_clipped_in_the_bam() {
    let tmp = tempfile::tempdir().unwrap();
    let recs = run(tmp.path(), "withad", true);
    assert!(
        !recs.is_empty(),
        "no BAM records — test would vacuously pass"
    );

    let clipped: Vec<_> = recs.iter().filter(|(c, _)| *c > 0).collect();
    assert!(
        !clipped.is_empty(),
        "adapters enabled with a 60 bp mean insert and 100 bp reads, but not one record \
         carries a soft clip. Adapter bases are being written as aligned sequence — the \
         char_to_cigar_kind defect."
    );

    // Known answer, independent of the implementation: read_len is fixed at 100, and
    // adapter fills exactly the bases past the insert, so every read is 100 query bases
    // and the clip is whatever is not genomic. A record clipped to its full length would
    // mean the genomic piece vanished.
    for (clip, query) in &recs {
        assert_eq!(
            *query, 100,
            "read_len is 100, so query-consuming ops must total 100; got {query}"
        );
        assert!(
            *clip < 100,
            "a read clipped over its entire length has no genomic anchor"
        );
    }
}

#[test]
fn without_adapters_nothing_is_soft_clipped() {
    // The must-not-fire case. `keep_short_fragments` gives the same short-insert library
    // WITHOUT adapter read-through, so the only difference is the feature under test.
    let tmp = tempfile::tempdir().unwrap();
    let recs = run(tmp.path(), "noad", false);
    assert!(
        !recs.is_empty(),
        "no BAM records — test would vacuously pass"
    );
    let total_clipped: usize = recs.iter().map(|(c, _)| *c).sum();
    assert_eq!(
        total_clipped, 0,
        "adapters are disabled, so no base should be soft-clipped; found {total_clipped}"
    );
}
