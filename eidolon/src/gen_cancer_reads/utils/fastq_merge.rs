//! Tag per-pass read names and concatenate FASTQs (replaces the `awk` step in
//! `tools/cancer_simulate.sh`).
//!
//! gen-reads names reads by genomic position (`EIDOLON_generated_<contig>_<start>_
//! <end>/<mate>`), so a normal and a tumor read sampled at the same coordinate
//! get IDENTICAL names. On chr22 at 30x combined we saw ~456k collisions out of
//! ~9.6M reads — which Picard MarkDuplicates silently drops (halving effective
//! coverage at those loci). Prefixing each pass's read names (`N_`/`T_`) fixes it.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::gen_cancer_reads::errors::GenCancerReadsError;
use eidolon_core::file_tools::block_gz::BlockGzWriter;

/// Stream each `(path, tag)` input, prefix the header line of every 4-line FASTQ
/// record with `<tag>_`, and write the concatenation to `out` (gzip).
pub fn tag_and_concat(inputs: &[(&Path, &str)], out: &Path) -> Result<(), GenCancerReadsError> {
    let f = File::create(out)?;
    let mut w = BlockGzWriter::new(BufWriter::new(f));
    for (path, tag) in inputs {
        let r = BufReader::new(MultiGzDecoder::new(File::open(path)?));
        for (i, line) in r.lines().enumerate() {
            let line = line?;
            // Record header is line 1 of every 4 (NR % 4 == 1). Quality lines
            // (line 4) can also start with '@', so key on position, not char.
            if i % 4 == 0 {
                match line.strip_prefix('@') {
                    Some(rest) => writeln!(w, "@{tag}_{rest}")?,
                    None => writeln!(w, "{line}")?, // defensive passthrough
                }
            } else {
                writeln!(w, "{line}")?;
            }
        }
    }
    w.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Read;

    fn write_gz(path: &Path, body: &str) {
        let f = File::create(path).unwrap();
        let mut w = GzEncoder::new(f, Compression::default());
        w.write_all(body.as_bytes()).unwrap();
        w.finish().unwrap();
    }

    fn read_gz(path: &Path) -> String {
        let mut s = String::new();
        MultiGzDecoder::new(File::open(path).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        s
    }

    #[test]
    fn tag_and_concat_prefixes_headers_and_preserves_records() {
        let dir = tempfile::tempdir().unwrap();
        let n = dir.path().join("n.fastq.gz");
        let t = dir.path().join("t.fastq.gz");
        let out = dir.path().join("merged.fastq.gz");

        // Note the second record's QUALITY line is "@@@@" (Phred '@' = Q31).
        // A naive "prefix every line starting with '@'" would corrupt it — the
        // function must key on record position (line 1 of every 4), not the
        // leading character.
        write_gz(
            &n,
            "@EIDOLON_generated_chr1_100_170/1\nACGT\n+\nIIII\n\
             @EIDOLON_generated_chr1_200_270/1\nTTTT\n+\n@@@@\n",
        );
        write_gz(&t, "@EIDOLON_generated_chr1_100_170/1\nGGGG\n+\nIIII\n");

        tag_and_concat(&[(&n as &Path, "N"), (&t as &Path, "T")], &out).unwrap();
        let merged = read_gz(&out);
        let lines: Vec<&str> = merged.lines().collect();

        // 2 normal records + 1 tumor record = 12 lines.
        assert_eq!(lines.len(), 12, "merged record count");

        // Headers (line 1 of each record) are tag-prefixed, normal block first.
        assert_eq!(lines[0], "@N_EIDOLON_generated_chr1_100_170/1");
        assert_eq!(lines[4], "@N_EIDOLON_generated_chr1_200_270/1");
        assert_eq!(lines[8], "@T_EIDOLON_generated_chr1_100_170/1");

        // Sequence / separator / quality lines pass through untouched — crucially
        // the "@@@@" quality line is NOT prefixed.
        assert_eq!(lines[1], "ACGT");
        assert_eq!(lines[6], "+");
        assert_eq!(
            lines[7], "@@@@",
            "quality line starting with '@' must not be tagged"
        );
        assert_eq!(lines[9], "GGGG");

        // The same coordinate (chr1_100_170) appears under both N_ and T_ — which
        // is the whole point: tagging prevents the QNAME collision.
        assert!(merged.contains("@N_EIDOLON_generated_chr1_100_170/1"));
        assert!(merged.contains("@T_EIDOLON_generated_chr1_100_170/1"));
    }
}
